#ifndef DUNE_COPASI_MODEL_DIFFUSION_REACTION_CC
#define DUNE_COPASI_MODEL_DIFFUSION_REACTION_CC

/**
 * This is the implementation of the ModelDiffusionReaction,
 * particularly, you want to notice that this is not an normal .cc
 * file but a header which has to be included when compiling.
 */

#include <dune/copasi/common/bit_flags.hh>
#include <dune/copasi/common/data_context.hh>
#include <dune/copasi/common/tiff_grayscale.hh>
#include <dune/copasi/common/parser_to_grid_function.hh>
// #include <dune/copasi/concepts/pdelab.hh>
// #include <dune/copasi/model/diffusion_reaction.hh>

#include <dune/assembler/common/trace.hh>
#include <dune/assembler/common/error_condition.hh>
#include <dune/assembler/concepts/multiindex.hh>
#include <dune/assembler/operator/inverse/differentiable_adapter.hh>
#include <dune/assembler/assembler/mass_stiffness/apply.hh>
#include <dune/assembler/assembler/mass_stiffness/jacobian.hh>
#include <dune/assembler/istl/linear_adapter.hh>
#include <dune/assembler/newton/newton.hh>
#include <dune/assembler/discrete_function_space/discrete_function.hh>

#include <dune/pdelab/common/function.hh>

#include <dune/istl/scalarproducts.hh>
#include <dune/istl/solvers.hh>
#include <dune/istl/umfpack.hh>
#include <dune/istl/preconditioners.hh>

#include <dune/localfunctions/common/interfaceswitch.hh>

#include <dune/grid/io/file/vtk.hh>

// #include <string>

// #include <sys/stat.h>

namespace Dune::Copasi {



  template<class M, class X, class Y, std::size_t l=1>
  class BlockJacobi : public Preconditioner<X,Y> {
  public:
    // //! \brief The matrix type the preconditioner is for.
    // typedef M matrix_type;
    // //! \brief The domain type of the preconditioner.
    // typedef X domain_type;
    // //! \brief The range type of the preconditioner.
    // typedef Y range_type;
    // //! \brief The field type of the preconditioner.
    // typedef typename X::field_type field_type;
    // //! \brief scalar type underlying the field_type
    // typedef Simd::Scalar<field_type> scalar_field_type;
    // //! \brief real scalar type underlying the field_type
    // typedef typename FieldTraits<scalar_field_type>::real_type real_field_type;

    /*! \brief Constructor.

       Constructor gets all parameters to operate the prec.
       \param A The matrix to operate on.
       \param n The number of iterations to perform.
       \param w The relaxation factor.
     */
    BlockJacobi (const M& A/*, int n*//*, real_field_type w*/)
      : _mat(A)/*, _n(n)*//*, _w(w)*/
    {
      // CheckIfDiagonalPresent<M,l>::check(_mat);
    }

    // BlockJacobi (const std::shared_ptr<const AssembledLinearOperator<M,X,Y>>& A, const ParameterTree& configuration)
    //   : BlockJacobi(A->getmat(), configuration)
    // {}

    // /*!
    //    \brief Constructor.

    //    \param A The matrix to operate on.
    //    \param configuration ParameterTree containing preconditioner parameters.

    //    ParameterTree Key | Meaning
    //    ------------------|------------
    //    iterations        | The number of iterations to perform. default=1
    //    relaxation        | The relaxation factor. default=1.0

    //    See \ref ISTL_Factory for the ParameterTree layout and examples.
    //  */
    // BlockJacobi (const M& A, const ParameterTree& configuration)
    //   : BlockJacobi(A, configuration.get<int>("iterations",1), configuration.get<real_field_type>("relaxation",1.0))
    // {}

    /*!
       \brief Prepare the preconditioner.

       \copydoc Preconditioner::pre(X&,Y&)
     */
    virtual void pre ([[maybe_unused]] X& x, [[maybe_unused]] Y& b)
    {}

    /*!
       \brief Apply the preconditioner.

       \copydoc Preconditioner::apply(X&,const Y&)
     */
    virtual void apply (X& x, const Y& b)
    {
      Dune::SeqJac<M,X,Y,l> bjac(_mat, 10, 0.5);
      bjac.apply(x,b);
      // auto rhs = b;
      // auto v = x;
      // solve(_mat, rhs, v, std::integral_constant<std::size_t, l>{});
      // x = v;
    }

    static void solve(const auto& mat, auto& x, auto& b, auto depth) {
      // if constexpr (depth == 0) {
      //   x = b/mat;
      //   // mat.solve(x,b);
      // } else {
      //   double weight = 1.;
      //   for (auto i = mat.begin(); i != mat.end(); ++i) {
      //     auto j = (*i).begin();

      //     auto& lx = x[j.index()];
      //     auto& lb = b[i.index()];

      //     using Assembler::linearTransformation;
      //     for (; j.index() != i.index(); ++j) {
      //       linearTransformation(lb, *j, lx, [](auto& y, auto a, auto x){y -= a*x;});
      //     }
      //     auto& diag=*j;
      //     for (; j != (*i).begin(); ++j)
      //       linearTransformation(lb, *j, lx, [](auto& y, auto a, auto x){y -= a*x;});

      //     solve(diag, lx, lb, std::integral_constant<std::size_t,depth-1>{});
      //   }
      //   using Assembler::axpy;
      //   axpy(x, weight, b);
      // }
    }


    /*!
       \brief Clean up.

       \copydoc Preconditioner::post(X&)
     */
    virtual void post ([[maybe_unused]] X& x)
    {}

    //! Category of the preconditioner (see SolverCategory::Category)
    virtual SolverCategory::Category category() const
    {
      return SolverCategory::sequential;
    }

  private:
    //! \brief The matrix we operate on.
    const M& _mat;
    //! \brief The number of steps to perform during apply.
    // int _n;
    // //! \brief The relaxation parameter to use.
    // real_field_type _w;
  };





template<class Traits>
ModelDiffusionReaction<Traits>::ModelDiffusionReaction(
  std::shared_ptr<Grid> grid,
  const Dune::ParameterTree& config,
  const CompartmentEntitySet& compartment_entity_set,
  BitFlags<ModelSetup::Stages> setup_policy)
  : _config(config)
  , _logger(Logging::Logging::componentLogger({}, "model"))
{
  if (_config.sub("compartments",true).getValueKeys().size() != 1)
    DUNE_THROW(IOError, "'compartments' section must contain one entry");

  _state.grid_storage() = std::move(grid);

  setup(setup_policy, compartment_entity_set);

  _logger.trace("ModelDiffusionReaction constructed"_fmt);
}

template<class Traits>
ModelDiffusionReaction<Traits>::~ModelDiffusionReaction()
{
  _logger.trace("ModelDiffusionReaction deconstructed"_fmt);
}

template<class Traits>
template<class GF>
void
ModelDiffusionReaction<Traits>::interpolate(State& state,
  const std::map<std::string, GF>& initial)
{
  _logger.debug("Set initial state from grid functions"_fmt);

  if (state.space().entitySet().size(0) == 0)
    return;

  Assembler::LocalContainerBuffer<Space, Coefficients> lcontainer{state.space(), state.coefficients()};

  std::vector<GF const *> initial_ptr(state.space().degree(), nullptr);

  std::size_t count = 0;
  for (std::size_t comp = 0; comp != state.space().degree(); ++comp) {
    auto space_comp = state.space().subSpace(Assembler::multiIndex(comp));
    auto it = initial.find(space_comp.name());
    if (it != end(initial)) {
      initial_ptr[comp] = &(it->second);
      ++count;
    }
  }
  if (count != initial.size())
    DUNE_THROW(IOError, "Map with initial conditions contains more components than there is functions");

  // loop once over the grid and interpolate
  auto lspace = state.space().localView();
  for (const auto& entity : elements(state.space().entitySet())) {
    // bind local function space to element
    lspace.bind(entity);
    lcontainer.clear(lspace);

    for (std::size_t comp = 0; comp != lspace.tree().degree(); ++comp) {
      if (initial_ptr[comp] == nullptr)
        continue;
      auto path = lspace.tree().child(comp).path();
      const auto& fe = lspace.tree().child(comp).finiteElement();
      using FiniteElement = std::remove_cvref_t<decltype(fe)>;
      using LocalBasis = typename FiniteElement::Traits::LocalBasisType;
      using FEDomain = typename LocalBasis::Traits::DomainType;
      using FERangeField = typename LocalBasis::Traits::RangeFieldType;
      auto lfunc = [&](const FEDomain& xlocal){
        FERangeField y;
        initial_ptr[comp]->evaluate(entity, xlocal, y);
        return SpeciesQuantity{y};
      };
      fe.localInterpolation().interpolate(lfunc, lcontainer[path]);
    }

    lcontainer.store(lspace);
    lspace.unbind();
  }
}

template<class Traits>
auto
ModelDiffusionReaction<Traits>::make_component_grid_function_space(
  const CompartmentEntitySet& entity_set,
  const std::string& name) const
{
  _logger.trace("Create a finite element map"_fmt);

  // create entity mapper context
  using Entity = typename Traits::Grid::LeafGridView::template Codim<0>::Entity;
  using Index = std::size_t;
  using EntityMapper = std::function<Index(Entity)>;
  EntityMapper&& em = [](const auto entity) {
    return entity.geometry().type().isCube() ? 0 : 1;
  };

  // get common geometry type on gridview
  if (not has_single_geometry_type(entity_set))
    DUNE_THROW(InvalidStateException,
               "Grid view has to have only one geometry type");
  GeometryType&& gt = entity_set.template begin<0>()->geometry().type();

  // create data context with entity mapper, geometry type and grid view
  auto&& ctx = Context::data_context(std::move(em), std::move(gt), entity_set);

  // create fem from factory
  std::shared_ptr<SpeciesFiniteElementMap> finite_element_map(Factory<SpeciesFiniteElementMap>::create(std::move(ctx)));

  auto constraints = std::make_shared<Dune::PDELab::NoConstraints>();

  _logger.trace("Setup grid function space for component {}"_fmt, name);
  auto comp_gfs = SpeciesSpace{SpeciesMergingStrategy{entity_set}, finite_element_map, constraints};
  comp_gfs.name(name);
  return comp_gfs;
}

template<class Traits>
auto
ModelDiffusionReaction<Traits>::make_compartment_function_space(
  const CompartmentEntitySet& entity_set
) const
{
  _logger.debug("Setup compartment grid function space"_fmt);
  const auto& compartment_names = _config.sub("compartments",true).getValueKeys();
  assert(compartment_names.size() == 1);
  auto components = _config.sub(compartment_names.front() + ".reaction").getValueKeys();

  std::vector<SpeciesSpace> comp_gfs_vec;
  for (const auto& name : components)
    comp_gfs_vec.push_back(make_component_grid_function_space(entity_set, name));

  CompartmentSpace compartment_space = makeCompositeDiscreteFunctionSpace(CompartmentMergingStrategy{entity_set}, comp_gfs_vec);
  compartment_space.name(compartment_names.front());

  _logger.info("No. of components on '{}' compartment: {}"_fmt, compartment_space.name(), compartment_space.degree());

  return compartment_space;
}

template<class Traits>
void
ModelDiffusionReaction<Traits>::setup_grid_function_space(
  State& state,
  const CompartmentEntitySet& entity_set
)
{
  TRACE_EVENT("dune", "Space::SetUp");
  if (not _state)
    state.time_point() = TimeQuantity{_config.get("time_stepping.begin", 0.)};
  auto comp_space = make_compartment_function_space(entity_set);
  state.set_space(Space{makeOrderedSpace(comp_space, entity_set)});
}

template<class Traits>
void
ModelDiffusionReaction<Traits>::setup_coefficient_vector(State& state)
{
  _logger.debug("Setup coefficient vector"_fmt);
  state.set_coefficients(state.space().makeContainer(CoefficientsBackend{}));
  state.space().resize(state.coefficients());
  Assembler::forEachContainerEntry(state.coefficients(), []<class T>(T& v, auto p){v = T{0.};});
}

template<class Traits>
void
ModelDiffusionReaction<Traits>::setup_initial_condition(State& state)
{
  TRACE_EVENT("dune", "InitialCondition");
  assert(_state);
  // get TIFF data if available
  const auto& tiff_config = _config.sub("data");
  std::map<std::string, std::shared_ptr<TIFFGrayscale<unsigned short>>> tiffs;
  for (const auto& tiff_key : tiff_config.getValueKeys())
    tiffs[tiff_key] = std::make_unique<TIFFGrayscale<unsigned short>>(tiff_config[tiff_key]);

  using GridFunction = ParserToGridFunctionAdapter<CompartmentEntitySet, TimeQuantity>;
  std::map<std::string, GridFunction> functions;

  for (std::size_t comp = 0; comp != state.space().degree(); ++comp) {
    auto space_comp = state.space().subSpace(Assembler::multiIndex(comp));
    auto [it, inserted] = functions.try_emplace(space_comp.name(), state.space().entitySet(), make_parser());
    if (not inserted)
      DUNE_THROW(IOError, "Two or more components share the same name");
    auto& function = it->second;
    std::string expr = _config[fmt::format("{}.initial.{}",  state.space().name(), space_comp.name())];
    function.parser().set_expression(expr);
    for (const auto& [name, tiff] : tiffs)
      function.parser().define_function(name,
        [=](const auto& x, const auto& y) {
          return std::invoke(*tiff, x, y);
        });
    function.setTime(state.time_point());
    function.parser().compile();
  }

  // evaluate compiled expressions as initial conditions
  interpolate(state, functions);
}

template<class Traits>
std::unique_ptr<typename ModelDiffusionReaction<Traits>::StepOperator>
ModelDiffusionReaction<Traits>::setup_step_operator(const ConstState& state) const
{
  TRACE_EVENT("dune", "StepOperator::SetUp");
  _logger.debug("Setup local operator"_fmt);
  const auto& domain_config = _config.sub(state.space().name(), true);

  //! Local operator
  using StiffnessLocalOperator = LocalOperatorDiffusionReactionCG<
    Space,
    typename SpeciesFiniteElementMap::Traits::FiniteElement::Traits::LocalBasisType::Traits>;

  _logger.trace("Create spatial local operator"_fmt);
  StiffnessLocalOperator stiff_lop{state.space(), domain_config};

  //! Temporal local operator
  using MassLocalOperator = TemporalLocalOperatorDiffusionReactionCG<
    Space,
    typename SpeciesFiniteElementMap::Traits::FiniteElement::Traits::LocalBasisType::Traits>;

  _logger.trace("Create temporal local operator"_fmt);
  MassLocalOperator mass_lop{state.space(), domain_config};

  _logger.debug("Create temporal operator"_fmt);

  using Weight = double;
  using RungeKutta = Assembler::RungeKuttaOperator<Coefficients, Residual, Weight, TimeQuantity, TimeQuantity>;
  auto rk_op = std::make_unique<RungeKutta>();

  using RKResidual = typename RungeKutta::InverseOperator::domain_type;
  using RKCoefficients = typename RungeKutta::InverseOperator::range_type;

  bool linear = true;
  // linear if all jacobians are zero
  const auto& jac_config = domain_config.sub("reaction.jacobian");
  for (const auto& jac_key : jac_config.getValueKeys()) {
    auto jac_expr = jac_config[jac_key];
    linear &= (jac_expr == "0") or (jac_expr == "0.0") or (jac_expr == ".0") or (jac_expr == "0.");
  }

  bool matrix_free = false;
  std::string solver = "RestartedGMRes";

  using Jacobian = decltype([](){
    if constexpr (CompartmentMergingStrategy::Blocked)
      return DynamicMatrix<BCRSMatrix<BCRSMatrix<double>>>{};
    else
      return DynamicMatrix<BCRSMatrix<double>>{};
  }());

  std::shared_ptr<Assembler::InverseDifferentiableOperator<RKCoefficients,RKResidual>> inverse_op;

  auto linear_solver_apply_op = [solver](Assembler::ForwardOperator<RKCoefficients,RKResidual>& derivative, double reduction, RKResidual& b, RKCoefficients& x) mutable {
    TRACE_EVENT("dune", "LinearSolver");
    static_assert(std::is_same_v<Coefficients,Residual>);
    auto istl_forward_op = std::make_shared<Dune::Assembler::ISTL::LinearAdapter<RKCoefficients,RKResidual>>([&derivative,&b](const RKCoefficients& coeff, RKResidual& residual){
      return derivative.apply(coeff, residual);
    });

    auto scalar_product_op = std::make_shared<Dune::ScalarProduct<RKResidual>>();
    std::shared_ptr<Preconditioner<RKCoefficients,RKResidual>> pre_op;
    pre_op = std::make_shared<Dune::Richardson<RKCoefficients,RKResidual>>(0.1);

    // try to get the jacobian...
    using JacobianOperator = Assembler::AssembledLinearOperator<Jacobian,RKCoefficients,RKResidual>;
    auto assembled_derivative = dynamic_cast<JacobianOperator const *>(&derivative);
    if (assembled_derivative) {
      // pre_op = std::make_shared<Dune::SeqSSOR<Jacobian,RKCoefficients,RKResidual,2>>(assembled_derivative->matrix(), 5, 1);
    }

    int verbosity = 4;
    int max_it = 100;
    Dune::InverseOperatorResult result;
    if (solver == "BiCGSTAB") {
      auto solver = Dune::BiCGSTABSolver<RKCoefficients>{istl_forward_op, scalar_product_op, pre_op, reduction, int(max_it), verbosity};
      solver.apply(x, b, result);
      DUNE_THROW(NotImplemented, "");
    } else if (solver == "RestartedGMRes") {
      int restart = 20;
      auto solver = Dune::RestartedGMResSolver<RKCoefficients, RKResidual>{istl_forward_op, scalar_product_op, pre_op, reduction, restart, int(max_it), verbosity};
      solver.apply(x, b, result);
    } else {
      DUNE_THROW(IOError, "Not known linear solver");
    }

    if (result.converged)
      return Dune::Assembler::ErrorCondition{};
    else
      return make_error_condition(Dune::Assembler::ConvergenceReason::DivergedNull);
  };


  if (linear) {
    auto linear_solve_op = std::make_unique<Assembler::InverseDifferentiableAdapter<RKCoefficients,RKResidual>>();

    std::optional<RKCoefficients> coeff_zero;
    linear_solve_op->setSolverApply([coeff_zero, linear_solver_apply_op](Assembler::DifferentiableOperator<RKCoefficients,RKResidual>& diff_op, RKResidual& b, RKCoefficients& x) mutable {
      TRACE_EVENT("dune", "LinearSolver::OneStepAdapter");
      static_assert(std::is_same_v<Coefficients,Residual>);
      auto derivative = diff_op.derivative(x, true, false);

      diff_op.apply(x, b).or_throw(); // residual is additive b += F(x)

      if (not coeff_zero) {
        coeff_zero.emplace(x);
        Assembler::forEachContainerEntry(*coeff_zero, []<class T>(T& v, auto){v = T{0};});
      }

      RKCoefficients z = *coeff_zero;

      double reduction{0.0001};
      // compute correction
      auto error_condition = linear_solver_apply_op(*derivative, reduction, b, z);
      if (error_condition)
        return error_condition;
      using Assembler::axpy;
      axpy(x,-1.,z);
      return Dune::Assembler::ErrorCondition{};
    });

    inverse_op = std::move(linear_solve_op);
  } else {
    auto newton_op = std::make_unique<Assembler::NewtonMethod<RKCoefficients,RKResidual, ResidualQuantity>>();

    newton_op->setLinearSolverApply(linear_solver_apply_op);
    newton_op->setNorm([](const auto& residual) -> ResidualQuantity {
      auto acc = ResidualQuantity{0}*ResidualQuantity{0};
      Assembler::forEachContainerEntry(residual, [&](auto& v, auto p){acc = v*v;});
      using std::sqrt;
      return sqrt(acc);
    });

    inverse_op = std::move(newton_op);
  }

  if (matrix_free) { // matrix free
    using MassStiffnesOperator = Assembler::MassStiffnessApplyOperator<RKCoefficients, RKResidual, Space, Space, MassLocalOperator, StiffnessLocalOperator, Weight, TimeQuantity, TimeQuantity, Assembler::DurationPosition::StiffnessNumerator>;
    auto diff_op = std::make_shared<MassStiffnesOperator>(state.space(), state.space(), mass_lop, stiff_lop);
    diff_op->setTimePoint(state.time_point());
    inverse_op->setDifferentiableOperator(diff_op);
  } else { // matrix based
    auto pattern_factory = [stiff_lop, mass_lop](const Space& test, const Space& trial, Jacobian& jac){
      auto& entry = jac[0][0];

      auto size_per_row = []<class MultiIndex>(MultiIndex i ,MultiIndex j) -> std::size_t {
        if constexpr (CompartmentMergingStrategy::Blocked) {
          assert(i.size() == j.size());
          assert(i.size() < 2);
          // TODO: put actual values from spaces
          if (i.size() == 0)
            return 5; // number of edges conecting to a vertex
          if (i.size() == 1)
            return 10; // number of component per vertex (P1)
          return 0;
        } else {
          return 5;
        }
      };

      auto pattern = [&](){
        if constexpr (CompartmentMergingStrategy::Blocked)
          return Assembler::BlockedSparsePattern<Assembler::LeafSparsePattern<Space, Space>>{test, {}, trial, {}, size_per_row};
        else
          return Assembler::LeafSparsePattern<Space, Space>{test, {}, trial, {}, size_per_row};
      }();

      spaceToPattern(stiff_lop, pattern);
      spaceToPattern(mass_lop, pattern);
      pattern.sort();
      patternToMatrix(pattern, entry);
      entry = 0.;

      // copy patterns into other runge-kutta jacobians
      for (std::size_t i = 0; i != jac.N(); ++i) {
        for (std::size_t j = 0; j != jac.M(); ++j) {
          if (i == 0 and j == 0) continue;
          jac[i][j] = entry;
        }
      }

      std::ofstream file{"mat-sd.svg"};
      writeSVGMatrix(jac, file);
    };

    using MassStiffnesOperator = Assembler::MassStiffnessOperator<RKCoefficients, RKResidual, Jacobian, Space, Space, MassLocalOperator, StiffnessLocalOperator, Weight, TimeQuantity, TimeQuantity, Assembler::DurationPosition::StiffnessNumerator>;
    auto diff_op = std::make_shared<MassStiffnesOperator>(state.space(), state.space(), mass_lop, stiff_lop, pattern_factory);
    diff_op->setTimePoint(state.time_point());
    inverse_op->setDifferentiableOperator(diff_op);
  }

  Assembler::ShuOsherTableau<Weight> rk_method = Assembler::alexander2Method();
  rk_op->setTableau(rk_method);
  rk_op->setInverseOperator(inverse_op);
  auto residual = state.space().makeContainer(CoefficientsBackend{});
  state.space().resize(residual);
  Assembler::forEachContainerEntry(residual, []<class T>(T& v, auto p){v = T{0.};});
  rk_op->setInitialResidual(std::move(residual));
  return rk_op;

}

template<class Traits>
void
ModelDiffusionReaction<Traits>::setup_vtk_writer(State& state)
{
  _logger.debug("Setup VTK writer"_fmt);
  // this map contains the written timesteps for a given path and is shared
  // among different `state`s copied from this one.
  auto mapu_out =
    std::make_shared<std::map<std::string, std::vector<double>>>();
  state.writer = [mapu_out](
                    const auto& state, const auto& path, bool append) mutable {
    TRACE_EVENT("dune", "Write::VTK");
    auto log_writer = Logging::Logging::componentLogger({}, "writer");

    // create directory if necessary
    auto path_entry = fs::directory_entry{ path };
    if (not path_entry.exists()) {
      log_writer.info("Creating output directory '{}'"_fmt,
                      path_entry.path().string());
      std::error_code ec{ 0, std::generic_category() };
      fs::create_directories(path_entry.path(), ec);
      if (ec)
        DUNE_THROW(IOError,
                   "\n Category: " << ec.category().name() << '\n'
                                   << "Value: " << ec.value() << '\n'
                                   << "Message: " << ec.message() << '\n');
    }

    // Recover old timestesps in case something was written before
    auto& timesteps = (*mapu_out)[path.string()];
    std::string name = fmt::format("{}-{}", path.filename().string(), state.space().name());
    if (not append) {
      timesteps.clear();
      log_writer.detail("Creating a time sequence file: '{}.pvd'"_fmt, name);
    } else {
      log_writer.trace("Overriding time sequence file: '{}.pvd'"_fmt, name);
    }

    // setup writer again with old timesteps if necessary
    auto writer = std::make_shared<VTKWriter<CompartmentEntitySet>>(state.space().entitySet(), Dune::VTK::conforming);
    auto sequential_writer = VTKSequenceWriter<CompartmentEntitySet>{ writer, name, path.string(), path.string() };
    sequential_writer.setTimeSteps(timesteps);

    for (std::size_t component = 0; component != state.space().degree(); ++component) {
      auto component_space = state.space().subSpace(Assembler::multiIndex(component));
      auto component_function = Assembler::ScalarDiscreteFunction{component_space, state.coefficients()};
      VTK::FieldInfo info{component_space.name(), VTK::FieldInfo::Type::scalar, 1};
      sequential_writer.vtkWriter()->addVertexData(component_function, info);
    }

#if HAVE_UNITS
    log_writer.detail("Writing solution for {:.2f}s time stamp"_fmt, state.time_point().number());
#else
    log_writer.detail("Writing solution for {:.2f}s time stamp"_fmt, state.time_point());
#endif
    log_writer.trace("Writing vtu file: '{0}/{0}-{1:0>5}.vtu'"_fmt,
                     name,
                     timesteps.size());
#if HAVE_UNITS
    auto time_point = state.time_point().number();
#else
    auto time_point = state.time_point();
#endif
    sequential_writer.write(time_point, Dune::VTK::base64);
    sequential_writer.vtkWriter()->clear();
    timesteps = sequential_writer.getTimeSteps();
  };
}

template<class Traits>
void
ModelDiffusionReaction<Traits>::setup(BitFlags<ModelSetup::Stages> setup_policy, const CompartmentEntitySet& entity_set)
{
  const auto& compartment_name = _config.sub("compartments",true).getValueKeys().front();
  _logger.detail("Setting up diffusion-reaction model for '{}' compartment"_fmt, compartment_name);

  try {
    if (setup_policy.test(ModelSetup::Stages::GridFunctionSpace))
      setup_grid_function_space(state(), entity_set);
    if (setup_policy.test(ModelSetup::Stages::CoefficientVector))
      setup_coefficient_vector(state());
    if (setup_policy.test(ModelSetup::Stages::InitialCondition))
      setup_initial_condition(state());
    if (setup_policy.test(ModelSetup::Stages::Writer))
      setup_vtk_writer(state());
  }
  catch(...)
  {
    // detailed report of the input paramter tree
    std::stringstream ss;
    _config.report(ss);
    _logger.error(2, "---- Diffusion-Reaction model parameter tree"_fmt);
    for (std::string line; std::getline(ss, line);)
      _logger.error(2, "{}"_fmt, line);
    _logger.error(2, "----"_fmt);
    throw;
  }
}

template<class Traits>
typename ModelDiffusionReaction<Traits>::SpeciesGridFunction
ModelDiffusionReaction<Traits>::as_function(const ConstState& state,
                                            Assembler::Concept::MultiIndex auto component_path)
{
  return Assembler::ScalarDiscreteFunction{state.grid_function_space->subSpace(component_path), state.coefficients};
}

template<class Traits>
typename ModelDiffusionReaction<Traits>::SpeciesGridFunction
ModelDiffusionReaction<Traits>::as_function(Assembler::Concept::MultiIndex auto component_path) const
{
  return ModelDiffusionReaction<Traits>::as_function(state(), component_path);
}

} // namespace Dune::Copasi

#endif // DUNE_COPASI_MODEL_DIFFUSION_REACTION_CC
