#ifndef DUNE_COPASI_MODEL_MULTIDOMAIN_DIFFUSION_REACTION_CC
#define DUNE_COPASI_MODEL_MULTIDOMAIN_DIFFUSION_REACTION_CC

/**
 * This is the implementation of the ModelMultiDomainDiffusionReaction,
 * particularly, you want to notice that this is not an normal .cc
 * file but a header which has to be included when compiling.
 */

#include <dune/copasi/model/multidomain_diffusion_reaction.hh>

#include <dune/assembler/common/trace.hh>

#include <random>

namespace Dune::Copasi {


// template<class X>
// class ScalarProduct : public Dune::ScalarProduct<X> {
// public:
//   //! export types, they come from the derived class
//   typedef X domain_type;
//   typedef typename X::field_type field_type;
//   typedef typename FieldTraits<field_type>::real_type real_type;

//   /*! \brief Dot product of two vectors.
//       It is assumed that the vectors are consistent on the interior+border
//       partition.
//     */
//   field_type dot (const X& x, const X& y) const override
//   {
//     return x[0][0].dot(y[0][0]) + x[0][1].dot(y[0][1])/0.001;
//   }

//   /*! \brief Norm of a right-hand side vector.
//       The vector must be consistent on the interior+border partition
//     */
//   real_type norm (const X& x) const override
//   {
//     return x[0][0].two_norm() + x[0][1].two_norm()/(0.001*0.001);
//   }

//   //! Category of the scalar product (see SolverCategory::Category)
//   virtual SolverCategory::Category category() const
//   {
//     return SolverCategory::sequential;
//   }

//   //! every abstract base class has a virtual destructor
//   virtual ~ScalarProduct () {}
// };





template<class Traits>
ModelMultiDomainDiffusionReaction<Traits>::ModelMultiDomainDiffusionReaction(
  std::shared_ptr<Grid> grid,
  const Dune::ParameterTree& config,
  BitFlags<ModelSetup::Stages> setup_policy)
  : _config(config)
  , _logger(Logging::Logging::componentLogger({}, "model"))
{
  _state.grid_storage() = std::move(grid);
  setup(setup_policy);
}

// template<class Traits>
// template<class GFGridView>
// auto
// ModelMultiDomainDiffusionReaction<Traits>::get_muparser_initial(
//   const ParameterTree& model_config,
//   const GFGridView& gf_grid_view,
//   bool compile)
// {
//   const auto& compartments = model_config.sub("compartments",true).getValueKeys();

//   using GridFunction = ExpressionToGridFunctionAdapter<EntitySet, SpeciesQuantity>;
//   std::vector<std::vector<std::shared_ptr<GridFunction>>> functions(compartments.size());

//   for (std::size_t i = 0; i < compartments.size(); i++) {
//     const std::string compartement = compartments[i];
//     std::size_t domain_id = model_config.sub("compartments",true).template get<std::size_t>(compartement);
//     assert(domain_id < compartments.size());
//     auto sub_model_config = model_config.sub(compartement,true);
//     functions[domain_id] = SubModel::get_muparser_initial(sub_model_config, gf_grid_view, compile);
//   }

//   return functions;
// }

template<class Traits>
template<class GF>
void
ModelMultiDomainDiffusionReaction<Traits>::interpolate(State& state,
  const std::map<std::array<std::string, 2>, GF>& initial)
{
  _logger.debug("Set initial state from grid functions"_fmt);

  if (state.space().entitySet().size(0) == 0)
    return;

  std::size_t count = 0;
  auto compartments_space = state.space().subSpace(compartments_path);
  std::vector<std::vector<GF const *>> domain_ptr(compartments_space.degree());
  for (std::size_t domain = 0; domain != compartments_space.degree(); ++domain) {
    auto domain_space = compartments_space.subSpace(Assembler::multiIndex(domain));
    domain_ptr[domain].assign(domain_space.degree(), nullptr);

    for (std::size_t comp = 0; comp != domain_space.degree(); ++comp) {
      auto space_comp = domain_space.subSpace(Assembler::multiIndex(comp));
      auto it = initial.find({domain_space.name(), space_comp.name()});
      if (it != end(initial)) {
        domain_ptr[domain][comp] = &(it->second);
        ++count;
      }
    }
  }

#if DUNE_COPASI_HAVE_MEMBRANE_SPACE
  std::vector<std::vector<GF const *>> membrane_ptr(_membrane_map.size());
  for (std::size_t membrane_id = 0; membrane_id != _membrane_map.size(); ++membrane_id) {
    auto [domain_i, domain_o] = _membrane_map[membrane_id];
    auto space_i = compartments_space.subSpace(Assembler::multiIndex(domain_i));
    auto space_o = compartments_space.subSpace(Assembler::multiIndex(domain_o));
    auto space_m = state.space().subSpace(Assembler::multiIndex(Indices::_1, membrane_id));
    membrane_ptr[membrane_id].assign(space_m.degree(), nullptr);

    for (std::size_t comp = 0; comp != space_m.degree(); ++comp) {
      auto space_comp = space_m.subSpace(Assembler::multiIndex(comp));
      auto it_io = initial.find({fmt::format("{}-{}", space_i.name(), space_o.name()), space_comp.name()});
      auto it_oi = initial.find({fmt::format("{}-{}", space_o.name(), space_i.name()), space_comp.name()});
      if ((domain_i != domain_o) and (it_io != end(initial)) and (it_oi != end(initial)))
        DUNE_THROW(IOError, fmt::format("Membrane component '{}' is encountred twice in initial conditions!", space_comp.name()));
      auto it = it_io != end(initial) ? it_io : it_oi;
      if (it != end(initial)) {
        membrane_ptr[membrane_id][comp] = &(it->second);
        ++count;
      }
    }
  }
  assert(std::is_sorted(begin(_membrane_map), end(_membrane_map)));
#endif
  if (count != initial.size())
    DUNE_THROW(IOError, "Map with initial conditions contains more components than there is functions");


  Assembler::LocalContainerBuffer<Space, Coefficients> lcontainer{state.space(), state.coefficients()};


  using SubDomainIndex = typename Grid::SubDomainIndex;
  // loop once over the grid and interpolate
  auto lspace = state.space().localView();
  const auto entity_set = state.space().entitySet();
  for (const auto& entity : elements(entity_set)) {
    // bind local function space to element
    lspace.bind(entity);
    lcontainer.clear(lspace);

    auto domain_set_i = entity_set.indexSet().subDomains(entity);
    assert(domain_set_i.size() == 1);
    SubDomainIndex domain_i = *(domain_set_i.begin());

    const auto& lspace_node_domain = TypeTree::child(lspace.tree(), join(compartments_path, Assembler::multiIndex(domain_i)));
    for (std::size_t comp = 0; comp != lspace_node_domain.degree(); ++comp) {
      if (domain_ptr[domain_i][comp] == nullptr)
        continue;
      auto path = lspace_node_domain.child(comp).path();
      const auto& fe = lspace_node_domain.child(comp).finiteElement();
      using FiniteElement = std::remove_cvref_t<decltype(fe)>;
      using LocalBasis = typename FiniteElement::Traits::LocalBasisType;
      using FEDomain = typename LocalBasis::Traits::DomainType;
      using FERange = typename LocalBasis::Traits::RangeType;
      auto lfunc = [&](const FEDomain& xlocal){
        FERange y;
        domain_ptr[domain_i][comp]->evaluate(entity, xlocal, y);
        return y;
      };
      fe.localInterpolation().interpolate(lfunc, lcontainer[path]);
    }

#if DUNE_COPASI_HAVE_MEMBRANE_SPACE
    for (const auto& is : intersections(entity_set, entity)) {
      SubDomainIndex domain_o;
      if (is.neighbor()) {
        auto domain_set_o = entity_set.indexSet().subDomains(is.outside());
        assert(domain_set_o.size() == 1);
        domain_o = *(domain_set_o.begin());
        if (domain_i == domain_o)
          continue;
      } else {
        domain_o = domain_i;
      }
      std::array<SubDomainIndex,2> domain_io = {std::min(domain_i, domain_o), std::max(domain_i, domain_o)};
      auto membrane_id_it = std::lower_bound(begin(_membrane_map), end(_membrane_map), domain_io);
      assert(membrane_id_it != end(_membrane_map));
      assert(domain_io == *membrane_id_it);
      std::size_t membrane_id = std::distance(begin(_membrane_map), membrane_id_it);

      const auto face_id = is.indexInInside();
      const auto& lspace_node_memb = lspace.tree().child(Indices::_1).child(membrane_id);
      for (std::size_t comp = 0; comp != lspace_node_memb.degree(); ++comp) {
        if (membrane_ptr[membrane_id][comp] == nullptr)
          continue;
        auto path = lspace_node_memb.child(comp).child(face_id).path();
        const auto& fe = lspace_node_memb.child(comp).child(face_id).finiteElement();
        using FiniteElement = std::remove_cvref_t<decltype(fe)>;
        using LocalBasis = typename FiniteElement::Traits::LocalBasisType;
        using FEDomain = typename LocalBasis::Traits::DomainType;
        using FERange = typename LocalBasis::Traits::RangeType;
        const auto& entity_m = entity.template subEntity<1>(face_id);
        const auto& grid_function = *membrane_ptr[membrane_id][comp];
        auto lfunc = [&](const FEDomain& xlocal){
          FERange y;
          grid_function.evaluate(entity_m, xlocal, y);
          return y;
        };
        fe.localInterpolation().interpolate(lfunc, lcontainer[path]);
      }
    }
#endif
    lcontainer.store(lspace);
    lspace.unbind();
  }
}

template<class Traits>
void
ModelMultiDomainDiffusionReaction<Traits>::setup_grid_function_space(State& state)
{
  _logger.debug("Setup grid function space"_fmt);

  const auto& compartments = _config.sub("compartments",true).getValueKeys();
  std::size_t domain_count = compartments.size();
  if (domain_count != state.grid().maxAssignedSubDomainIndex() + 1)
    DUNE_THROW(IOError, "Number of compartments do not match with sub-domains in the grid");
  std::vector<std::optional<CompartmentSpace>> gfs_map(domain_count);

  if (not state) {
    state.time_point() = _config.get("time_stepping.begin",0.);
  }

  for (std::size_t domain = 0; domain != domain_count; ++domain) {
    const std::string compartement = compartments[domain];
    std::size_t domain_id = _config.sub("compartments",true).template get<std::size_t>(compartement);
    if (domain_id >= domain_count)
      DUNE_THROW(IOError, "Compartment index is larger than number of sub-domains in the grid");

    // remove other compartments from sub model config, this will make the sub
    // model to only see one compartment tree
    auto sub_model_config = _config;
    sub_model_config.sub("compartments") = ParameterTree{};
    auto compartment_sect = "compartments." + compartement;
    sub_model_config[compartment_sect] = _config[compartment_sect];
    CompartmentEntitySet sub_grid_view = state.grid().subDomain(domain_id).leafGridView();

    _logger.trace("Create a sub model ({}) for compartment '{}'"_fmt, domain_id, compartement);
    auto sub_model = SubModel{state.grid_storage(), sub_model_config, sub_grid_view, ModelSetup::Stages::None};
    // extract grid function space from the sub model
    if (gfs_map[domain_id].has_value())
      DUNE_THROW(IOError, "Space for domain " << domain_id << " is already assigned");
    gfs_map[domain_id].emplace(sub_model.make_compartment_function_space(sub_grid_view));
  }

  std::vector<CompartmentSpace> gfs_vec;
  for (auto& gfs : gfs_map)
    gfs_vec.emplace_back(gfs.value());

  _logger.trace("Setup of multi-compartment grid function space"_fmt);

  MultiCompartmentSpace multi_compartment_space = makeCompositeDiscreteFunctionSpace(MultiCompartmentMergingStrategy{}, gfs_vec);

#if DUNE_COPASI_HAVE_MEMBRANE_SPACE
  _logger.debug("Setup membranes grid function space"_fmt);

  auto subDomain = [&](const auto& entity){
    auto domain_set = state.grid().leafGridView().indexSet().subDomains(entity);
    assert(domain_set.size() == 1);
    return *(domain_set.begin());
  };

  std::set<std::array<SubDomainIndex,2>> membrane_set;
  for (const auto& entity : elements(state.grid().leafGridView())) {
    for (const auto& intersection : intersections(state.grid().leafGridView(), entity)) {
      auto domain_i = subDomain(intersection.inside());
      
      if (intersection.neighbor()) {
        auto domain_o = subDomain(intersection.outside());
        membrane_set.insert({std::min(domain_i, domain_o), std::max(domain_i, domain_o)});
      } else if (intersection.boundary()) {
        membrane_set.insert({domain_i, domain_i});
      }
    }
  }

  _membrane_map.assign(begin(membrane_set), end(membrane_set));

  std::vector<MembraneSpace> mem_space;
  auto mem_entity_set = EntitySet{state.grid().leafGridView()};
  MembraneSpeciesMergingStrategy mem_species_ms{mem_entity_set};
  MembraneMergingStrategy        mem_ms{mem_entity_set};


  for (auto [domain_i, domain_o] : _membrane_map) {

    using FiniteElement = MembraneSpeciesFiniteElementMap::BaseFiniteElement;
    auto finite_element_map = std::make_shared<MembraneSpeciesFiniteElementMap>(
      mem_entity_set,
      SkeletonFiniteElementMap{mem_entity_set}, // TODO
      domain_i, domain_o,
      FiniteElement{} // TODO
    );

    auto constraints = std::make_shared<Dune::PDELab::NoConstraints>();

    const auto& space_i = multi_compartment_space.child(domain_i);
    const auto& space_o = multi_compartment_space.child(domain_o);

    std::vector<MembraneSpeciesSpace> mem_spe_gfs_vec;

    const auto& config_io = _config.sub(space_i.name() + ".boundary." + space_o.name());
    const auto& config_oi = _config.sub(space_o.name() + ".boundary." + space_i.name());

    const auto& mem_component_io = config_io.sub("reaction").getValueKeys();
    const auto& mem_component_oi = config_oi.sub("reaction").getValueKeys();

    std::set<std::string> mem_component(begin(mem_component_io), end(mem_component_io));
    for (const auto& component : mem_component_oi)
      if (mem_component.insert(component).second)
        DUNE_THROW(IOError, "Membrane component " << component << " is listed on both sides of the membrane");

    _logger.info("No. of components on '{}-{}' membreane: {}"_fmt, space_i.name(), space_o.name() , mem_component.size());

    for (const auto& name : mem_component) {
      mem_spe_gfs_vec.emplace_back(MembraneSpeciesSpace{mem_species_ms, finite_element_map, constraints});
      mem_spe_gfs_vec.back().name(name);
    }
    mem_space.emplace_back(makeCompositeDiscreteFunctionSpace(mem_ms, mem_spe_gfs_vec));
    mem_space.back().name(fmt::format("{}-{}", space_i.name(), space_o.name()));
  }

  MultiMembraneSpace mult_mem_space = makeCompositeDiscreteFunctionSpace(MultiMembraneMergingStrategy{}, mem_space);

  MultiCompartmentMembraneSpace _space = makeCompositeDiscreteFunctionSpace(MultiCompartmentMembraneMergingStrategy{}, std::tuple{multi_compartment_space, mult_mem_space});

  Space space = makeOrderedSpace(_space, EntitySet{state.grid().leafGridView()});
  state.set_space(std::move(space));
#else
  Space space = makeOrderedSpace(multi_compartment_space, EntitySet{state.grid().leafGridView()});
  state.set_space(std::move(space));
#endif
}

template<class Traits>
void
ModelMultiDomainDiffusionReaction<Traits>::setup_coefficient_vector(State& state)
{
  _logger.debug("Setup coefficient vector"_fmt);
  state.set_coefficients(state.space().makeContainer(CoefficientsBackend{}));
  state.space().resize(state.coefficients());
  Assembler::forEachContainerEntry(state.coefficients(), []<class T>(T& v, auto p){v = T{0.};});
}

template<class Traits>
void
ModelMultiDomainDiffusionReaction<Traits>::setup_initial_condition(State& state)
{
  TRACE_EVENT("dune", "InitialCondition");
  assert(state);
  // get TIFF data if available
  const auto& tiff_config = _config.sub("data");
  std::map<std::string, std::shared_ptr<TIFFGrayscale<unsigned short>>> tiffs;
  for (const auto& tiff_key : tiff_config.getValueKeys())
    tiffs[tiff_key] = std::make_unique<TIFFGrayscale<unsigned short>>(tiff_config[tiff_key]);

  using RF = double;
  EntitySet entity_set = state.space().entitySet();
  using GridFunction = ParserToGridFunctionAdapter<EntitySet, RF>;
  std::map<std::array<std::string, 2>, GridFunction> functions;

  auto compartments_space = state.space().subSpace(compartments_path);
  for (std::size_t domain_i = 0; domain_i != compartments_space.degree(); ++domain_i) {
    // domain functions
    auto space_i = compartments_space.subSpace(Assembler::multiIndex(domain_i));
    for (std::size_t comp = 0; comp != space_i.degree(); ++comp) {
      auto space_comp = space_i.subSpace(Assembler::multiIndex(comp));
      auto [it, inserted] = functions.try_emplace({space_i.name(), space_comp.name()}, entity_set, make_parser());
      if (not inserted)
        DUNE_THROW(IOError, "Two or more components share the same name");
      auto& function = it->second;
      std::string expr = _config[fmt::format("{}.initial.{}", space_i.name(), space_comp.name())];
      function.parser().set_expression(expr);
      for (const auto& [name, tiff] : tiffs)
        function.parser().define_function(name,
          [=](const RF& x, const RF& y) {
            return std::invoke(*tiff, x, y);
          });
      function.setTime(state.time_point());
      function.parser().compile();
    }
  }
#if DUNE_COPASI_HAVE_MEMBRANE_SPACE
    // membrane functions
  for (std::size_t membrane_id = 0; membrane_id != _membrane_map.size(); ++membrane_id) {
    auto [domain_i, domain_o] = _membrane_map[membrane_id];
    auto space_i = compartments_space.subSpace(Assembler::multiIndex(domain_i));
    auto space_o = compartments_space.subSpace(Assembler::multiIndex(domain_o));
    auto space_m = state.space().subSpace(Assembler::multiIndex(Indices::_1, membrane_id));

    auto config_io = _config.sub(fmt::format("{}.boundary.{}.initial", space_i.name(), space_o.name()));
    auto config_oi = _config.sub(fmt::format("{}.boundary.{}.initial", space_o.name(), space_i.name()));
   
    for (std::size_t comp = 0; comp != space_m.degree(); ++comp) {
      auto space_comp = space_m.subSpace(Assembler::multiIndex(comp));
      auto [it, inserted] = functions.try_emplace({space_m.name(), space_comp.name()}, entity_set, make_parser());
      if (not inserted)
        DUNE_THROW(IOError, "Two or more components share the same name");
      auto& function = it->second;
      std::string expr;
      if (config_io.hasKey(space_comp.name()))
        expr = config_io[space_comp.name()];
      else if (config_oi.hasKey(space_comp.name()))
        expr = config_oi[space_comp.name()];
      else
        DUNE_THROW(IOError, fmt::format("Expression for initial value of component '{}' in membrane '{}' does not exist on either boundary", space_comp.name(), space_m.name()));
      function.parser().set_expression(expr);
      for (const auto& [name, tiff] : tiffs)
        function.parser().define_function(name,
          [=](const RF& x, const RF& y) {
            return std::invoke(*tiff, x, y);
          });
      function.setTime(state.time_point());
      function.parser().compile();
    }
  }
#endif

  // std::uniform_real_distribution<double> unif(0.0, 0.5);
  // std::default_random_engine re;
  // re.seed(0);
  // Dune::Assembler::forEachContainerEntry(
  //   state.coefficients(), [&]<class V>(V& value, auto) { value = V{unif(re)}; });

  // state.coefficients() = 1.; // TODO REMOVE!!!!

  // evaluate compiled expressions as initial conditions
  interpolate(state, functions);
}

template<class Traits>
std::unique_ptr<typename ModelMultiDomainDiffusionReaction<Traits>::StepOperator>
ModelMultiDomainDiffusionReaction<Traits>::setup_step_operator(const ConstState& state) const
{
  _logger.debug("Setup local operator"_fmt);

  auto compartment_path = join(compartments_path, Assembler::multiIndex(std::size_t{0}));
  using CompartmentSpace = decltype(state.space().subSpace(compartment_path));

  _logger.trace("Create spatial local operator"_fmt);
  using CompartmentStiffnessLocalOperator = LocalOperatorDiffusionReactionCG<
    CompartmentSpace,
    typename Traits::SpeciesFiniteElementMap::Traits::FiniteElement::Traits::LocalBasisType::Traits>;

  using StiffnessLocalOperator = LocalOperatorMultiDomainDiffusionReaction<Space,CompartmentStiffnessLocalOperator>;
  StiffnessLocalOperator stiff_lop{state.space(), _config};

  _logger.trace("Create temporal local operator"_fmt);
  using CompartmentMassLocalOperator = TemporalLocalOperatorDiffusionReactionCG<
    CompartmentSpace,
    typename Traits::SpeciesFiniteElementMap::Traits::FiniteElement::Traits::LocalBasisType::Traits>;

  using MassLocalOperator = TemporalLocalOperatorMultiDomainDiffusionReaction<Space,CompartmentMassLocalOperator>;
  MassLocalOperator mass_lop{state.space(), _config};

  _logger.debug("Create temporal operator"_fmt);

  using Weight = double;
  using RungeKutta = Assembler::RungeKuttaOperator<Coefficients, Residual, Weight, Time, Time>;
  auto rk_op = std::make_unique<RungeKutta>();

  using RKResidual = typename RungeKutta::InverseOperator::domain_type;
  using RKCoefficients = typename RungeKutta::InverseOperator::range_type;

  bool linear = false;
  // // linear if all jacobians are zero
  // auto compartments_space = state.space().subSpace(compartments_path);
  // for (std::size_t domain_i = 0; domain_i != compartments_space.degree(); ++domain_i) {
  //   auto domain_space_i = compartments_space.subSpace(Assembler::multiIndex(domain_i));
  //   const auto& domain_config_i = _config.sub(domain_space_i.name(), true);
  //   const auto& jac_config = domain_config_i.sub("reaction.jacobian", true);
  //   for (const auto& jac_key : jac_config.getValueKeys()) {
  //     auto jac_expr = jac_config[jac_key];
  //     linear &= (jac_expr == "0") or (jac_expr == "0.0") or (jac_expr == ".0") or (jac_expr == "0.");
  //   }
  //   for (std::size_t domain_o = 0; domain_o != compartments_space.degree(); ++domain_o) {
  //     auto domain_space_o = compartments_space.subSpace(Assembler::multiIndex(domain_o));
  //     const auto& mem_config = domain_config_i.sub(fmt::format("boundary.{}", domain_space_o.name()));
  //     const auto& mem_jac_config = mem_config.sub("reaction.jacobian");
  //     for (const auto& jac_key : mem_jac_config.getValueKeys()) {
  //       auto jac_expr = mem_jac_config[jac_key];
  //       linear &= (jac_expr == "0") or (jac_expr == "0.0") or (jac_expr == ".0") or (jac_expr == "0.");
  //     }
  //     const auto& out_jac_config = mem_config.sub("outflow.jacobian");
  //     for (const auto& jac_key : out_jac_config.getValueKeys()) {
  //       auto jac_expr = out_jac_config[jac_key];
  //       linear &= (jac_expr == "0") or (jac_expr == "0.0") or (jac_expr == ".0") or (jac_expr == "0.");
  //     }
  //   }
  // }

  bool matrix_free = false;
  std::string solver = "BiCGSTAB";

  using Jacobian = decltype([](){
    constexpr std::size_t depth = MultiCompartmentMergingStrategy::Blocked
#if DUNE_COPASI_HAVE_MEMBRANE_SPACE
                                + MultiCompartmentMembraneMergingStrategy::Blocked
#endif
                                + CompartmentMergingStrategy::Blocked;
    if constexpr (depth == 0)
      return DynamicMatrix<BCRSMatrix<double>>{};
    else if constexpr (depth == 1)
      return DynamicMatrix<BCRSMatrix<BCRSMatrix<double>>>{};
    else if constexpr (depth == 2)
      return DynamicMatrix<BCRSMatrix<BCRSMatrix<BCRSMatrix<double>>>>{};
    else if constexpr (depth == 3)
      return DynamicMatrix<BCRSMatrix<BCRSMatrix<BCRSMatrix<BCRSMatrix<double>>>>>{};
    else
      static_assert(AlwaysFalse<Space>{});
  }());

  std::shared_ptr<Assembler::InverseDifferentiableOperator<RKCoefficients,RKResidual>> inverse_op;

  auto linear_solver_apply_op = [solver](Assembler::ForwardOperator<RKCoefficients,RKResidual>& derivative, double reduction, RKResidual& b, RKCoefficients& x) mutable {
    uint64_t solver_timestamp = perfetto::TrackEvent::GetTraceTimeNs();
    TRACE_EVENT("dune", "LinearSolver", solver_timestamp);

    static_assert(std::is_same_v<Coefficients,Residual>);
    auto istl_forward_op = std::make_shared<Dune::Assembler::ISTL::LinearAdapter<RKCoefficients,RKResidual>>([&derivative,&b](const RKCoefficients& coeff, RKResidual& residual){
      derivative.apply(coeff, residual).or_throw();
    });

    auto scalar_product_op = std::make_shared<ScalarProduct<RKCoefficients>>();
    std::shared_ptr<Preconditioner<RKCoefficients,RKResidual>> pre_op;

    pre_op = std::make_shared<Dune::Richardson<RKCoefficients,RKResidual>>(0.1);

    // try to get the jacobian...
    // using JacobianOperator = Assembler::AssembledLinearOperator<Jacobian,RKCoefficients,RKResidual>;
    // auto assembled_derivative = dynamic_cast<JacobianOperator const *>(&derivative);
    // if (assembled_derivative) {
    //   // pre_op = std::make_shared<Dune::SeqSSOR<Jacobian,RKCoefficients,RKResidual,2>>(assembled_derivative->matrix(), 5, 1);
    // }

    int verbosity = 4;
    int max_it = 100;
    
    Dune::InverseOperatorResult result;
    if (solver == "BiCGSTAB") {
      auto solver = Dune::BiCGSTABSolver<RKCoefficients>{istl_forward_op, scalar_product_op, pre_op, reduction, int(max_it), verbosity};
      solver.apply(x, b, result);
    } else if (solver == "RestartedGMRes") {
      int restart = 20;
      auto solver = Dune::RestartedGMResSolver<RKCoefficients>{istl_forward_op, scalar_product_op, pre_op, reduction, restart, int(max_it), verbosity};
      solver.apply(x, b, result);
    } else {
      DUNE_THROW(IOError, "Not known linear solver");
    }

    TRACE_COUNTER("dune", "LinearSolver::Iterations",         solver_timestamp, result.iterations);
    TRACE_COUNTER("dune", "LinearSolver::Reduction",          solver_timestamp, result.reduction);
    TRACE_COUNTER("dune", "LinearSolver::Converged",          solver_timestamp, result.converged);
    TRACE_COUNTER("dune", "LinearSolver::ConvergenceRate",    solver_timestamp, result.conv_rate);
    TRACE_COUNTER("dune", "LinearSolver::ConditionEstimate",  solver_timestamp, result.condition_estimate);

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

      double reduction = 0.1;
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
    auto newton_op = std::make_unique<Assembler::NewtonMethod<RKCoefficients,RKResidual>>();

    newton_op->setLinearSolverApply(linear_solver_apply_op);
    newton_op->setNorm([](const auto& residual){
      TRACE_EVENT("dune", "Norm");
      using Dune::dot;
      using std::sqrt;
      ScalarProduct<RKCoefficients> scalar_product_op;
      return scalar_product_op.norm(residual);
    });

    inverse_op = std::move(newton_op);
  }

  if (matrix_free) { // matrix free
    using MassStiffnesOperator = Assembler::MassStiffnessApplyOperator<RKCoefficients, RKResidual, Space, Space, MassLocalOperator, StiffnessLocalOperator, Weight, Time, Time, Assembler::DurationPosition::StiffnessNumerator>;
    auto diff_op = std::make_shared<MassStiffnesOperator>(state.space(), state.space(), mass_lop, stiff_lop);
    diff_op->setTimePoint(state.time_point());
    inverse_op->setDifferentiableOperator(diff_op);
  } else { // matrix based
    auto pattern_factory = [stiff_lop, mass_lop](const Space& test, const Space& trial, Jacobian& jac){
      auto& entry = jac[0][0];


      auto size_per_row = []<class MultiIndex>(MultiIndex i ,MultiIndex j) -> std::size_t {
        // if constexpr (MultiCompartmentMergingStrategy::Blocked and CompartmentMergingStrategy::Blocked) {
        //   assert(i.size() == j.size());
        //   assert(i.size() < 3);
        //   // TODO: put actual values from spaces
        //   if (i.size() == 0)
        //     return 5; // number of neigboring compartments
        //   if (i.size() == 1)
        //     return 5; // number of edges conecting to a vertex
        //   if (i.size() == 2)
        //     return 10; // number of component per vertex (P1)
        //   return 0;
        // } else if constexpr (not MultiCompartmentMergingStrategy::Blocked and not CompartmentMergingStrategy::Blocked) {
          return 100; /////!!!!!!!!!!!!!!!!!!!!!!! pattern depends on this???? maybe buffer is not working...
        // } else {
        //   static_assert(AlwaysFalse<MultiIndex>{});
        // }
      };

      auto pattern = [&](){
        if constexpr (blockLevel<Jacobian>() == 2) {
          return Assembler::LeafSparsePattern<Space, Space>{test, {}, trial, {}, size_per_row};
        } else if constexpr (blockLevel<Jacobian>() == 3) {
          return  Assembler::BlockedSparsePattern<Assembler::LeafSparsePattern<Space, Space>>{test, {}, trial, {}, size_per_row};
        } else if constexpr (blockLevel<Jacobian>() == 4) {
          return  Assembler::BlockedSparsePattern<Assembler::BlockedSparsePattern<Assembler::LeafSparsePattern<Space, Space>>>{test, {}, trial, {}, size_per_row};
        } else if constexpr (blockLevel<Jacobian>() == 5) {
          return  Assembler::BlockedSparsePattern<Assembler::BlockedSparsePattern<Assembler::BlockedSparsePattern<Assembler::LeafSparsePattern<Space, Space>>>>{test, {}, trial, {}, size_per_row};
        } else {
          static_assert(AlwaysFalse<Space>{});
        }
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

      std::ofstream file{"mat-md.svg"};
      writeSVGMatrix(jac, file);
    };

    using MassStiffnesOperator = Assembler::MassStiffnessOperator<RKCoefficients, RKResidual, Jacobian, Space, Space, MassLocalOperator, StiffnessLocalOperator, Weight, Time, Time, Assembler::DurationPosition::StiffnessNumerator>;
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
ModelMultiDomainDiffusionReaction<Traits>::setup_vtk_writer(State& state)
{
  _logger.debug("Setup VTK writer"_fmt);
  // this map contains the written timesteps for a given path and is shared
  // among different `state`s copied from this one.
  auto mapu_out =
    std::make_shared<std::map<std::string, std::vector<double>>>();

  state.writer = [
    mapu_out
#if DUNE_COPASI_HAVE_MEMBRANE_SPACE
    , membrane_map = this->_membrane_map
#endif
  ](
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

    // let's avoid throubles and write a different file for each sub-domain
    // https://public.kitware.com/pipermail/paraview/2014-July/031732.html

    log_writer.detail("Writing solution for {:.2f}s time stamp"_fmt, state.time_point());
    auto compartments_space = state.space().subSpace(compartments_path);

    for (std::size_t compartment = 0; compartment != compartments_space.degree(); ++compartment) {

      // setup writer again with old timesteps if necessary
      using namespace Dune::Indices;
      auto compartment_space = compartments_space.subSpace(Assembler::multiIndex(compartment));

      using CES = typename decltype(compartment_space)::EntitySet;
      auto writer = std::make_shared<VTKWriter<CES>>(compartment_space.entitySet(), Dune::VTK::conforming);

      std::string name = fmt::format("{}-{}", path.filename().string(), compartment_space.name());
      if (not append) {
        timesteps.clear();
        log_writer.detail("Creating a time sequence file: '{}.pvd'"_fmt, name);
      } else {
        log_writer.trace("Overriding time sequence file: '{}.pvd'"_fmt, name);
      }

      auto sequential_writer = VTKSequenceWriter{ writer, name, path.string(), path.string()};
      sequential_writer.setTimeSteps(timesteps);

      for (std::size_t component = 0; component != compartment_space.degree(); ++component) {
        const auto& component_space = compartment_space.subSpace(Assembler::multiIndex(component));
        auto component_function = Assembler::ScalarDiscreteFunction{component_space, state.coefficients()};
        VTK::FieldInfo info{component_space.name(), VTK::FieldInfo::Type::scalar, 1};
        sequential_writer.vtkWriter()->addVertexData(component_function, info);
      }

      log_writer.trace("Writing vtu file: '{0}/{0}-{1:0>5}.vtu'"_fmt, name, timesteps.size());

      sequential_writer.write(state.time_point(), Dune::VTK::base64);
      sequential_writer.vtkWriter()->clear();
    }

#if DUNE_COPASI_HAVE_MEMBRANE_SPACE
    // we have to write data manually for membranes

    using Geometry = typename Grid::LeafIntersection::Geometry;
    using GlobalCoordinate = Geometry::GlobalCoordinate;
    constexpr auto dim = Geometry::coorddimension;
    static_assert(dim == 2 or dim == 3);

    std::vector< std::vector<GlobalCoordinate> > membrane_points(membrane_map.size());
    std::vector< std::vector<double> > membrane_point_data(membrane_map.size());

    auto multi_mem_space = state.space().subSpace(Assembler::multiIndex(Indices::_1));

    Assembler::LocalContainerBuffer<decltype(multi_mem_space),const Coefficients> lcontainer(multi_mem_space, state.coefficients());

    using Range = typename MembraneSpeciesFiniteElementMap::Traits::FiniteElement::Traits::LocalBasisType::Traits::RangeType;
    std::vector<Range> range;

    Dune::MultipleCodimMultipleGeomTypeMapper mapper{ state.space().entitySet(), Dune::mcmgElementLayout() };

    auto uniqueIndex = [&]<class Entity>(const Entity& entity) {
      return mapper.index(entity);
    };
    auto subDomain = [&](const auto& entity){
      auto domain_set = state.grid().leafGridView().indexSet().subDomains(entity);
      assert(domain_set.size() == 1);
      return *(domain_set.begin());
    };


    auto lspace = multi_mem_space.localView();
    for (const auto& entity : elements(state.grid().leafGridView())) {

      auto id_i = uniqueIndex(entity);
      auto domain_i = subDomain(entity);

      for (const auto& intersection : intersections(state.grid().leafGridView(), entity)) {
        auto membrane_it = end(membrane_map);
        if (intersection.neighbor()) {
          auto id_o = uniqueIndex(intersection.outside());
          if (id_i > id_o) continue;
          auto domain_o = subDomain(intersection.outside());
          auto domains = std::array{std::min(domain_i, domain_o), std::max(domain_i, domain_o)};
          membrane_it = std::find(begin(membrane_map), end(membrane_map), domains);
        } else if (intersection.boundary()) {
          membrane_it = std::find(begin(membrane_map), end(membrane_map), std::array{domain_i, domain_i});
        } else {
          continue;
        }

        assert(membrane_it != end(membrane_map));
        std::size_t membrane_id = std::distance(begin(membrane_map), membrane_it);

        auto geo = intersection.geometry();
        if (not intersection.conforming())
          DUNE_THROW(RangeError, "Membrane intersections shall be conforming");
        if (dim == 2 and not geo.type().isLine())
          DUNE_THROW(RangeError, "Membrane intersections in 2D shall be lines");
        if (dim == 3 and not geo.type().isTriangle())
          DUNE_THROW(RangeError, "Membrane intersections in 3D shall be triangles");

        int sub_index = intersection.indexInInside();
        const auto& entity = intersection.inside();

        lspace.bind(entity);
        lcontainer.load(lspace);

        // evaluate function in geometry corners
        const auto& node = lspace.tree().child(membrane_id);
        auto geo_f = entity.template subEntity<1>(sub_index).geometry();
        for (int c = 0; c != geo_f.corners(); ++c) {
        // store corner coordinates of intersection
          auto corner = geo_f.corner(c);
          membrane_points[membrane_id].push_back(corner);
          auto xlocal = geo_f.local(corner);

          for (std::size_t comp = 0; comp != node.degree(); ++comp){
            const auto& node_comp = node.child(comp).child(sub_index);

            node_comp.finiteElement().localBasis().evaluateFunction(xlocal, range);
            double value = 0;
            for (std::size_t dof = 0; dof < node_comp.finiteElement().size(); ++dof) {
              value += lcontainer(node_comp, dof) * range[dof];
            }
            membrane_point_data[membrane_id].push_back(value);
          }
        }

        lcontainer.clear(lspace);
        lspace.unbind();
      }
    }

    for (std::size_t membrane_id = 0; membrane_id != membrane_map.size(); ++membrane_id) {
      assert(not empty(membrane_points[membrane_id]));
      auto mem_space = multi_mem_space.subSpace(Assembler::multiIndex(membrane_id));
      std::string filename = fmt::format("{0}/{0}-{1}-{2:0>5}.vtu", path.filename().string(), mem_space.name(), timesteps.size());
      log_writer.trace("Writing vtu file: '{}'"_fmt, filename);
      std::ofstream file{ filename };

      std::size_t cells = (dim == 2)
        ? membrane_points[membrane_id].size()/2
        : membrane_points[membrane_id].size()/3;
      file <<
        "<VTKFile type=\"UnstructuredGrid\" version=\"0.1\" byte_order=\"BigEndian\">\n"
        "  <UnstructuredGrid>\n"
        "    <Piece NumberOfPoints=\"" << membrane_points[membrane_id].size() << "\" NumberOfCells=\""<< cells << "\">\n"
        "      <Points>\n"
        "        <DataArray type=\"Float32\" NumberOfComponents=\"3\" Format=\"ascii\">\n";

      for (auto point : membrane_points[membrane_id]) {
        for (std::size_t i = 0; i != dim; ++i)
          file << point[i] << " ";
        if constexpr (dim == 2)
          file << 0;
        file << "\n";
      }

      file <<
        "        </DataArray>\n"
        "      </Points>\n"
        "      <Cells>\n"
        "        <DataArray type=\"Int32\" Name=\"connectivity\" Format=\"ascii\">\n";

      for (std::size_t i = 0; i != membrane_points[membrane_id].size(); ++i)
        file << i << " ";

      file <<
        "\n        </DataArray>\n"
        "        <DataArray type=\"Int32\" Name=\"offsets\" Format=\"ascii\">\n";

      if constexpr (dim == 2) {// line
        for (std::size_t i = 0; i != cells; ++i)
          file << 2+i*2 << " ";
        file <<
          "\n        </DataArray>\n"
          "        <DataArray type=\"Int32\" Name=\"types\" Format=\"ascii\">\n";
        for (std::size_t i = 0; i != cells; ++i)
          file << 3 << " ";
      } else { // triangle
        for (std::size_t i = 0; i != cells; ++i)
          file << 3+i*3 << " ";
        file <<
          "\n        </DataArray>\n"
          "        <DataArray type=\"Int32\" Name=\"types\" Format=\"ascii\">\n";
        for (std::size_t i = 0; i != cells; ++i)
          file << 5 << " ";
      }

      file <<
        "\n        </DataArray>\n"
        "      </Cells>\n"
        "      <PointData Scalars=\"scalars\">\n";

      for (std::size_t comp = 0 ; comp != mem_space.degree(); ++comp) {
        auto comp_space = mem_space.subSpace(Assembler::multiIndex(comp));
        file <<
          "      <DataArray type=\"Float32\" Name=\"" << comp_space.name() << "\" Format=\"ascii\">\n";
        assert(membrane_point_data[membrane_id].size() == membrane_points[membrane_id].size() * mem_space.degree());
        for (std::size_t i = 0; i != membrane_points[membrane_id].size(); ++i)
          file << membrane_point_data[membrane_id][i*mem_space.degree()+comp] << " ";

        file <<
          "      </DataArray>\n";
      }

      file <<
        "      </PointData>\n"
        "    </Piece>\n"
        "  </UnstructuredGrid>\n"
        "</VTKFile>\n";
      file.close();
    }

    // TODO: write pvd file
#endif

    timesteps.push_back(state.time_point());
  };
}

template<class Traits>
void
ModelMultiDomainDiffusionReaction<Traits>::setup(
  BitFlags<ModelSetup::Stages> setup_policy)
{
  _logger.detail("Setting up multi-compartment diffusion-reaction model"_fmt);

  try
  {
    if (setup_policy.test(ModelSetup::Stages::GridFunctionSpace))
      setup_grid_function_space(state());
    if (setup_policy.test(ModelSetup::Stages::CoefficientVector))
      setup_coefficient_vector(state());
    if (setup_policy.test(ModelSetup::Stages::InitialCondition))
      setup_initial_condition(state());
    if (setup_policy.test(ModelSetup::Stages::Writer))
      setup_vtk_writer(state());
  }
  catch(...)
  {
    // detailed report of the configuration paramter tree
    std::stringstream ss;
    _config.report(ss);
    _logger.error(
      2, "---- Multidomain Diffusion-Reaction model parameter tree"_fmt);
    for (std::string line; std::getline(ss, line);)
      _logger.error(2, "{}"_fmt, line);
    _logger.error(2, "----"_fmt);
    throw;
  }
}

// template<class Traits>
// auto
// ModelMultiDomainDiffusionReaction<Traits>::get_data_handler(
//   const ConstState& state) const
// {
//   std::vector<std::shared_ptr<DataHandler>> data(_domains);

//   const auto& compartments = _config.sub("compartments",true).getValueKeys();
//   EntityTransformation et(_grid);
//   for (std::size_t i = 0; i < _domains; ++i) {
//     const std::string compartement = compartments[i];
//     int sub_domain_id =
//       _config.sub("compartments",true).template get<int>(compartement);
//     SubDomainGridView sub_grid_view =
//       _grid->subDomain(sub_domain_id).leafGridView();
//     data[i] = std::make_shared<DataHandler>(
//       *state.grid_function_space, *state.coefficients, sub_grid_view, et);
//   }
//   return data;
// }

// template<class Traits>
// auto
// ModelMultiDomainDiffusionReaction<Traits>::get_grid_function(
//   const ConstState& state,
//   std::size_t domain,
//   std::size_t comp) const -> std::shared_ptr<ComponentGridFunction>
// {
//   auto data = get_data_handler(state);
//   return std::make_shared<ComponentGridFunction>(
//     data[domain]->_lfs.child(domain).child(comp), data[domain]);
// }

// template<class Traits>
// auto
// ModelMultiDomainDiffusionReaction<Traits>::get_grid_function(
//   std::size_t domain,
//   std::size_t comp) const -> std::shared_ptr<ComponentGridFunction>
// {
//   return get_grid_function(const_state(), domain, comp);
// }

// template<class Traits>
// auto
// ModelMultiDomainDiffusionReaction<Traits>::get_grid_functions(
//   const ConstState& state) const
//   -> std::vector<std::vector<std::shared_ptr<ComponentGridFunction>>>
// {
//   const auto& compartments = _config.sub("compartments",true).getValueKeys();
//   std::vector<std::vector<std::shared_ptr<ComponentGridFunction>>>
//     grid_functions(_domains);

//   for (std::size_t domain_i = 0; domain_i < _domains; domain_i++) {
//     const std::string compartement = compartments[domain_i];
//     std::size_t domain =
//       _config.sub("compartments",true).template get<std::size_t>(compartement);
//     const auto& model_config = _config.sub(compartments[domain_i],true);
//     const auto& vars = model_config.sub("diffusion",true);
//     grid_functions[domain_i].resize(vars.getValueKeys().size());
//     for (std::size_t var_i = 0; var_i < vars.getValueKeys().size(); var_i++) {
//       grid_functions[domain_i][var_i] =
//         get_grid_function(state, domain, var_i);
//     }
//   }
//   return grid_functions;
// }

// template<class Traits>
// auto
// ModelMultiDomainDiffusionReaction<Traits>::get_grid_functions() const
//   -> std::vector<std::vector<std::shared_ptr<ComponentGridFunction>>>
// {
//   return get_grid_functions(const_state());
// }

} // namespace Dune::Copasi

#endif // DUNE_COPASI_MODEL_MULTIDOMAIN_DIFFUSION_REACTION_CC
