#ifndef DUNE_COPASI_MODEL_MAKE_STEP_OPERATOR_HH
#define DUNE_COPASI_MODEL_MAKE_STEP_OPERATOR_HH

#include <dune/pdelab/common/algebra.hh>
#include <dune/pdelab/common/container_entry.hh>
#include <dune/pdelab/common/convergence/reason.hh>
#include <dune/pdelab/common/error_condition.hh>
#include <dune/pdelab/common/trace.hh>
#include <dune/pdelab/common/execution.hh>
#include <dune/pdelab/concepts/basis.hh>
#include <dune/pdelab/operator/adapter.hh>
#include <dune/pdelab/operator/forward/instationary/assembler.hh>
#include <dune/pdelab/operator/forward/instationary/coefficients.hh>
#include <dune/pdelab/operator/forward/runge_kutta.hh>
#include <dune/pdelab/operator/inverse/newton.hh>
#include <dune/pdelab/operator/operator.hh>
#include <dune/pdelab/pattern/basis_to_pattern.hh>
#include <dune/pdelab/pattern/pattern_to_matrix.hh>
#include <dune/pdelab/pattern/sparsity_pattern.hh>
#include <dune/copasi/solver/istl/factory/inverse.hh>

#include <dune/istl/io.hh>
#include <dune/istl/operators.hh>
#include <dune/istl/preconditioners.hh>
#include <dune/istl/scalarproducts.hh>
#include <dune/istl/solvers.hh>
#include <dune/istl/superlu.hh>
#include <dune/istl/umfpack.hh>

#include <dune/common/overloadset.hh>

#include <dune/istl/bvector.hh>

#include <spdlog/spdlog.h>

#include <algorithm>
#include <concepts>
#include <functional>
#include <memory>

// memory_resource on MacOS has a minimum runtime target of MacOS >= 14 (see https://developer.apple.com/xcode/cpp/)
#if defined(__cpp_lib_memory_resource) && ((defined(__MAC_OS_X_VERSION_MIN_REQUIRED)  && __MAC_OS_X_VERSION_MIN_REQUIRED  < 140000) || (defined(__IPHONE_OS_VERSION_MIN_REQUIRED) && __IPHONE_OS_VERSION_MIN_REQUIRED < 170000))
#undef __cpp_lib_memory_resource
#endif

#if (__cpp_lib_memory_resource >= 201603L) && (__cpp_lib_polymorphic_allocator >= 201902L)
#include <memory_resource>
#endif

namespace Dune::Copasi {

namespace Impl {
inline auto jacobian_selector =
  overload([](BlockVector<BlockVector<double>>) -> BCRSMatrix<BCRSMatrix<double>> { return {}; },
           [](BlockVector<BlockVector<BlockVector<double>>>)
             -> BCRSMatrix<BCRSMatrix<BCRSMatrix<double>>> { return {}; },
           [](BlockVector<BlockVector<BlockVector<BlockVector<double>>>>)
             -> BCRSMatrix<BCRSMatrix<BCRSMatrix<BCRSMatrix<double>>>> { return {}; });
}

template<class Domain, std::copy_constructible Range>
class LinearSolver : public PDELab::Operator<Range, Domain>
{

  using Jacobian = decltype(Impl::jacobian_selector(Range{}));

  // converts a forward operator into a dune-istl linear operator
  class MatrixFreeAdapter : public Dune::LinearOperator<Domain, Range>
  {
  public:
    MatrixFreeAdapter(PDELab::Operator<Domain, Range>& forward_op)
      : _forward{ forward_op }
    {
    }

    void apply(const Domain& x, Range& y) const override
    {
      PDELab::forEachContainerEntry(
        PDELab::Execution::par_unseq, y, []<class T>(T& v) { v = T{ 0 }; });
      _forward.apply(x, y).or_throw();
    }

    using typename Dune::LinearOperator<Domain, Range>::field_type;

    //! apply operator to x, scale and add:  \f$ y = y + \alpha A(x) \f$
    void applyscaleadd(field_type alpha, const Domain& x, Range& y) const override
    {
      tmp = y;
      this->apply(x, tmp);
      Dune::PDELab::axpy(PDELab::Execution::par_unseq, y, alpha, tmp);
    }

    SolverCategory::Category category() const override
    {
      return SolverCategory::Category::sequential;
    }

  private:
    PDELab::Operator<Domain, Range>& _forward;
    mutable Range tmp;
  };

public:
  LinearSolver(const ParameterTree& config)
    : _config{ config }
  {}

  virtual PDELab::ErrorCondition apply(Range& range, Domain& domain) override
  {
    [[maybe_unused]] uint64_t solver_timestamp = perfetto::TrackEvent::GetTraceTimeNs();
    TRACE_EVENT("dune", "LinearSolver", solver_timestamp);
    static_assert(std::is_same_v<Range, Domain>);
    auto& forward = this->template get<PDELab::Operator<Domain, Range>>("forward");

    Jacobian const * jac = nullptr;
    if (forward.hasKey("container"))
      jac = &forward.template get<const Jacobian>("container");

    // choose a solver from the solver factory based on whether the operator is matrix-free
#if (__cpp_lib_memory_resource >= 201603L) && (__cpp_lib_polymorphic_allocator >= 201902L)
    auto mr = std::pmr::monotonic_buffer_resource{};
    std::pmr::polymorphic_allocator<std::byte> alloc{&mr};
#else
    auto alloc = std::allocator<void>{};
#endif
    std::shared_ptr<InverseOperator<Domain, Range>> solver;
    if (jac) {
      using Op = Dune::MatrixAdapter<Jacobian, Domain, Range>;
      solver = ISTL::makeInverseOperator(std::make_shared<Op>(*jac), _config, alloc);
    } else {
      using Op = MatrixFreeAdapter;
      solver = ISTL::makeInverseOperator(std::make_shared<Op>(forward), _config, alloc);
    }
    if (not solver)
      throw format_exception(InvalidStateException{}, "Solver was configured incorrectly!");

    // relative tolerance may be dynamically changed by the newton solver
    auto rel_tol = this->template get<double>("convergence_condition.relative_tolerance");
    // do solve the linear system
    Dune::InverseOperatorResult result;
    solver->apply(domain, range, rel_tol, result);

    TRACE_COUNTER("dune", "LinearSolver::Iterations", solver_timestamp, result.iterations);
    TRACE_COUNTER("dune", "LinearSolver::Reduction",  solver_timestamp, result.reduction);
    TRACE_COUNTER("dune", "LinearSolver::Converged",  solver_timestamp, result.converged);

    if (result.converged) {
      return PDELab::ErrorCondition{};
    } else {
      return make_error_condition(PDELab::Convergence::Reason::DivergedByDivergenceTolarance);
    }
  }

  virtual PDELab::ErrorCondition apply(const Range& range, Domain& domain) override
  {
    Domain tmp_range = range;
    return this->apply(tmp_range, domain);
  }

private:
  std::shared_ptr<Preconditioner<Domain, Range>> _preconditioner;
  ParameterTree _config;
};

template<class Coefficients,
         class Residual,
         class ResidualQuantity,
         class TimeQuantity,
         PDELab::Concept::Basis Basis>
[[nodiscard]] inline static std::unique_ptr<PDELab::OneStep<Coefficients>>
make_step_operator(const ParameterTree& config,
                   const Basis& basis,
                   const auto& mass_local_operator,
                   const auto& stiffness_local_operator)
{
  TRACE_EVENT("dune", "StepOperator::SetUp");

  using RKCoefficients = Dune::BlockVector<Coefficients>;
  using RKResidual = Dune::BlockVector<Residual>;

  if (basis.dimension() == 0) {
    throw format_exception(InvalidStateException{},
                           "Basis has dimension 0, make sure to have at least one 'scalar_field' "
                           "with a non-empty 'compartment'");
  }

  std::shared_ptr<PDELab::Operator<RKResidual, RKCoefficients>> runge_kutta_inverse;

  bool const is_linear = mass_local_operator.localAssembleIsLinear() and
                         stiffness_local_operator.localAssembleIsLinear();
  if (is_linear) {
    spdlog::info("Local operator is linear");
  } else {
    spdlog::info("Local operator is non-linear");
  }

  static_assert(std::same_as<Coefficients, Residual>);

  using RKJacobian = decltype(Impl::jacobian_selector(RKCoefficients{}));

  // configure linear solver
  const auto& lsover_config = config.sub("linear_solver");
  std::shared_ptr<PDELab::Operator<RKResidual, RKCoefficients>> const linear_solver =
    std::make_shared<LinearSolver<RKResidual, RKCoefficients>>(lsover_config);

  auto svg_path = lsover_config.get("layout.writer.svg.path", "");
  const bool matrix_free = lsover_config.get("matrix_free", false);
  const auto lin_solver = lsover_config.get("type", std::string{ "UMFPack" });
  auto rel_tol = lsover_config.get("convergence_condition.relative_tolerance", 1e-4);
  linear_solver->get("convergence_condition.relative_tolerance") = rel_tol;

  if (lin_solver == "CG" and not is_linear) {
    spdlog::warn("Using conjugate gradient in a non-linear solver is not recommended");
  }

  if (is_linear) {
    std::optional<RKCoefficients> coeff_zero;
    std::shared_ptr<PDELab::Operator<RKCoefficients, RKResidual>> derivative;
    auto linear_runge_kutta_inverse =
      std::make_unique<PDELab::OperatorAdapter<RKResidual, RKCoefficients>>(
        [coeff_zero, linear_solver, derivative](PDELab::Operator<RKResidual, RKCoefficients>& self,
                                    RKResidual& b,
                                    RKCoefficients& x) mutable {
          TRACE_EVENT("dune", "LinearSolver::DefectCorrection");
          static_assert(std::is_same_v<Coefficients, Residual>);
          auto& forward =
            self.template get<PDELab::Operator<RKCoefficients, RKResidual>>("forward");

          forward.apply(x, b).or_throw(); // residual is additive b += F(x)

          if (not coeff_zero) {
            coeff_zero.emplace(x);
            PDELab::forEachContainerEntry(*coeff_zero, []<class T>(T& v) { v = T{ 0 }; });
          }

          RKCoefficients z = *coeff_zero;

          // set jacobian to the linear solver (note that we reuse the same matrix between timesteps)
          linear_solver->get("forward") = derivative = forward.derivative(x, derivative);

          // compute correction
          auto error_condition = linear_solver->apply(b, z);

          if (error_condition) {
            return error_condition;
          }
          PDELab::axpy(x, -1., z);
          return Dune::PDELab::ErrorCondition{};
        });

    runge_kutta_inverse = std::move(linear_runge_kutta_inverse);
  } else {
    spdlog::info("Creating non-linear solver with 'Newton' type");
    auto newton_op =
      std::make_unique<PDELab::NewtonOperator<RKCoefficients, RKResidual, ResidualQuantity>>();

    // configure non-linear solver
    const auto& nlsover_config = config.sub("nonlinear_solver");
    auto nrel_tol = nlsover_config.get("convergence_condition.relative_tolerance", 1e-4);
    auto nit_range = nlsover_config.get("convergence_condition.iteration_range",
                                        std::array<std::size_t, 2>{ 0, 40 });
    auto lin_thld = nlsover_config.get("linearization_threshold", 0.);
    auto dxinv_fixed_tol = nlsover_config.get("dx_inverse_fixed_tolerance", false);
    auto dxinv_min_rel_tol = nlsover_config.get("dx_inverse_min_relative_tolerance", 0.1);
    auto nnorm = nlsover_config.get("norm", "l_2");
    auto verbosity = nlsover_config.get("verbosity", 1);

    newton_op->get("verbosity") = verbosity;
    newton_op->get("convergence_condition.relative_tolerance") = nrel_tol;
    newton_op->get("convergence_condition.iteration_range") = { nit_range[0], nit_range[1] };
    if (nlsover_config.hasKey("convergence_condition.absolute_tolerance")) {
      newton_op->get("convergence_condition.absolute_tolerance") =
        nlsover_config.template get<double>("convergence_condition.absolute_tolerance");
    }
    newton_op->get("linearization_threshold") = lin_thld;
    newton_op->get("dx_inverse_fixed_tolerance") = dxinv_fixed_tol;
    newton_op->get("dx_inverse_min_relative_tolerance") = dxinv_min_rel_tol;

    newton_op->get("dx_inverse") = linear_solver;
    if (nnorm == "l_2") {
      newton_op->setNorm(
        [](const auto& residual) -> ResidualQuantity { return residual.two_norm2(); });
    } else if (nnorm == "l_inf") {
      newton_op->setNorm(
        [](const auto& residual) -> ResidualQuantity { return residual.infinity_norm(); });
    } else if (nnorm == "l_1") {
      newton_op->setNorm(
        [](const auto& residual) -> ResidualQuantity { return residual.one_norm(); });
    } else {
      throw format_exception(IOError{}, "Not known nonlinear norm solver of type '{}'", nnorm);
    }

    runge_kutta_inverse = std::move(newton_op);
  }

  std::shared_ptr<PDELab::Operator<RKCoefficients, RKResidual>> instationary_op;
  if (matrix_free) { // matrix free
    instationary_op = PDELab::makeInstationaryMatrixFreeAssembler<RKCoefficients, RKResidual>(
      basis, basis, mass_local_operator, stiffness_local_operator);
  } else { // matrix based
    auto pattern_factory = [basis, stiffness_local_operator, mass_local_operator, svg_path](
                             const auto& /*op*/, RKJacobian& jac) {
      // TODO(sospinar): resize outer jacobian wrt instationary coefficients
      if (jac.buildStage() == RKJacobian::BuildStage::notAllocated) {
        jac.setBuildMode(RKJacobian::implicit);
        jac.setImplicitBuildModeParameters(1,1.);
        jac.setSize(1,1);
        jac.entry(0,0);
        jac.compress();
        // jac = RKJacobian(1,1,1,RKJacobian::row_wise);
        // for(auto row=jac.createbegin(); row!=jac.createend(); ++row)
        //   for (std::size_t col = 0; col != row.size(); ++col)
        //     row.insert(col);
      }

      // all jacobian entries of the RK jacobian have the same pattern
      auto& entry = jac[0][0];

      // number of components per compartment
      std::vector<std::size_t> comp_size(basis.degree(), 0);
      for (std::size_t i = 0; i != comp_size.size(); ++i)
        comp_size[i] = basis.subSpace(TypeTree::treePath(i )).degree();

      using SizePrefix = typename Basis::SizePrefix;

      // blocked case
      auto comp_max = std::ranges::max(comp_size);
      std::function<std::size_t(SizePrefix, SizePrefix)> block_size;
      if (basis.degree() == basis.size()) { // compartment is blokced
        block_size = [=](SizePrefix, SizePrefix j) -> std::size_t {
          switch (j.size()) {
            case 0:
              return comp_size.size(); // number of compartmetns
            case 1:
              return comp_size[j[0]] * 6; // aprox number of neighboring entities time the number of
                                          // components (assume 2D & P1)
            default:
              std::terminate();
          }
        };
      } else { // component is blokced
        block_size = [=](SizePrefix, SizePrefix j) -> std::size_t {
          switch (j.size()) {
            case 0:
              return comp_max * 6; // max number of components times aprox number of neighboring
                                   // entities (assume 2D)
            case 1:
              return comp_max; // number of components (assume P1)
            default:
              std::terminate();
          }
        };
      }

      // blocked-blocked case
      auto block_block_size = [=]<class MultiIndex>(MultiIndex, MultiIndex j) -> std::size_t {
        switch (j.size()) {
          case 0:
            return comp_size.size(); // number of compartmetns
          case 1:
            return 6; // aprox number of neighboring entities (assume 2D)
          case 2:
            return comp_size[j[1]]; // number of components (assume P1)
          default:
            std::terminate();
        }
        std::terminate();
      };

      auto pattern_selector = overload(
        [&](BCRSMatrix<double>&) -> PDELab::LeafSparsePattern<Basis, Basis> {
          return { basis, {}, basis, {}, comp_max * 6 };
        },
        [&](BCRSMatrix<BCRSMatrix<double>>&)
          -> PDELab::BlockedSparsePattern<PDELab::LeafSparsePattern<Basis, Basis>> {
          return { basis, {}, basis, {}, block_size };
        },
        [&](BCRSMatrix<BCRSMatrix<BCRSMatrix<double>>>&)
          -> PDELab::BlockedSparsePattern<
            PDELab::BlockedSparsePattern<PDELab::LeafSparsePattern<Basis, Basis>>> {
          return { basis, {}, basis, {}, block_block_size };
        });

      auto pattern = pattern_selector(entry);

      PDELab::basisToPattern(stiffness_local_operator, pattern);
      PDELab::basisToPattern(mass_local_operator, pattern);
      pattern.sort();
      PDELab::patternToMatrix(pattern, entry);
      entry = 0.;

      // copy patterns into other runge-kutta jacobians
      for (std::size_t i = 0; i != jac.N(); ++i) {
        for (std::size_t j = 0; j != jac.M(); ++j) {
          if (i == 0 and j == 0)
            continue;
          jac[i][j] = entry;
        }
      }

      if (not svg_path.empty()) {
        auto path = std::filesystem::path{ svg_path }.replace_extension("svg");
        spdlog::info("Writing matrix pattern in svg file: '{}'", path.string());
        std::ofstream file{ path.string() };
        writeSVGMatrix(file, jac);
      }
    };

    instationary_op =
      PDELab::makeInstationaryMatrixBasedAssembler<RKCoefficients, RKResidual, RKJacobian>(
        basis, basis, mass_local_operator, stiffness_local_operator, pattern_factory);
  }

  using RungeKutta = PDELab::RungeKutta<RKCoefficients, RKResidual, TimeQuantity, TimeQuantity>;
  auto runge_kutta = std::make_unique<RungeKutta>();

  runge_kutta->get("inverse") = runge_kutta_inverse;
  runge_kutta->get("inverse.forward") = instationary_op;
  const auto& type = config.get("type", "Alexander2");
  std::shared_ptr<PDELab::InstationaryCoefficients> inst_coeff;

  auto pdelab2coeff = [](auto pdlab_param) {
    return std::make_unique<PDELab::InstationaryCoefficients>(pdlab_param);
  };

  spdlog::info("Creating time-stepping solver with '{}' type", type);

  if (type == "ExplicitEuler") {
    inst_coeff = pdelab2coeff(Dune::PDELab::ExplicitEulerParameter<double>{});
  } else if (type == "ImplicitEuler") {
    inst_coeff = pdelab2coeff(Dune::PDELab::ImplicitEulerParameter<double>{});
  } else if (type == "Heun") {
    inst_coeff = pdelab2coeff(Dune::PDELab::HeunParameter<double>{});
  } else if (type == "Shu3") {
    inst_coeff = pdelab2coeff(Dune::PDELab::Shu3Parameter<double>{});
  } else if (type == "RungeKutta4") {
    inst_coeff = pdelab2coeff(Dune::PDELab::RK4Parameter<double>{});
  } else if (type == "Alexander2") {
    inst_coeff = pdelab2coeff(Dune::PDELab::Alexander2Parameter<double>{});
  } else if (type == "FractionalStepTheta") {
    inst_coeff = pdelab2coeff(Dune::PDELab::FractionalStepParameter<double>{});
  } else if (type == "Alexander3") {
    inst_coeff = pdelab2coeff(Dune::PDELab::Alexander3Parameter<double>{});
  } else {
    throw format_exception(IOError{}, "Not known linear solver of type '{}'", type);
  }

  runge_kutta->get("instationary_coefficients") = inst_coeff;
  return runge_kutta;
}

} // namespace Dune::Copasi

#endif // DUNE_COPASI_MODEL_MAKE_STEP_OPERATOR_HH
