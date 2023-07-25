#ifndef DUNE_COPASI_MODEL_MAKE_STEP_OPERATOR_HH
#define DUNE_COPASI_MODEL_MAKE_STEP_OPERATOR_HH

#include <dune/pdelab/common/algebra.hh>
#include <dune/pdelab/common/container_entry.hh>
#include <dune/pdelab/common/convergence/reason.hh>
#include <dune/pdelab/common/error_condition.hh>
#include <dune/pdelab/common/trace.hh>
#include <dune/pdelab/concepts/basis.hh>
#include <dune/pdelab/operator/adapter.hh>
#include <dune/pdelab/operator/forward/instationary/assembler.hh>
#include <dune/pdelab/operator/forward/instationary/coefficients.hh>
#include <dune/pdelab/operator/forward/runge_kutta.hh>
#include <dune/pdelab/operator/inverse/newton.hh>
#include <dune/pdelab/operator/operator.hh>
#include <dune/pdelab/pattern/basis_to_pattern.hh>
#include <dune/pdelab/pattern/sparsity_pattern.hh>
#include <dune/pdelab/pattern/pattern_to_matrix.hh>

#include <dune/istl/operators.hh>
#include <dune/istl/preconditioners.hh>
#include <dune/istl/scalarproducts.hh>
#include <dune/istl/solvers.hh>
#include <dune/istl/superlu.hh>
#include <dune/istl/umfpack.hh>
#include <dune/istl/io.hh>

#include <dune/common/overloadset.hh>

#include <concepts>
#include <functional>

#include <dune/istl/bvector.hh>

#include <spdlog/spdlog.h>

#include <memory>

namespace Dune::Copasi {

template<class Domain, std::copy_constructible Range>
class LinearSolver : public PDELab::Operator<Range, Domain>
{

  // converts a forward operator into a dune-istl linear operator
  class ISTLLinearAdapter : public Dune::LinearOperator<Domain, Range>
  {
  public:
    ISTLLinearAdapter(PDELab::Operator<Domain, Range>& forward_op)
      : _forward{ forward_op }
    {
    }

    void apply(const Domain& x, Range& y) const override
    {
      PDELab::forEachContainerEntry(
        std::execution::par_unseq, y, []<class T>(T& v) { v = T{ 0 }; });
      _forward.apply(x, y).or_throw();
    }

    using typename Dune::LinearOperator<Domain, Range>::field_type;

    //! apply operator to x, scale and add:  \f$ y = y + \alpha A(x) \f$
    void applyscaleadd(field_type alpha, const Domain& x, Range& y) const override
    {
      tmp = y;
      this->apply(x, tmp);
      Dune::PDELab::axpy(std::execution::par_unseq, y, alpha, tmp);
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
  LinearSolver() {}

  virtual PDELab::ErrorCondition apply(Range& range, Domain& domain) override
  {
    uint64_t solver_timestamp = perfetto::TrackEvent::GetTraceTimeNs();
    TRACE_EVENT("dune", "LinearSolver", solver_timestamp);
    static_assert(std::is_same_v<Range, Domain>);
    auto& forward = this->template get<PDELab::Operator<Domain, Range>>("forward");

    auto istl_forward_op = std::make_shared<ISTLLinearAdapter>(forward);

    auto scalar_product_op = std::make_shared<Dune::ScalarProduct<Range>>();
    std::shared_ptr<Preconditioner<Domain, Range>> pre_op;
    pre_op = std::make_shared<Dune::Richardson<Domain, Range>>(0.1);

    auto solver = this->template get<std::string>("type");
    auto verbosity = this->template get<int>("verbosity");
    auto rel_tol = this->template get<double>("convergence_condition.relative_tolerance");
    const auto& it_range = this->get("convergence_condition.iteration_range").as_vector();
    std::size_t min_iterations = unwrap_property_ref<const std::size_t>(it_range[0]);
    std::size_t max_iterations = unwrap_property_ref<const std::size_t>(it_range[1]);

    Dune::InverseOperatorResult result;
    PDELab::ErrorCondition ec{};
    if (solver == "BiCGSTAB") {
      auto istl_solver = Dune::BiCGSTABSolver<Domain>(
        istl_forward_op, scalar_product_op, pre_op, rel_tol, int(max_iterations), verbosity);
      istl_solver.apply(domain, range, result);
    } else if (solver == "CG") {
      auto istl_solver = CGSolver<Domain>(
        istl_forward_op, scalar_product_op, pre_op, rel_tol, max_iterations, verbosity);
      istl_solver.apply(domain, range, result);
    } else if (solver == "RestartedGMRes") {
      auto restart = this->get("restart", std::size_t{ 40 });
      auto istl_solver = RestartedGMResSolver<Domain, Range>(istl_forward_op,
                                                                          scalar_product_op,
                                                                          pre_op,
                                                                          rel_tol,
                                                                          restart,
                                                                          int(max_iterations),
                                                                          verbosity);
      istl_solver.apply(domain, range, result);
    } else if (solver == "SuperLU") {
#if HAVE_SUPERLU
      if constexpr (std::same_as<Domain,BlockVector<BlockVector<double>>>) {
        if (domain.size() != 1 or range.size() != 1)
          throw format_exception(IOError{}, "SuperLU solver can only be used on semi-implicit or explicit time stepping methods");
        if (not forward.hasKey("container"))
          throw format_exception(IOError{}, "SuperLU solver cannot be a matrix-free operator");
        auto& matrix = forward.template get<DynamicMatrix<BCRSMatrix<double>>>("container");
        auto istl_solver = Dune::SuperLU<BCRSMatrix<double>>(matrix[0][0]);
        istl_solver.apply(domain[0], range[0], result);
      } else {
        throw format_exception(IOError{}, "SuperLU solver can only be used on flat matrices for semi-implicit or explicit time stepping methods");
      }
#else
      throw format_exception(IOError{}, "SuperLU solver is not available");
#endif
    } else if (solver == "UMFPack") {
#if HAVE_SUITESPARSE_UMFPACK
      if constexpr (std::same_as<Domain,BlockVector<BlockVector<double>>>) {
        if (domain.size() != 1 or range.size() != 1)
          throw format_exception(IOError{}, "UMFPack solver can only be used on semi-implicit or explicit time stepping methods");
        if (not forward.hasKey("container"))
          throw format_exception(IOError{}, "UMFPack solver cannot be a matrix-free operator");
        auto& matrix = forward.template get<DynamicMatrix<BCRSMatrix<double>>>("container");
        auto istl_solver = Dune::UMFPack<BCRSMatrix<double>>(matrix[0][0]);
        istl_solver.apply(domain[0], range[0], result);
      } else {
        throw format_exception(IOError{}, "UMFPack solver can only be used on flat matrices for semi-implicit or explicit time stepping methods");
      }
#else
      throw format_exception(IOError{}, "UMFPack solver is not available");
#endif
    } else {
      throw format_exception(IOError{}, "Not known linear solver of type '{}'", solver);
    }

    TRACE_COUNTER("dune", "LinearSolver::Iterations", solver_timestamp, result.iterations);
    TRACE_COUNTER("dune", "LinearSolver::Reduction", solver_timestamp, result.reduction);
    TRACE_COUNTER("dune", "LinearSolver::Converged", solver_timestamp, result.converged);

    if (result.iterations <= static_cast<int>(min_iterations)) {
      spdlog::warn("Minimum of {} iteration requested but {} were performed",
                   min_iterations,
                   result.iterations);
    }

    if (result.converged) {
      return ec;
    } else {
      return make_error_condition(PDELab::Convergence::Reason::DivergedByDivergenceTolarance);
    }
  }

  virtual PDELab::ErrorCondition apply(const Range& range, Domain& domain) override
  {
    Domain tmp_range = range;
    return this->apply(tmp_range, domain);
  }
};

template<class Coefficients, class Residual, class ResidualQuantity, class TimeQuantity, PDELab::Concept::Basis Basis>
[[nodiscard]] inline static std::unique_ptr<PDELab::OneStep<Coefficients>>
make_step_operator(const ParameterTree& config,
                   const Basis& basis,
                   const auto& mass_local_operator,
                   const auto& stiffness_local_operator)
{
  TRACE_EVENT("dune", "StepOperator::SetUp");

  using RKCoefficients = Dune::BlockVector<Coefficients>;
  using RKResidual = Dune::BlockVector<Residual>;

  std::shared_ptr<PDELab::Operator<RKResidual, RKCoefficients>> runge_kutta_inverse;

  bool const is_linear = mass_local_operator.localAssembleIsLinear() and
                         stiffness_local_operator.localAssembleIsLinear();
  if (is_linear) {
    spdlog::info("Local operator is linear");
  } else {
    spdlog::info("Local operator is non-linear");
  }

  auto jacobian_selector = overload(
    [](BlockVector<double>) -> BCRSMatrix<double> {return {};},
    [](BlockVector<BlockVector<double>>) -> BCRSMatrix<BCRSMatrix<double>> { return {}; },
    [](BlockVector<BlockVector<BlockVector<double>>>) -> BCRSMatrix<BCRSMatrix<BCRSMatrix<double>>> { return {}; }
  );

  static_assert(std::same_as<Coefficients, Residual>);

  using Jacobian = decltype(jacobian_selector(Coefficients{}));
  using RKJacobian = DynamicMatrix<Jacobian>;

  std::shared_ptr<PDELab::Operator<RKResidual, RKCoefficients>> const linear_solver =
    std::make_shared<LinearSolver<RKResidual, RKCoefficients>>();

  // configure linear solver
  const auto& lsover_config = config.sub("linear_solver");
  auto svg_path = lsover_config.get("layout.writer.svg.path", "");
  const bool matrix_free = lsover_config.get("matrix_free", false);
  const auto lin_solver = lsover_config.get("type", std::string{ "UMFPack" });
  auto rel_tol = lsover_config.get("convergence_condition.relative_tolerance", 1e-4);
  auto it_range =
    lsover_config.get("convergence_condition.iteration_range", std::array<std::size_t, 2>{ 0, 40 });
  linear_solver->get("type") = lin_solver;
  linear_solver->get("convergence_condition.relative_tolerance") = rel_tol;
  linear_solver->get("convergence_condition.iteration_range") = { it_range[0], it_range[1] };
  linear_solver->get("verbosity") = lsover_config.get("verbosity", int{ 0 });
  linear_solver->get("restart") = lsover_config.get("restart", std::size_t{ 40 });
  spdlog::info("Creating linear solver with '{}' type", lin_solver);

  if (lin_solver == "CG" and not is_linear) {
    spdlog::warn("Using conjugate gradient in a non-linear solver is not recommended");
  }

  if (is_linear) {
    std::optional<RKCoefficients> coeff_zero;
    auto linear_runge_kutta_inverse =
      std::make_unique<PDELab::OperatorAdapter<RKResidual, RKCoefficients>>(
        [coeff_zero, linear_solver](PDELab::Operator<RKResidual, RKCoefficients>& self,
                                    RKResidual& b,
                                    RKCoefficients& x) mutable {
          TRACE_EVENT("dune", "LinearSolver::DefectCorrection");
          static_assert(std::is_same_v<Coefficients, Residual>);
          auto& forward =
            self.template get<PDELab::Operator<RKCoefficients, RKResidual>>("forward");

          forward.apply(x, b).or_throw(); // residual is additive b += F(x)

          if (not coeff_zero) {
            coeff_zero.emplace(x);
            PDELab::forEachContainerEntry(*coeff_zero, []<class T>(T& v, auto) { v = T{ 0 }; });
          }

          RKCoefficients z = *coeff_zero;

          if (not linear_solver->hasKey("forward"))
            linear_solver->get("forward") = forward.derivative(x);
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
                             const auto& op, RKJacobian& jac) {
      // TODO(sospinar): resize outer jacobian wrt instationary coefficients
      jac.resize(1, 1);

      // all jacobian entries of the RK jacobian have the same pattern
      auto& entry = jac[0][0];

      // TODO(sospinar): fix pattern in PDELab!
      auto size_per_row = []<class MultiIndex>(MultiIndex i, MultiIndex j) -> std::size_t {
        //       if constexpr (CompartmentMergingStrategy::Blocked) {
        //         assert(i.size() == j.size());
        //         assert(i.size() < 2);
        //         // TODO: put actual values from spaces
        //         if (i.size() == 0)
        //           return 5; // number of edges connecting to a vertex
        //         if (i.size() == 1)
        //           return 10; // number of component per vertex (P1)
        //         return 0;
        //       } else {
        //          return 5;
        //       }
        return 100;
      };

      auto pattern_selector = overload(
        [&](BCRSMatrix<double>&) -> PDELab::LeafSparsePattern<Basis, Basis> {
          return { basis, {}, basis, {}, size_per_row };
        },
        [&](BCRSMatrix<BCRSMatrix<double>>&)
          -> PDELab::BlockedSparsePattern<PDELab::LeafSparsePattern<Basis, Basis>> {
          return { basis, {}, basis, {}, size_per_row };
        },
        [&](BCRSMatrix<BCRSMatrix<BCRSMatrix<double>>>&)
          -> PDELab::BlockedSparsePattern<
            PDELab::BlockedSparsePattern<PDELab::LeafSparsePattern<Basis, Basis>>> {
          return { basis, {}, basis, {}, size_per_row };
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
          if (i == 0 and j == 0) continue;
          jac[i][j] = entry;
        }
      }

      if (not svg_path.empty()) {
        auto path = fs::path{svg_path}.replace_extension("svg");
        spdlog::info("Writing matrix pattern in svg file: '{}'", path.string());
        std::ofstream file{ path.string() };
        writeSVGMatrix(jac, file);
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
      return std::make_shared<PDELab::InstationaryCoefficients>(pdlab_param);
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
