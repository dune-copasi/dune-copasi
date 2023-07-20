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

#include <dune/istl/operators.hh>
#include <dune/istl/preconditioners.hh>
#include <dune/istl/scalarproducts.hh>
#include <dune/istl/solvers.hh>
#include <dune/istl/superlu.hh>
#include <dune/istl/umfpack.hh>

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

    // // try to get the jacobian...
    // using JacobianOperator = PDELab::AssembledLinearOperator<Jacobian,Domain,Range>;
    // auto assembled_derivative = dynamic_cast<JacobianOperator const *>(&derivative);
    // if (assembled_derivative) {
    //   // pre_op =
    //   std::make_shared<Dune::SeqSSOR<Jacobian,Domain,Range,2>>(assembled_derivative->matrix(), 5,
    //   1);
    // }

    auto solver = this->template get<std::string>("type");
    auto verbosity = this->template get<u_char>("verbosity");
    auto rel_tol = this->template get<double>("convergence_condition.relative_tolerance");
    const auto& it_range = this->get("convergence_condition.iteration_range").as_vector();
    std::size_t min_iterations = unwrap_property_ref<const std::size_t>(it_range[0]);
    std::size_t max_iterations = unwrap_property_ref<const std::size_t>(it_range[1]);

    Dune::InverseOperatorResult result;
    PDELab::ErrorCondition ec{};
    std::unique_ptr<IterativeSolver<Domain, Range>> istl_solver;
    if (solver == "BiCGSTAB") {
      istl_solver = std::make_unique<Dune::BiCGSTABSolver<Domain>>(
        istl_forward_op, scalar_product_op, pre_op, rel_tol, int(max_iterations), verbosity);
    } else if (solver == "CG") {
      istl_solver = std::make_unique<CGSolver<Domain>>(
        istl_forward_op, scalar_product_op, pre_op, rel_tol, max_iterations, verbosity);
    } else if (solver == "RestartedGMRes") {
      auto restart = this->get("restart", std::size_t{ 40 });
      istl_solver = std::make_unique<RestartedGMResSolver<Domain, Range>>(istl_forward_op,
                                                                          scalar_product_op,
                                                                          pre_op,
                                                                          rel_tol,
                                                                          restart,
                                                                          int(max_iterations),
                                                                          verbosity);

    } else {
      throw format_exception(IOError{}, "Not known linear solver of type '{}'", solver);
    }
    istl_solver->apply(domain, range, result);

    TRACE_COUNTER("dune", "LinearSolver::Iterations", solver_timestamp, result.iterations);
    TRACE_COUNTER("dune", "LinearSolver::Reduction", solver_timestamp, result.reduction);
    TRACE_COUNTER("dune", "LinearSolver::Converged", solver_timestamp, result.converged);

    if (result.iterations <= min_iterations) {
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

template<class Coefficients, class Residual, class ResidualQuantity, class TimeQuantity>
[[nodiscard]] inline static std::unique_ptr<PDELab::OneStep<Coefficients>>
make_step_operator(const ParameterTree& config,
                   const PDELab::Concept::Basis auto& basis,
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

  // using Jacobian = decltype([](){
  //   if constexpr (CompartmentMergingStrategy::Blocked)
  //     return DynamicMatrix<BCRSMatrix<BCRSMatrix<double>>>{};
  //   else
  //     return DynamicMatrix<BCRSMatrix<double>>{};
  // }());

  std::shared_ptr<PDELab::Operator<RKResidual, RKCoefficients>> const linear_solver =
    std::make_shared<LinearSolver<RKResidual, RKCoefficients>>();

  // configure linear solver
  const auto& lsover_config = config.sub("linear_solver");
  const bool matrix_free = lsover_config.get("matrix_free", true);
  const auto lin_solver = lsover_config.get("type", std::string{ "RestartedGMRes" });
  auto rel_tol = lsover_config.get("convergence_condition.relative_tolerance", 1e-4);
  auto it_range =
    lsover_config.get("convergence_condition.iteration_range", std::array<std::size_t, 2>{ 0, 40 });
  linear_solver->get("type") = lin_solver;
  linear_solver->get("convergence_condition.relative_tolerance") = rel_tol;
  linear_solver->get("convergence_condition.iteration_range") = { it_range[0], it_range[1] };
  linear_solver->get("verbosity") = lsover_config.get("verbosity", u_char{ 1 });
  linear_solver->get("restart") = lsover_config.get("restart", std::size_t{ 40 });
  spdlog::info("Creating linear solver with '{}' method", lin_solver);

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

          linear_solver->get("forward") = forward.derivative(x);
          // compute correction
          auto error_condition = linear_solver->apply(b, z);

          if (error_condition) {
            return error_condition;
          }
          PDELab::axpy(x, -1., z);
          linear_solver->get("forward") = nullptr;
          return Dune::PDELab::ErrorCondition{};
        });

    runge_kutta_inverse = std::move(linear_runge_kutta_inverse);
  } else {
    spdlog::info("Creating non-linear solver with 'Newton' method");
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
    DUNE_THROW(NotImplemented, "");
    //   //   auto pattern_factory = [stiffness_local_operator, mass_local_operator](const Basis&
    //   test, const Basis& trial, Jacobian& jac){
    //   //     auto& entry = jac[0][0];

    //   //     auto size_per_row = []<class MultiIndex>(MultiIndex i ,MultiIndex j) -> std::size_t
    //   {
    //   //       if constexpr (CompartmentMergingStrategy::Blocked) {
    //   //         assert(i.size() == j.size());
    //   //         assert(i.size() < 2);
    //   //         // TODO: put actual values from spaces
    //   //         if (i.size() == 0)
    //   //           return 5; // number of edges connecting to a vertex
    //   //         if (i.size() == 1)
    //   //           return 10; // number of component per vertex (P1)
    //   //         return 0;
    //   //       } else {
    //   //         return 5;
    //   //       }
    //   //     };

    //   //     auto pattern = [&](){
    //   //       if constexpr (CompartmentMergingStrategy::Blocked)
    //   //         return PDELab::BlockedSparsePattern<PDELab::LeafSparsePattern<Basis,
    //   Basis>>{test, {}, trial, {}, size_per_row};
    //   //       else
    //   //         return PDELab::LeafSparsePattern<Basis, Basis>{test, {}, trial, {},
    //   size_per_row};
    //   //     }();

    //   //     spaceToPattern(stiffness_local_operator, pattern);
    //   //     spaceToPattern(mass_local_operator, pattern);
    //   //     pattern.sort();
    //   //     patternToMatrix(pattern, entry);
    //   //     entry = 0.;

    //   //     // copy patterns into other runge-kutta jacobians
    //   //     for (std::size_t i = 0; i != jac.N(); ++i) {
    //   //       for (std::size_t j = 0; j != jac.M(); ++j) {
    //   //         if (i == 0 and j == 0) continue;
    //   //         jac[i][j] = entry;
    //   //       }
    //   //     }

    //   //     std::ofstream file{"mat-sd.svg"};
    //   //     writeSVGMatrix(jac, file);
    //   //   };

    //   //   using MassStiffnesOperator = PDELab::MassStiffnessOperator<RKCoefficients, RKResidual,
    //   Jacobian, Basis, Basis, MassLocalOperator, StiffnessLocalOperator, Weight, TimeQuantity,
    //   TimeQuantity, PDELab::DurationPosition::StiffnessNumerator>;
    //   //   auto diff_op = std::make_shared<MassStiffnesOperator>(state.basis(), state.basis(),
    //   mass_local_operator, stiffness_local_operator, pattern_factory);
    //   //   diff_op->setTimePoint(state.time());
    //   //   rk_invesrse->setDifferentiableOperator(diff_op);
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

  spdlog::info("Creating time-stepping solver with '{}' method", type);

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
