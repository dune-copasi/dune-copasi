#ifndef DUNE_COPASI_MODEL_MAKE_STEP_OPERATOR_HH
#define DUNE_COPASI_MODEL_MAKE_STEP_OPERATOR_HH

#include <dune/pdelab/common/convergence/reason.hh>
#include <dune/pdelab/common/error_condition.hh>
#include <dune/pdelab/common/trace.hh>
#include <dune/pdelab/concepts/basis.hh>
#include <dune/pdelab/operator/adapter.hh>
#include <dune/pdelab/operator/forward/instationary/assembler.hh>
#include <dune/pdelab/operator/forward/instationary/coefficients.hh>
#include <dune/pdelab/operator/forward/runge_kutta.hh>
#include <dune/pdelab/operator/inverse/istl_adapter.hh>
#include <dune/pdelab/operator/inverse/newton.hh>

#include <dune/istl/bvector.hh>

#include <spdlog/spdlog.h>

#include <memory>

namespace Dune::Copasi {

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

  bool const matrix_free = config.get("matrix_free", true);
  std::string const solver = "RestartedGMRes";

  // using Jacobian = decltype([](){
  //   if constexpr (CompartmentMergingStrategy::Blocked)
  //     return DynamicMatrix<BCRSMatrix<BCRSMatrix<double>>>{};
  //   else
  //     return DynamicMatrix<BCRSMatrix<double>>{};
  // }());

  std::shared_ptr<PDELab::Operator<RKResidual, RKCoefficients>> const linear_solver =
    std::make_shared<Dune::PDELab::ISTL::LinearSolver<RKResidual, RKCoefficients>>();

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
          linear_solver->get("convergence_condition.relative_tolerance") = 0.1;
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
    auto newton_op =
      std::make_unique<PDELab::NewtonOperator<RKCoefficients, RKResidual, ResidualQuantity>>();

    newton_op->get("dx_inverse") = linear_solver;
    newton_op->setNorm(
      [](const auto& residual) -> ResidualQuantity { return residual.two_norm2(); });

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
  runge_kutta->get("instationary_coefficients") =
    std::make_shared<PDELab::InstationaryCoefficients>(Dune::PDELab::Alexander2Parameter<double>{});
  return runge_kutta;
}

} // namespace Dune::Copasi

#endif // DUNE_COPASI_MODEL_MAKE_STEP_OPERATOR_HH
