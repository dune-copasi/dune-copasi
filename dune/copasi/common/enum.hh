#ifndef DUNE_COPASI_ENUM_HH
#define DUNE_COPASI_ENUM_HH

namespace Dune::Copasi {

/**
 * @brief      This class describes a model setup policy.
 */
enum class ModelSetupPolicy
{
  None,
  GridFunctionSpace,
  CoefficientVector,
  Constraints,
  LocalOperator,
  GridOperator,
  Solver,
  Writer,
  All
};

/**
 * @brief      This class describes an adaptivity policy.
 */
enum class AdaptivityPolicy
{
  None
};

/**
 * @brief      This class describes a jacobian method.
 */
enum class JacobianMethod
{
  Analytical,
  Numerical
};

} // namespace Dune::Copasi

#endif // DUNE_COPASI_ENUM_HH