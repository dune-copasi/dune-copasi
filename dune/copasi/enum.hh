#ifndef DUNE_COPASI_ENUM_HH
#define DUNE_COPASI_ENUM_HH

namespace Dune::Copasi {

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

enum class AdaptivityPolicy
{
  None
};

enum class JacobianMethod
{
  Analytical,
  Numerical
};

} // namespace Dune::Copasi

#endif // DUNE_COPASI_ENUM_HH