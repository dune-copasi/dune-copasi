#ifndef DUNE_COPASI_ENUM_HH
#define DUNE_COPASI_ENUM_HH

#include <dune/copasi/common/bit_flags.hh>

namespace Dune::Copasi {

//! Helper objects to define how a model should be setup
namespace ModelSetup {

/**
 * @brief     List of stages in a model setup
 * @details   A setup policy can be made by combining several stages.
 * @see       BitFlags
 */
enum class Stages
{
  None                = 1 << 0,
  GridFunctionSpace   = 1 << 1,
  CoefficientVector   = 1 << 2,
  InitialCondition    = 1 << 3,
  Constraints         = 1 << 4,
  Writer              = 1 << 5
};

//! Policy that setup the grid function spaces in a model
static constexpr Stages setup_grid_function_space = Stages::GridFunctionSpace;

//! Policy that setup the coefficient vectors in a model
static constexpr Stages setup_coefficient_vector =
  setup_grid_function_space | Stages::CoefficientVector;

//! Policy that setup the initial condition in a model
static constexpr Stages setup_initial_condition =
  setup_coefficient_vector | Stages::InitialCondition;

//! Policy that setup the constraints in a model
static constexpr Stages setup_constraints =
  setup_coefficient_vector | Stages::Constraints;

//! Policy that setup the writer in a model
static constexpr Stages setup_writer =
  setup_coefficient_vector | Stages::Writer;

} // namespace ModelSetup

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
