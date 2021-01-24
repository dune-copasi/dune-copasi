#ifndef DUNE_COPASI_MODEL_BASE_HH
#define DUNE_COPASI_MODEL_BASE_HH

#include <dune/copasi/common/enum.hh>

#include <dune/logging/logging.hh>

#include <dune/common/parametertree.hh>

namespace Dune::Copasi {

/**
 * @brief      Class for model base.
 * @todo       Add Dorie to the credits.
 */
class ModelBase
{
public:
  /**
   * @brief      Constructs the model
   *
   * @param[in]  config  The configuration file for the model.
   */
  ModelBase(const Dune::ParameterTree& config);

  /**
   * @brief      Destroys the model
   */
  virtual ~ModelBase();

  /**
   * @brief      Sets the adaptivity policy.
   *
   * @param[in]  adapt_policy  The adaptivity policy.
   */
  void set_policy(AdaptivityPolicy adapt_policy);

  /**
   * @brief      Returns the current adaptivity policy.
   *
   * @return     The current adaptivity policy.
   */
  AdaptivityPolicy adaptivity_policy() const;

  /**
   * @brief      Mark the grid for adaptivity.
   */
  virtual void mark_grid();

  /**
   * @brief      Operations before adaptation of the grid.
   */
  virtual void pre_adapt_grid();

  /**
   * @brief      Adapt the grid together it every dependency of the grid
   * @details   (e.g. solution vector and grid function spaces).
   */
  virtual void adapt_grid();

  /**
   * @brief      Operations after adaptation of the grid.
   */
  virtual void post_adapt_grid();

protected:
  Logging::Logger _logger;

private:
  AdaptivityPolicy _adapt_policy;
};

} // namespace Dune::Copasi

#endif // DUNE_COPASI_MODEL_BASE_HH
