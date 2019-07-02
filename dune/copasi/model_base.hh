#ifndef DUNE_COPASI_MODEL_BASE_HH
#define DUNE_COPASI_MODEL_BASE_HH

#include <dune/copasi/enum.hh>

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
  ~ModelBase();

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

  /**
   * @brief      Method that provides the begin time of the model.
   *
   * @return     Begin time of the model.
   */
  double& begin_time();

  //! @copydoc ModelBase::begin_time()
  const double& begin_time() const;

  /**
   * @brief      Method that provides the end time of the model.
   *
   * @return     End time of the model.
   */
  double& end_time();

  //! @copydoc ModelBase::end_time()
  const double& end_time() const;

  /**
   * @brief      Method that provides the current time of the model.
   *
   * @return     Current time of the model.
   */
  double& current_time();

  //! @copydoc ModelBase::current_time()
  const double& current_time() const;

  /**
   * @brief      Suggest a time step to the model.
   *
   * @param[in]  dt    Suggestion for the internal time step of the model.
   */
  virtual void suggest_timestep(double dt) = 0;

  /**
   * @brief      Performs one steps in direction to end_time().
   * @details    The time-step should never result on a bigger step than the one
   *             suggested in suggest_timestep().
   */
  virtual void step() = 0;

  /**
   * @brief      Runs the model performing steps until current_time() equals
   *             end_time().
   */
  virtual void run();

protected:
  Logging::Logger _logger;

private:
  AdaptivityPolicy _adapt_policy;
  double _begin_time, _end_time, _current_time;
};

} // namespace Dune::Copasi

#endif // DUNE_COPASI_MODEL_BASE_HH