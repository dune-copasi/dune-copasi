#ifndef DUNE_COPASI_MODEL_BASE_HH
#define DUNE_COPASI_MODEL_BASE_HH

#include <dune/copasi/enum.hh>

#include <dune/logging/logging.hh>

#include <dune/common/exceptions.hh>
#include <dune/common/float_cmp.hh>

namespace Dune::Copasi {

/**
 * @brief      Class for model base.
 * @todo       Add Dorie to the credits.
 */
class ModelBase {

public:

  /**
   * @brief      Sets the adaptivity policy.
   *
   * @param[in]  adapt_policy  The adaptivity policy.
   */
  void set_policy(AdaptivityPolicy adapt_policy)
  {
    _adapt_policy = adapt_policy;
  }

  /**
   * @brief      Returns the current adaptivity policy.
   *
   * @return     The current adaptivity policy.
   */
  AdaptivityPolicy adaptivity_policy() const {return _adapt_policy;}

  /**
   * @brief      Mark the grid for adaptivity.
   */
  virtual void mark_grid()
  {
    DUNE_THROW(Dune::NotImplemented, "'mark_grid' not implemented!");
  }


  /**
   * @brief      Operations before adaptation of the grid.
   */
  virtual void pre_adapt_grid() {};

  /**
   * @brief      Adapt the grid together it every dependency of the grid (e.g.
   *             solution vector and grid function spaces).
   */
  virtual void adapt_grid()
  {
    if (_adapt_policy != AdaptivityPolicy::None) {
      DUNE_THROW(Dune::NotImplemented, "'adapt_grid' not implemented");
    }
    else {
      DUNE_THROW(Dune::InvalidStateException, "Invalid adaptation policy");
    }
  }

  /**
   * @brief      Operations after adaptation of the grid.
   */
  virtual void post_adapt_grid() {};

  /**
   * @brief      Method that provides the begin time of the model.
   *
   * @return     Begin time of the model.
   */
  virtual double begin_time() const = 0;

  /**
  * @brief      Method that provides the end time of the model.
  *
  * @return     End time of the model.
  */
  virtual double end_time() const = 0;

  /**
   * @brief      Method that provides the current time of the model.
   *
   * @return     Current time of the model.
   */
  virtual double current_time() const = 0;

  /*-----------------------------------------------------------------------*//**
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
  virtual void run()
  {
    auto do_step = [&](const auto& time)
    {
      return Dune::FloatCmp::lt(time, end_time());
    };


    while( do_step(current_time()) )
    {
      step();

      if (adaptivity_policy() != AdaptivityPolicy::None)
        if ( do_step(current_time()) )
        {
          mark_grid();
          pre_adapt_grid();
          adapt_grid();
          post_adapt_grid();
        }
    }
  }

private:
  AdaptivityPolicy _adapt_policy;
};

} // Dune::Copasi namespace

#endif // DUNE_COPASI_MODEL_BASE_HH