#ifndef DUNE_COPASI_COMMON_STEPPERS_HH
#define DUNE_COPASI_COMMON_STEPPERS_HH

#include <dune/copasi/common/exceptions.hh>
#include <dune/copasi/common/fmt_style.hh>

#include <dune/pdelab/common/convergence/reason.hh>
#include <dune/pdelab/common/trace.hh>
#include <dune/pdelab/operator/operator.hh>

#include <dune/common/float_cmp.hh>
#include <dune/common/parametertree.hh>

#include <fmt/core.h>

#include <spdlog/spdlog.h>

#include <any>
#include <concepts>
#include <cstddef>
#include <memory>
#include <tuple>

namespace Dune::Copasi {

/**
 * @brief Simple time stepper of one-step operators
 * @details The polymorphic method `do_step` advances one time step from
 * input to output state with a given duration. Then, this class is able to
 * evolve a system until a final time.
 *
 * @tparam TimeQuantity       Type representing time quantities
 * @tparam DurationQuantity   Type representing duration quantities
 * @tparam State              Type representing states
 */
template<std::copyable State, class TimeQuantity, class DurationQuantity>
class TimeStepper
{
  std::function<TimeQuantity&(State&)> _time_mut;
  std::function<const TimeQuantity&(const State&)> _time_const;

public:
  /**
   * @brief Construct a new time stepper
   *
   * @param time_mut    Functor that obtains a mutable time quantity from a
   * given state
   * @param time_const  Functor that obtains a constant time quantity from a
   * given constant state
   */
  TimeStepper(std::function<TimeQuantity&(State&)> time_mut,
              std::function<const TimeQuantity&(const State&)> time_const)
    : _time_mut(std::move(time_mut))
    , _time_const(std::move(time_const))
  {
  }

  //! Construct a new time stepper
  template<std::same_as<void> = void>
    requires requires(const State& cstate, State& mstate) {
      {
        cstate.time
      } -> std::convertible_to<const TimeQuantity&>;
      {
        mstate.time
      } -> std::convertible_to<TimeQuantity&>;
    }
  TimeStepper()
    : TimeStepper{ &State::time, &State::time }
  {
  }

  TimeStepper(const TimeStepper&) = delete;
  TimeStepper(TimeStepper&&) = delete;

  TimeStepper& operator=(const TimeStepper&) = delete;
  TimeStepper& operator=(TimeStepper&&) = delete;

  //! Time stepper destructor
  virtual ~TimeStepper() = default;

  /**
   * @brief Perform one step of the system
   *
   * @details The input state is advanced a delta time `dt` and placed in `out`
   * state.
   *
   * @param one_step  System that contain suitable operators to advance in time
   * @param in        Input state to advance time from
   * @param out       Output state where result will be placed
   * @param dt        Mutable delta time to perform step
   * @return error condition describing the reason of (not) convergence of the
   * timestep
   */
  virtual PDELab::ErrorCondition do_step(PDELab::OneStep<State>& one_step,
                                         const State& in,
                                         State& out,
                                         DurationQuantity& dt) const
  {
    [[maybe_unused]] uint64_t timestamp = perfetto::TrackEvent::GetTraceTimeNs();
    TRACE_EVENT("dune", "TimeStep", timestamp);
    TRACE_COUNTER("dune", "Stepper::ΔTime", timestamp, dt);

    TimeQuantity time = _time_const(in);
    TimeQuantity time_next = time + dt;
    spdlog::info("Evaluating time step: {:.2e}s + {:.2e}s -> {:.2e}s",
                 DUNE_COPASI_FMT_STYLED_BOLD(time),
                 dt,
                 DUNE_COPASI_FMT_STYLED_BOLD(time_next));

    // set stepper with time-stepping parameters
    one_step["duration"] = dt;
    one_step["time"] = time;

    // advance one step from `in` and store result in `out`
    auto error = one_step.apply(in, out);

    if (not error) {
      _time_mut(out) = time_next; // set new time in resulting state
    } else {
      spdlog::warn("Time step failed: {}", error.message());
    }
    TRACE_COUNTER("dune", "Stepper::Converged", timestamp, bool(not error));
    return error;
  }

  /**
   * @brief Evolve the system until `end_time`
   *
   * @details The input state is advanced by time steps `dt` to approach
   * `end_time`. The application of each timestep is performed using
   * the polymorphic `do_step` method, which may adapt the timestep `dt`
   * dynamically.
   *
   * @param one_step  System that contain suitable operators to advance in time
   * @param in        Input state to advance time from
   * @param out       Output state where result will be placed.
   * @param dt        Mutable delta time to perform step
   * @param end_time  Final time that `out` state must reach
   * @param callable  A function called with an state at the end of each
   * successful step
   * @param snap_to_end_time True if `end_time` must be reached exactly,
   * otherwise it will be reached from below with the provided dt
   */
  PDELab::ErrorCondition evolve(
    PDELab::OneStep<State>& one_step,
    const State& in,
    State& out,
    DurationQuantity& dt,
    TimeQuantity end_time,
    std::function<void(const State&)> callable = [](const State&) {},
    bool snap_to_end_time = true) const
  {
    TimeQuantity time = _time_const(in);
    spdlog::info("Evolving system: {:.2e}s -> {:.2e}s", time, end_time);
    out = in;
    auto prev_out = in;
    // if snap to end time, stop loop 2*dt before final time, otherwise 1*dt
    const double stop_dt = snap_to_end_time ? 2. : 1.;
    PDELab::ErrorCondition error;
    while (FloatCmp::le<TimeQuantity>(_time_const(out) + stop_dt * dt, end_time)) {
      std::swap(prev_out, out);
      error = do_step(one_step, prev_out, out, dt);
      if (error) {
        spdlog::warn("Evolving system could not approach final time");
        return error;
      }
      callable(out);
    }

    if (snap_to_end_time) {
      std::swap(prev_out, out);
      error = snap_to_time(one_step, prev_out, out, dt, end_time, callable);
    }
    return error;
  }

  /**
   * @brief Snap the output to a specific time
   * @details This method is used to advance towards a snap_time
   * avoiding very small timesteps. It work best if the snap time is a couple
   * of timesteps ahead of the suggested dt. It differs from evolve in that the
   * steps might be done with a suitable stepper for wider range of timesteps.
   *
   * @param one_step  System that contain suitable operators to advance in time.
   * @param dt        Suggested delta time to reach snap time. It will be
   *                  modified to reach snap time exactly.
   * @param snap_time Final time that `out` state must reach
   * @param callable  A function called with an state at the end of each
   * successful step
   */
  PDELab::ErrorCondition snap_to_time(
    PDELab::OneStep<State>& one_step,
    const State& in,
    State& out,
    DurationQuantity& dt,
    TimeQuantity snap_time,
    std::function<void(const State&)> callable = [](const State&) {}) const
  {
    spdlog::info("Snapping current time to {:.5e}", snap_time);

    // let's avoid getting into an infinite loop
    // If max_snap_count is reached something is wrong
    const std::size_t max_snap_count = 100;
    std::size_t snap_count = 0;

    // reduce last timestep adaptively until time is exactly reached
    out = in;
    auto prev_out = in;
    TimeQuantity time;
    PDELab::ErrorCondition error;
    while (FloatCmp::lt(time = _time_const(out), snap_time)) {
      // find number of time steps to fit in the (time,snap_time) range
      auto timesteps_ratio = ((snap_time - time) / dt);
      using std::ceil;
      const auto timesteps_until_end = static_cast<int>(ceil(timesteps_ratio));
      if (timesteps_until_end <= 0) {
        throw format_exception(MathError{}, "Timestep doesn't make advances towards snap step");
      }
      DurationQuantity new_dt = (snap_time - time) / timesteps_until_end;
      spdlog::info("Snapping step size {:.5e}s -> {:.5e}s", dt, new_dt);
      dt = new_dt;

      std::swap(prev_out, out);
      error = do_step(one_step, prev_out, out, dt);

      if (not error) {
        callable(out);
      } else {
        out = prev_out;
        if (max_snap_count == snap_count++) {
          throw format_exception(RangeError{}, "Snapping time exceeded maximum iteration count");
        }
        spdlog::info("Reducing step size: {:.5e}s -> {:.5e}s", dt, dt * 0.5);
        dt *= 0.5;
      }
    }
    return error;
  }
};

/**
 * @brief Adapt an stepper into a simple adaptive time stepper
 * @details When a time step is not successful, the delta time will be decreased
 * until step is successful or the minimum time step is reached. When step is
 * successful, delta time is increased up to a maximum value
 *
 * @tparam TimeQuantity       Type representing time quantities
 * @tparam DurationQuantity   Type representing duration quantities
 * @tparam State              Type representing states
 */
template<class State, class TimeQuantity, class DurationQuantity>
class SimpleAdaptiveStepper : public TimeStepper<State, TimeQuantity, DurationQuantity>
{
  using Base = TimeStepper<State, TimeQuantity, DurationQuantity>;

public:
  /**
   * @brief Construct a new simple adaptive stepper
   *
   * @param time_mut        Functor that obtains a mutable time quantity from a
   * given state
   * @param time_const      Functor that obtains a constant time quantity from a
   * given constant state
   * @param decrease_factor Factor to decrease time step on failure
   * @param increase_factor Factor to increase time step on success
   * @param min_step        Minimum step size
   * @param max_step        Maximum step size
   */
  SimpleAdaptiveStepper(std::function<TimeQuantity&(State&)> time_mut,
                        std::function<const TimeQuantity&(const State&)> time_const,
                        double decrease_factor = 0.5,
                        double increase_factor = 1.5,
                        std::optional<DurationQuantity> min_step = {},
                        std::optional<DurationQuantity> max_step = {})
    : Base{ std::move(time_mut), std::move(time_const) }
    , _min_step{ min_step }
    , _max_step{ max_step }
    , _decrease_factor{ decrease_factor }
    , _increase_factor{ increase_factor }
  {
    if (FloatCmp::gt(_decrease_factor, 1.)) {
      throw format_exception(
        IOError{}, "Decrease factor {} must be in the range of (0,1]", _decrease_factor);
    }
    if (FloatCmp::lt(_decrease_factor, 0.)) {
      throw format_exception(
        IOError{}, "Decrease factor {} must be in the range of (0,1]", _decrease_factor);
    }
    if (FloatCmp::lt(_increase_factor, 1.)) {
      throw format_exception(
        IOError{}, "Increase factor {} must be in the range of (1,∞)", _increase_factor);
    }
  }

  /**
   * @brief Construct a new simple adaptive stepper
   *
   * @param decrease_factor Factor to decrease time step on failure
   * @param increase_factor Factor to increase time step on success
   * @param min_step        Minimum step size
   * @param max_step        Maximum step size
   */
  template<std::same_as<void> = void>
    requires requires(const State& cstate, State& mstate) {
      {
        cstate.time
      } -> std::convertible_to<const TimeQuantity&>;
      {
        mstate.time
      } -> std::convertible_to<TimeQuantity&>;
    }
  explicit SimpleAdaptiveStepper(double decrease_factor = 0.5,
                                 double increase_factor = 1.5,
                                 std::optional<DurationQuantity> min_step = {},
                                 std::optional<DurationQuantity> max_step = {})
    : SimpleAdaptiveStepper{
      &State::time, &State::time, decrease_factor, increase_factor, min_step, max_step,
    }
  {
  }

  /**
   * @brief Perform one simple adaptive time step of the system
   *
   * @details When a time step is not successful, the delta time will be
   * decreased until step is successful or the minimum time step is reached.
   * When step is successful, delta time is increased up to a maximum value
   *
   * @param one_step  System that contain suitable operators to advance in time
   * @param in        Input state to advance time from
   * @param out       Output state where result will be placed
   * @param dt        Mutable delta time to perform step
   * @return error condition describing the reason of (not) convergence of the
   * timestep
   */
  PDELab::ErrorCondition do_step(PDELab::OneStep<State>& one_step,
                                 const State& in,
                                 State& out,
                                 DurationQuantity& dt) const override
  {
    PDELab::ErrorCondition error = check_dt(dt);
    if (error) {
      return error;
    }
    // try to do first step
    error = Base::do_step(one_step, in, out, dt);

    // iterate until we succeed or reach a lower time step limit
    while (error) {
      spdlog::warn("Reducing step size: {:.2e}s -> {:.2e}s", dt, dt * _decrease_factor);
      error = check_dt(dt = dt * _decrease_factor);
      if (error) {
        return error;
      }
      error = Base::do_step(one_step, in, out, dt);
    }

    // incerease time step
    auto new_step = dt * _increase_factor;
    if (_max_step) {
      new_step = std::clamp(new_step, -std::abs(_max_step.value()), std::abs(_max_step.value()));
    }
    spdlog::info("Increasing step size: {:.2e}s -> {:.2e}s", dt, new_step);

    dt = new_step;
    return error;
  }

private:
  PDELab::ErrorCondition check_dt(DurationQuantity dt) const
  {
    if (_min_step and Dune::FloatCmp::lt(std::abs(dt), std::abs(*_min_step))) {
      spdlog::warn("Time step '{}' lower limit '{}' was reached ", dt, *_min_step);
      return make_error_condition(Dune::PDELab::Convergence::Reason::DivergedNull);
    }
    if (_max_step and Dune::FloatCmp::gt(std::abs(dt), std::abs(*_max_step))) {
      spdlog::warn("Time step '{}' upper limit '{}' was reached ", dt, *_max_step);
      return make_error_condition(Dune::PDELab::Convergence::Reason::DivergedNull);
    }
    return {};
  }

  std::optional<DurationQuantity> _min_step, _max_step;
  double _decrease_factor{}, _increase_factor{};
};

} // namespace Dune::Copasi

#endif // DUNE_COPASI_COMMON_STEPPERS_HH
