#ifndef DUNE_COPASI_COMMON_STEPPERS_HH
#define DUNE_COPASI_COMMON_STEPPERS_HH

// this is because without it name lookup fails to find PDELab::native(...)
#include <dune/pdelab/backend/istl.hh>

#include <dune/pdelab/instationary/onestep.hh>
#include <dune/pdelab/solver/newton.hh>
#include <dune/pdelab/function/callableadapter.hh>

#include <dune/logging.hh>

#include <dune/common/float_cmp.hh>
#include <dune/common/parametertree.hh>

#include <any>
#include <memory>
#include <tuple>

namespace Dune::Copasi {

/**
 * @brief Base stepper for static polymorphic steppers
 *
 * @details  The setpper implementation must contain the method `do_step` with
 * input and output states. Then, this class will be able to instantiate
 * `do_step` without in/out states and the methods to evolve a system until a
 * final time.
 *
 * @tparam Impl    Stepper implementation
 */
template<class Impl>
class BaseStepper
{
public:
  /**
   * @brief Perform one step in the system
   *
   * @details The input state is advanced a delta time `dt` and placed in `out`
   * state. Signature of `do_step` method that `Impl` must implement
   *
   * @tparam System   System that contain suitable operators to advance in time
   * @tparam State    State to perform time-advance with
   * @tparam Time     Valid time type to operate with
   * @param system    System that contain suitable operators to advance in time
   * @param in        Input state to advance time from
   * @param out       Output state where result will be placed. If step is not
   * possible, out is guaranteed to be false when casted to bool
   * @param dt        Delta time to perform step. If the step is not possible,
   * stepper may modify it to a suitable value
   */
  template<class System, class State, class Time>
  void do_step(System& system, const State& in, State& out, Time& dt) const
  {
    DUNE_THROW(NotImplemented,
               "The derived type '" << className<Impl>()
                                    << "' lacks implementation of do_step");
  }

  /**
   * @brief Perform one step in the system
   *
   * @details The system state is advanced a delta time `dt`
   *
   * @tparam System   System that contain suitable operators to advance in time
   * @tparam Time     Valid time type to operate with
   * @param system    System that contain suitable operators to advance in time.
   * If the step is not possible, the system stays with its initial state
   * @param dt        Delta time to perform step. If the step is not possible,
   * stepper may modify it to a suitable value
   */
  template<class System, class Time>
  void do_step(System& system, Time& dt) const
  {
    auto in = system.state();
    auto out = in;
    asImpl().do_step(system, in, out, dt);
    if (out)
      system.set_state(out);
  }

  /**
   * @brief Evolve the system until `end_time`
   *
   * @details The input state is advanced by time steps `dt` to approach
   * `end_time`. The application of each timestep is performed using
   * the `do_step` method, which may adapt the timestep `dt` adaptively.
   *
   * @tparam System   System that contain suitable operators to advance in time
   * @tparam State    State to perform time-advance with
   * @tparam Time     Valid time type to operate with
   * @param system    System that contain suitable operators to advance in time
   * @param in        Input state to advance time from
   * @param out       Output state where result will be placed. If step is not
   * possible, out is guaranteed to be false when casted to bool
   * @param dt        Delta time to perform first step. If the step is not
   * possible, stepper may modify it to a suitable value
   * @param end_time  Final time that `out` state must reach
   * @param callable  A function called with an state at the end of each
   * succesful step
   * @param snap_to_end_time True if `end_time` must be reached exactly,
   * otherwise it will be reached from below with the provided dt.
   */
  template<class System, class State, class Time, class Callable>
  void evolve(
    System& system,
    const State& in,
    State& out,
    Time& dt,
    const Time& end_time,
    Callable&& callable = [](const auto& state) {},
    bool snap_to_end_time = true) const
  {
    const auto& logger = asImpl().logger();
    logger.notice("Evolving system: {:.2e}s -> {:.2e}s"_fmt, in.time, end_time);
    out = in;
    auto prev_out = in;
    // if snap to end time, stop loop 2 dt before final time, otherwise 1 dt
    double stop_dt = snap_to_end_time ? 2. : 1.;
    while (FloatCmp::le<Time>(out.time+stop_dt*dt, end_time)) {
      std::swap(prev_out, out);
      asImpl().do_step(system, prev_out, out, dt);
      if (not out) {
        logger.warn("Evolving system could not approach final time"_fmt);
        break;
      }
      callable(out);
    }

    if (out and snap_to_end_time)
    {
      std::swap(prev_out, out);
      asImpl().snap_to_time(system,prev_out,out,dt,end_time,callable);
    }
  }

  /**
   * @brief Evolve the system until `end_time`
   *
   * @details The input state is advanced by time steps `dt` to approach
   * `end_time`. The application of each timestep is performed using
   * the `do_step` method, which may adapt the timestep `dt` adaptively.
   *
   * @tparam System   System that contain suitable operators to advance in time
   * @tparam Time     Valid time type to operate with
   * @param system    System that contain suitable operators to advance in time.
   * If the step is not possible, the system stays with its initial state
   * @param dt        Delta time to perform first step. If the step is not
   * possible, stepper may modify it to a suitable value
   * @param end_time  Final time that `out` state must reach
   * @param callable  A function called with an state at the end of each
   * succesful step
   * @param snap_to_end_time True if `end_time` must be reached exactly,
   * otherwise it will be reached from below with the provided dt.
   */
  template<class System, class Time, class Callable>
  void evolve(
    System& system,
    Time& dt,
    const Time& end_time,
    Callable&& callable = [](const auto& state) {},
    bool snap_to_end_time = true) const
  {
    auto in = system.state();
    auto out = in;
    asImpl().evolve(system, in, out, dt, end_time, callable, snap_to_end_time);
    if (out)
      system.set_state(out);
  }

  /**
   * @brief Snap the output to a specific time
   * @details This method is used to advance towards a snap_time
   * avoiding very small timesteps. It work best if the snap time is a couple
   * of timesteps ahead of the suggested dt. It differs from evolve in that the
   * steps might be done with a suitable stepper for wider range of timesteps.
   *
   * @tparam System   System that contain suitable operators to advance in time
   * @tparam Time     Valid time type to operate with
   * @param system    System that contain suitable operators to advance in time.
   * @param dt        Suggested delta time to to reach snap time. It will be
   *                  modified to reach snap time exactly.
   * @param snap_time  Final time that `out` state must reach
   * @param callable  A function called with an state at the end of each
   *                  succesful step
   */
  template<class System, class State, class Time, class Callable>
  void snap_to_time(System& system,
                    const State& in,
                    State& out,
                    Time& dt,
                    const Time& snap_time,
                    Callable&& callable) const
  {
    // use self implementation as default stepper
    snap_to_time(asImpl(), system, in, out, dt, snap_time, callable);
  }

protected:
  /**
   * @copydoc BaseStepper::snap_to_time()
   * @details Default implementation
   *
   * @tparam Stepper Class that fullfils the stepper interface.
   * @param stepper object to perfom `do_step` with.
   */
  template<class Stepper, class System, class State, class Time, class Callable>
  static void snap_to_time(const Stepper& stepper,
                           System& system,
                           const State& in,
                           State& out,
                           Time& dt,
                           const Time& snap_time,
                           Callable&& callable)
  {
    auto& logger = stepper.logger();

    logger.info(2, "Snapping current time to {:.5e}"_fmt, snap_time);

    // let's avoid getting into an infinite loop
    // If max_it is reached something is wrong
    std::size_t max_it = 100, it = 0;

    // reduce last timestep adaptively until time is exactly reached
    out = in;
    auto prev_out = in;
    while (FloatCmp::lt<Time>(out.time, snap_time)) {
      // find number of time steps to fit in the (out.time,snap_time) range
      int timesteps_until_end = std::ceil((snap_time - out.time) / dt);
      if (timesteps_until_end <= 0)
        DUNE_THROW(MathError,
                   "Timestep doesn't make advances towards snap step");
      Time new_dt = (snap_time - out.time) / timesteps_until_end;
      logger.detail(2, "Snapping step size {:.5e}s -> {:.5e}s"_fmt, dt, new_dt);
      dt = new_dt;

      std::swap(prev_out, out);
      stepper.do_step(system, prev_out, out, dt);

      if (out) {
        callable(out);
      } else {
        out = prev_out;
        if (max_it == it++) {
          DUNE_THROW(RangeError,
                    "Snapping time exceeded maximum interation count");
          return;
        }
        logger.detail(
          2, "Reducing step size: {:.5e}s -> {:.5e}s"_fmt, dt, dt * 0.5);
        dt *= 0.5;
      }
    }
  }

private:
  //! Cast to implementation (Barton–Nackman trick)
  const Impl& asImpl() const { return static_cast<const Impl&>(*this); }
};

/**
 * @brief Runge-Kutta stepper for PDEs
 *
 * @details This class is able to advance in time PDE models from PDELab
 * implemented in dune-copasi using [Runge-Kutta
 * methods](https://en.wikipedia.org/wiki/List_of_Runge%E2%80%93Kutta_methods).
 * Fully implicit Runge-Kutta are not available.
 *
 * @warning Due the PDELab objects lifetime, advancing in time may modify system
 * local and grid operators. Thus, when applying steps to different states, be
 * sure that system operators are semantically const correct
 *
 * @tparam T    Time type to step systems with
 */
template<class T = double>
class RKStepper : public BaseStepper<RKStepper<T>>
{
public:
  using Time = T;

  /**
   * @brief Construct a new object
   *
   * @param rk_method   PDELab time stepping Runge-Kutta paramameters
   * @param solver_parameters parameter options for the newton solver
   */
  RKStepper(
    std::unique_ptr<PDELab::TimeSteppingParameterInterface<Time>> rk_method,
    const ParameterTree& solver_parameters)
    : _rk_method{ std::move(rk_method) }
    , _solver_parameters{ solver_parameters }
    , _logger{ Logging::Logging::componentLogger({}, "stepper") }
  {
    _logger.detail("Setting up time stepper"_fmt);
    _logger.trace("Stepper methd: {}"_fmt, _rk_method->name());
  }

  /**
   * @brief Construct a new object
   *
   * @param rk_method   Keyword with the RK method type
   * @param solver_parameters parameter options for the newton solver
   */
  RKStepper(std::string rk_method, const ParameterTree& solver_parameters)
    : RKStepper(make_rk_method(rk_method), solver_parameters)
  {}

  /**
   * @brief Make RK parameters
   *
   * @param type  a Keyword with the RK method type
   * @return auto  a pointer to PDELab time stepping Runge-Kutta paramameters
   */
  static std::unique_ptr<PDELab::TimeSteppingParameterInterface<Time>>
  make_rk_method(const std::string type)
  {
    if (type == "explicit_euler")
      return std::make_unique<PDELab::ExplicitEulerParameter<Time>>();
    if (type == "implicit_euler")
      return std::make_unique<PDELab::ImplicitEulerParameter<Time>>();
    if (type == "heun")
      return std::make_unique<PDELab::HeunParameter<Time>>();
    if (type == "shu_3")
      return std::make_unique<PDELab::Shu3Parameter<Time>>();
    if (type == "runge_kutta_4")
      return std::make_unique<PDELab::RK4Parameter<Time>>();
    if (type == "alexander_2")
      return std::make_unique<PDELab::Alexander2Parameter<Time>>();
    if (type == "fractional_step_theta")
      return std::make_unique<PDELab::FractionalStepParameter<Time>>();
    if (type == "alexander_3")
      return std::make_unique<PDELab::Alexander3Parameter<Time>>();

    DUNE_THROW(NotImplemented, "Not known '" << type << "' Runge Kutta method");
  }

  //! @copydoc BaseStepper::do_step()
  template<class System>
  void do_step(System& system, Time& dt) const
  {
    BaseStepper<RKStepper<T>>::do_step(system, dt);
  }

  //! @copydoc BaseStepper::do_step()
  template<class System, class State>
  void do_step(const System& system,
               const State& in,
               State& out,
               Time& dt) const
  {
    assert(out.grid_function_space and out.grid);

    using Coefficients = typename State::Coefficients;
    const bool cout_redirected = Dune::Logging::Logging::isCoutRedirected();
    if (not cout_redirected)
      Dune::Logging::Logging::redirectCout(_logger.name(),
                                           Dune::Logging::LogLevel::detail);

    _logger.detail("Trying step: {:.2e}s + {:.2e}s -> {:.2e}s"_fmt,
                   in.time,
                   dt,
                   in.time + dt);

  auto& gfs = *in.grid_function_space;
  auto b0lambda = [&](const auto& i, const auto& x) { return false; };
  auto b0 = PDELab::makeBoundaryConditionFromCallable(gfs.gridView(), b0lambda);

  _logger.trace("Assemble constraints"_fmt);
  auto& cc = *system.get_constraints();
  Dune::PDELab::constraints(b0, gfs, cc);

  _logger.detail("Constrained dofs: {} of {}"_fmt, cc.size(), gfs.globalSize());

    auto& solver = get_solver(system);

    const auto& x_in = in.coefficients;
    auto& x_out = out.coefficients;
    if (not x_out or x_out == x_in)
      x_out = std::make_shared<Coefficients>(*x_in);

    // try & catch is the only way we have to check if the solver was successful
    try {
      solver.apply(in.time, dt, *x_in, *x_out);
      solver.result(); // triggers an exception if solver failed

      _logger.notice("Time Step: {:.2e}s + {:.2e}s -> {:.2e}s"_fmt,
                     in.time,
                     dt,
                     in.time + dt);
      out.time = in.time + dt;

    } catch (const Dune::PDELab::NewtonError& e) {
      out.coefficients.reset(); // set output to an invalid state
      _logger.warn(2, "Step failed (NewtonError)"_fmt);
      _logger.trace("{}"_fmt, e.what());
    } catch (const Dune::PDELab::TerminateError& e) {
      out.coefficients.reset(); // set output to an invalid state
      _logger.warn(2, "Step failed (NewtonError::TerminateError)"_fmt);
      _logger.trace("{}"_fmt, e.what());
    } catch (const Dune::PDELab::LineSearchError& e) {
      out.coefficients.reset(); // set output to an invalid state
      _logger.warn(2, "Step failed (NewtonError::LineSearchError)"_fmt);
      _logger.trace("{}"_fmt, e.what());
    } catch (const std::exception& e) {
      // other types are just logged an thrown
      out.coefficients.reset(); // set output to an invalid state
      _logger.warn(2, "Step failed (MathError)"_fmt);
      _logger.trace("{}"_fmt, e.what());
      throw;
    }

    if (not cout_redirected)
      Dune::Logging::Logging::restoreCout();
  }

  //! Return stepper logger
  const Logging::Logger& logger() const { return _logger; }

private:
  //! Setup and cache solver for a given system
  template<class System>
  auto& get_solver(const System& system) const
  {
    using Coefficients = typename System::State::Coefficients;
    using InstationaryGridOperator = typename System::InstationaryGridOperator;
    using LinearSolver = Dune::PDELab::ISTLBackend_NOVLP_BCGS_SSORk<InstationaryGridOperator>;
    using NonLinearOperator =
      PDELab::NewtonMethod<InstationaryGridOperator, LinearSolver>;
    using OneStepOperator = PDELab::
      OneStepMethod<Time, InstationaryGridOperator, NonLinearOperator, Coefficients>;

    using InternalState = std::tuple<System const*,
                                     InstationaryGridOperator const*,
                                     std::shared_ptr<LinearSolver>,
                                     std::shared_ptr<NonLinearOperator>,
                                     std::shared_ptr<OneStepOperator>>;

    auto& grid_operator = *system.get_instationary_grid_operator();

    // If internal data correspond to input, return cached one-step-operator.
    if (_internal_state.has_value() and
        _internal_state.type() == typeid(InternalState)) {
      auto& internal_state = *std::any_cast<InternalState>(&_internal_state);
      if (std::get<0>(internal_state) == &system and
          std::get<1>(internal_state) == &grid_operator)
        return *std::get<4>(internal_state);
    }

    auto linear_solver = std::make_unique<LinearSolver>(grid_operator);

    _logger.trace("Get non-linear operator"_fmt);
    auto non_linear_operator =
      std::make_unique<NonLinearOperator>(grid_operator, *linear_solver);

    // Add settings to the newton solver
    non_linear_operator->setParameters(_solver_parameters);

    _logger.trace("Get one step operator"_fmt);
    auto one_step_operator = std::make_unique<OneStepOperator>(
      *_rk_method, grid_operator, *non_linear_operator);
    one_step_operator->setVerbosityLevel(5);

    _internal_state =
      std::make_any<InternalState>(&system,
                                   &grid_operator,
                                   std::move(linear_solver),
                                   std::move(non_linear_operator),
                                   std::move(one_step_operator));
    return *std::get<4>(*std::any_cast<InternalState>(&_internal_state));
  }

private:
  std::unique_ptr<const PDELab::TimeSteppingParameterInterface<Time>>
    _rk_method;
  const ParameterTree _solver_parameters;
  const Logging::Logger _logger;
  mutable std::any _internal_state;
};

/**
 * @brief Adapt an stepper into a simple adaptive time stepper
 *
 * @details   When an step is not successful, the delta time will be decreased
 * until step is successful or minimum time step is reached. When step is
 * successful, delta time is increased up to a maximum value
 *
 * @tparam SimpleStepper  Stepper that does not adapt its delta time on failure
 */
template<class SimpleStepper>
class SimpleAdaptiveStepper
  : public BaseStepper<SimpleAdaptiveStepper<SimpleStepper>>
{
  using Base = BaseStepper<SimpleAdaptiveStepper<SimpleStepper>>;
public:
  using Time = typename SimpleStepper::Time;

  /**
   * @brief Construct a new Simple Adaptive Stepper object
   *
   * @param stepper Stepper to repeat steps with
   * @param min_step Minimum step size
   * @param max_step Maximum step size
   * @param decrease_factor Factor to decrease time step on failure
   * @param increase_factor Factor to increase time step on success
   */
  SimpleAdaptiveStepper(SimpleStepper&& stepper,
                        Time min_step,
                        Time max_step,
                        double decrease_factor = 0.5,
                        double increase_factor = 1.5)
    : _stepper{ std::move(stepper) }
    , _min_step{ min_step }
    , _max_step{ max_step }
    , _decrease_factor{ decrease_factor }
    , _increase_factor{ increase_factor }
  {
    assert(FloatCmp::lt(decrease_factor, 1.));
    assert(FloatCmp::gt(increase_factor, 1.));
  }

  //! @copydoc BaseStepper::do_step()
  template<class System>
  void do_step(System& system, Time& dt) const
  {
    BaseStepper<SimpleAdaptiveStepper<SimpleStepper>>::do_step(system, dt);
  }

  //! @copydoc BaseStepper::do_step()
  template<class System, class State>
  void do_step(const System& system,
               const State& in,
               State& out,
               Time& dt) const
  {
    if (Dune::FloatCmp::lt<Time>(dt, _min_step))
      DUNE_THROW(InvalidStateException,
                 "time-step '" << dt
                               << "' is less than the minimum allowed step '"
                               << _min_step << "'");

    if (Dune::FloatCmp::gt<Time>(dt, _max_step))
      DUNE_THROW(InvalidStateException,
                 "time-step '" << dt
                               << "' is greater than the maximum allowed step '"
                               << _max_step << "'");

    _stepper.do_step(system, in, out, dt);
    while (not out) {
      logger().detail(2,
        "Reducing step size: {:.2e}s -> {:.2e}s"_fmt, dt, dt * _decrease_factor);
      dt = dt * _decrease_factor;
      if (FloatCmp::lt<Time>(dt, _min_step))
        DUNE_THROW(MathError,
                   "Time step '" << dt
                                 << "' is less than the minimun allowed step '"
                                 << _min_step << "'");
      _stepper.do_step(system, in, out, dt);
    }

    auto new_step = std::min<Time>(_max_step, dt * _increase_factor);
    logger().detail(2,
      "Increasing step size: {:.2e}s -> {:.2e}s"_fmt, dt, new_step);
    dt = new_step;
  }


  /**
   * @copydoc BaseStepper::snap_to_time()
   * @details This class overloads the default implementation in oder to make
   * steps with the simple stepper. That is, without restrictions on maximum and
   * minimum timesteps.
   */
  template<class System, class State, class Time, class Callable>
  void snap_to_time(System& system,
                    const State& in,
                    State& out,
                    Time& dt,
                    const Time& time,
                    Callable&& callable) const
  {
    // use snap_to_time on the wrapped stepper
    _stepper.snap_to_time(system, in, out, dt, time, callable);
  }

  //! Return stepper logger
  const Logging::Logger& logger() const { return _stepper.logger(); }

private:
  const SimpleStepper _stepper;
  const Time _min_step, _max_step;
  const double _decrease_factor, _increase_factor;
};

/**
 * @brief Adapt an stepper into an event time stepper
 *
 * @details   For a given predicate, this stepper finds the moment (or event)
 *            in time when the predicate changes its result. If the predicate
 *            does not change, it simply applyies the do_step of the underlying
 *            stepper.
 *
 * @tparam Stepper  Underlying stepper that will advance the system on time
 * @tparam Predicate Functor that accepts system states and signals event
 * changes
 */
template<class Stepper, class Predicate>
class EventStepper : public BaseStepper<EventStepper<Stepper, Predicate>>
{
  using Base = BaseStepper<EventStepper<Stepper, Predicate>>;

public:
  using Time = typename Stepper::Time;

  EventStepper(Stepper&& stepper,
               Predicate&& predicate,
               double threshold = 1e-4,
               std::size_t max_it = 20)
    : _stepper(std::move(stepper))
    , _predicate(std::move(predicate))
    , _threshold(threshold)
    , _max_it(max_it)
  {}

  /**
   * @copydoc BaseStepper::do_step()
   *
   * @param out       Output state where result will be placed. If step is not
   * possible, out is guaranteed to be false when casted to bool. If the
   * predicate signals an event change, out will hold the change before the
   * event
   */
  template<class System, class State>
  void do_step(const System& system,
               const State& in,
               State& out,
               Time& dt) const
  {
    std::size_t it = _max_it;
    _stepper.do_step(system, in, out, dt);
    if (not out)
      return;
    // let's make an binary search of the approiated time step
    Time start = 0.;
    Time end = dt;
    const bool pre_in = _predicate(in);
    bool pre_out = _predicate(std::as_const(out));
    State test_out = out;
    while ((pre_in != pre_out) and (0 != it--) and
           ((end - start) / dt > _threshold)) {
      Time test = start + (end - start) / 2.;
      _stepper.do_step(system, in, test_out, test);
      bool pre_out = _predicate(std::as_const(test_out));
      if (pre_in == pre_out) {
        out = test_out;
        start = test;
      } else {
        end = test;
      }
    }
  }

  //! Return stepper logger
  const Logging::Logger& logger() const { return _stepper.logger(); }

private:
  const Stepper _stepper;
  Predicate _predicate;
  const double _threshold;
  const std::size_t _max_it;
};

/**
 * @brief Make a default adaptive simple stepper
 *
 * @tparam Time type of time
 * @param config  Configuration file for the time stepping scheme
 * @return auto  An stepper
 */
template<class Time = double>
auto
make_default_stepper(const ParameterTree& config)
{
  auto log = Logging::Logging::componentLogger({}, "stepper");
  auto rk_type = config.get("rk_method","alexander_2");
  auto min_step = config.template get<double>("min_step");
  auto max_step = config.template get<double>("max_step");
  auto decrease_factor = config.get("decrease_factor", 0.9);
  auto increase_factor = config.get("increase_factor", 1.1);

  log.trace("Increase factor: {}"_fmt,increase_factor);
  log.trace("Decrease factor: {}"_fmt,decrease_factor);
  log.trace("Runge-Kutta method: {}"_fmt,rk_type);

  const auto& np = config.sub("newton", true);

  // make parameters compatible with dune-copasi 1.0
  ParameterTree param;
  param["Reduction"] = np["reduction"];
  param["UseMaxNorm"] = np.get("use_max_norm", "false");
  param["MinLinearReduction"] = np["min_linear_reduction"];
  param["FixedLinearReduction"] = np["fixed_linear_reduction"];
  param["AbsoluteLimit"] = np["absolute_limit"];
  param["ReassembleThreshold"] = np["reassemble_threshold"];
  param["KeepMatrix"] = np["keep_matrix"];
  param["Terminate.MaxIterations"] = np["max_iterations"];
  param["Terminate.ForceIteration"] = np["force_iteration"];
  param["LineSearchStrategy"] = np["linear_search.strategy"];
  param["LineSearch.MaxIterations"] = np["linear_search.max_iterations"];
  param["LineSearch.DampingFactor"] = np["linear_search.damping_factor"];
  param["LineSearch.AcceptBest"] = np.get("linear_search.accept_best", "false");

  using RKStepper = Dune::Copasi::RKStepper<double>;
  auto rk_stepper = RKStepper{ rk_type, param };
  using Stepper = Dune::Copasi::SimpleAdaptiveStepper<RKStepper>;
  return Stepper{
    std::move(rk_stepper), min_step, max_step, decrease_factor, increase_factor
  };
}

} // namespace Dune::Copasi

#endif // DUNE_COPASI_COMMON_STEPPERS_HH
