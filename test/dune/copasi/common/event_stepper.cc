#ifdef HAVE_DUNE_COPASI_CONFIG_H
#include <dune/copasi/config.h>
#endif

#include <dune/copasi/common/stepper.hh>

#include <dune/logging.hh>

#include <dune/common/exceptions.hh>

#include <cassert>
#include <iostream>

struct SimpleState
{
  double time;
  double value;
  operator bool() { return true; }
};

struct ExplicitEuler : public Dune::Copasi::BaseStepper<ExplicitEuler>
{
  using Time = double;

  template<class System, class State>
  void do_step(const System& system,
               const State& in,
               State& out,
               Time& dt) const
  {
    out = in;
    out.time += dt;
    out.value = in.value + dt * system(out);
  }

  //! Return stepper logger
  const Dune::Logging::Logger& logger() const { return _logger; }

  Dune::Logging::Logger _logger =
    Dune::Logging::Logging::componentLogger({}, "stepper");
};

int
main(int argc, char** argv)
{
  bool failed = false;

  // initialize mpi
  auto& mpi_helper = Dune::MPIHelper::instance(argc, argv);
  auto comm = mpi_helper.getCollectiveCommunication();
  Dune::Logging::Logging::init(comm, {});
  try {

    // let's solve a very simple explicit euler steps
    // https://en.wikipedia.org/wiki/Euler_method#Using_step_size_equal_to_1_(h_=_1)
    auto f = [](const auto& state) { return state.value; };

    ExplicitEuler explicit_steper{};

    auto find_value = [](const auto& x) { return x.value > 10; };

    auto event_stepper = Dune::Copasi::EventStepper{
      std::move(explicit_steper), std::move(find_value), 1e-7, 200
    };

    SimpleState in{ 0., 1. }, out;
    double dt = 1.;
    event_stepper.do_step(f, in, out, dt);
    if (Dune::FloatCmp::ne(out.value, 2.))
      DUNE_THROW(Dune::MathError, "Not expected value");
    event_stepper.do_step(f, out, in, dt);
    if (Dune::FloatCmp::ne(in.value, 4.))
      DUNE_THROW(Dune::MathError, "Not expected value");
    event_stepper.do_step(f, in, out, dt);
    if (Dune::FloatCmp::ne(out.value, 8.))
      DUNE_THROW(Dune::MathError, "Not expected value");
    event_stepper.do_step(f, out, in, dt);
    if (Dune::FloatCmp::ne(in.value, 10.))
      DUNE_THROW(Dune::MathError, "Not expected value");
    if (Dune::FloatCmp::ne(in.time, 3.25))
      DUNE_THROW(Dune::MathError, "Not expected time");

  } catch (Dune::Exception& e) {
    std::cerr << "Dune reported error: " << e << std::endl;
    failed = true;
  } catch (...) {
    std::cerr << "Unknown exception thrown!" << std::endl;
    failed = true;
  }
  return failed;
}
