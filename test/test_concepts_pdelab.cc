#ifdef HAVE_CONFIG_H
#include "config.h"
#endif

#include <dune/copasi/concepts/pdelab.hh>

#include <dune/pdelab/common/function.hh>
#include <dune/pdelab/function/callableadapter.hh>

#include <dune/grid/yaspgrid.hh>

#include <dune/common/exceptions.hh>
#include <dune/common/math.hh>

#include <cassert>
#include <complex>

// class wuth the grid function interface
template<typename T>
class F
  : public Dune::PDELab::FunctionInterface<
      Dune::PDELab::FunctionTraits<T,
                                   2,
                                   Dune::FieldVector<T, 2>,
                                   T,
                                   1,
                                   Dune::FieldVector<T, 1>>,
      F<T>>
{
public:
  inline void evaluate(const Dune::FieldVector<T, 2>& x,
                       Dune::FieldVector<T, 1>& y) const
  {
    y = sin(3 * Dune::StandardMathematicalConstants<T>::pi() * x[0]) *
        cos(7 * Dune::StandardMathematicalConstants<T>::pi() * x[1]);
  }
};

int
main(int argc, char** argv)
{
  try {
    bool passed = true;

    using namespace Dune::Copasi::Concept;
    [[maybe_unused]] auto& mpi_helper = Dune::MPIHelper::instance(argc, argv);

    // check functions
    passed &= isPDELabFunction<F<double>>();
    if (not passed)
      DUNE_THROW(Dune::Exception, "");
    passed &= isPDELabFunction<F<int>>();
    if (not passed)
      DUNE_THROW(Dune::Exception, "");
    passed &= isPDELabFunction<F<std::complex<double>>>();
    if (not passed)
      DUNE_THROW(Dune::Exception, "");

    passed &= not isPDELabFunction<double>();
    if (not passed)
      DUNE_THROW(Dune::Exception, "");
    passed &= not isPDELabFunction<int>();
    if (not passed)
      DUNE_THROW(Dune::Exception, "");
    passed &= not isPDELabFunction<std::complex<double>>();
    if (not passed)
      DUNE_THROW(Dune::Exception, "");

    // we need a grid for the following checks
    int constexpr dimDomain = 3;
    using Grid = Dune::YaspGrid<dimDomain>;
    using GridView = typename Grid::LeafGridView;
    using DF = typename Grid::ctype;

    Dune::FieldVector<DF, 3> L3(1.0);
    L3[0] = 4;
    Grid grid(L3, { { 4, 4 } });
    GridView grid_view = grid.leafGridView();

    // check local callable
    auto pdelab_local_callable = [](const auto& e, const auto& x) {
      return e.geometry().global(x);
    };

    using PDELabLocalCallable = decltype(pdelab_local_callable);
    passed &= isPDELabLocalCallable<GridView, PDELabLocalCallable>();
    if (not passed)
      DUNE_THROW(Dune::Exception, "");
    passed &= isPDELabCallable<GridView, PDELabLocalCallable>();
    if (not passed)
      DUNE_THROW(Dune::Exception, "");

    passed &= not isPDELabLocalCallable<GridView, F<double>>();
    if (not passed)
      DUNE_THROW(Dune::Exception, "");
    passed &= not isPDELabLocalCallable<GridView, double>();
    if (not passed)
      DUNE_THROW(Dune::Exception, "");

    // check global callable
    auto pdelab_global_callable = [](const auto& x) { return x; };

    using PDELabGlobalCallable = decltype(pdelab_global_callable);
    passed &= isPDELabGlobalCallable<GridView, PDELabGlobalCallable>();
    if (not passed)
      DUNE_THROW(Dune::Exception, "");
    passed &= isPDELabCallable<GridView, PDELabGlobalCallable>();
    if (not passed)
      DUNE_THROW(Dune::Exception, "");

    passed &= not isPDELabGlobalCallable<GridView, F<double>>();
    if (not passed)
      DUNE_THROW(Dune::Exception, "");
    passed &= not isPDELabGlobalCallable<GridView, double>();
    if (not passed)
      DUNE_THROW(Dune::Exception, "");

    passed &= not isPDELabGlobalCallable<GridView, PDELabLocalCallable>();
    if (not passed)
      DUNE_THROW(Dune::Exception, "");
    passed &= not isPDELabLocalCallable<GridView, PDELabGlobalCallable>();
    if (not passed)
      DUNE_THROW(Dune::Exception, "");

    // check grid function from the pdelab creator from callables
    [[maybe_unused]] auto gf_from_local =
      Dune::PDELab::makeGridFunctionFromCallable(grid_view,
                                                 pdelab_local_callable);
    [[maybe_unused]] auto gf_from_global =
      Dune::PDELab::makeGridFunctionFromCallable(grid_view,
                                                 pdelab_global_callable);

    // TODO: For some reason, they are not recoginized as grid functions
    // passed &= isPDELabGridFunction<decltype(gf_from_local)>();
    // if (not passed) DUNE_THROW(Dune::Exception, "");
    // passed &= isPDELabGridFunction<decltype(gf_from_global)>();
    // if (not passed) DUNE_THROW(Dune::Exception, "");

    return not passed;
  } catch (Dune::Exception& e) {
    std::cerr << "Dune reported error: " << e << std::endl;
    return 1;
  } catch (...) {
    std::cerr << "Unknown exception thrown!" << std::endl;
    return 1;
  }
}