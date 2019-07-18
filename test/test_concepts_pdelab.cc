#ifdef HAVE_CONFIG_H
# include "config.h"
#endif

#include <dune/copasi/concepts/pdelab.hh>

#include <dune/pdelab/common/function.hh>
#include <dune/pdelab/function/callableadapter.hh>

#include <dune/grid/yaspgrid.hh>

#include <dune/common/math.hh>

#include <cassert>
#include <complex>

// class wuth the grid function interface
template<typename T>
class F : public Dune::PDELab::FunctionInterface<
  Dune::PDELab::FunctionTraits<T,2,Dune::FieldVector<T,2>,
                               T,1,Dune::FieldVector<T,1> >,
  F<T> >
{
public:
  inline void evaluate (const Dune::FieldVector<T,2>& x,
                        Dune::FieldVector<T,1>& y) const
  {
    y = sin(3*Dune::StandardMathematicalConstants<T>::pi()*x[0])
      * cos(7*Dune::StandardMathematicalConstants<T>::pi()*x[1]);
  }
};


int main(int argc, char** argv)
{
  try{
    bool passed = true;

    using namespace Dune::Copasi::Concept;

    // check functions
    passed &= isPDELabFunction<F<double>>();
    passed &= isPDELabFunction<F<int>>();
    passed &= isPDELabFunction<F<std::complex<double>>>();

    passed &= not isPDELabFunction<double>();
    passed &= not isPDELabFunction<int>();
    passed &= not isPDELabFunction<std::complex<double>>();

    // we need a grid for the following checks
    int constexpr dimDomain = 3;
    using Grid = Dune::YaspGrid<dimDomain>;
    using GridView = typename Grid::LeafGridView;
    using DF = typename Grid::ctype;

    Dune::FieldVector<DF,3> L3(1.0); L3[0] = 4;
    Grid grid(L3,{{4,4}});
    GridView grid_view = grid.leafGridView();

    // check local callable
    auto pdelab_local_callable = [](const auto& e, const auto& x)
    {
      return e.geometry().global(x);
    };

    using PDELabLocalCallable = decltype(pdelab_local_callable);
    passed &= isPDELabLocalCallable<GridView,PDELabLocalCallable>();
    passed &= isPDELabCallable<GridView,PDELabLocalCallable>();

    passed &= not isPDELabLocalCallable<GridView,F<double>>();
    passed &= not isPDELabLocalCallable<GridView,double>();

    // check global callable
    auto pdelab_global_callable = [](const auto& x)
    {
      return x;
    };

    using PDELabGlobalCallable = decltype(pdelab_global_callable);
    passed &= isPDELabGlobalCallable<GridView,PDELabGlobalCallable>();
    passed &= isPDELabCallable<GridView,PDELabGlobalCallable>();

    passed &= not isPDELabGlobalCallable<GridView,F<double>>();
    passed &= not isPDELabGlobalCallable<GridView,double>();

    passed &= not isPDELabGlobalCallable<GridView,PDELabLocalCallable>();
    passed &= not isPDELabLocalCallable<GridView,PDELabGlobalCallable>();

    // check grid function from the pdelab creator from callables
    auto gf_from_local = Dune::PDELab::makeGridFunctionFromCallable(grid_view,pdelab_local_callable);
    auto gf_from_global = Dune::PDELab::makeGridFunctionFromCallable(grid_view,pdelab_global_callable);

    passed &= isPDELabGridFunction<decltype(gf_from_local)>();
    passed &= isPDELabGridFunction<decltype(gf_from_global)>();

    return not passed;
  }
  catch (Dune::Exception &e){
    std::cerr << "Dune reported error: " << e << std::endl;
  }
  catch (...){
    std::cerr << "Unknown exception thrown!" << std::endl;
  }
}