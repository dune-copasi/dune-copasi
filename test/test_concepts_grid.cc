#ifdef HAVE_CONFIG_H
#include "config.h"
#endif

#include <dune/copasi/concepts/grid.hh>

#include <dune/grid/multidomaingrid.hh>
#include <dune/grid/uggrid.hh>
#include <dune/grid/yaspgrid.hh>

#include <dune/common/exceptions.hh>

#include <cassert>
#include <complex>

template<class HG, int i = 2>
bool
test_multidomaingrid()
{
  bool passed = true;

  using namespace Dune::Copasi::Concept;

  // create the multidomain grid from the host grid
  constexpr int dim = HG::dimension;
  using MDGTraits = Dune::mdgrid::DynamicSubDomainCountTraits<dim, 1>;
  using Grid = Dune::mdgrid::MultiDomainGrid<HG, MDGTraits>;

  // check that it is a multidomain grid with its respective subdomain grids
  passed &= isGrid<Grid>();
  if (not passed)
    DUNE_THROW(Dune::Exception, "");
  passed &= isMultiDomainGrid<Grid>();
  if (not passed)
    DUNE_THROW(Dune::Exception, "");
  passed &= isSubDomainGrid<typename Grid::SubDomainGrid>();
  if (not passed)
    DUNE_THROW(Dune::Exception, "");

  if constexpr (i > 0)
    // test for nested multidomains
    return test_multidomaingrid<Grid, i - 1>();
  else
    return passed;
}

int
main(int argc, char** argv)
{
  [[maybe_unused]] auto& mpi_helper = Dune::MPIHelper::instance(argc, argv);

  try {
    bool passed = true;

    using namespace Dune::Copasi::Concept;

    // check grids yaspgrid dim 1
    passed &= isGrid<Dune::YaspGrid<1>>();
    if (not passed)
      DUNE_THROW(Dune::Exception, "");
    passed &= not isMultiDomainGrid<Dune::YaspGrid<1>>();
    if (not passed)
      DUNE_THROW(Dune::Exception, "");
    passed &= not isSubDomainGrid<Dune::YaspGrid<1>>();
    if (not passed)
      DUNE_THROW(Dune::Exception, "");
    passed &= test_multidomaingrid<Dune::YaspGrid<1>>();
    if (not passed)
      DUNE_THROW(Dune::Exception, "");

    // check grids yaspgrid dim 2
    passed &= isGrid<Dune::YaspGrid<2>>();
    if (not passed)
      DUNE_THROW(Dune::Exception, "");
    passed &= not isMultiDomainGrid<Dune::YaspGrid<2>>();
    if (not passed)
      DUNE_THROW(Dune::Exception, "");
    passed &= not isSubDomainGrid<Dune::YaspGrid<2>>();
    if (not passed)
      DUNE_THROW(Dune::Exception, "");
    passed &= test_multidomaingrid<Dune::YaspGrid<2>>();
    if (not passed)
      DUNE_THROW(Dune::Exception, "");

    // check grids yaspgrid dim 5
    passed &= isGrid<Dune::YaspGrid<5>>();
    if (not passed)
      DUNE_THROW(Dune::Exception, "");
    passed &= not isMultiDomainGrid<Dune::YaspGrid<5>>();
    if (not passed)
      DUNE_THROW(Dune::Exception, "");
    passed &= not isSubDomainGrid<Dune::YaspGrid<5>>();
    if (not passed)
      DUNE_THROW(Dune::Exception, "");
    passed &= test_multidomaingrid<Dune::YaspGrid<5>>();
    if (not passed)
      DUNE_THROW(Dune::Exception, "");

    // check grids uggrid dim 2
    passed &= isGrid<Dune::UGGrid<2>>();
    if (not passed)
      DUNE_THROW(Dune::Exception, "");
    passed &= not isMultiDomainGrid<Dune::UGGrid<2>>();
    if (not passed)
      DUNE_THROW(Dune::Exception, "");
    passed &= not isSubDomainGrid<Dune::UGGrid<2>>();
    if (not passed)
      DUNE_THROW(Dune::Exception, "");
    passed &= test_multidomaingrid<Dune::UGGrid<2>>();
    if (not passed)
      DUNE_THROW(Dune::Exception, "");

    // check grids uggrid dim 3
    passed &= isGrid<Dune::UGGrid<3>>();
    if (not passed)
      DUNE_THROW(Dune::Exception, "");
    passed &= not isMultiDomainGrid<Dune::UGGrid<3>>();
    if (not passed)
      DUNE_THROW(Dune::Exception, "");
    passed &= not isSubDomainGrid<Dune::UGGrid<3>>();
    if (not passed)
      DUNE_THROW(Dune::Exception, "");
    passed &= test_multidomaingrid<Dune::UGGrid<3>>();
    if (not passed)
      DUNE_THROW(Dune::Exception, "");

    return not passed;
  } catch (Dune::Exception& e) {
    std::cerr << "Dune reported error: " << e << std::endl;
    return 1;
  } catch (...) {
    std::cerr << "Unknown exception thrown!" << std::endl;
    return 1;
  }
}