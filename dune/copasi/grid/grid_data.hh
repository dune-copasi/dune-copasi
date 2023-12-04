#ifndef DUNE_COPASI_GRID_GRID_DATA_HH
#define DUNE_COPASI_GRID_GRID_DATA_HH

#include <dune/copasi/common/exceptions.hh>
#include <dune/copasi/common/filesystem.hh>

#include <dune/pdelab/operator/operator.hh>

#include <dune/functions/gridfunctions/gridviewfunction.hh>

#include <dune/grid/concepts/grid.hh>
#include <dune/grid/concepts/gridview.hh>

#include <dune/common/parametertree.hh>

#include <any>
#include <map>
#include <string>
#include <memory>
#include <unordered_map>

namespace Dune::Copasi {

template<Dune::Concept::Grid Grid_,
         Dune::Concept::GridView GridView_ = typename Grid_::LeafGridView>
struct GridData
{

  GridData(){};

  // Make the gmsh id available --> testing grid data (TO DO: @Dylan)
  //std::function<double(HostEntity)> _gmsh_id_map;

};

};

#endif // DUNE_COPASI_GRID_GRID_DATA_HH
