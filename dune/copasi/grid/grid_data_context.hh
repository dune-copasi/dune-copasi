#ifndef DUNE_COPASI_GRID_DATA_CONTEXT_HH
#define DUNE_COPASI_GRID_DATA_CONTEXT_HH

#include <dune/copasi/common/exceptions.hh>
#include <dune/copasi/parser/factory.hh>
#include <dune/copasi/grid/cell_data.hh>
#include <dune/copasi/concepts/grid.hh>

#include <dune/grid/common/exceptions.hh>
#include <dune/grid/common/gridinfo.hh>
#include <dune/grid/common/mcmgmapper.hh>
#include <dune/grid/common/rangegenerators.hh>
#include <dune/grid/io/file/gmshreader.hh>
#include <dune/grid/uggrid/uggridfactory.hh>
#include <dune/grid/utility/structuredgridfactory.hh>

#include <spdlog/spdlog.h>

#include <fmt/core.h>

#include <dune/common/fvector.hh>
#include <dune/common/exceptions.hh>
#include <dune/common/parametertree.hh>
#include <dune/common/overloadset.hh>

#include <unordered_map>
#include <vector>
#include <any>

namespace Dune::Copasi {

// Class that manages the grid data
template<Dune::Concept::GridView GV>
class GridDataContext
{

  constexpr static auto dim = GV::dimension;

  // The necessary types
  using IndexSet = typename GV::IndexSet;
  using DataType = double;

public:

  // This is before the grid is constructed ---> only generic parameters
  explicit GridDataContext(const ParameterTree& config, GV grid_view)
    : _cell_data( grid_view.grid().levelGridView(0) )
  {
    auto grid_config = config.sub("grid");

    if (grid_config.hasKey("path")) {
        has_griddata = true;
        auto grid_path = grid_config.template get<std::string>("path");
        // load grid data
        load_grid_data(grid_config);
    }
  };

  void update_cell_data(const Concept::IndexableEntity<IndexSet> auto& entity,
                        std::vector<DataType>& cell_data,
                        std::vector<bool>& cell_mask) const {
    // Check whether grid is created from Gmsh path -> if not return
    // (cannot assign data to a rule-based grid created by refine!)
    if( not has_griddata )
      return;

    _cell_data.getData(entity, cell_data, cell_mask);
  }

  //! Number of elements per entity in the container
  std::size_t cell_data_size() const{
    return _cell_data.size();
  }

  //! Return the keys of cell data
  std::vector<std::string> keys() const{
    return _cell_data.keys();
  }

private:

  /**
  * @brief Load the griddata
  *
  * The parameter tree is unwinded and used to load and create the needed variables
  *
  * @param[in]  ParameterTree     The parameter tree containing the grid configuration properties
  */
  void load_grid_data(Dune::ParameterTree& grid_config){

    auto& griddata_config = grid_config.sub("griddata");

    // reserve the data needed to store the griddata
    // --> Every element can contain griddata for every key
    // --> Fixed look up table adresses providing (potentially) highest performance
    // --> First spot is always reserved for the gmsh_id (notice +1)
    _cell_data.reserve(griddata_config.getSubKeys().size() + 1 );

    add_gmsh_id(grid_config);

    for (const auto& sub : griddata_config.getSubKeys()) {
      const auto& sub_config = griddata_config.sub(sub, true);
      const std::string datafile_path = sub_config["path"];
      const std::string type = sub_config["type"];
      if (type == "scalar") {
        add_scalar_griddata(datafile_path, sub);
      } else if (type == "vector") {
        throw format_exception(NotImplemented{}, "vector griddata is not implemented");
      } else if (type == "tensor") {
        throw format_exception(NotImplemented{}, "tensor griddata is not implemented");
      } else {
        throw format_exception(IOError{}, "Unknown type {}", type);
      }
    }
  };


private:

  // Read the gmsh file & add 'gmsh_id' to cell_data
  void add_gmsh_id(Dune::ParameterTree& grid_config){

    auto grid_parser = [](){
      if constexpr (Concept::MultiDomainGrid<typename GV::Grid>) {
        using Grid = typename GV::Grid::HostGrid;
        auto grid_factory = Dune::GridFactory<Grid>{};
        return Dune::GmshReaderParser<Grid>{ grid_factory, false, false };
      }else{
        using Grid = typename GV::Grid;
        auto grid_factory = Dune::GridFactory<Grid>{};
        return Dune::GmshReaderParser<Grid>{ grid_factory, false, false };
      }
    }();

    auto grid_path = grid_config.template get<std::string>("path");

    //auto grid_factory = Dune::GridFactory<Grid>{};
    //auto grid_parser = Dune::GmshReaderParser<Grid>{ grid_factory, false, false };

    grid_parser.read(grid_path);
    auto& gmsh_physical_entity = grid_parser.elementIndexMap();

    std::map<std::size_t, DataType> data;
    for(std::size_t i = 0; i < gmsh_physical_entity.size(); i++)
      data[i] = gmsh_physical_entity[i];

    _cell_data.addData("gmsh_id", data);

  };

  // Read a datafile for the scalar griddata type & add it to cell_data
  void add_scalar_griddata(std::string datafile_path, std::string key){
    std::cout << "Adding datafile: " << datafile_path << std::endl;

    // read the data file:
    std::ifstream gridDataFile;
    gridDataFile.open(datafile_path);

    if(!gridDataFile.is_open()){
      throw format_exception(IOError{}, "Could not open datafile");
    }

    std::string line;
    int numDataLines = 0;
    bool foundNumDataLines = false;

    // Find first non comment line and format it to an integer which should
    // provide the number of data lines in the file.
    while (std::getline(gridDataFile, line)) {
        if (line.empty() || line[0] == '#') {
            continue; // Ignore comment lines
        }
        // Check if this line contains the number of data lines
        std::istringstream iss(line);
        if (!(iss >> numDataLines)) {
            throw format_exception(IOError{}, "IO format error -> format of number of data lines not correct!");
        }
        foundNumDataLines = true;
        break; // Found the number of data lines, exit the loop
    }

    if (!foundNumDataLines) {
        throw format_exception(IOError{}, "IO format error -> number of data lines not found!");
    }

    std::map<std::size_t, DataType> data;

    // Read the data lines (integer followed by double)
    int intValue;
    double doubleValue;
    for (int i = 0; i < numDataLines; ++i) {
        if (!(gridDataFile >> intValue >> doubleValue)) {
            throw format_exception(IOError{}, "IO format error -> error reading the dataformat");
        }

        data[intValue-1] = doubleValue; // Gmsh starts at 1 -> The current index starts at 0 !!!!
    }
    gridDataFile.close();
    // add the data to the cell data manager
    _cell_data.addData(key, data);
  };

  // ---------------------------------------------------------------------------
  // disable functionality if grid is not read from Gmsh path
  bool has_griddata = false;
  // the underlying cell data
  CellData<typename GV::Grid::LevelGridView, DataType> _cell_data;
};

} // namespace Dune::Copasi

#endif // DUNE_COPASI_GRID_DATA_CONTEXT_HH
