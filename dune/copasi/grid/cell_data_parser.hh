#ifndef DUNE_COPASI_GRID_CELL_DATA_PARSER_HH
#define DUNE_COPASI_GRID_CELL_DATA_PARSER_HH

#include <dune/copasi/common/exceptions.hh>
#include <dune/copasi/grid/cell_data.hh>
#include <dune/copasi/concepts/grid.hh>

#include <dune/grid/common/rangegenerators.hh>

#include <spdlog/spdlog.h>

#include <dune/common/exceptions.hh>
#include <dune/common/parametertree.hh>

#include <vector>
#include <fstream>

namespace Dune::Copasi {

namespace Impl {

  // Read a datafile for the scalar griddata type & add it to cell_data
  template<Dune::Concept::GridView GV, class T>
  void cell_data_scalar_parser(CellData<GV,T>& cell_data, std::string datafile_path, std::string key) {
    spdlog::info("Reading grid cell data '{}'", datafile_path);

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

    std::map<std::size_t, T> data;

    // Read the data lines (integer followed by double)
    int intValue;
    double doubleValue;
    for (int i = 0; i < numDataLines; ++i) {
        if (!(gridDataFile >> intValue >> doubleValue)) {
            throw format_exception(IOError{}, "IO format error -> error reading the dataformat");
        }

        // Note that this follows 0-base index.
        // If user wants gmsh compatible indices, (s)he needs to subtract 1.
        data[intValue] = doubleValue;
    }
    gridDataFile.close();
    // add the data to the cell data manager
    cell_data.addData(key, data);
  }
}

/**
 * @brief Parse data file with cell data and add its values to cell data
 * @details The format of the file is documented in the `config_opts.json`.
 * In particular, the cell index set in the file will be assigned to the
 * index set of the grid view used to construct the cell data object. This means
 * that this is mostly useful for leaf grid views and level 0 grid vies.
 * Be careful to choose the correct one!
 * 
 * @author Dylan Vermoortele
 * 
 * @tparam GV  The grid view type
 * @tparam T   The type of the cell data values
 * @param grid_config  A config containing sub-sections with the paths for each cell data file 
 * @param cell_data    A cell data object to store the parsed cell data
 */
template<Dune::Concept::GridView GV, class T>
void cell_data_parser(const ParameterTree& grid_config, CellData<GV,T>& cell_data) {

  auto& cell_data_config = grid_config.sub("cell_data");

  // reserve the data needed to store the griddata
  // --> Every element can contain cell data for every key
  // --> Fixed look up table addresses providing (potentially) highest performance
  // --> Cell data may already contain data (e.g. gmsh_id)
  cell_data.reserve(cell_data.size() + cell_data_config.getSubKeys().size());

  for (const auto& sub : cell_data_config.getSubKeys()) {
    const auto& sub_config = cell_data_config.sub(sub, true);
    const std::string datafile_path = sub_config["path"];
    const std::string type = sub_config["type"];
    if (type == "scalar") {
      Impl::cell_data_scalar_parser(cell_data, datafile_path, sub);
    } else if (type == "vector") {
      throw format_exception(NotImplemented{}, "vector griddata is not implemented");
    } else if (type == "tensor") {
      throw format_exception(NotImplemented{}, "tensor griddata is not implemented");
    } else {
      throw format_exception(IOError{}, "Unknown type {}", type);
    }
  }
}

} // namespace Dune::Copasi

#endif // DUNE_COPASI_GRID_CELL_DATA_PARSER_HH
