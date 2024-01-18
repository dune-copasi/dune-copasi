#ifndef DUNE_COPASI_PARSER_GRID_CONTEXT_HH
#define DUNE_COPASI_PARSER_GRID_CONTEXT_HH

#include <dune/copasi/common/exceptions.hh>
#include <dune/copasi/parser/factory.hh>

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

  // helper class to obtain HostEntity
  template<Dune::Concept::Grid MDGrid>
  class ParserGridMapper
  {

  public:

    using SDGrid = typename MDGrid::SubDomainGrid;
    using HostGrid = typename MDGrid::HostGrid;

    using HostEntity = typename HostGrid::template Codim<0>::Entity;
    using SDEntity = typename SDGrid::template Codim<0>::Entity;
    using MDEntity = typename MDGrid::template Codim<0>::Entity;

    ParserGridMapper(std::shared_ptr<MDGrid> md_grid_ptr): _md_grid_ptr(md_grid_ptr) {};

    const HostEntity& getHostEntity(const SDEntity& entity) {
        return _md_grid_ptr->hostEntity(entity);
    };

    const HostEntity& getHostEntity(const MDEntity&  entity) {
        return _md_grid_ptr->hostEntity(entity);
    };

    const HostEntity& getHostEntity(const HostEntity& entity) {
        return entity;
    };

  private:
    std::shared_ptr<const MDGrid> _md_grid_ptr = nullptr;

  };

// Class that manages the grid data
template<Dune::Concept::Grid MDGrid>
class ParserGridContext
{

public:

  using SDGrid = typename MDGrid::SubDomainGrid;
  using HostGrid = typename MDGrid::HostGrid;

  using HostEntity = typename HostGrid::template Codim<0>::Entity;
  using SDEntity = typename SDGrid::template Codim<0>::Entity;
  using MDEntity = typename MDGrid::template Codim<0>::Entity;

  // This is before the grid is constructed ---> only generic parameters
  explicit ParserGridContext(const ParameterTree& config = {}) {};

  void configure(Dune::ParameterTree& grid_config, std::shared_ptr<MDGrid> md_grid_ptr, std::shared_ptr<HostGrid> host_grid_ptr){

    _host_grid_ptr = host_grid_ptr;
    _md_grid_ptr = md_grid_ptr;

    _parser_grid_mapper = std::make_shared<ParserGridMapper<MDGrid>>(md_grid_ptr);

    if (grid_config.hasKey("path")) {
        auto grid_path = grid_config.template get<std::string>("path");

        auto host_grid_factory = Dune::GridFactory<HostGrid>{};
        auto host_grid_parser = Dune::GmshReaderParser<HostGrid>{ host_grid_factory, false, false };

        host_grid_parser.read(grid_path);
        auto& gmsh_physical_entity = host_grid_parser.elementIndexMap();

        auto mcmg = MultipleCodimMultipleGeomTypeMapper<typename HostGrid::LeafGridView>(_host_grid_ptr->leafGridView(), mcmgElementLayout());

        _update_gmsh_id = [gmsh_physical_entity, mcmg](const HostEntity& entity) {
            HostEntity _entity = entity;
            while (_entity.hasFather()) {
              _entity = _entity.father();
            }
            assert(_entity.level() == 0);
            return gmsh_physical_entity[mcmg.index(_entity)];
          };

        _update_index = [gmsh_physical_entity, mcmg](const HostEntity& entity) {
            HostEntity _entity = entity;
            while (_entity.hasFather()) {
              _entity = _entity.father();
            }
            assert(_entity.level() == 0);
            return mcmg.index(_entity);
          };

        // load grid data
        load_grid_data(grid_config);

    }
  };

  void add_grid_context(Parser& parser) const;

  // update the griddata
  void update_grid_data(std::unordered_map<std::string, double>& cell_data, const HostEntity& entity) const{

    for (auto& node : cell_data){
      // obtain the functor to update value
      auto it = _cell_data.find(node.first);
      if(it != _cell_data.end()){
        std::size_t index = _update_index(entity);
        const double* ptr = it->second(index);
        cell_data[node.first] = (ptr!=nullptr) ? *(ptr) : 0.0;
      }else{
        cell_data[node.first] = 0.0;
      }
    }
    
  }


  std::function<double (const HostEntity& entity)> get_gmsh_id() const{
    return _update_gmsh_id;
  }

  std::shared_ptr<ParserGridMapper<MDGrid>> get_parser_grid_mapper() const {
    return _parser_grid_mapper;
  }

  const std::unordered_map<std::string, std::function<double*(std::size_t)>>& get_cell_data() const{
    return _cell_data;
  }

private:

  // Load the grid data files based on the config file parameter tree
  void load_grid_data(Dune::ParameterTree& grid_config){

    auto& griddata_config = grid_config.sub("griddata");

    for (const auto& sub : griddata_config.getSubKeys()) {
      const auto& sub_config = griddata_config.sub(sub, true);
      const std::string type = sub_config["type"];
      if (type == "scalar") {

        const std::string datafile_path = sub_config["path"];
        std::cout << "datafile path: " << datafile_path << std::endl;

        // unorder map uses a hashtable to search for data element
        std::unordered_map<int, double> data;

        // ---> create separate function ---<
        // ======================= FILE IO =====================================
        // read the data file:
        std::ifstream gridDataFile;
        gridDataFile.open(datafile_path);

        if(!gridDataFile.is_open()){
          throw format_exception(IOError{}, "Could not open datafile");
        }

        std::string line;
        int numDataLines = 0;
        bool foundNumDataLines = false;

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

        data.reserve(numDataLines);

        // Read the data lines (integer followed by double)
        int intValue;
        double doubleValue;
        for (int i = 0; i < numDataLines; ++i) {
            if (!(gridDataFile >> intValue >> doubleValue)) {
                throw format_exception(IOError{}, "IO format error -> error readig the dataformat");
            }

            data[intValue-1] = doubleValue; // Gmsh starts at 1 -> The current index starts at 0 !!!!
        }

        gridDataFile.close();
        // ======================= END FILE IO =================================

        // create a functor to access the unordered map with hashtable to gaurantee O(1) look-up time
        std::function<double*(std::size_t)> cell_data_functor = [data = std::move(data)](std::size_t index) mutable  -> double*{
          if( auto it = data.find(index); it != data.end())
            return &(it->second);
          else
            return nullptr;
        };

        std::cout << "size of data: " << data.size() << std::endl;

        _cell_data[sub] = cell_data_functor;

      } else if (type == "vector") {
        throw format_exception(NotImplemented{}, "vector griddata is not implemented");
      } else if (type == "tensor") {
        throw format_exception(NotImplemented{}, "tensor griddata is not implemented");
      } else {
        throw format_exception(IOError{}, "Unknown type {}", type);
      }
    }
  };

  std::shared_ptr<const MDGrid> _md_grid_ptr = nullptr;
  std::shared_ptr<const HostGrid> _host_grid_ptr = nullptr;
  std::shared_ptr<ParserGridMapper<MDGrid>> _parser_grid_mapper;
  std::function<double (const HostEntity& entity)> _update_gmsh_id;
  // a lambda giving the index based on the entity
  std::function<std::size_t (const HostEntity& entity)> _update_index;
  // unordered map giving a lambda pointing to the cell data
  std::unordered_map<std::string, std::function<double*(std::size_t)>> _cell_data;
};

} // namespace Dune::Copasi

#endif // DUNE_COPASI_PARSER_GRID_CONTEXT_HH
