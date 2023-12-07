#ifndef DUNE_COPASI_PARSER_GRID_CONTEXT_HH
#define DUNE_COPASI_PARSER_GRID_CONTEXT_HH

#include <dune/copasi/parser/factory.hh>

#include <dune/grid/common/exceptions.hh>
#include <dune/grid/common/gridinfo.hh>
#include <dune/grid/common/mcmgmapper.hh>
#include <dune/grid/common/rangegenerators.hh>
#include <dune/grid/io/file/gmshreader.hh>
#include <dune/grid/uggrid/uggridfactory.hh>
#include <dune/grid/utility/structuredgridfactory.hh>

#include <dune/common/fvector.hh>
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
  explicit ParserGridContext(const ParameterTree& config = {}) {

  };

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

    }
  };

  void add_grid_context(Parser& parser) const;

  std::function<double (const HostEntity& entity)> get_gmsh_id() const{
    return _update_gmsh_id;
  }

  std::shared_ptr<ParserGridMapper<MDGrid>> get_parser_grid_mapper() const {
    return _parser_grid_mapper;
  }

private:

  std::shared_ptr<const MDGrid> _md_grid_ptr = nullptr;
  std::shared_ptr<const HostGrid> _host_grid_ptr = nullptr;
  std::shared_ptr<ParserGridMapper<MDGrid>> _parser_grid_mapper;
  std::function<double (const HostEntity& entity)> _update_gmsh_id;
};

} // namespace Dune::Copasi

#endif // DUNE_COPASI_PARSER_GRID_CONTEXT_HH
