#ifndef DUNE_COPASI_GRID_MAP_DATA_HANDLE_HH
#define DUNE_COPASI_GRID_MAP_DATA_HANDLE_HH

#include <dune/grid/common/datahandleif.hh>

namespace Dune::Copasi {

template<class IdSet, class Map>
class MapDataHandle
  : public Dune::CommDataHandleIF<MapDataHandle<IdSet,Map>, typename Map::mapped_type>
{
public:
  using DataType = typename Map::mapped_type;

public:

  bool contains (int dim, int codim) const
  {
    return _codim == codim;
  }

  bool fixedSize (int dim, int codim) const
  {
    return true;
  }

  template<class Entity>
  size_t size (Entity& entity) const
  {
    return 1;
  }

  template<class MessageBuffer, class Entity>
  void gather (MessageBuffer& buff, const Entity& entity) const
  {
    const auto& id = _id_set.id(entity);
    buff.write(_map.at(id));
  }

  template<class MessageBuffer, class Entity>
  void scatter (MessageBuffer& buff, const Entity& entity, size_t)
  {
    const auto& id = _id_set.id(entity);
    buff.read(_map[id]);
  }

  MapDataHandle (const IdSet& id_set, Map& map, const std::size_t& codim = 0)
    : _id_set(id_set)
    , _map(map)
    , _codim(codim)
  {}

private:
  const IdSet& _id_set;
  Map& _map;
  const std::size_t _codim;
};

} // namespace Dune::Copasi

#endif // DUNE_COPASI_GRID_MAP_DATA_HANDLE_HH
