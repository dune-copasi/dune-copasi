#ifndef DUNE_COPASI_MULTIDOMAIN_LOCAL_FINITE_ELEMENT_MAP_HH
#define DUNE_COPASI_MULTIDOMAIN_LOCAL_FINITE_ELEMENT_MAP_HH

#include <dune/pdelab/finiteelementmap/finiteelementmap.hh>

#include <dune/copasi/dynamic_local_finite_element.hh>

namespace Dune::Copasi {

template<class FiniteElementMap, class GridView>
class MultiDomainLocalFiniteElementMap
  : public PDELab::LocalFiniteElementMapInterface<
      PDELab::LocalFiniteElementMapTraits<DynamicPowerLocalFiniteElement<
        typename FiniteElementMap::Traits::FiniteElement>>,
      MultiDomainLocalFiniteElementMap<FiniteElementMap, GridView>>
{
  using BaseFiniteElement = typename FiniteElementMap::Traits::FiniteElement;
  using FiniteElement = DynamicPowerLocalFiniteElement<BaseFiniteElement>;

public:
  static const int dimension = FiniteElementMap::dimension;

  MultiDomainLocalFiniteElementMap(GridView grid_view,
                                   FiniteElementMap fem,
                                   BaseFiniteElement base_finite_element,
                                   std::size_t power_size = 1)
    : _grid_view(grid_view)
    , _power_size(power_size)
    , _fem(fem)
    , _fe_cache(NULL)
    , _fe_null(base_finite_element, 0)
  {}

  template<
    typename =
      std::enable_if_t<std::is_default_constructible_v<BaseFiniteElement>, int>>
  MultiDomainLocalFiniteElementMap(GridView grid_view,
                                   FiniteElementMap fem,
                                   std::size_t power_size = 1)
    : MultiDomainLocalFiniteElementMap(grid_view,
                                       fem,
                                       BaseFiniteElement{},
                                       power_size)
  {}

  ~MultiDomainLocalFiniteElementMap() { delete _fe_cache; }

  template<class EntityType>
  const FiniteElement& find(const EntityType& e) const
  {
    if constexpr (Concept::isSubDomainGrid<typename GridView::Grid>()) {
      using SubDomainGrid = typename GridView::Grid;
      using MultiDomainGrid = typename SubDomainGrid::MultiDomainGrid;

      constexpr int codim = EntityType::codimension;

      auto sub_domain = _grid_view.grid().domain();

      using MultiDomainEntity =
        typename MultiDomainGrid::template Codim<codim>::Entity;
      if constexpr (std::is_same_v<EntityType, MultiDomainEntity>) {
        bool in_grid_view = _grid_view.grid()
                              .multiDomainGrid()
                              .leafIndexSet()
                              .subDomains(e)
                              .contains(sub_domain);
        if (not in_grid_view)
          return _fe_null;
      }
    }
    auto base_fe = _fem.find(e);
    // cache the last used base finite element
    if (_base_fe_cache != &base_fe) // TODO: Cache more than one fe
    {
      _base_fe_cache = &base_fe;
      if (_fe_cache != NULL)
        delete _fe_cache;
      _fe_cache = new FiniteElement(base_fe, _power_size);
    }
    return *_fe_cache;
  }

  bool fixedSize() const { return false; }

  std::size_t size(GeometryType gt) const
  {
    return _power_size * _fem.size(gt);
  }

  bool hasDOFs(int codim) const
  {
    return (_power_size != 0) && _fem.hasDOFs(codim);
  }

  std::size_t maxLocalSize() const { return _power_size * _fem.maxLocalSize(); }

private:
  GridView _grid_view;
  std::size_t _power_size;
  FiniteElementMap _fem;
  mutable BaseFiniteElement* _base_fe_cache;
  mutable FiniteElement* _fe_cache;
  const FiniteElement _fe_null;
};

} // namespace Dune::Copasi

#endif // DUNE_COPASI_MULTIDOMAIN_LOCAL_FINITE_ELEMENT_MAP_HH
