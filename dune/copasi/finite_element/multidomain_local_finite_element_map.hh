#ifndef DUNE_COPASI_MULTIDOMAIN_LOCAL_FINITE_ELEMENT_MAP_HH
#define DUNE_COPASI_MULTIDOMAIN_LOCAL_FINITE_ELEMENT_MAP_HH

#include <dune/copasi/finite_element/dynamic_local_finite_element.hh>

#include <dune/pdelab/finiteelementmap/finiteelementmap.hh>

namespace Dune::Copasi {

/**
 * @brief      This class describes a multi domain local finite element map.
 * @details    This class wrapps a usual PDELab finite element map into a
 *             dynamic power finite element map. If the entity to be map does
 *             not belong to the grid view, it will return a 0 power element map
 *             that turns to return a finite element with no degrees of freedom.
 *             This behaviour is useful when the pdelab machinary is operating
 *             in the whole grid (e.g. multidomain grid) but you want to have
 *             different finite elements per sub domain.
 *
 * @tparam     FiniteElementMap  The original finite element map to wrap
 * @tparam     GridView          The grid view where the finite element will be
 *                               not zero
 */
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

  /**
   * @brief      Constructs a new instance.
   *
   * @param[in]  grid_view            The grid view (may be a sub domain)
   * @param[in]  fem                  The finite element map to wrap
   * @param[in]  base_finite_element  The base finite element of the finite
   *                                  element map
   * @param[in]  power_size           The power size for the sub domain part
   *                                  (this is always 0 outside of the
   *                                  subdomain)
   */
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

  /**
   * @brief      Constructs a new instance.
   * @details    This constructor is only available if the base finite element
   *             is default constructiible
   *
   * @param[in]  grid_view   The grid view
   * @param[in]  fem         The finite element map to wrap
   * @param[in]  power_size  The power size for the sub domain part (this is
   *                         always 0 outside of the subdomain)
   *
   * @tparam     <unnamed>   Template helper to allow default construction of
   *                         base finite element constructor. Internal use only
   */
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

  /**
   * @brief      Destroys the object.
   */
  ~MultiDomainLocalFiniteElementMap() { delete _fe_cache; }

  /**
   * @brief      Searches for the finite element for entity e.
   * @todo       Cache more than one finite element
   *
   * @param[in]  e           The entity
   *
   * @tparam     EntityType  The entity
   *
   * @return     A dynamic power local finite element
   *             @DynamicPowerLocalFiniteElement
   */
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
      } else {
        DUNE_THROW(NotImplemented,
                   "Method not implemented for subdomain entites");
      }
    }
    auto base_fe = _fem.find(e);

    // cache the last used base finite element
    if (_base_fe_cache != &base_fe) {
      _base_fe_cache = &base_fe;
      if (_fe_cache != NULL)
        delete _fe_cache;
      _fe_cache = new FiniteElement(base_fe, _power_size);
    }
    return *_fe_cache;
  }

  /**
   * @brief      Returns true if this finite element map has a fixed size
   *
   * @return     Always false for this type
   */
  bool fixedSize() const { return false; }

  /**
   * @brief      Size for a given geometry type
   *
   * @param[in]  gt    The geometry type
   *
   * @return     The size for the given geometry type
   */
  std::size_t size(GeometryType gt) const
  {
    return _power_size * _fem.size(gt);
  }

  /**
   * @brief      Describes wheter a codim has associated degrees of freedom
   *
   * @param[in]  codim  The codim
   *
   * @return     True if codim has dregees of freedom
   */
  bool hasDOFs(int codim) const
  {
    return (_power_size != 0) && _fem.hasDOFs(codim);
  }

  /**
   * @brief      Maximum local size for all finite elements
   *
   * @return     The maximim size
   */
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
