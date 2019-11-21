#ifndef DUNE_COPASI_MULTIDOMAIN_LOCAL_FINITE_ELEMENT_MAP_HH
#define DUNE_COPASI_MULTIDOMAIN_LOCAL_FINITE_ELEMENT_MAP_HH

#include <dune/copasi/common/factory.hh>
#include <dune/copasi/finite_element_map/dynamic_power.hh>
#include <dune/copasi/grid/has_single_geometry_type.hh>

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
  : public DynamicPowerLocalFiniteElementMap<FiniteElementMap>
{
  using Base = DynamicPowerLocalFiniteElementMap<FiniteElementMap>;

  using BaseFiniteElement = typename FiniteElementMap::Traits::FiniteElement;
  using FiniteElement = DynamicPowerLocalFiniteElement<BaseFiniteElement>;

  using SubDomainGrid = typename GridView::Grid;
  using MultiDomainGrid = typename SubDomainGrid::MultiDomainGrid;

  static_assert(Concept::isSubDomainGrid<typename GridView::Grid>(),
                "This class is only meant to be used for multidomain grids and "
                "its subdomains");

public:
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
  MultiDomainLocalFiniteElementMap(const GridView& grid_view,
                                   FiniteElementMap fem,
                                   BaseFiniteElement base_finite_element,
                                   std::size_t power_size = 1)
    : Base(fem, power_size)
    , _grid_view(grid_view)
    , _fe_null(base_finite_element, 0)
  {}

  /**
   * @brief      Constructs a new instance.
   * @details    This constructor is only available if the base finite element
   *             is default constructible
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
    bool default_constructible = std::is_default_constructible_v<BaseFiniteElement>,
    class = std::enable_if_t<default_constructible>>
  MultiDomainLocalFiniteElementMap(const GridView& grid_view,
                                   FiniteElementMap fem,
                                   std::size_t power_size = 1)
    : MultiDomainLocalFiniteElementMap(grid_view,
                                       fem,
                                       BaseFiniteElement{},
                                       power_size)
  {}

  /**
   * @brief      Searches for the finite element for entity e.
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
    constexpr int codim = EntityType::codimension;

    using MultiDomainEntity =
      typename MultiDomainGrid::template Codim<codim>::Entity;

    using SubDomainEntity =
      typename SubDomainGrid::template Codim<codim>::Entity;

    auto sub_domain = _grid_view.grid().domain();
    const auto& md_grid = _grid_view.grid().multiDomainGrid();

    bool in_grid_view = true;

    if constexpr (std::is_same_v<EntityType, MultiDomainEntity>) {
      in_grid_view = md_grid.leafIndexSet().subDomains(e).contains(sub_domain);
    } else if constexpr (std::is_same_v<EntityType, SubDomainEntity>) {
      in_grid_view = md_grid.leafIndexSet()
                       .subDomains(md_grid.multiDomainEntity(e))
                       .contains(sub_domain);
    }

    if (in_grid_view)
      return Base::find(e);
    else
      return _fe_null;
  }

  /**
   * @brief      Returns true if this finite element map has a fixed size
   *
   * @return     Always false for this type
   */
  bool fixedSize() const { return false; }

private:
  GridView _grid_view;
  const FiniteElement _fe_null;
};

template<class BaseLocalFiniteElementMap, class GridView>
struct Factory<MultiDomainLocalFiniteElementMap<BaseLocalFiniteElementMap,GridView>>
{
public:
  template<class Ctx>
  static auto create(const Ctx& ctx)
  {
    static_assert(Concept::has_method_grid_view<Ctx>());
    static_assert(Concept::isSubDomainGrid<typename Ctx::GridView::Grid>());

    auto sub_domain_gv = ctx.grid_view().grid().multiDomainGrid().leafGridView();
    auto sub_somain_ctx = Context::make_grid_view(sub_domain_gv);
    auto base_fem = Factory<BaseLocalFiniteElementMap>::create(sub_somain_ctx);
    using FEM = MultiDomainLocalFiniteElementMap<BaseLocalFiniteElementMap,GridView>;

    using BaseFE = typename BaseLocalFiniteElementMap::Traits::FiniteElement;

    if (not has_single_geometry_type(ctx.grid_view()))
      DUNE_THROW(InvalidStateException,"Grid view has to have only one grid view");

    GeometryType gt = ctx.grid_view().template begin<0>()->geometry().type();

    auto gt_ctx = Context::inject_geometry_type(Ctx{ctx},gt);
    auto base_fe = Factory<BaseFE>::create(gt_ctx);

    if constexpr (Concept::has_method_power_size<Ctx>())
      return std::make_unique<FEM>(gt_ctx.grid_view(),*base_fem,*base_fe,gt_ctx.power_size());
    else
      return std::make_unique<FEM>(gt_ctx.grid_view(),*base_fem,*base_fe);
  }
};

} // namespace Dune::Copasi

#endif // DUNE_COPASI_MULTIDOMAIN_LOCAL_FINITE_ELEMENT_MAP_HH