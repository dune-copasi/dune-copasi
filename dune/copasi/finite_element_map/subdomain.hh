#ifndef DUNE_COPASI_MULTIDOMAIN_LOCAL_FINITE_ELEMENT_MAP_HH
#define DUNE_COPASI_MULTIDOMAIN_LOCAL_FINITE_ELEMENT_MAP_HH

#include <dune/copasi/common/factory.hh>
#include <dune/copasi/common/data_context.hh>
#include <dune/copasi/finite_element_map/dynamic_power.hh>
#include <dune/copasi/grid/has_single_geometry_type.hh>

namespace Dune::Copasi {

/**
 * @brief      This class describes a sub domain local finite element map.
 * @details    This class wrapps a usual PDELab finite element map into a
 *             dynamic power finite element map. If the entity to be map does
 *             not belong to the grid view, it will return a 0 power element map
 *             that turns to return a finite element with no degrees of freedom.
 *             This behaviour is useful when the pdelab machinary is operating
 *             in the whole grid (e.g. multidomain grid) but you want to have
 *             different finite elements per sub domain.
 * @ingroup    FiniteElementMap
 * 
 * @tparam     FiniteElementMap  The original finite element map to wrap
 * @tparam     SubGridView       The grid view where the finite element will be
 *                               not zero
 */
template<class FiniteElementMap, class SubGridView>
class SubDomainLocalFiniteElementMap
  : public DynamicPowerLocalFiniteElementMap<FiniteElementMap>
{
  using Base = DynamicPowerLocalFiniteElementMap<FiniteElementMap>;

  using BaseFiniteElement = typename FiniteElementMap::Traits::FiniteElement;
  using FiniteElement = DynamicPowerLocalFiniteElement<BaseFiniteElement>;

  static_assert(Concept::isSubDomainGrid<typename SubGridView::Grid>(),
                "This class is only meant to be used for multidomain grids and "
                "its subdomains");

  using SubDomainGrid = typename SubGridView::Grid;
  using MultiDomainGrid = typename SubDomainGrid::MultiDomainGrid;

public:
  /**
   * @brief      Constructs a new instance.
   *
   * @param[in]  grid_view            The sub grid view
   * @param[in]  fem                  The finite element map to wrap
   * @param[in]  base_finite_element  The base finite element of the finite
   *                                  element map
   * @param[in]  power_size           The power size for the sub domain part
   *                                  (this is always 0 outside of the
   *                                  subdomain)
   */
  SubDomainLocalFiniteElementMap(const SubGridView& grid_view,
                                   const FiniteElementMap& fem,
                                   const BaseFiniteElement& base_finite_element,
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
  SubDomainLocalFiniteElementMap(const SubGridView& grid_view,
                                   const FiniteElementMap& fem,
                                   std::size_t power_size = 1)
    : SubDomainLocalFiniteElementMap(grid_view,
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
    } else {
      static_assert(Dune::AlwaysFalse<EntityType>::value,
                    "Onle mult/sub-domain entities are allowed");
    }

    return in_grid_view ? Base::find(e) : _fe_null;
  }

  /**
   * @brief      Returns true if this finite element map has a fixed size
   *
   * @return     Always false for this type
   */
  static bool constexpr fixedSize() { return false; }

private:
  SubGridView _grid_view;
  const FiniteElement _fe_null;
};

/**
 * @brief      Factory for SubDomainLocalFiniteElementMap instances
 * @ingroup    Factory, FiniteElementMap
 * @tparam     <unnamed>  Template paramenters of the SubDomainLocalFiniteElementMap
 */
template<class BaseLocalFiniteElementMap, class SubGridView>
struct Factory<SubDomainLocalFiniteElementMap<BaseLocalFiniteElementMap,SubGridView>>
{
  /**
   * @brief      Create method
   *
   * @param      ctx   @ref DataContext containing a grid view of the type SubGridView and
   *                   sufficient data to create the base finite element map and its base finite element
   *                   from another factory
   *
   * @tparam     Ctx   Universal reference to the @ref DataContext
   *
   * @return     Instance of SubDomainLocalFiniteElementMap
   */
  template<class Ctx>
  static auto create(Ctx&& ctx)
  {
    using dCtx = std::decay_t<Ctx>;
    static_assert(dCtx::has( Context::Tag<SubGridView>{} ));

    const auto& sub_domain_gv = ctx.view( Context::Tag<SubGridView>{} );

    // create a multi domain context
    auto multi_domain_gv = sub_domain_gv.grid().multiDomainGrid().leafGridView();
    using MultiDomainGridView = std::decay_t<decltype(multi_domain_gv)>;

    auto multi_domain_ctx = Context::DataContext<MultiDomainGridView,dCtx>(multi_domain_gv,std::forward<Ctx>(ctx));

    // base finite element and finite element map are created with multidomain contex
    using BaseFE = typename BaseLocalFiniteElementMap::Traits::FiniteElement;
    auto base_fe = Factory<BaseFE>::create(multi_domain_ctx);

    auto base_fem = Factory<BaseLocalFiniteElementMap>::create(std::move(multi_domain_ctx));
    using FEM = SubDomainLocalFiniteElementMap<BaseLocalFiniteElementMap,SubGridView>;

    // final finite element map is created with sub domain context
    return std::make_unique<FEM>(sub_domain_gv,*base_fem,*base_fe);
  }
};

} // namespace Dune::Copasi

#endif // DUNE_COPASI_MULTIDOMAIN_LOCAL_FINITE_ELEMENT_MAP_HH
