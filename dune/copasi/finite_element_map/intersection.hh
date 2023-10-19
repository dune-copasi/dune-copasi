#ifndef DUNE_COPASI_INTERSECTION_LOCAL_FINITE_ELEMENT_MAP_HH
#define DUNE_COPASI_INTERSECTION_LOCAL_FINITE_ELEMENT_MAP_HH

#include <dune/copasi/common/factory.hh>
#include <dune/copasi/common/data_context.hh>
#include <dune/copasi/finite_element_map/dynamic_power.hh>
#include <dune/copasi/grid/has_single_geometry_type.hh>

#include <dune/grid/common/mcmgmapper.hh>

namespace Dune::Copasi {

template<class FiniteElementMap, class MultiDomainGridView>
class IntersectionLocalFiniteElementMap
  : public DynamicPowerLocalFiniteElementMap<FiniteElementMap>
{
  using Base = DynamicPowerLocalFiniteElementMap<FiniteElementMap>;
  using MultiDomainGrid = typename MultiDomainGridView::Grid;
  using SubDomainGrid = typename MultiDomainGrid::SubDomainGrid;
  using Mapper = Dune::MultipleCodimMultipleGeomTypeMapper<MultiDomainGridView>;

public:
  using BaseFiniteElement = typename FiniteElementMap::Traits::FiniteElement;
  using FiniteElement = DynamicPowerLocalFiniteElement<BaseFiniteElement>;

  static_assert(Concept::isMultiDomainGrid<typename MultiDomainGridView::Grid>(),
                "This class is only meant to be used for multidomain grids and "
                "its subdomains");


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
  IntersectionLocalFiniteElementMap(const MultiDomainGridView& grid_view,
                                   const FiniteElementMap& fem,
                                   std::size_t domain_i, std::size_t domain_o,
                                   const BaseFiniteElement& base_finite_element,
                                   std::size_t power_size = 1)
    : Base(fem, power_size)
    , _domain_i{domain_i}
    , _domain_o{domain_o}
    , _grid_view(grid_view)
    , _fe_null(base_finite_element, 0)
    , _mapper{ _grid_view, Dune::mcmgLayout(Dune::Codim<1>{}) }
  {
    _boundary.assign(_mapper.size(), false);
    for (const auto& element : elements(_grid_view))
      for (const auto& intersection : intersections(_grid_view, element)) {
        auto subEntity = element.template subEntity<1>(intersection.indexInInside());
        _boundary[_mapper.index(subEntity)] = intersection.boundary();
      }
  }

  /**
   * @brief      Constructs a new instance.
   * @details    This constructor is only available if the base finite element
   *             is default constructible
   *
   * @param[in]  grid_view   The grid view
   * @param[in]  fem         The finite element map to wrap
   * @param[in]  power_size  The power size for the sub domain part (this is
   *                         always 0 outside of the intersection)
   *
   * @tparam     <unnamed>   Template helper to allow default construction of
   *                         base finite element constructor. Internal use only
   */
  template<
    bool default_constructible = std::is_default_constructible_v<BaseFiniteElement>,
    class = std::enable_if_t<default_constructible>>
  IntersectionLocalFiniteElementMap(const MultiDomainGridView& grid_view,
                                   const FiniteElementMap& fem,
                                   std::size_t domain_i, std::size_t domain_o,
                                   std::size_t power_size = 1)
    : IntersectionLocalFiniteElementMap(grid_view,
                                       fem,
                                       domain_i, domain_o,
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
    static_assert(EntityType::codimension == 1);
    using MultiDomainEntity = typename MultiDomainGrid::template Codim<1>::Entity;
    using SubDomainEntity = typename SubDomainGrid::template Codim<1>::Entity;

    const auto& md_entity = [&](){
      if constexpr (std::is_same_v<EntityType, MultiDomainEntity>)
        return e;
      else if constexpr (std::is_same_v<EntityType, SubDomainEntity>)
        return _grid_view.grid().multiDomainEntity(e);
      else
        static_assert(AlwaysFalse<EntityType>{}, "not known type");
    }();

    if (_domain_i == _domain_o)
      return _boundary[_mapper.index(md_entity)] ? Base::find(md_entity) : _fe_null;
    else {
      const auto& domain_set = _grid_view.indexSet().subDomains(md_entity);
      bool in_membrane = domain_set.contains(_domain_i) and domain_set.contains(_domain_o);
      return in_membrane ? Base::find(md_entity) : _fe_null;
    }
  }

  /**
   * @brief      Returns true if this finite element map has a fixed size
   *
   * @return     Always false for this type
   */
  static bool constexpr fixedSize() { return false; }

private:
  std::size_t _domain_i, _domain_o;
  MultiDomainGridView _grid_view;
  const FiniteElement _fe_null;
  Mapper _mapper;
  std::vector<bool> _boundary;
};

// /**
//  * @brief      Factory for IntersectionLocalFiniteElementMap instances
//  * @ingroup    Factory, FiniteElementMap
//  * @tparam     <unnamed>  Template paramenters of the IntersectionLocalFiniteElementMap
//  */
// template<class BaseLocalFiniteElementMap, class MultiDomainGridView>
// struct Factory<IntersectionLocalFiniteElementMap<BaseLocalFiniteElementMap,MultiDomainGridView>>
// {
//   /**
//    * @brief      Create method
//    *
//    * @param      ctx   @ref DataContext containing a grid view of the type MultiDomainGridView and
//    *                   sufficient data to create the base finite element map and its base finite element
//    *                   from another factory
//    *
//    * @tparam     Ctx   Universal reference to the @ref DataContext
//    *
//    * @return     Instance of IntersectionLocalFiniteElementMap
//    */
//   template<class Ctx>
//   static auto create(Ctx&& ctx)
//   {
//     using dCtx = std::decay_t<Ctx>;
//     static_assert(dCtx::has( Context::Tag<MultiDomainGridView>{} ));

//     const auto& sub_domain_gv = ctx.view( Context::Tag<MultiDomainGridView>{} );

//     // create a multi domain context
//     auto multi_domain_gv = sub_domain_gv.grid().multiDomainGrid().leafGridView();
//     using MultiDomainGridView = std::decay_t<decltype(multi_domain_gv)>;

//     auto multi_domain_ctx = Context::DataContext<MultiDomainGridView,dCtx>(multi_domain_gv,std::forward<Ctx>(ctx));

//     // base finite element and finite element map are created with multidomain contex
//     using BaseFE = typename BaseLocalFiniteElementMap::Traits::FiniteElement;
//     auto base_fe = Factory<BaseFE>::create(multi_domain_ctx);

//     auto base_fem = Factory<BaseLocalFiniteElementMap>::create(std::move(multi_domain_ctx));
//     using FEM = IntersectionLocalFiniteElementMap<BaseLocalFiniteElementMap,MultiDomainGridView>;

//     // final finite element map is created with sub domain context
//     return std::make_unique<FEM>(sub_domain_gv,*base_fem,*base_fe);
//   }
// };

} // namespace Dune::Copasi

#endif // DUNE_COPASI_INTERSECTION_LOCAL_FINITE_ELEMENT_MAP_HH
