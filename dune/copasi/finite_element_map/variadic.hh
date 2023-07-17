#ifndef DUNE_COPASI_VARIADIC_LOCAL_FINITE_ELEMENT_HH
#define DUNE_COPASI_VARIADIC_LOCAL_FINITE_ELEMENT_HH

#include <dune/copasi/common/factory.hh>
#include <dune/copasi/common/data_context.hh>
#include <dune/copasi/finite_element_map/virtual.hh>

#include <dune/localfunctions/common/virtualinterface.hh>
#include <dune/localfunctions/common/virtualwrappers.hh>

#include <array>
#include <memory>
#include <type_traits>

namespace Dune::Copasi {

/**
 * @brief      This class describes a variadic local finite element map.
 * @details    The idea is that this class models a local finite element map
 *             and returns one of the variadic types depending on the requested
 *             entity. For doing so, this class virtualize every argument of the
 *             variadic local finite element map and returns a (polymorphic)
 *             reference of the base interface depending on the entity mapper.
 *             The entity mapper is a function that receives an entity and returns an
 *             index that correspond to the i-th variadic local finite element map.
 *             For example, this class could be used to return Pk finite elements on simplices
 *             and Qk for cubes. Another example would be to return P0 on cubes and Pk on simplices.
 * @ingroup    FiniteElementMap
 *
 * @tparam     Entity                  Entity type
 * @tparam     LocalFiniteElementMaps  Local finite element maps
 */
template<class Entity, class... LocalFiniteElementMaps>
class VariadicLocalFiniteElementMap
{
  //! Base interface to return
  using Interface = std::common_type_t<typename VirtualLocalFiniteElementMapWrapper<LocalFiniteElementMaps, Entity>::Interface...>;

  //! Number of local finite elements
  static constexpr std::size_t _size = sizeof...(LocalFiniteElementMaps);

  //! size for static loops
  using _integral_size = std::integral_constant<std::size_t,_size>;

  /**
   * @brief      Makes an array of virtualized finite element map pointers from a tuple of finite element maps.
   *
   * @param      finite_element_maps  Tuple of unique pointers to finite element maps
   *
   * @tparam     FEMs                 Finite element map types
   *
   * @return     Array of unique pointers to virtualized finite element maps
   */
  template<class... FEMs>
  static std::array<std::unique_ptr<const Interface>,_size>
  make_unique_array_from_tuple(std::tuple<std::unique_ptr<FEMs>...>&& finite_element_maps)
  {
    using VariadicFEMs = std::tuple<std::decay_t<LocalFiniteElementMaps>...>;
    std::array<std::unique_ptr<const Interface>,_size>  fem_unique_array;
    Dune::Hybrid::forEach(Dune::range(_integral_size{}), [&](auto i) {
      auto& fem = std::get<i>(finite_element_maps);
      using FEM = std::decay_t<decltype(*fem)>;
      using VariadicFEM = std::tuple_element_t<i,VariadicFEMs>;
      static_assert(std::is_same_v<FEM,VariadicFEM>,
        "Arguments for the VariadicLocalFiniteElementMap constructor should be the same as the variadic types");
      static_assert(not std::is_polymorphic_v<FEM>);
      fem_unique_array[i] = std::make_unique<VirtualLocalFiniteElementMapWrapper<FEM, Entity>>(std::move(fem));
    });
    return fem_unique_array;
  }

  /**
   * @brief      Makes an array of virtualized finite element map pointers from finite element maps.
   *
   * @param      finite_element_maps  The finite element maps
   *
   * @tparam     FEMs                 The finite element map types
   *
   * @return     Array of unique pointers to virtualized finite element maps
   */
  template<class... FEMs>
  static std::array<std::unique_ptr<const Interface>,_size>
  make_unique_array_from_variadic(FEMs&&... finite_element_maps)
  {
    auto tuple = std::make_tuple(std::make_unique<std::decay_t<FEMs>>(std::forward<FEMs>(finite_element_maps))...);
    return make_unique_array_from_tuple(std::move(tuple));
  }

public:
  using Traits = typename Interface::Traits;

  //! dimension of finite element map
  static constexpr int dimension = std::common_type_t<std::integral_constant<int,LocalFiniteElementMaps::dimension>...>::value;

  /**
   * @brief      Constructs a new instance.
   *
   * @param[in]  entity_mapper        The entity mapper
   * @param      finite_element_maps  Array of virtualized finite element maps unique pointers
   *
   * @tparam     EntityMapper         A functor which when evaluated with an entity returns an (std::size_t) index
   */
  template<class EntityMapper>
  explicit VariadicLocalFiniteElementMap(const EntityMapper& entity_mapper, std::array<std::unique_ptr<const Interface>,_size>&& finite_element_maps)
    : _entity_mapper(entity_mapper)
    , _finite_element_maps(std::move(finite_element_maps))
  {}

  /**
   * @brief      Constructs a new instance.
   *
   * @param[in]  entity_mapper        The entity mapper
   * @param      finite_element_maps  Tuple of finite element maps unique pointers
   *
   * @tparam     EntityMapper         A functor which when evaluated with an entity returns an (std::size_t) index
   */
  template<class EntityMapper>
  explicit VariadicLocalFiniteElementMap(const EntityMapper& entity_mapper, std::tuple<std::unique_ptr<LocalFiniteElementMaps>...>&& finite_element_maps)
    : VariadicLocalFiniteElementMap(entity_mapper, make_unique_array_from_tuple(std::move(finite_element_maps)))
  {}

  /**
   * @brief      Constructs a new instance.
   *
   * @param[in]  entity_mapper        The entity mapper
   * @param      finite_element_maps  The finite element maps
   *
   * @tparam     EntityMapper         A functor which when evaluated with an entity returns an (std::size_t) index
   * @tparam     FEMs                 Finite element map types
   */
  template<class EntityMapper, class... FEMs>
  explicit VariadicLocalFiniteElementMap(const EntityMapper& entity_mapper, FEMs&&... finite_element_maps)
    : VariadicLocalFiniteElementMap(entity_mapper, make_unique_array_from_variadic(std::forward<FEMs>(finite_element_maps)...))
  {}

  //! Default destructor
  ~VariadicLocalFiniteElementMap() = default;

  /**
   * @brief      Searches for a finite element for entity e
   *
   * @param[in]  e     Entity
   *
   * @return     A finite element for i-th finite element map according to the entity mapper
   */
  const typename Traits::FiniteElementType&
  find (const typename Traits::EntityType& e) const
  {
    std::size_t index = _entity_mapper(e);
    assert(index < _size);
    assert(index >= 0);
    return _finite_element_maps[index]->find(e);
  }

  /**
   * @brief      Fixed size of the degrees of freedom wrt entities
   *
   * @return     False if there is more than one varaidic fem, otherwise theresult of the wrapped fem
   */
  bool fixedSize() const
  {
    return _size == 1 ? _finite_element_maps[0]->fixedSize() : false;
  }

  /**
   * @brief      Size for the geomertry type gt
   *
   * @param[in]  gt    Geometry type
   *
   * @return     The size if there is only one variadic fem, otherwise throw
   */
  std::size_t size(GeometryType gt) const
  {
    if constexpr (_size == 1)
      return _finite_element_maps[0]->size(gt);
    else
      DUNE_THROW(PDELab::FiniteElementMapError,
         "\tThis method should not be called for not fixed sizes");
  }

  /**
   * @brief      Determines if codim has degrees of freedom
   *
   * @param[in]  codim  The codim
   *
   * @return     true if codim has degrees of freedom for any variadic fem
   */
  bool hasDOFs(int codim) const
  {
    bool has_dofs = false;
    for (auto&& fem : _finite_element_maps)
      has_dofs |= fem->hasDOFs(codim);
    return has_dofs;
  }

  /**
   * @brief      Max local size for this finite element
   *
   * @return     Max local size of the variadic maxLocalSizes
   */
  std::size_t maxLocalSize() const
  {
    std::size_t max_local_size(0);
    for (auto&& fem : _finite_element_maps)
      max_local_size = std::max(max_local_size,fem->maxLocalSize());
    return max_local_size;
  }

private:
  const std::function<std::size_t(const typename Traits::EntityType&)> _entity_mapper;
  std::array<std::unique_ptr<const Interface>,_size> _finite_element_maps;
};

/**
 * @brief      Factory for SubDomainLocalFiniteElementMap instances
 * @ingroup    Factory, FiniteElementMap
 * @tparam     <unnamed>  Template parameters of the SubDomainLocalFiniteElementMap
 */
template<class Entity, class... LocalFiniteElementMaps>
struct Factory<VariadicLocalFiniteElementMap<Entity,LocalFiniteElementMaps...>>
{
  /**
   * @brief      Create method
   *
   * @param      ctx   @ref DataContext containing a entity mapper and sufficient data to
   *                   create finite element maps of the type LocalFiniteElementMaps... from
   *                   another factory.
   *
   * @tparam     Ctx   Universal reference to the @ref DataContext
   *
   * @return     Instance of SubDomainLocalFiniteElementMap
   */
  template<class Ctx>
  static auto create(Ctx&& ctx)
  {
    using FEM = VariadicLocalFiniteElementMap<Entity,LocalFiniteElementMaps...>;

    using Index = std::size_t;
    using EntityMapper = std::function<Index(Entity)>;
    using dCtx = std::decay_t<Ctx>;
    static_assert(dCtx::has( Context::Tag<EntityMapper>{} ), "Context does not contain a EntityMapper");

    return std::make_unique<FEM>(ctx.view( Context::Tag<EntityMapper>{} ),
                                 std::make_tuple(Factory<LocalFiniteElementMaps>::create(ctx)...));
  }
};


} // namespace Dune::Copasi

#endif // DUNE_COPASI_VARIADIC_LOCAL_FINITE_ELEMENT_HH
