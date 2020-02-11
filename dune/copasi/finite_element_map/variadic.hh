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

template<class GV, class... LocalFiniteElementMaps>
class VariadicLocalFiniteElementMap
{
  using Interface = std::common_type_t<typename VirtualLocalFiniteElementMapWrapper<LocalFiniteElementMaps, GV>::Interface...>;
  static constexpr std::size_t _size = sizeof...(LocalFiniteElementMaps);
  using _integral_size = std::integral_constant<std::size_t,_size>;

  template<class... FEMs>
  static std::array<std::unique_ptr<const Interface>,_size>
  make_unique_array_from_tuple(std::tuple<std::unique_ptr<const FEMs>...>&& finite_element_maps)
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
      fem_unique_array[i] = std::make_unique<VirtualLocalFiniteElementMapWrapper<FEM, GV>>(std::move(fem));
    });
    return fem_unique_array;
  }

  // template<class... FEMs>
  // static std::array<std::unique_ptr<const Interface>,_size>
  // make_unique_array_from_variadic(FEMs&&... finite_element_maps)
  // {
  //   return make_unique_array_from_tuple(std::forward_as_tuple(std::make_unique<const std::decay_t<FEMs>>(std::forward<FEMs>(finite_element_maps))...));
  // }

public:
  using Traits = typename Interface::Traits;

  // if this fails, finite element maps do not have common dimension!
  static constexpr int dimension = std::common_type_t<std::integral_constant<int,LocalFiniteElementMaps::dimension>...>::value;

  template<class EntityMapper>
  explicit VariadicLocalFiniteElementMap(const EntityMapper& entity_mapper, std::array<std::unique_ptr<const Interface>,_size>&& finite_element_maps)
    : _entity_mapper(entity_mapper)
    , _finite_element_maps(std::move(finite_element_maps))
  {}

  template<class EntityMapper>
  explicit VariadicLocalFiniteElementMap(const EntityMapper& entity_mapper, std::tuple<std::unique_ptr<const LocalFiniteElementMaps>...>&& finite_element_maps)
    : VariadicLocalFiniteElementMap(entity_mapper, make_unique_array_from_tuple(std::move(finite_element_maps)))
  {}

  // template<class EntityMapper, class... FEMs>
  // explicit VariadicLocalFiniteElementMap(const EntityMapper& entity_mapper, FEMs&&... finite_element_maps)
  //   : VariadicLocalFiniteElementMap(entity_mapper, make_unique_array_from_variadic(std::forward<FEMs>(finite_element_maps)...))
  // {}

  ~VariadicLocalFiniteElementMap()
  {}

  const typename Traits::FiniteElementType&
  find (const typename Traits::EntityType& e) const
  {
    assert(_entity_mapper(e) < _size);
    assert(_entity_mapper(e) >= 0);
    return _finite_element_maps[_entity_mapper(e)]->find(e);
  }

  bool fixedSize() const
  {
    return _size == 1 ? _finite_element_maps[0]->fixedSize() : false;
  }

  std::size_t size(GeometryType gt) const
  {
    if constexpr (_size == 1)
      return _finite_element_maps[0]->size(gt);
    else
      DUNE_THROW(PDELab::FiniteElementMapError,
        "this method should not be called for not fixed sizes");
  }

  bool hasDOFs(int codim) const
  {
    bool has_dofs = false;
    for (auto&& fem : _finite_element_maps)
      has_dofs |= fem->hasDOFs(codim);
    return has_dofs;
  }

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

template<class GV, class... LocalFiniteElementMaps>
struct Factory<VariadicLocalFiniteElementMap<GV,LocalFiniteElementMaps...>>
{
  template<class Ctx>
  static auto create(Ctx&& ctx)
  {
    using FEM = VariadicLocalFiniteElementMap<GV,LocalFiniteElementMaps...>;

    using Entity = typename GV::template Codim<0>::Entity;
    using Index = int; // TODO
    using EntityMapper = std::function<Index(Entity)>;
    using dCtx = std::decay_t<Ctx>;
    static_assert(dCtx::has( Context::Tag<EntityMapper>{} ), "Invalid provided context");

    return std::make_unique<FEM>(ctx.view( Context::Tag<EntityMapper>{} ),
                                 std::forward_as_tuple(Factory<LocalFiniteElementMaps>::create(ctx)...));
  }
};


} // namespace Dune::Copasi

#endif // DUNE_COPASI_VARIADIC_LOCAL_FINITE_ELEMENT_HH