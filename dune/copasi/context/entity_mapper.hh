#ifndef DUNE_COPASI_CONTEXT_ENTITY_MAPPER_HH
#define DUNE_COPASI_CONTEXT_ENTITY_MAPPER_HH

#include <dune/copasi/context/base.hh>
#include <dune/copasi/context/entity.hh>

#include <type_traits>
#include <memory>

namespace Dune::Copasi::Signature {

struct EntityMapper {};
constexpr EntityMapper entity_mapper;

}

namespace Dune::Copasi::Context {

template<class Ctx, class E>
struct EntityMapperCtx : public Ctx
{
  using Entity = E;
  using Index = std::size_t;
  using EntityMapper = std::function<Index(const Entity&)>;

  using Ctx::has;
  using Ctx::get;

  static_assert(not Ctx::has(Signature::entity_mapper));

  template<class EM>
  EntityMapperCtx(const Ctx& ctx, EM&& entity_mapper)
    : Ctx(ctx)
    , _entity_mapper(std::forward<EntityMapper>(entity_mapper))
  {}

  template<class EM>
  EntityMapperCtx(Ctx&& ctx, EM&& entity_mapper)
    : Ctx(std::move(ctx))
    , _entity_mapper(std::forward<EntityMapper>(entity_mapper))
  {}

  static bool constexpr has(Signature::EntityMapper)
  {
    return true;
  }

  inline const EntityMapper& get(Signature::EntityMapper) const
  {
    return _entity_mapper;
  }

  inline const EntityMapper& entity_mapper() const
  {
    return get(Signature::entity_mapper);
  }

private:
  EntityMapper _entity_mapper;
};

} // namespace Dune::Copasi::Context

#endif // DUNE_COPASI_CONTEXT_ENTITY_MAPPER_HH
