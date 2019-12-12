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

template<class E, class Ctx>
struct EntityMapperCtx : public Ctx
{
  using Entity = E;
  using Index = std::size_t;
  using EntityMapper = std::function<Index(const Entity&)>;

  using Ctx::has;
  using Ctx::get;
  using Ctx::_bound;

  EntityMapperCtx(Ctx&& ctx)
    : Ctx(std::move(ctx))
  {
    _bound = false;
  }

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

  const Index& index() const
  {
    return _index;
  }

  template<class BindCtx>
  void bind(const BindCtx& bind_ctx)
  {
    _index = _entity_mapper(bind_ctx.get(Signature::entity));
    Ctx::bind(bind_ctx);
  }

  void unbind()
  {
#ifndef NDEBUG
    _index = std::numeric_limits<Index>::max();
#endif
    Ctx::unbind();
  }

private:
  EntityMapper _entity_mapper;
  Index _index;
};

template<class Entity, class EntityMapper>
auto make(Signature::EntityMapper, EntityMapper&& entity_mapper)
{
  EntityCtx<Entity,BaseCtx> et_ctx{BaseCtx{}};
  return EntityMapperCtx<Entity,EntityCtx<Entity,BaseCtx>>{std::move(et_ctx),entity_mapper};
}

} // namespace Dune::Copasi::Context

#endif // DUNE_COPASI_CONTEXT_ENTITY_MAPPER_HH
