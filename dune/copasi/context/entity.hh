#ifndef DUNE_COPASI_CONTEXT_ENTITY_HH
#define DUNE_COPASI_CONTEXT_ENTITY_HH

#include <dune/copasi/context/base.hh>

#include <type_traits>
#include <memory>

namespace Dune::Copasi::Signature {

struct Entity {};
constexpr Entity entity;

}

namespace Dune::Copasi::Context {

template<class E, class Ctx>
struct EntityCtx : public Ctx
{
  using Entity = E;

  using Ctx::has;
  using Ctx::get;
  using Ctx::_bound;

  EntityCtx(Ctx&& ctx)
    : Ctx(std::move(ctx))
  {
    _bound = false;
  }

  EntityCtx(Ctx&& ctx, const Entity& entity)
    : Ctx(std::move(ctx))
    , _entity(entity)
  {}

  const Entity& get(Signature::Entity) const
  {
    return _entity;
  }

  const Entity& entity() const
  {
    return get(Signature::entity);
  }

  template<class BindCtx>
  void bind(const BindCtx& bind_ctx)
  {
    _entity = bind_ctx.get(Signature::entity);
    Ctx::bind(bind_ctx);
  }

  void unbind()
  {
#ifndef NDEBUG
    _entity = Entity{};
#endif
    Ctx::unbind();
  }

private:
  Entity _entity;
};

template<class Entity>
auto make(Signature::Entity, Entity&& entity)
{
  return EntityCtx<Entity,BaseCtx>{BaseCtx{},std::forward<Entity>(entity)};
}

} // namespace Dune::Copasi::Context

#endif // DUNE_COPASI_CONTEXT_ENTITY_HH
