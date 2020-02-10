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

template<class Ctx, class E>
struct EntityCtx : public Ctx
{
  using Entity = E;

  using Ctx::has;
  using Ctx::get;

  static_assert(not Ctx::has(Signature::entity));

  EntityCtx(const Ctx& ctx, const Entity& entity)
    : Ctx(ctx)
    , _entity(entity)
  {}

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

private:
  Entity _entity;
};

} // namespace Dune::Copasi::Context

#endif // DUNE_COPASI_CONTEXT_ENTITY_HH
