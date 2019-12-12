#ifndef DUNE_COPASI_CONTEXT_GEOMETRY_TYPE_HH
#define DUNE_COPASI_CONTEXT_GEOMETRY_TYPE_HH

#include <dune/copasi/context/entity.hh>
#include <dune/copasi/context/base.hh>

#include <dune/geometry/type.hh>

#include <dune/common/typetraits.hh>

#include <type_traits>

namespace Dune::Copasi::Signature {

struct GeometryType {};
constexpr GeometryType geometry_type;

}

namespace Dune::Copasi::Context {

template<class Ctx>
struct GeometryTypeCtx : public Ctx
{
  using GeometryType = Dune::GeometryType;

  using Ctx::has;
  using Ctx::get;
  using Ctx::_bound;

  GeometryTypeCtx(Ctx&& ctx)
    : Ctx(std::move(ctx))
  {
    _bound = false;
  }

  GeometryTypeCtx(Ctx&& ctx, const Dune::GeometryType& geometry_type)
    : Ctx(std::move(ctx))
    , _gt(geometry_type)
  {}

  static bool constexpr has(Signature::GeometryType)
  {
    return true;
  }

  inline const Dune::GeometryType& get(Signature::GeometryType) const
  {
    return _gt;
  }

  template<class BindCtx>
  void inline bind(const BindCtx& bind_ctx)
  {
    if constexpr (Ctx::has(Signature::geometry_type))
      _gt = bind_ctx.get(Signature::geometry_type);
    else if constexpr (Ctx::has(Signature::entity))
      _gt = bind_ctx.get(Signature::entity).geometry().type();
    else
      static_assert(Dune::AlwaysFalse<Ctx>::value);
    Ctx::bind(bind_ctx);
  }

private:
  Dune::GeometryType _gt;
};

template<class GeometryType>
auto make(Signature::GeometryType, GeometryType&& geometry_type)
{
  return GeometryTypeCtx<BaseCtx>{BaseCtx{},std::forward<GeometryType>(geometry_type)};
}

} // namespace Dune::Copasi::Context

#endif // DUNE_COPASI_CONTEXT_GEOMETRY_TYPE_HH
