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

  static_assert(not Ctx::has(Signature::geometry_type));

  GeometryTypeCtx(const Ctx& ctx, const Dune::GeometryType& geometry_type)
    : Ctx(ctx)
    , _gt(geometry_type)
  {}

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


  const Dune::GeometryType& geometry_type() const
  {
    return get(Signature::geometry_type);
  }

private:
  Dune::GeometryType _gt;
};

} // namespace Dune::Copasi::Context

#endif // DUNE_COPASI_CONTEXT_GEOMETRY_TYPE_HH
