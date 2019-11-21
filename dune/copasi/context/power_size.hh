#ifndef DUNE_COPASI_CONTEXT_POWER_SIZE_HH
#define DUNE_COPASI_CONTEXT_POWER_SIZE_HH

#include <dune/copasi/concepts/has_method.hh>
#include <dune/copasi/context/base.hh>

#include <type_traits>

namespace Dune::Copasi::Context {

template<class Ctx, class = std::void_t<>>
struct PowerSizeCtx : public Ctx
{
  using PowerSize = std::size_t;

  PowerSizeCtx(Ctx&& ctx)
    : Ctx(std::move(ctx))
  {}

  std::size_t power_size() const
  {
    return _power_size;
  }

  void set_power_size(std::size_t power_size)
  {
    _power_size = power_size;
  }

  template<class BindCtx>
  std::enable_if_t<Concept::has_method_power_size<BindCtx>()>
  bind(const BindCtx& bind_ctx)
  {
    set_power_size(bind_ctx.power_size());
    Ctx::bind(bind_ctx);
  }

private:
  std::size_t _power_size;
};

template<class Ctx>
struct PowerSizeCtx<Ctx,std::enable_if_t<Concept::has_method_geometry_type<Ctx>()>> : public Ctx
{
  PowerSizeCtx(Ctx&& ctx)
    : Ctx(std::move(ctx))
  {}
};

template<class Ctx>
auto inject_power_size(Ctx&& ctx, std::size_t power_size)
{
  PowerSizeCtx<Ctx> ps_ctx{std::move(ctx)};
  ps_ctx.set_power_size(power_size);
  return ps_ctx;
}

auto make_power_size(std::size_t power_size)
{
  return inject_power_size(BaseCtx{},power_size);
}

} // namespace Dune::Copasi::Context

#endif // DUNE_COPASI_CONTEXT_POWER_SIZE_HH
