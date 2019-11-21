#ifndef DUNE_COPASI_CONTEXT_GRID_VIEW_HH
#define DUNE_COPASI_CONTEXT_GRID_VIEW_HH

#include <dune/copasi/concepts/has_method.hh>
#include <dune/copasi/context/base.hh>

#include <type_traits>

namespace Dune::Copasi::Context {

template<class GV, class Ctx, class = std::void_t<>>
struct GridViewCtx : public Ctx
{
  using GridView = GV;

  GridViewCtx(Ctx&& ctx)
    : Ctx(std::move(ctx))
  {}

  const GV& grid_view() const
  {
    return *_grid_view;
  }

  void set_grid_view(const GV& grid_view)
  {
    _grid_view = &grid_view;
  }

  template<class BindCtx>
  std::enable_if_t<Concept::has_method_grid_view<BindCtx>()>
  bind(const BindCtx& bind_ctx)
  {
    set_grid_view(bind_ctx.grid_view());
    Ctx::bind(bind_ctx);
  }

private:
  const GV* _grid_view;
};

template<class Ctx>
struct GridViewCtx<Ctx,std::enable_if_t<Concept::has_method_grid_view<Ctx>()>> : public Ctx
{
  GridViewCtx(Ctx&& ctx)
    : Ctx(std::move(ctx))
  {}
};

template<class Ctx, class GV>
auto inject_grid_view(Ctx&& ctx, const GV& grid_view)
{
  Dune::Copasi::Context::GridViewCtx<GV,Ctx> gv_ctx{std::move(ctx)};
  gv_ctx.set_grid_view(grid_view);
  return gv_ctx;
}

template<class GV>
auto make_grid_view(const GV& grid_view)
{
  return inject_grid_view(BaseCtx{},grid_view);
}

} // namespace Dune::Copasi::Context

#endif // DUNE_COPASI_CONTEXT_GRID_VIEW_HH
