#ifndef DUNE_COPASI_CONTEXT_GRID_VIEW_HH
#define DUNE_COPASI_CONTEXT_GRID_VIEW_HH

#include <dune/copasi/context/base.hh>

#include <type_traits>
#include <memory>

namespace Dune::Copasi::Signature {

struct GridView {};
constexpr GridView grid_view;

}

namespace Dune::Copasi::Context {

template<class GV, class Ctx>
struct GridViewCtx : public Ctx
{
  using GridView = GV;

  using Ctx::has;
  using Ctx::get;
  using Ctx::_bound;

  GridViewCtx(Ctx&& ctx)
    : Ctx(std::move(ctx))
  {
    _bound = false;
  }

  GridViewCtx(const GridViewCtx& other_ctx)
    : Ctx(other_ctx)
  {
    _grid_view = std::make_unique<GridView>(*other_ctx._grid_view);
  }

  GridViewCtx(Ctx&& ctx, const GridView& grid_view)
    : Ctx(std::move(ctx))
    , _grid_view(std::make_unique<GridView>(grid_view))
  {}

  ~GridViewCtx()
  {
    unbind();
  }

  static bool constexpr has(Signature::GridView)
  {
    return true;
  }

  inline const GV& get(Signature::GridView) const
  {
    return *_grid_view;
  }

  inline void unbind()
  {
    _grid_view.reset();
    Ctx::unbind();
  }

  template<class BindCtx>
  inline void bind(const BindCtx& bind_ctx)
  {
    _grid_view = std::make_unique<GridView>(bind_ctx.get(Signature::grid_view));
    Ctx::bind(bind_ctx);
  }

private:
  std::unique_ptr<const GridView> _grid_view;
};

template<class GridView>
auto make(Signature::GridView, const GridView& grid_view)
{
  return GridViewCtx<GridView,BaseCtx>{BaseCtx{},grid_view};
}

} // namespace Dune::Copasi::Context

#endif // DUNE_COPASI_CONTEXT_GRID_VIEW_HH
