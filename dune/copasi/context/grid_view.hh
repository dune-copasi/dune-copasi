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

template<class Ctx, class GV>
struct GridViewCtx : public Ctx
{
  using GridView = GV;

  using Ctx::has;
  using Ctx::get;

  static_assert(not Ctx::has(Signature::grid_view));

  GridViewCtx(const Ctx& ctx, const GridView& grid_view)
    : Ctx(ctx)
    , _grid_view(std::make_unique<GridView>(grid_view))
  {}

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
  {}

  static bool constexpr has(Signature::GridView)
  {
    return true;
  }

  inline const GV& get(Signature::GridView) const
  {
    return *_grid_view;
  }

  inline const GV& grid_view() const
  {
    return get(Signature::grid_view);
  }

private:
  std::unique_ptr<const GridView> _grid_view;
};

} // namespace Dune::Copasi::Context

#endif // DUNE_COPASI_CONTEXT_GRID_VIEW_HH
