#ifndef DUNE_COPASI_CONTEXT_BASE_HH
#define DUNE_COPASI_CONTEXT_BASE_HH

#include <type_traits>

namespace Dune::Copasi::Context {

struct BaseCtx
{
  BaseCtx() {}

  template<class Ctx>
  BaseCtx(Ctx&& ctx)
  {}


  template<class BindCtx>
  void bind(const BindCtx& bind_ctx)
  {}
};

} // namespace Dune::Copasi::Context

#endif // DUNE_COPASI_CONTEXT_BASE_HH
