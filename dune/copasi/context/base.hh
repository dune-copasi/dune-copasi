#ifndef DUNE_COPASI_CONTEXT_BASE_HH
#define DUNE_COPASI_CONTEXT_BASE_HH

#include <type_traits>

namespace Dune::Copasi::Context {

struct BaseCtx
{
  BaseCtx()
    : _bound(false)
  {}

  template<class Ctx>
  BaseCtx(Ctx&& ctx)
    : _bound(true)
  {}

  // does this context has the method get for this signature?
  template<class Signature>
  static bool constexpr has(Signature)
  {
    return false;
  }

  template<class Signature>
  void get(Signature) const
  {
    static_error(Dune::AlwaysFalse<Signature>::value);
  }

  template<class BindCtx>
  void bind(const BindCtx& bind_ctx)
  {
    assert(not _bound);
    _bound = true;
  }

  inline void unbind()
  {
    _bound = false;
  }

  inline bool bound() const
  {
    return _bound;
  }

protected:
  bool _bound;
};

} // namespace Dune::Copasi::Context

#endif // DUNE_COPASI_CONTEXT_BASE_HH
