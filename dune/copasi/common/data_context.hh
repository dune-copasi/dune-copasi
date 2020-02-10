#ifndef DUNE_COPASI_CONTEXT_BASE_HH
#define DUNE_COPASI_CONTEXT_BASE_HH

#include <type_traits>

namespace Dune::Copasi::Context {

template<class T>
struct Tag {};

struct DataContextBase
{
  DataContextBase()
  {}

  template<class Ctx>
  DataContextBase(Ctx&& ctx)
  {}

  // does this context has the method get for this signature?
  template<class T>
  static bool constexpr has(Tag<T>)
  {
    return false;
  }

  template<class T>
  const T& view(Tag<T>) const
  {
    static_error(Dune::AlwaysFalse<T>::value, "Data context does not contain type T");
  }

  template<class T>
  T&& get(T)
  {
    static_error(Dune::AlwaysFalse<T>::value, "Data context does not contain type T");
  }

  friend std::ostream& operator<<(std::ostream& os, const DataContextBase& ctx)
  {
    return os;
  }
};

template<class T, class Ctx = DataContextBase>
struct DataContext : public Ctx
{
  static_assert(std::is_copy_constructible_v<T> or std::is_move_constructible_v<T>);
  static_assert(not Ctx::has( Tag<T>{} ), "Data context already contains type T");

  template<class _T>
  static constexpr bool is_valid_t = std::is_same_v<std::decay_t<T>,std::decay_t<_T>>;
  template<class _Ctx>
  static constexpr bool is_valid_ctx = std::is_same_v<std::decay_t<Ctx>,std::decay_t<_Ctx>>;

  using Ctx::get;
  using Ctx::has;
  using Ctx::view;

  template<class _T, class _Ctx,
    class = std::enable_if_t<is_valid_t<_T>&&is_valid_ctx<_Ctx>>>
  DataContext(_T&& value, _Ctx&& ctx)
    : Ctx(std::forward<_Ctx>(ctx))
    , _value(std::forward<_T>(value))
  {}

  template<class _T,
    bool is_valid_t = std::is_same_v<std::decay_t<T>,std::decay_t<_T>>,
    bool is_ctx_default_ctble = std::is_default_constructible_v<Ctx>,
    class = std::enable_if_t<is_valid_t&&is_ctx_default_ctble>>
  DataContext(_T&& value)
    : _value(std::forward<_T>(value))
  {}

  // does this context has the method get for this signature?
  static bool constexpr has(Tag<T>)
  {
    return true;
  }

  const T& view(Tag<T>) const
  {
    return _value;
  }

  T&& get(Tag<T>)
  {
    return std::move(_value);
  }

  // Print all contained types in os
  friend std::ostream& operator<<(std::ostream& os, const DataContext<T,Ctx>& ctx)
  {
    std::string ending = std::is_same_v<Ctx,DataContextBase> ? ";" : ", ";
    os << "\t" << Dune::className<T>() << ending;
    os << *static_cast<const Ctx *>(&ctx);
    return os;
  }
private:
  T _value;
};

} // namespace Dune::Copasi::Context

#endif // DUNE_COPASI_CONTEXT_BASE_HH
