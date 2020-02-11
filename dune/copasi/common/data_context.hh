#ifndef DUNE_COPASI_CONTEXT_BASE_HH
#define DUNE_COPASI_CONTEXT_BASE_HH

#include <type_traits>

namespace Dune::Copasi::Context {

/**
 * @brief      Tags any type with a unique default constructible struct
 * @ingroup    DataContext
 */
template<class T>
struct Tag {};

/**
 * @brief      Base data context context
 * @ingroup    DataContext
 */
struct DataContextBase
{
  //! Default constructor
  DataContextBase()
  {}

  /**
   * @brief      Method to check if context contains a type
   *
   * @param[in]  <unnamed>  Tagged version of the type to check
   *
   * @tparam     T          The type to check
   */
  template<class T>
  static bool constexpr has(Tag<T>)
  {
    return false;
  }

  /**
   * @brief      Gets a view on the stored value for type T
   *
   * @param[in]  <unnamed>  Tagged version of the type to view
   *
   * @tparam     T          The type to view
   *
   * @return     A const reference to the stored value
   */
  template<class T>
  const T& view(Tag<T>) const
  {
    static_error(Dune::AlwaysFalse<T>::value, "Data context does not contain type T");
  }

  /**
   * @brief      Gets the ownership of the stored value for type T
   *
   * @param[in]  <unnamed>  Tagged version of the type to get
   *
   * @tparam     T          The type to get
   *
   * @return     A rvalue to the stored value
   */
  template<class T>
  T&& get(Tag<T>)
  {
    static_error(Dune::AlwaysFalse<T>::value, "Data context does not contain type T");
  }

  /**
   * @brief      Writes type names to the output stream
   *
   * @param      os    The output stream
   * @param[in]  ctx   The context to print
   *
   * @return     The output stream with the type names
   */
  friend std::ostream& operator<<(std::ostream& os, const DataContextBase& ctx)
  {
    return os;
  }
};

/**
 * @brief      This class describes a data context.
 * @ingroup    DataContext
 *
 * @tparam     T     { description }
 * @tparam     Ctx   { description }
 */
template<class T, class Ctx = DataContextBase>
class DataContext : public Ctx
{
  static_assert(std::is_copy_constructible_v<T> or std::is_move_constructible_v<T>);
  static_assert(not Ctx::has( Tag<T>{} ), "Data context already contains type T");

  template<class _T>
  static constexpr bool is_valid_t = std::is_same_v<std::decay_t<T>,std::decay_t<_T>>;
  template<class _Ctx>
  static constexpr bool is_valid_ctx = std::is_same_v<std::decay_t<Ctx>,std::decay_t<_Ctx>>;

public:
  using Ctx::get;
  using Ctx::has;
  using Ctx::view;

  template<class _T, class _Ctx,
  class = std::enable_if_t<is_valid_t<_T> and is_valid_ctx<_Ctx>>>
  DataContext(_T&& value, _Ctx&& ctx)
    : Ctx(std::forward<_Ctx>(ctx))
    , _value(std::forward<_T>(value))
  {}

  template<class _T,
    bool is_ctx_default_ctble = std::is_default_constructible_v<Ctx>,
    class = std::enable_if_t<is_valid_t<_T> and is_ctx_default_ctble>>
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

template<class T>
auto data_context(T&& value)
{
  return Context::DataContext<std::decay_t<T>>(std::forward<T>(value));
}

template<class Arg0, class... Args>
auto data_context(Arg0&& arg, Args&&... args)
{
  auto base_ctx = data_context(std::forward<Args>(args)...);
  return Context::DataContext<std::decay_t<Arg0>, decltype(base_ctx)>(arg, std::move(base_ctx));
}

} // namespace Dune::Copasi::Context

#endif // DUNE_COPASI_CONTEXT_BASE_HH
