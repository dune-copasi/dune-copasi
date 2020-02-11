#ifndef DUNE_COPASI_CONTEXT_BASE_HH
#define DUNE_COPASI_CONTEXT_BASE_HH

#include <type_traits>

/**
@defgroup DataContext Data Context
@{
  @brief Generalized and extensible static polymorphic data holder.

  A data context is a generic form to pass around arbitrary unique amount of data following
  a defined interface. Moreover, a data context is always extensible to hold more data
  at compile time. It has a very similar concept and use cases as `**kwargs` in python.
  It is called context because the content of the object depends on the context where it was called.

  A data context knows at comile-time whether it contains a type identified
  with certain signature.

  @code{.cpp}
  template<class Ctx>
  void foo(Ctx&& ctx)
  {
    using GV = ...;
    // check whether context contains the type GV
    constexpr bool result = Ctx::has( Context::Tag<GV>{} );

    // Get a view on the value GV contained in the context
    const auto& view = ctx.view( Context::Tag<GV>{} );

    // Create a new context containing a value 10 for the type `int`
    // Notice that the old context can be moved into the new one if possible
    auto new_ctx = Context::DataContext<Ctx,int>{10,std::forward<Ctx>(ctx))};

    // move the value contained in `int` outside of the context
    auto value = new_ctx.get( Context::Tag<int>{} );

    // If not sure what types a data context contain, one can print them with
    // Notice that it does not give information whether the data is valid or not
    std::cout << "Output: " << new_ctx << std::endl
    // Output: <data contained in ctx>, int;
  }
  @endcode

  It's important to notice that this is also possible with variadic templates or
  tuples, however, they usually involve a lot of template metaprogramming to identify
  types and to add more types to it.

  An easy way to create a data context with different values is using the method Context::data_context

  @code{.cpp}
  auto ctx = Context::data_context(int(1), double(2.0), std::string("I'm stored in a context"));
  @endcode
@}
*/

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
 * @tparam     T     Type to store in the data context
 * @tparam     Ctx   Context to be extended
 */
template<class T, class Ctx = DataContextBase>
class DataContext : public Ctx
{
  // ensure that we have means to store T
  static_assert(std::is_copy_constructible_v<T> or std::is_move_constructible_v<T>);

  // ensure that each type is associated with a value only once
  static_assert(not Ctx::has( Tag<T>{} ), "Data context already contains type T");

  // helper function to assert a universal reference of T
  template<class _T>
  static constexpr bool is_valid_t = std::is_same_v<std::decay_t<T>,std::decay_t<_T>>;

  // helper function to assert a universal reference of Ctx
  template<class _Ctx>
  static constexpr bool is_valid_ctx = std::is_same_v<std::decay_t<Ctx>,std::decay_t<_Ctx>>;

public:
  // Inherit methods from the base context
  using Ctx::get;
  using Ctx::has;
  using Ctx::view;

  /**
   * @brief      Constructs a new instance
   *
   * @param      value      The value
   * @param      ctx        The base context
   *
   * @tparam     _T         Universal reference of T
   * @tparam     _Ctx       Universal reference of Ctx
   * @tparam     <unnamed>  Helper SFINAE to disable constructor when _T and _Ctx are not valid
   */
  template<class _T, class _Ctx,
  class = std::enable_if_t<is_valid_t<_T> and is_valid_ctx<_Ctx>>>
  DataContext(_T&& value, _Ctx&& ctx)
    : Ctx(std::forward<_Ctx>(ctx))
    , _value(std::forward<_T>(value))
  {}

  /**
   * @brief      Constructs a new instance.
   *
   * @param      value                 The value
   *
   * @tparam     _T                    Universal reference of T
   * @tparam     is_ctx_default_ctble  Helper SFINAE bool for default constructible base context
   * @tparam     <unnamed>             Helper SFINAE to disable constructor when _T is not valid and Ctx is not default construcible
   */
  template<class _T,
    bool is_ctx_default_ctble = std::is_default_constructible_v<Ctx>,
    class = std::enable_if_t<is_valid_t<_T> and is_ctx_default_ctble>>
  DataContext(_T&& value)
    : _value(std::forward<_T>(value))
  {}

  /**
   * @brief      Method to check if context contains a type
   *
   * @tparam     T          The type to check
   */
  static bool constexpr has(Tag<T>)
  {
    return true;
  }

  /**
   * @brief      Gets a view on the stored value for type T
   *
   * @tparam     T          The type to view
   *
   * @return     A const reference to the stored value
   */
  const T& view(Tag<T>) const
  {
    return _value;
  }

  /**
   * @brief      Gets the ownership of the stored value for type T
   *
   * @tparam     T          The type to get
   *
   * @return     A rvalue to the stored value
   */
  T&& get(Tag<T>)
  {
    return std::move(_value);
  }

  /**
   * @brief      Writes type names to the output stream
   *
   * @param      os    The output stream
   * @param[in]  ctx   The context to print
   *
   * @return     The output stream with the type names
   */
  friend std::ostream& operator<<(std::ostream& os, const DataContext<T,Ctx>& ctx)
  {
    std::string ending = std::is_same_v<Ctx,DataContextBase> ? ";" : ", ";
    os << "\t" << Dune::className<T>() << ending;
    os << *static_cast<const Ctx *>(&ctx);
    return os;
  }
private:
  //! Actual storage data of type T
  T _value;
};

/**
 * @brief      Create a data context for one value
 *
 * @param      value  The value to store in the data context
 *
 * @tparam     T      Universal reference of the value to store
 *
 * @return     A context storing the value 
 */
template<class T>
auto data_context(T&& value)
{
  return Context::DataContext<std::decay_t<T>>(std::forward<T>(value));
}

/**
 * @brief      Create a data context for several values
 *
 * @param      arg   The first argument in the data context
 * @param      args  The rest of arguments in the data context
 *
 * @tparam     Arg0  Universal reference of the first value to store
 * @tparam     Args  Universal references of the rest of values to store
 *
 * @return     A context storing all the values
 */
template<class Arg0, class... Args>
auto data_context(Arg0&& arg, Args&&... args)
{
  auto base_ctx = data_context(std::forward<Args>(args)...);
  return Context::DataContext<std::decay_t<Arg0>, decltype(base_ctx)>(arg, std::move(base_ctx));
}

} // namespace Dune::Copasi::Context

#endif // DUNE_COPASI_CONTEXT_BASE_HH
