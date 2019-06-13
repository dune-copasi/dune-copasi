#ifndef DUNE_COPASI_UTIL_HH
#define DUNE_COPASI_UTIL_HH

namespace Dune {
namespace Copasi {

/**
 * @brief      Tag a type with an ID for a given set of types.
 * @details    Only the first appearance of a type is tagged, hence, duplicated
 *             types should be avoided.
 *
 * @tparam     Type        Type to identify
 * @tparam     id          index to the first type on the set
 * @tparam     SetOfTypes  A variadic set of types where Type should be
 *                         contained
 */
template<class Type, int id, class... SetOfTypes>
struct argument_id;

#ifndef DOXYGEN

template<class Type, int id>
struct argument_id<Type,id>
{
  static_assert(true,
      "SetOfTypes... do not contain Type");
};

template<class Type, int id, class Head, class... SetOfTypes>
struct argument_id<Type,id,Head,SetOfTypes...> 
  : public std::conditional_t<
            std::is_same_v<Type,Head>,
              std::integral_constant<int,id>,
              argument_id<Type,id+1,SetOfTypes...>> 
{};
#endif // DOXYGEN

} // Copasi namespace
} // Dune namespace

#endif // DUNE_COPASI_UTIL_HH