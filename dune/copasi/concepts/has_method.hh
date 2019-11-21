#ifndef DUNE_COPASI_CONCEPTS_HAS_METHOD_HH
#define DUNE_COPASI_CONCEPTS_HAS_METHOD_HH

#include <dune/functions/common/functionconcepts.hh>

/**
 * @ingroup Concepts
 */
namespace Dune::Copasi::Concept {

using namespace Dune::Concept;

//! Concept for methods (with no arguments)
#define CONCEPT_FOR_METHOD(X)                                  \
struct HasMethod_##X                                           \
{                                                              \
private:                                                       \
  template<class T, class = std::void_t<>>                     \
  struct has_method : public std::false_type {};               \
                                                               \
  template<class T>                                            \
  struct has_method<T,                                         \
      std::void_t<decltype(std::declval<T>().X())>>            \
    : public std::true_type {};                                \
                                                               \
public:                                                        \
  template<class T>                                            \
  auto require(T&& t) ->                                       \
    decltype(requireTrue<has_method<T>::value>());             \
};                                                             \
                                                               \
/**                                                            \
 * @brief     Check if a type has a method X                   \
 * @tparam G        The type to check                          \
 * @return true     if the type has method X                   \
 * @return false    if the type does not have method X         \
 */                                                            \
template<class T>                                              \
static constexpr bool                                          \
has_method_##X()                                               \
{                                                              \
  return models<Concept::HasMethod_##X, T>();                  \
}                                                              \

CONCEPT_FOR_METHOD(power_size);
CONCEPT_FOR_METHOD(geometry_type);
CONCEPT_FOR_METHOD(order);
CONCEPT_FOR_METHOD(grid_view);

} // namespace Dune::Copasi::Concept

#endif // DUNE_COPASI_CONCEPTS_HAS_METHOD_HH