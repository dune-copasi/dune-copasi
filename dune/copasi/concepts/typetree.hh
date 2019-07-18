#ifndef DUNE_COPASI_CONCEPTS_TYPETREE_HH
#define DUNE_COPASI_CONCEPTS_TYPETREE_HH

#include <dune/functions/common/functionconcepts.hh>

/**
 * @ingroup Concepts
 */
namespace Dune::Copasi::Concept {

using namespace Dune::Concept;

//! Concept for dune type tree nodes
struct TypeTreeNode
{
  template<class T>
  auto require(T&& t) -> decltype(
    requireConvertible<bool>(T::isLeaf),
    requireConvertible<bool>(T::isPower),
    requireConvertible<bool>(T::isComposite),
    requireConvertible<std::size_t>(T::CHILDREN),
    requireType<typename T::NodeTag>(),
    requireConvertible<std::size_t>(T::degree())
  );
};

/**
 * @brief     Check if a type is dune type tree node
 * @tparam G        The type to check
 * @return true     if the type is a dune type tree node
 * @return false    if the type is not a dune type tree node
 */
template<class T>
static constexpr bool
isTypeTreeNode()
{
  return models<Concept::TypeTreeNode, T>();
}

} // namespace Dune::Copasi::Concept

#endif // DUNE_COPASI_CONCEPTS_TYPETREE_HH