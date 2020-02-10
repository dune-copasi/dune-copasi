#ifndef DUNE_COPASI_MODEL_PARAMETER_TREE_HH
#define DUNE_COPASI_MODEL_PARAMETER_TREE_HH

#include <dune/typetree/typetree.hh>

namespace Dune::Copasi {
namespace ModelParameterTree {

struct Equation
  : public Dune::TypeTree::LeafNode
{
  std::string name;
  std::string reaction;
  std::string diffusion;

  Dune::TypeTree::HybridTreePath<std::size_t,std::size_t,std::size_t,std::size_t> path;
  Dune::TypeTree::HybridTreePath<std::size_t,std::size_t,std::size_t> user_path;
};

struct SubDomain
  : public Dune::TypeTree::DynamicPowerNode<Equation>
{
  inline Equation& equation(std::size_t i)
  {
    return this->child(i);
  }

  inline const Equation& equation(std::size_t i) const
  {
    return this->child(i);
  }
};

struct Method
  : public Dune::TypeTree::DynamicPowerNode<SubDomain>
{
  inline SubDomain& sub_domain(std::size_t i)
  {
    return this->child(i);
  }

  inline const SubDomain& sub_domain(std::size_t i) const
  {
    return this->child(i);
  }
};

struct Operator
  : public Dune::TypeTree::DynamicPowerNode<Method>
{
  inline Method& method(std::size_t i)
  {
    return this->child(i);
  }

  inline const Method& method(std::size_t i) const
  {
    return this->child(i);
  }
};

} // namespace ModelParameterTree
} // namespace Dune::Copasi

#endif // DUNE_COPASI_MODEL_PARAMETER_TREE_HH