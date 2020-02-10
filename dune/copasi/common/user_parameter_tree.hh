#ifndef DUNE_COPASI_USER_PARAMETER_TREE_HH
#define DUNE_COPASI_USER_PARAMETER_TREE_HH

#include <dune/typetree/typetree.hh>

namespace Dune::Copasi {

enum class Methods {
  ConformingGalerkinSimplex,
  FiniteVolume
  };

namespace UserParameter {

struct Equation
  : public Dune::TypeTree::LeafNode
{
  Dune::Copasi::Methods method;
  std::size_t op;
  std::string variable;
  std::string reaction_part;
  std::string diffusion_part;
  Dune::TypeTree::HybridTreePath<std::size_t,std::size_t,std::size_t> path;
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

struct MultiDomain
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

} // namespace UserParameter
} // namespace Dune::Copasi

#endif // DUNE_COPASI_USER_PARAMETER_TREE_HH