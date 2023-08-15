#ifndef DUNE_COPASI_CONCEPTS_COMPARTMENT_TREE_HH
#define DUNE_COPASI_CONCEPTS_COMPARTMENT_TREE_HH

#include <dune-copasi-config.hh>

#include <dune/pdelab/concepts/local_basis.hh>
#include <dune/pdelab/concepts/tree.hh>

/**
 * @ingroup Concepts
 */
namespace Dune::Copasi::Concept {

template<class Tree, std::size_t dim>
concept ScalarLocalBasisNode =
  PDELab::Concept::LocalBasisTree<Tree> &&
  PDELab::Concept::LeafTreeNode<Tree> && (Tree::dimDomain == dim);

template<class Tree>
concept CompartmentScalarLocalBasisNode =
  PDELab::Concept::LocalBasisTree<Tree> &&
  ScalarLocalBasisNode<Tree, Tree::Element::dimension>;

template<class Tree>
concept MembraneScalarLocalBasisNode =
  PDELab::Concept::LocalBasisTree<Tree> &&
  ScalarLocalBasisNode<Tree, Tree::Element::dimension - 1>;

template<class Tree>
concept CompartmentLocalBasisNode =
  PDELab::Concept::LocalBasisTree<Tree> &&
  PDELab::Concept::VectorTreeNode<Tree> &&
  CompartmentScalarLocalBasisNode<typename Tree::ChildType>;

template<class Tree>
concept MembraneLocalBasisNode =
  PDELab::Concept::LocalBasisTree<Tree> &&
  PDELab::Concept::VectorTreeNode<Tree> &&
  MembraneScalarLocalBasisNode<typename Tree::ChildType>;

template<class Tree>
concept MembraneSubEntitiesLocalBasisNode =
  PDELab::Concept::LocalBasisTree<Tree> &&
  PDELab::Concept::VectorTreeNode<Tree> &&
  MembraneLocalBasisNode<typename Tree::ChildType>;

template<class Tree>
concept CompartmentMembraneLocalBasisNode =
  PDELab::Concept::LocalBasisTree<Tree> &&
  PDELab::Concept::TupleTreeNode<Tree> && (Tree::degree() == 2) &&
  requires(Tree tree) {
    {
      tree.child(Indices::_0)
    } -> CompartmentLocalBasisNode;
    {
      tree.child(Indices::_1)
    } -> MembraneSubEntitiesLocalBasisNode;
  };

template<class Tree>
concept MultiCompartmentLocalBasisNode =
  PDELab::Concept::LocalBasisTree<Tree> &&
  PDELab::Concept::VectorTreeNode<Tree> &&
  CompartmentLocalBasisNode<typename Tree::ChildType>;

template<class Tree>
concept MultiMembraneLocalBasisNode =
  PDELab::Concept::LocalBasisTree<Tree> &&
  PDELab::Concept::VectorTreeNode<Tree> &&
  MembraneSubEntitiesLocalBasisNode<typename Tree::ChildType>;

template<class Tree>
concept MultiCompartmentMembraneLocalBasisNode =
  PDELab::Concept::LocalBasisTree<Tree> &&
  PDELab::Concept::TupleTreeNode<Tree> && (Tree::degree() == 2) &&
  requires(Tree tree) {
    {
      tree.child(Indices::_0)
    } -> MultiCompartmentLocalBasisNode;
    {
      tree.child(Indices::_1)
    } -> MembraneSubEntitiesLocalBasisNode;
  };

template<class Tree>
concept CompartmentLocalBasisTree =
  Concept::CompartmentScalarLocalBasisNode<Tree> ||
  Concept::CompartmentLocalBasisNode<Tree> ||
  Concept::MultiCompartmentLocalBasisNode<Tree>;

template<class Tree>
concept MembraneLocalBasisTree =
  Concept::MembraneScalarLocalBasisNode<Tree> ||
  Concept::MembraneSubEntitiesLocalBasisNode<Tree> ||
  Concept::MembraneLocalBasisNode<Tree> ||
  Concept::MultiMembraneLocalBasisNode<Tree>;

template<class Tree>
concept LocalBasisTree =
  Concept::MultiCompartmentMembraneLocalBasisNode<Tree> ||
  Concept::CompartmentLocalBasisTree<Tree>;

} // namespace Dune::Copasi::Concept

#endif // DUNE_COPASI_CONCEPTS_COMPARTMENT_TREE_HH
