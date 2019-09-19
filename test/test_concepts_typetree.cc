#ifdef HAVE_CONFIG_H
#include "config.h"
#endif

#include <dune/copasi/concepts/typetree.hh>

#include <dune/typetree/typetree.hh>

#include <dune/common/exceptions.hh>

#include <cassert>
#include <complex>

int
main(int argc, char** argv)
{
  try {
    bool passed = true;

    using namespace Dune::Copasi::Concept;

    using Node = Dune::TypeTree::LeafNode;
    passed &= isTypeTreeNode<Node>();

    using PowerNode = Dune::TypeTree::PowerNode<Node, 5>;
    passed &= isTypeTreeNode<PowerNode>();

    using CompositeNode = Dune::TypeTree::CompositeNode<Node, PowerNode>;
    passed &= isTypeTreeNode<CompositeNode>();

    return not passed;
  } catch (Dune::Exception& e) {
    std::cerr << "Dune reported error: " << e << std::endl;
  } catch (...) {
    std::cerr << "Unknown exception thrown!" << std::endl;
  }
}