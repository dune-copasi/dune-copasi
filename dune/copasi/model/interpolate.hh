#ifndef DUNE_COPASI_MODEL_INTERPOLATE_HH
#define DUNE_COPASI_MODEL_INTERPOLATE_HH

#include <dune/copasi/concepts/compartment_tree.hh>

#include <dune/pdelab/common/local_container.hh>
#include <dune/pdelab/concepts/basis.hh>

#include <dune/typetree/treecontainer.hh>

#include <spdlog/spdlog.h>

#include <optional>
#include <string>
#include <unordered_map>

namespace Dune::Copasi {

template<PDELab::Concept::Basis Basis, class GridFunction>
  requires Concept::LocalBasisTree<typename Basis::LocalView::Tree>
inline static void
interpolate(const Basis& basis,
            auto& coefficients,
            const std::unordered_map<std::string, GridFunction>& initial)
{
  spdlog::info("Set initial state from grid functions");

  if (basis.entitySet().size(0) == 0) {
    return;
  }

  auto lspace = basis.localView();

  PDELab::LocalContainerBuffer lcontainer{ basis, &coefficients };
  std::vector<double> buff;

  // create a tree container with optional local functions
  lspace.bind(*basis.entitySet().template begin<0>());
  using LocalFunction = decltype(localFunction(std::declval<GridFunction>()));
  auto lfuncs = TypeTree::makeTreeContainer(
    lspace.tree(), [&](const auto& lbasis_node) -> std::optional<LocalFunction> {
      const auto& compartment_basis = basis.subSpace(pop_back(lbasis_node.orderingViewPath()));
      const auto& component_basis = basis.subSpace(lbasis_node.orderingViewPath());
      auto it =
        initial.find(fmt::format("{}.{}", compartment_basis.name(), component_basis.name()));
      if (it != end(initial)) {
        return localFunction(it->second);
      }
      return std::nullopt;
    });
  lspace.unbind();

  // loop once over the grid and interpolate
  for (const auto& entity : elements(basis.entitySet())) {
    // bind local views to element
    lspace.bind(entity);
    lcontainer.clear(lspace);
    buff.clear();

    forEachLeafNode(lspace.tree(), [&](const auto& lbasis_node, auto path) {
      auto& lfunc = lfuncs[path];
      if (lbasis_node.size() == 0 or not lfunc) {
        return;
      }
      lfunc->bind(lbasis_node.element());
      const auto& fe = lbasis_node.finiteElement();
      buff.assign(fe.size(), 0.);
      fe.localInterpolation().interpolate(*lfunc, buff);
      for (std::size_t dof = 0; dof != lbasis_node.size(); ++dof) {
        lcontainer.accumulate(lbasis_node, dof, buff[dof]);
      }
      lfunc->unbind();
    });

    lcontainer.store(lspace);
    lspace.unbind();
  }
}

} // namespace Dune::Copasi

#endif // DUNE_COPASI_MODEL_INTERPOLATE_HH
