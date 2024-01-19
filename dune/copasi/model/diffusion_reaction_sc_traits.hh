#ifndef DUNE_COPASI_MODEL_DIFFUSION_REACTION_SINGLE_COMAPARTMENT_TRAITS_HH
#define DUNE_COPASI_MODEL_DIFFUSION_REACTION_SINGLE_COMAPARTMENT_TRAITS_HH

#include <dune/pdelab/basis/merging_strategy.hh>
#include <dune/pdelab/finiteelementmap/pkfem.hh>

#include <dune/grid/concepts/grid.hh>
#include <dune/grid/concepts/gridview.hh>

namespace Dune::Copasi {

template<Dune::Concept::Grid G,
         Dune::Concept::GridView ES,
         std::size_t Order = 1,
         class RQ = double,
         class TQ = double,
         bool ScalarBlocking = false>
struct ModelDiffusionPkReactionTraits
{
  using Grid = G;
  using CompartmentEntitySet = ES;
  static_assert(Order != 0);
  using ScalarFiniteElementMap =
    PDELab::PkLocalFiniteElementMap<CompartmentEntitySet, double, double, Order>;

  using ScalarMergingStrategy = PDELab::EntityGrouping<CompartmentEntitySet, false>;
  using CompartmentMergingStrategy = PDELab::EntityGrouping<CompartmentEntitySet, ScalarBlocking>;
  using RangeQuatinty = RQ;
  using TimeQuantity = TQ;
};

} // namespace Dune::Copasi

#endif // DUNE_COPASI_MODEL_DIFFUSION_REACTION_SINGLE_COMAPARTMENT_TRAITS_HH
