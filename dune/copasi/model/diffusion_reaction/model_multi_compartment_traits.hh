#ifndef DUNE_COPASI_MODEL_DIFFUSION_REACTION_MULTI_COMPARTMENT_TRAITS_HH
#define DUNE_COPASI_MODEL_DIFFUSION_REACTION_MULTI_COMPARTMENT_TRAITS_HH

// file: diffusion reaction for multi compartment models

#include <dune/copasi/concepts/grid.hh>

#include <dune/pdelab/basis/merging_strategy.hh>
#include <dune/pdelab/finiteelementmap/pkfem.hh>

#include <dune/grid/concepts/grid.hh>
#include <dune/grid/concepts/gridview.hh>

namespace Dune::Copasi {

template<class BaseTraits, bool CompartmentBlocking = false>
  requires Concept::MultiDomainGrid<typename BaseTraits::Grid> &&
           Concept::SubDomainGrid<typename BaseTraits::CompartmentEntitySet::Grid>
struct ModelMultiCompartmentDiffusionReactionPkTraits : public BaseTraits
{
  using MultiCompartmentEntitySet = typename BaseTraits::Grid::LeafGridView;
  using MultiCompartmentMergingStrategy = PDELab::Lexicographic<CompartmentBlocking>;
};

} // namespace Dune::Copasi

#endif // DUNE_COPASI_MODEL_DIFFUSION_REACTION_MULTI_COMPARTMENT_TRAITS_HH
