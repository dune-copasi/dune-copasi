#ifndef DUNE_COPASI_MODEL_DIFFUSION_REACTION_MULTI_COMPARTMENT_HH
#define DUNE_COPASI_MODEL_DIFFUSION_REACTION_MULTI_COMPARTMENT_HH

// file: diffusion reaction for multi compartment models

#include <dune/copasi/concepts/grid.hh>
#include <dune/copasi/model/local_equations/functor_factory.hh>
#include <dune/copasi/model/model.hh>
#include <dune/copasi/model/constraints.hh>

#include <dune/pdelab/basis/merging_strategy.hh>
#include <dune/pdelab/basis/prebasis/composite.hh>
#include <dune/pdelab/basis/prebasis/leaf.hh>
#include <dune/pdelab/finiteelementmap/pkfem.hh>

#include <dune/grid/concepts/grid.hh>
#include <dune/grid/concepts/gridview.hh>

namespace Dune::Copasi {

template<class Traits>
class ModelMultiCompartmentDiffusionReaction
  : public Model<typename Traits::Grid,
                 typename Traits::CompartmentEntitySet,
                 typename Traits::RangeQuatinty,
                 typename Traits::TimeQuantity>
{
  using Base = Model<typename Traits::Grid,
                     typename Traits::CompartmentEntitySet,
                     typename Traits::RangeQuatinty,
                     typename Traits::TimeQuantity>;

public:
  using State = typename Base::State;
  using Grid = typename Traits::Grid;
  using TimeQuantity = typename Traits::TimeQuantity;
  using ScalarQuantity = typename Traits::RangeQuatinty;
  using CompartmentEntitySet = typename Traits::CompartmentEntitySet;
  using MultiCompartmentEntitySet = typename Traits::MultiCompartmentEntitySet;

  using ScalarFiniteElementMap = typename Traits::ScalarFiniteElementMap;

  using ScalarMergingStrategy = typename Traits::ScalarMergingStrategy;
  using CompartmentMergingStrategy = typename Traits::CompartmentMergingStrategy;
  using MultiCompartmentMergingStrategy = typename Traits::MultiCompartmentMergingStrategy;

  using ScalarPreBasis =
    PDELab::PreBasis<ScalarMergingStrategy, ScalarFiniteElementMap, Constraints<CompartmentEntitySet>>;
  using CompartmentPreBasis = PDELab::PreBasisVector<CompartmentMergingStrategy, ScalarPreBasis>;
  using MultiCompartmentPreBasis =
    PDELab::PreBasisVector<MultiCompartmentMergingStrategy, CompartmentPreBasis>;

  using ResidualQuantity = ScalarQuantity;

  using GridFunction = typename Base::GridFunction;

  explicit ModelMultiCompartmentDiffusionReaction(
    std::shared_ptr<const FunctorFactory<Grid::dimensionworld>> functor_factory,
    std::shared_ptr<const CellData<typename Grid::LeafGridView, ScalarQuantity>> cell_data = nullptr)
    : _functor_factory{ std::move(functor_factory) }
    , _cell_data{ std::move(cell_data) }
  {
    assert(_functor_factory);
  }

  std::unique_ptr<State> make_state(const std::shared_ptr<const Grid>&,
                                    const ParameterTree&) const override;

  void interpolate(State&, const std::unordered_map<std::string, GridFunction>&) const override;

  std::unordered_map<std::string, GridFunction> make_initial(const Grid&,
                                                             const ParameterTree&) const override;

  GridFunction make_compartment_function(const std::shared_ptr<const State>&,
                                         std::string_view) const override;

  std::unique_ptr<PDELab::OneStep<State>> make_step_operator(const State&,
                                                             const ParameterTree&) const override;

  void write_vtk(const State&, const std::filesystem::path&, bool) const override;

  std::map<std::string, double> reduce(const State&, const ParameterTree&) const override;

private:
  static MultiCompartmentPreBasis make_multi_compartment_pre_basis(const Grid&,
                                                                   const ParameterTree&,
                                                                   std::shared_ptr<const FunctorFactory<Grid::dimensionworld>>);

  static void setup_coefficient_vector(State&);
  static CompartmentEntitySet get_entity_set(const Grid&, std::size_t);

  mutable std::unordered_map<std::string, std::vector<double>> _writer_timesteps;
  std::shared_ptr<const FunctorFactory<Grid::dimensionworld>> _functor_factory;
  std::shared_ptr<const CellData<typename Grid::LeafGridView, ScalarQuantity>> _cell_data;
};

} // namespace Dune::Copasi

#ifndef DUNE_COPASI_PRECOMPILED_MODE
#include <dune/copasi/model/diffusion_reaction_mc.impl.hh>
#endif

#endif // DUNE_COPASI_MODEL_DIFFUSION_REACTION_MULTI_COMPARTMENT_HH
