#ifndef DUNE_COPASI_MODEL_DIFFUSION_REACTION_SINGLE_COMAPARTMENT_HH
#define DUNE_COPASI_MODEL_DIFFUSION_REACTION_SINGLE_COMAPARTMENT_HH

// file: diffusion reaction for single compartment models

#include <dune/copasi/model/constraints.hh>
#include <dune/copasi/model/functor_factory.hh>
#include <dune/copasi/model/model.hh>
#include <dune/copasi/grid/boundary_entity_mapper.hh>
#include <dune/copasi/grid/cell_data.hh>

#include <dune/pdelab/basis/prebasis/composite.hh>
#include <dune/pdelab/basis/prebasis/leaf.hh>

#include <dune/grid/concepts/grid.hh>
#include <dune/grid/concepts/gridview.hh>

namespace Dune::Copasi::DiffusionReaction {

template<class Traits>
class ModelSingleCompartment
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
  using ScalarFiniteElementMap = typename Traits::ScalarFiniteElementMap;
  using ScalarMergingStrategy = typename Traits::ScalarMergingStrategy;
  using ScalarPreBasis =
    PDELab::PreBasis<ScalarMergingStrategy, ScalarFiniteElementMap, Constraints<CompartmentEntitySet>>;
  using CompartmentMergingStrategy = typename Traits::CompartmentMergingStrategy;
  using CompartmentPreBasis = PDELab::PreBasisVector<CompartmentMergingStrategy, ScalarPreBasis>;
  using ResidualQuantity = ScalarQuantity;

  using GridFunction = typename Base::GridFunction;

  explicit ModelSingleCompartment(
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

  void write_vtk(const State&, const std::filesystem::path&, bool = true) const override;

  std::map<std::string, double> reduce(const State&, const ParameterTree&) const override;

  static CompartmentPreBasis make_compartment_pre_basis(const CompartmentEntitySet&,
                                                        std::string_view,
                                                        const std::vector<std::string>&,
                                                        const ParameterTree& = {},
                                                        std::shared_ptr<const FunctorFactory<Grid::dimensionworld>> = nullptr);

private:
  static ScalarPreBasis make_scalar_field_pre_basis(std::shared_ptr<BoundaryEntityMapper<CompartmentEntitySet>>, const CompartmentEntitySet&, std::string_view, const ParameterTree&, std::shared_ptr<const FunctorFactory<Grid::dimensionworld>>);
  static void setup_basis(State&, const Grid&, const ParameterTree&, std::shared_ptr<const FunctorFactory<Grid::dimensionworld>>);
  static void setup_coefficient_vector(State&);
  static CompartmentEntitySet get_entity_set(const Grid&, std::size_t);

  mutable std::unordered_map<std::string, std::vector<double>> _writer_timesteps;
  std::shared_ptr<const FunctorFactory<Grid::dimensionworld>> _functor_factory;
  std::shared_ptr<const CellData<typename Grid::LeafGridView, ScalarQuantity>> _cell_data;
};

} // namespace Dune::Copasi::DiffusionReaction

#ifndef DUNE_COPASI_PRECOMPILED_MODE
#include <dune/copasi/model/diffusion_reaction/model_single_compartment.impl.hh>
#endif

#endif // DUNE_COPASI_MODEL_DIFFUSION_REACTION_SINGLE_COMAPARTMENT_HH
