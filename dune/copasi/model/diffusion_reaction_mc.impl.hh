#ifndef DUNE_COPASI_MODEL_DIFFUSION_REACTION_MULTI_COMPARTMENT_IMPL_HH
#define DUNE_COPASI_MODEL_DIFFUSION_REACTION_MULTI_COMPARTMENT_IMPL_HH

#include <dune/copasi/common/ostream_to_spdlog.hh>
#include <dune/copasi/concepts/grid.hh>
#include <dune/copasi/local_operator/diffusion_reaction/continuous_galerkin.hh>
#include <dune/copasi/model/diffusion_reaction_mc.hh>
#include <dune/copasi/model/diffusion_reaction_sc.hh>
#include <dune/copasi/model/reduce.hh>
#include <dune/copasi/model/interpolate.hh>
#include <dune/copasi/model/make_initial.hh>
#include <dune/copasi/model/make_step_operator.hh>

#include <dune/pdelab/basis/backend/istl.hh>
#include <dune/pdelab/basis/basis.hh>
#include <dune/pdelab/basis/discrete_global_function.hh>
#include <dune/pdelab/common/trace.hh>
#include <dune/pdelab/concepts/multiindex.hh>

#include <dune/grid/io/file/vtk.hh>

#include <spdlog/spdlog.h>

#include <fmt/core.h>

#ifdef DUNE_COPASI_PRECOMPILED_MODE
#warning "Including this file in pre-compiled mode may defeat the purpose of pre-compilation"
#endif

namespace Dune::Copasi {

template<class Traits>
void
ModelMultiCompartmentDiffusionReaction<Traits>::interpolate(
  State& state,
  const std::unordered_map<std::string, GridFunction>& initial) const
{
  using MultiCompartmentBasis =
    PDELab::Basis<MultiCompartmentEntitySet, MultiCompartmentPreBasis, TypeTree::HybridTreePath<>>;

  using CoefficientsBackend = PDELab::ISTLUniformBackend<ScalarQuantity>;
  using Coefficients = typename MultiCompartmentBasis::template Container<CoefficientsBackend>;
  using MultiCompartmentBasis =
    PDELab::Basis<MultiCompartmentEntitySet, MultiCompartmentPreBasis, TypeTree::HybridTreePath<>>;

  using CoefficientsBackend = PDELab::ISTLUniformBackend<ScalarQuantity>;
  using Coefficients = typename MultiCompartmentBasis::template Container<CoefficientsBackend>;
  const auto& basis = any_cast<const MultiCompartmentBasis&>(state.basis);
  auto& coefficients = any_cast<Coefficients&>(state.coefficients);

  Dune::Copasi::interpolate(basis, coefficients, initial);
}

template<class Traits>
auto
ModelMultiCompartmentDiffusionReaction<Traits>::make_multi_compartment_pre_basis(
  const Grid& grid,
  const ParameterTree& config,
  std::shared_ptr<const FunctorFactory<Grid::dimensionworld>> functor_factory) -> MultiCompartmentPreBasis
{
  TRACE_EVENT("dune", "Basis::SetUp");
  spdlog::info("Setup grid function space");

  const auto& compartments_config = config.sub("compartments", true);

  std::vector<CompartmentPreBasis> compartment_pre_basis_vec;

  std::map<std::string, std::vector<std::string>> compartment2componets;
  const auto& scalar_fields_config = config.sub("scalar_field", true);
  for (const auto& component : scalar_fields_config.getSubKeys())
    compartment2componets[config[fmt::format("scalar_field.{}.compartment", component)]].push_back(component);

  for (const auto& compartment : compartments_config.getSubKeys()) {
    using SubDomainIndex = typename Grid::SubDomainIndex;
    const auto& components = compartment2componets[compartment];
    SubDomainIndex domain_id =
      compartments_config.sub(compartment, true).template get<SubDomainIndex>("id");

    CompartmentEntitySet sub_grid_view = grid.subDomain(domain_id).leafGridView();

    CompartmentPreBasis compartment_pre_basis =
      ModelDiffusionReaction<Traits>::make_compartment_pre_basis(
        sub_grid_view, compartment, components, scalar_fields_config, functor_factory);
    compartment_pre_basis_vec.emplace_back(compartment_pre_basis);
  }

  spdlog::info("Setup of multi-compartment grid function space");

  MultiCompartmentPreBasis multi_compartment_pre_basis =
    composite(MultiCompartmentMergingStrategy{}, compartment_pre_basis_vec);

  multi_compartment_pre_basis.name("compartments");
  return multi_compartment_pre_basis;
}

template<class Traits>
void
ModelMultiCompartmentDiffusionReaction<Traits>::setup_coefficient_vector(State& state)
{
  using MultiCompartmentBasis =
    PDELab::Basis<MultiCompartmentEntitySet, MultiCompartmentPreBasis, TypeTree::HybridTreePath<>>;

  using CoefficientsBackend = PDELab::ISTLUniformBackend<ScalarQuantity>;
  using Coefficients = typename MultiCompartmentBasis::template Container<CoefficientsBackend>;

  using MultiCompartmentBasis =
    PDELab::Basis<MultiCompartmentEntitySet, MultiCompartmentPreBasis, TypeTree::HybridTreePath<>>;

  using CoefficientsBackend = PDELab::ISTLUniformBackend<ScalarQuantity>;
  using Coefficients = typename MultiCompartmentBasis::template Container<CoefficientsBackend>;

  spdlog::info("Setup coefficient vector");
  const auto& basis = any_cast<const MultiCompartmentBasis&>(state.basis);
  state.coefficients = Coefficients{ basis.makeContainer(CoefficientsBackend{}) };
}

template<class Traits>
auto
ModelMultiCompartmentDiffusionReaction<Traits>::make_state(const std::shared_ptr<const Grid>& grid,
                                                           const ParameterTree& config) const
  -> std::unique_ptr<State>
{
  auto state_ptr = std::make_unique<State>();
  state_ptr->basis = makeBasis(MultiCompartmentEntitySet{ grid->leafGridView() },
                               make_multi_compartment_pre_basis(*grid, config, _functor_factory));
  setup_coefficient_vector(*state_ptr);
  state_ptr->grid = grid;
  state_ptr->time = TimeQuantity{ 0. };
  return state_ptr;
}

template<class Traits>
auto
ModelMultiCompartmentDiffusionReaction<Traits>::make_compartment_function(
  const std::shared_ptr<const State>& state,
  std::string_view name) const -> GridFunction
{
  using ScalarBasis = PDELab::Basis<CompartmentEntitySet,
                                    MultiCompartmentPreBasis,
                                    TypeTree::HybridTreePath<std::size_t, size_t>>;
  using MultiCompartmentBasis =
    PDELab::Basis<MultiCompartmentEntitySet, MultiCompartmentPreBasis, TypeTree::HybridTreePath<>>;
  using CompartmentBasis = PDELab::
    Basis<CompartmentEntitySet, MultiCompartmentPreBasis, TypeTree::HybridTreePath<std::size_t>>;

  using CoefficientsBackend = PDELab::ISTLUniformBackend<ScalarQuantity>;
  using Coefficients = typename MultiCompartmentBasis::template Container<CoefficientsBackend>;

  using ScalarBasis = PDELab::Basis<CompartmentEntitySet,
                                    MultiCompartmentPreBasis,
                                    TypeTree::HybridTreePath<std::size_t, size_t>>;
  using MultiCompartmentBasis =
    PDELab::Basis<MultiCompartmentEntitySet, MultiCompartmentPreBasis, TypeTree::HybridTreePath<>>;
  using CompartmentBasis = PDELab::
    Basis<CompartmentEntitySet, MultiCompartmentPreBasis, TypeTree::HybridTreePath<std::size_t>>;

  using CoefficientsBackend = PDELab::ISTLUniformBackend<ScalarQuantity>;
  using Coefficients = typename MultiCompartmentBasis::template Container<CoefficientsBackend>;

  const auto& basis = any_cast<const MultiCompartmentBasis&>(state->basis);
  const auto& coeff = any_cast<const Coefficients&>(state->coefficients);
  std::shared_ptr<const Coefficients> coeff_ptr(state, &coeff);

  MultiCompartmentBasis multi_compartment_basis = basis.subSpace(TypeTree::treePath());
  for (std::size_t compartment = 0; compartment != multi_compartment_basis.degree();
       ++compartment) {
    CompartmentBasis compartment_basis =
      multi_compartment_basis.subSpace(TypeTree::treePath(compartment));
    for (std::size_t species = 0; species != compartment_basis.degree(); ++species) {
      ScalarBasis species_basis = compartment_basis.subSpace(TypeTree::treePath(species));
      if (species_basis.name() == name)
        return makeDiscreteGlobalBasisFunction(species_basis, coeff_ptr);
    }
  }
  throw format_exception(RangeError{}, "State doesn't contain any function with name: {}", name);
}

template<class Traits>
auto
ModelMultiCompartmentDiffusionReaction<Traits>::make_initial(const Grid& grid,
                                                             const ParameterTree& config) const
  -> std::unordered_map<std::string, GridFunction>
{
  return Dune::Copasi::make_initial<GridFunction>(grid, config, *_functor_factory);
}

template<class Traits>
auto
ModelMultiCompartmentDiffusionReaction<Traits>::reduce(const State& state, const ParameterTree& config) const
  -> std::map<std::string, double>
{
  using MultiCompartmentBasis =
    PDELab::Basis<MultiCompartmentEntitySet, MultiCompartmentPreBasis, TypeTree::HybridTreePath<>>;
  using CoefficientsBackend = PDELab::ISTLUniformBackend<ScalarQuantity>;
  using Coefficients = typename MultiCompartmentBasis::template Container<CoefficientsBackend>;

  const auto& basis = any_cast<const MultiCompartmentBasis&>(state.basis);
  const auto& coeff = any_cast<const Coefficients&>(state.coefficients);

  return Dune::Copasi::reduce(basis, coeff, state.time, config, _functor_factory);
}

template<class Traits>
auto
ModelMultiCompartmentDiffusionReaction<Traits>::make_step_operator(
  const State& state,
  const ParameterTree& config) const -> std::unique_ptr<PDELab::OneStep<State>>
{
  using MultiCompartmentBasis =
    PDELab::Basis<MultiCompartmentEntitySet, MultiCompartmentPreBasis, TypeTree::HybridTreePath<>>;

  using CoefficientsBackend = PDELab::ISTLUniformBackend<ScalarQuantity>;
  using Coefficients = typename MultiCompartmentBasis::template Container<CoefficientsBackend>;
  using ResidualBackend = PDELab::ISTLUniformBackend<ResidualQuantity>;
  using Residual = typename MultiCompartmentBasis::template Container<ResidualBackend>;

  using MultiCompartmentBasis =
    PDELab::Basis<MultiCompartmentEntitySet, MultiCompartmentPreBasis, TypeTree::HybridTreePath<>>;

  using CoefficientsBackend = PDELab::ISTLUniformBackend<ScalarQuantity>;
  using Coefficients = typename MultiCompartmentBasis::template Container<CoefficientsBackend>;
  using ResidualBackend = PDELab::ISTLUniformBackend<ResidualQuantity>;
  using Residual = typename MultiCompartmentBasis::template Container<ResidualBackend>;

  const auto& basis = any_cast<const MultiCompartmentBasis&>(state.basis);

  using LocalOperator = LocalOperatorDiffusionReactionCG<
    MultiCompartmentBasis,
    typename ScalarFiniteElementMap::Traits::FiniteElement::Traits::LocalBasisType::Traits>;

  spdlog::info("Creating mass/stiffness local operator");
  LocalOperator stiff_lop{ basis,
                           LocalOperatorType::Stiffness,
                           config.get("is_linear", false),
                           config.sub("scalar_field"),
                           _functor_factory };
  LocalOperator mass_lop{ basis,
                          LocalOperatorType::Mass,
                          config.get("is_linear", false),
                          config.sub("scalar_field"),
                          _functor_factory };

  std::shared_ptr one_step =
    Dune::Copasi::make_step_operator<Coefficients, Residual, ResidualQuantity, TimeQuantity>(
      config.sub("time_step_operator"), basis, mass_lop, stiff_lop);

  // type erase the original runge kutta operator
  auto type_erased_one_step = std::make_unique<PDELab::OperatorAdapter<State, State>>(
    [one_step](PDELab::Operator<State, State>& self, State& domain, State& range) mutable {
      auto log_guard = ostream2spdlog();
      // copy contents of this operator into runge-kutta operator
      static_cast<PDELab::PropertyTree&>(*one_step) =
        static_cast<const PDELab::PropertyTree&>(self);
      Coefficients& domain_coeff = any_cast<Coefficients&>(domain.coefficients);
      Residual& range_coeff = any_cast<Residual&>(range.coefficients);
      return one_step->apply(domain_coeff, range_coeff);
    });

  // assign properties of the original one step to the type erased one
  static_cast<PDELab::PropertyTree&>(*type_erased_one_step) =
    static_cast<const PDELab::PropertyTree&>(*one_step);

  auto residual_ptr = std::make_shared<Residual>(basis.makeContainer(ResidualBackend{}));
  type_erased_one_step->get("initial_residual") = residual_ptr;
  type_erased_one_step->get("time") = state.time;
  return type_erased_one_step;
}

template<class Traits>
void
ModelMultiCompartmentDiffusionReaction<Traits>::write_vtk(const State& state,
                                                      const fs::path& path,
                                                      bool append) const
{
  using CompartmentBasis = PDELab::
    Basis<CompartmentEntitySet, MultiCompartmentPreBasis, TypeTree::HybridTreePath<std::size_t>>;
  using MultiCompartmentBasis =
    PDELab::Basis<MultiCompartmentEntitySet, MultiCompartmentPreBasis, TypeTree::HybridTreePath<>>;
  using ScalarBasis = PDELab::Basis<CompartmentEntitySet,
                                    MultiCompartmentPreBasis,
                                    TypeTree::HybridTreePath<std::size_t, size_t>>;
  using CoefficientsBackend = PDELab::ISTLUniformBackend<ScalarQuantity>;
  using Coefficients = typename MultiCompartmentBasis::template Container<CoefficientsBackend>;

  const auto& basis = any_cast<const MultiCompartmentBasis&>(state.basis);
  const auto& coeff = any_cast<const Coefficients&>(state.coefficients);
  // warning: we use Dune::stackobject_to_shared_ptr to avoid copying the coefficients vector,
  // be we need to check that no one else took ownership of this pointer when we leave this function!
  std::shared_ptr<const Coefficients> coeff_ptr = Dune::stackobject_to_shared_ptr(coeff);

  // create directory if necessary
  auto path_entry = fs::directory_entry{ path };
  if (not path_entry.exists()) {
    spdlog::info("Creating output directory '{}'", path_entry.path().string());
    std::error_code ec{ 0, std::generic_category() };
    fs::create_directories(path_entry.path(), ec);
    if (ec)
      throw format_exception(IOError{},
                             "Category: {}\n"
                             "Value: {}\n"
                             "Message: {}",
                             ec.category().name(),
                             ec.value(),
                             ec.message());
  }

  spdlog::info("Writing solution for {:.2f}s time stamp", state.time);
  MultiCompartmentBasis multi_compartment_basis = basis.subSpace(TypeTree::treePath());

  // Recover old timestesps in case something was written before
  auto& timesteps = _writer_timesteps[path.string()];
  std::vector<double> tmp_timesteps;

  for (std::size_t compartment = 0; compartment != multi_compartment_basis.degree();
       ++compartment) {
    CompartmentBasis compartment_basis =
      multi_compartment_basis.subSpace(TypeTree::treePath(compartment));
    std::string name = fmt::format("{}-{}", path.filename().string(), compartment_basis.name());
    if (not append) {
      timesteps.clear();
      spdlog::info("Creating a time sequence file: '{}.pvd'", name);
    } else {
      spdlog::info("Overriding time sequence file: '{}.pvd'", name);
    }

    // setup writer again with old timesteps if necessary
    auto writer = std::make_shared<VTKWriter<CompartmentEntitySet>>(compartment_basis.entitySet(),
                                                                    Dune::VTK::conforming);
    auto sequential_writer =
      VTKSequenceWriter<CompartmentEntitySet>{ writer, name, path.string(), path.string() };
    sequential_writer.setTimeSteps(timesteps);

    for (std::size_t species = 0; species != compartment_basis.degree(); ++species) {
      ScalarBasis species_basis =
        compartment_basis.subSpace(TypeTree::treePath(std::size_t{ species }));
      sequential_writer.vtkWriter()->addVertexData(
        makeDiscreteGlobalBasisFunction(species_basis, coeff_ptr),
        VTK::FieldInfo{ species_basis.name(), VTK::FieldInfo::Type::scalar, 1 });
    }
    spdlog::info("Writing vtu file: '{0}/{0}-{1:0>5}.vtu'", name, timesteps.size());
    sequential_writer.write(state.time, Dune::VTK::base64);
    sequential_writer.vtkWriter()->clear();

    tmp_timesteps = sequential_writer.getTimeSteps();
  }

  timesteps = tmp_timesteps;

  if (coeff_ptr.use_count() != 1)
    throw format_exception(
      InvalidStateException{},
      "Fake shared pointer from coefficient vector may have been leaked outsie of this function!");
}

} // namespace Dune::Copasi

#endif // DUNE_COPASI_MODEL_DIFFUSION_REACTION_MULTI_COMPARTMENT_IMPL_HH
