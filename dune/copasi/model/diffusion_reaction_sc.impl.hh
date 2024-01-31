#ifndef DUNE_COPASI_MODEL_DIFFUSION_REACTION_SINGLE_COMPARTMENT_IMPL_HH
#define DUNE_COPASI_MODEL_DIFFUSION_REACTION_SINGLE_COMPARTMENT_IMPL_HH

#include <dune/copasi/model/diffusion_reaction_sc.hh>
#include <dune/copasi/model/reduce.hh>
#include <dune/copasi/model/interpolate.hh>
#include <dune/copasi/model/make_initial.hh>
#include <dune/copasi/model/make_diffusion_reaction_step_operator.hh>
#include <dune/copasi/model/local_equations/functor_factory_parser.hh>

#include <dune/copasi/common/exceptions.hh>
#include <dune/copasi/common/ostream_to_spdlog.hh>
#include <dune/copasi/concepts/grid.hh>
#include <dune/copasi/grid/has_single_geometry_type.hh>
#include <dune/copasi/grid/boundary_entity_mapper.hh>
#include <dune/copasi/parser/factory.hh>

#include <dune/pdelab/basis/backend/istl.hh>
#include <dune/pdelab/basis/basis.hh>
#include <dune/pdelab/basis/discrete_global_function.hh>
#include <dune/pdelab/common/partition/identity.hh>
#include <dune/pdelab/common/trace.hh>
#include <dune/pdelab/concepts/multiindex.hh>

#include <dune/grid/io/file/vtk.hh>

#include <spdlog/spdlog.h>

#ifdef DUNE_COPASI_PRECOMPILED_MODE
#warning "Including this file in pre-compiled mode may defeat the purpose of pre-compilation"
#endif

namespace Dune::Copasi {

template<class Traits, Dune::Concept::Grid MDGrid>
auto
ModelDiffusionReaction<Traits, MDGrid>::get_entity_set(const Grid& grid, std::size_t subdomain)
  -> CompartmentEntitySet
{
  if constexpr (std::same_as<typename Grid::LeafGridView, CompartmentEntitySet>) {
    return grid.leafGridView();
  } else if constexpr (Concept::SubDomainGrid<typename CompartmentEntitySet::Grid>) {
    static_assert(std::same_as<typename Grid::SubDomainGrid::LeafGridView, CompartmentEntitySet>);
    return grid.subDomain(subdomain).leafGridView();
  }
  throw format_exception(NotImplemented{}, "Not known mapping from Grid to CompartmentEntitySet");
}

template<class Traits, Dune::Concept::Grid MDGrid>
void
ModelDiffusionReaction<Traits, MDGrid>::interpolate(
  State& state,
  const std::unordered_map<std::string, GridFunction>& initial) const
{
  using CompartmentBasis = PDELab::Basis<PDELab::EntitySetPartitioner::Identity<CompartmentEntitySet>, CompartmentPreBasis>;
  using CoefficientsBackend = PDELab::ISTLUniformBackend<ScalarQuantity>;
  using Coefficients = typename CompartmentBasis::template Container<CoefficientsBackend>;
  const auto& basis = any_cast<const CompartmentBasis&>(state.basis);
  auto& coefficients = any_cast<Coefficients&>(state.coefficients);

  Dune::Copasi::interpolate(basis, coefficients, initial);
}

template<class Traits, Dune::Concept::Grid MDGrid>
auto
ModelDiffusionReaction<Traits, MDGrid>::make_scalar_field_pre_basis(std::shared_ptr<BoundaryEntityMapper<CompartmentEntitySet>> boundary_mapper,
                                                            const CompartmentEntitySet& entity_set,
                                                            std::string_view name,
                                                            const ParameterTree& scalar_field_config,
                                                            std::shared_ptr<const FunctorFactory< MDGrid >> functor_factory) -> ScalarPreBasis
{
  spdlog::info("Setup basis functions for component '{}'", name);
  auto scalar_field_pre_basis =
    ScalarPreBasis{ ScalarMergingStrategy{ entity_set },
                    std::make_shared<ScalarFiniteElementMap>(entity_set),
                    Constraints{boundary_mapper, scalar_field_config.sub("constrain"), functor_factory} };
  scalar_field_pre_basis.name(name);
  return scalar_field_pre_basis;
}

template<class Traits, Dune::Concept::Grid MDGrid>
auto
ModelDiffusionReaction<Traits, MDGrid>::make_compartment_pre_basis(
  const CompartmentEntitySet& entity_set,
  std::string_view compartment_name,
  const std::vector<std::string>& scalar_field_names,
  const ParameterTree& scalar_fields_config,
  std::shared_ptr<const FunctorFactory< MDGrid >> functor_factory) -> CompartmentPreBasis
{
  spdlog::info("Setup compartment basis functions for compartment '{}'", compartment_name);

  if (entity_set.size(0) == 0)
    spdlog::warn("Compartment '{}' is empty", compartment_name);

  auto boundary_mapper = std::make_shared<BoundaryEntityMapper<CompartmentEntitySet>>(entity_set);
  std::vector<ScalarPreBasis> scalar_field_pre_basis;
  for (const auto& name : scalar_field_names) {
    scalar_field_pre_basis.push_back(make_scalar_field_pre_basis(boundary_mapper, entity_set, name, scalar_fields_config.sub(name), functor_factory));
  }

  CompartmentPreBasis compartment_pre_basis =
    composite(CompartmentMergingStrategy{ entity_set }, scalar_field_pre_basis );
  compartment_pre_basis.name(compartment_name);

  spdlog::info("No. of components on '{}' compartment: {}",
               compartment_pre_basis.name(),
               compartment_pre_basis.degree());

  return compartment_pre_basis;
}

template<class Traits, Dune::Concept::Grid MDGrid>
void
ModelDiffusionReaction<Traits, MDGrid>::setup_basis(State& state,
                                            const Grid& grid,
                                            const ParameterTree& config,
                                            std::shared_ptr<const FunctorFactory< MDGrid >> functor_factory)
{
  TRACE_EVENT("dune", "Basis::SetUp");
  using CompartmentBasis = PDELab::Basis<PDELab::EntitySetPartitioner::Identity<CompartmentEntitySet>, CompartmentPreBasis>;
  const auto& compartments_config = config.sub("compartments", true);
  const auto& compartments = compartments_config.getSubKeys();
  if (compartments.size() != 1) {
    throw format_exception(InvalidStateException{}, "Config file should only have one compartment");
  }
  auto compartment = compartments.front();
  auto entity_set = get_entity_set(
    grid, compartments_config.sub(compartment, true).template get<std::size_t>("id"));

  std::vector<std::string> components;
  const auto& scalar_fields_config = config.sub("scalar_field", true);
  for (const auto& scalar_field : scalar_fields_config.getSubKeys()) {
    if (scalar_fields_config.sub(scalar_field)["compartment"] == compartment) {
      components.push_back(scalar_field);
    }
  }
  auto comp_space = make_compartment_pre_basis(entity_set, compartment, components, scalar_fields_config, functor_factory);
  state.basis = CompartmentBasis{ makeBasis(entity_set, comp_space) };
}

template<class Traits, Dune::Concept::Grid MDGrid>
void
ModelDiffusionReaction<Traits, MDGrid>::setup_coefficient_vector(State& state)
{
  spdlog::info("Setup coefficient vector");
  using CompartmentBasis = PDELab::Basis<PDELab::EntitySetPartitioner::Identity<CompartmentEntitySet>, CompartmentPreBasis>;
  using CoefficientsBackend = PDELab::ISTLUniformBackend<ScalarQuantity>;
  using Coefficients = typename CompartmentBasis::template Container<CoefficientsBackend>;
  const auto& basis = any_cast<const CompartmentBasis&>(state.basis);
  state.coefficients = Coefficients{ basis.makeContainer(CoefficientsBackend{}) };
}

template<class Traits, Dune::Concept::Grid MDGrid>
auto
ModelDiffusionReaction<Traits, MDGrid>::make_state(const std::shared_ptr<const Grid>& grid,
                                           const ParameterTree& config) const
  -> std::unique_ptr<State>
{
  auto state_ptr = std::make_unique<State>();
  setup_basis(*state_ptr, *grid, config, _functor_factory);
  setup_coefficient_vector(*state_ptr);
  state_ptr->grid = grid;
  state_ptr->time = TimeQuantity{ 0. };
  return state_ptr;
}

template<class Traits, Dune::Concept::Grid MDGrid>
typename ModelDiffusionReaction<Traits, MDGrid>::GridFunction
ModelDiffusionReaction<Traits, MDGrid>::make_compartment_function(const std::shared_ptr<const State>& state,
                                                          std::string_view name) const
{
  using CompartmentBasis = PDELab::Basis<PDELab::EntitySetPartitioner::Identity<CompartmentEntitySet>, CompartmentPreBasis>;
  using CoefficientsBackend = PDELab::ISTLUniformBackend<ScalarQuantity>;
  using Coefficients = typename CompartmentBasis::template Container<CoefficientsBackend>;
  const auto& basis = any_cast<const CompartmentBasis&>(state->basis);
  const auto& coeff = any_cast<const Coefficients&>(state->coefficients);
  std::shared_ptr<const Coefficients> const coeff_ptr(state, &coeff);

  for (std::size_t component = 0; component != basis.degree(); ++component) {
    auto leaf_space = basis.subSpace(TypeTree::treePath(std::size_t{ component }));
    if (leaf_space.name() == name) {
      return makeDiscreteGlobalBasisFunction(leaf_space, coeff_ptr);
    }
  }
  throw format_exception(RangeError{}, "State doesn't contain any function with name: {}", name);
}

template<class Traits, Dune::Concept::Grid MDGrid>
auto
ModelDiffusionReaction<Traits, MDGrid>::make_initial(const Grid& grid, const ParameterTree& config) const
  -> std::unordered_map<std::string, GridFunction>
{
  return Dune::Copasi::make_initial<GridFunction>(grid, config, *_functor_factory);
}

template<class Traits, Dune::Concept::Grid MDGrid>
auto
ModelDiffusionReaction<Traits, MDGrid>::make_step_operator(const State& state,
                                                   const ParameterTree& config) const
  -> std::unique_ptr<PDELab::OneStep<State>>
{
  using CompartmentBasis = PDELab::Basis<PDELab::EntitySetPartitioner::Identity<CompartmentEntitySet>, CompartmentPreBasis>;
  using CoefficientsBackend = PDELab::ISTLUniformBackend<ScalarQuantity>;
  using Coefficients = typename CompartmentBasis::template Container<CoefficientsBackend>;
  using ResidualBackend = PDELab::ISTLUniformBackend<ResidualQuantity>;
  using Residual = typename CompartmentBasis::template Container<ResidualBackend>;
  using LocalBasisTraits = typename ScalarFiniteElementMap::Traits::FiniteElement::Traits::LocalBasisType::Traits;

  const auto& basis = any_cast<const CompartmentBasis&>(state.basis);

  std::shared_ptr one_step = make_diffusion_reaction_step_operator<LocalBasisTraits, Coefficients, Residual, ResidualQuantity, TimeQuantity, MDGrid>(config, basis, 1, _functor_factory);

  // type erase the original runge kutta operator
  auto type_erased_one_step = std::make_unique<PDELab::OperatorAdapter<State, State>>(
    [one_step](PDELab::Operator<State, State>& self, State& domain, State& range) mutable {
      auto log_guard = ostream2spdlog();
      // copy contents of this operator into runge-kutta operator
      static_cast<PDELab::PropertyTree&>(*one_step) =
        static_cast<const PDELab::PropertyTree&>(self);
      auto& domain_coeff = any_cast<Coefficients&>(domain.coefficients);
      auto& range_coeff = any_cast<Residual&>(range.coefficients);
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

template<class Traits, Dune::Concept::Grid MDGrid>
auto
ModelDiffusionReaction<Traits, MDGrid>::reduce(const State& state, const ParameterTree& config) const
  -> std::map<std::string, double>
{
  using CompartmentBasis = PDELab::Basis<PDELab::EntitySetPartitioner::Identity<CompartmentEntitySet>, CompartmentPreBasis>;
  using CoefficientsBackend = PDELab::ISTLUniformBackend<ScalarQuantity>;
  using Coefficients = typename CompartmentBasis::template Container<CoefficientsBackend>;
  const auto& basis = any_cast<const CompartmentBasis&>(state.basis);
  const auto& coeff = any_cast<const Coefficients&>(state.coefficients);

  return Dune::Copasi::reduce(basis, coeff, state.time, config, _functor_factory);
}

template<class Traits, Dune::Concept::Grid MDGrid>
void
ModelDiffusionReaction<Traits, MDGrid>::write_vtk(const State& state,
                                          const std::filesystem::path& path,
                                          bool append) const
{
  using CompartmentBasis = PDELab::Basis<PDELab::EntitySetPartitioner::Identity<CompartmentEntitySet>, CompartmentPreBasis>;
  using CoefficientsBackend = PDELab::ISTLUniformBackend<ScalarQuantity>;
  using Coefficients = typename CompartmentBasis::template Container<CoefficientsBackend>;
  const auto& basis = any_cast<const CompartmentBasis&>(state.basis);
  const auto& coeff = any_cast<const Coefficients&>(state.coefficients);
  // warning: we use Dune::stackobject_to_shared_ptr to avoid copying the coefficients vector,
  // be we need to check that no one else took ownership of this pointer when we leave this
  // function!
  std::shared_ptr<const Coefficients> const coeff_ptr = Dune::stackobject_to_shared_ptr(coeff);

  // create directory if necessary
  auto path_entry = std::filesystem::directory_entry{ path };
  if (not path_entry.exists()) {
    spdlog::info("Creating output directory '{}'", path_entry.path().string());
    std::error_code ec{ 0, std::generic_category() };
    std::filesystem::create_directories(path_entry.path(), ec);
    if (ec) {
      throw format_exception(IOError{},
                             "\n Category: {}\nValue: {}\nMessage: {}\n",
                             ec.category().name(),
                             ec.value(),
                             ec.message());
    }
  }

  // Recover old timestesps in case something was written before
  auto& timesteps = _writer_timesteps[path.string()];
  std::string name = fmt::format("{}-{}", path.filename().string(), basis.name());
  if (not append) {
    timesteps.clear();
    spdlog::info("Creating a time sequence file: '{}.pvd'", name);
  } else {
    spdlog::info("Overriding time sequence file: '{}.pvd'", name);
  }

  // setup writer again with old timesteps if necessary
  auto writer =
    std::make_shared<VTKWriter<CompartmentEntitySet>>(basis.entitySet(), Dune::VTK::conforming);
  auto sequential_writer =
    VTKSequenceWriter<CompartmentEntitySet>{ writer, name, path.string(), path.string() };
  sequential_writer.setTimeSteps(timesteps);

  for (std::size_t component = 0; component != basis.degree(); ++component) {
    auto const component_basis =
      basis.subSpace(TypeTree::treePath(std::size_t{ component }));
    sequential_writer.vtkWriter()->addVertexData(
      makeDiscreteGlobalBasisFunction(component_basis, coeff_ptr),
      VTK::FieldInfo{ component_basis.name(), VTK::FieldInfo::Type::scalar, 1 });
  }

  spdlog::info("Writing solution for {:.2f}s time stamp", state.time);
  spdlog::info("Writing vtu file: '{0}/{0}-{1:0>5}.vtu'", name, timesteps.size());
  TRACE_EVENT("dune", "WriteVTK");
  sequential_writer.write(state.time, Dune::VTK::base64);
  sequential_writer.vtkWriter()->clear();
  timesteps = sequential_writer.getTimeSteps();

  if (coeff_ptr.use_count() != 1) {
    throw format_exception(
      InvalidStateException{},
      "Fake shared pointer from coefficient vector may have been leaked outside of this"
      "function!");
  }
}

} // namespace Dune::Copasi

#endif // DUNE_COPASI_MODEL_DIFFUSION_REACTION_SINGLE_COMPARTMENT_IMPL_HH
