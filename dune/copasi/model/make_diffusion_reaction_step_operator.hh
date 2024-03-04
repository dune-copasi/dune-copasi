#ifndef DUNE_COPASI_MODEL_MAKE_DIFFUSION_REACTION_STEP_OPERATOR_HH
#define DUNE_COPASI_MODEL_MAKE_DIFFUSION_REACTION_STEP_OPERATOR_HH

#include <dune/copasi/local_operator/diffusion_reaction/continuous_galerkin.hh>
#include <dune/copasi/model/make_step_operator.hh>

#if HAVE_METIS && DUNE_COPASI_CONCURRENT_ASSEMBLY
#include <dune/pdelab/common/partition/metis.hh>
#endif

#include <dune/pdelab/common/partition/simple.hh>
#include <dune/pdelab/common/execution.hh>
#include <dune/pdelab/operator/operator.hh>

#include <spdlog/spdlog.h>

#include <memory>

namespace Dune::Copasi {

template<class LocalBasisTraits,
         class Coefficients,
         class Residual,
         class ResidualQuantity,
         class TimeQuantity,
         PDELab::Concept::Basis Basis>
[[nodiscard]] inline static std::unique_ptr<PDELab::Operator<Coefficients, Coefficients>>
make_diffusion_reaction_step_operator(const ParameterTree& config,
                                      const Basis& basis,
                                      std::size_t halo,
                                      std::shared_ptr<const FunctorFactory<Basis::EntitySet::GridView::dimension>> functor_factory,
                                      std::shared_ptr<const GridDataContext<typename Basis::EntitySet::GridView>> grid_data_context )
{

  std::unique_ptr<PDELab::Operator<Coefficients, Coefficients>> one_step;
  const auto& assembly_cfg = config.sub("assembly");

  auto make_one_step_op = [&, functor_factory, grid_data_context]<class ExecutionPolicy, PDELab::Concept::Basis OperatorBasis>(
                            ExecutionPolicy execution_policy, const OperatorBasis& operator_basis) {
    spdlog::info("Creating mass/stiffness local operator");
    const auto& time_step_cfg = config.sub("time_step_operator");
    const auto& scalar_field_cfg = config.sub("scalar_field");
    bool is_linear = config.get("is_linear", false);
    using LocalOperator =
      LocalOperatorDiffusionReactionCG<OperatorBasis, LocalBasisTraits, ExecutionPolicy>;
    LocalOperator const stiff_lop(operator_basis,
                                  LocalOperatorType::Stiffness,
                                  is_linear,
                                  scalar_field_cfg,
                                  functor_factory,
                                  grid_data_context,
                                  execution_policy);
    LocalOperator const mass_lop(operator_basis,
                                 LocalOperatorType::Mass,
                                 is_linear,
                                 scalar_field_cfg,
                                 functor_factory,
                                 grid_data_context,
                                 execution_policy);
    return Dune::Copasi::make_step_operator<Coefficients, Residual, ResidualQuantity, TimeQuantity>(
      time_step_cfg, operator_basis, mass_lop, stiff_lop);
  };

  auto entity_set = basis.entitySet();
  std::size_t part_patches =
    assembly_cfg.get("partition.patches", static_cast<std::size_t>(entity_set.size(0) / 40));

#if DUNE_COPASI_CONCURRENT_ASSEMBLY
  auto make_concurrent_one_step_op =
    [&]<class ExecutionPolicy>(const ExecutionPolicy execution_policy) {
      auto entity_set = basis.entitySet();
      auto default_part_type =
#if HAVE_METIS
      "METIS";
#else
      "Simple";
#endif // HAVE_METIS
      auto part_type = assembly_cfg.get("partition.type", default_part_type);
      auto part_coloring = assembly_cfg.get("partition.coloring", "None");
      spdlog::info("  Concurrent Partitioning:");
      spdlog::info("    Type:      {}", part_type);
      spdlog::info("    Patches:   {}", part_patches);
      spdlog::info("    Halo:      {}", halo);
      spdlog::info("    Coloring:  {}", part_coloring);
      if (part_type == "METIS") {
#if HAVE_METIS
        if (part_coloring == "None") {
          return make_one_step_op(
            execution_policy,
            basis.subSpace(
              PDELab::EntitySetPartitioner::MetisColored{ entity_set, part_patches, halo },
              TypeTree::treePath()));
        } else if (part_coloring == "DSatur") {
          return make_one_step_op(
            execution_policy,
            basis.subSpace(PDELab::EntitySetPartitioner::Metis{ entity_set, part_patches, halo },
                           TypeTree::treePath()));
        } else {
          throw format_exception(IOError{}, "Not known coloring algorithm '{}' known", part_coloring);
        }
#else
        throw format_exception(IOError{}, "'METIS' partitioner is not available on this executable");
#endif // HAVE_METIS
      } else if (part_type == "Simple") {
        if (part_coloring == "None") {
          return make_one_step_op(
            execution_policy,
            basis.subSpace(
              PDELab::EntitySetPartitioner::Simple{ entity_set, part_patches, halo },
              TypeTree::treePath()));
        } else if (part_coloring == "DSatur") {
          return make_one_step_op(
            execution_policy,
            basis.subSpace(PDELab::EntitySetPartitioner::SimpleColored{ entity_set, part_patches, halo },
                           TypeTree::treePath()));
        } else {
          throw format_exception(IOError{}, "Not known coloring algorithm '{}' known", part_coloring);
        }
      } else {
        throw format_exception(IOError{}, "Not known partition algorithm '{}' known", part_type);
      }
    };

  auto exec_policy = assembly_cfg.get("type", "concurrent");
#else
  auto exec_policy = assembly_cfg.get("type", "sequential");
#endif // DUNE_COPASI_CONCURRENT_ASSEMBLY

  spdlog::info("Creating time-step operator with '{}' execution policy", exec_policy);
  if (exec_policy != "sequential" and part_patches < std::thread::hardware_concurrency()) {
    exec_policy = "sequential";
    spdlog::warn("Too few patches '{}', execution policy of assembler is switched to 'sequential'", part_patches);
  }

  if (exec_policy == "sequential") {
    one_step = make_one_step_op(PDELab::Execution::seq, basis);
#if DUNE_COPASI_CONCURRENT_ASSEMBLY
  } else if (exec_policy == "concurrent") {
    one_step = make_concurrent_one_step_op(PDELab::Execution::par);
#endif // DUNE_COPASI_CONCURRENT_ASSEMBLY
  } else {
    throw format_exception(IOError{}, "Not known execution '{}' policy known", exec_policy);
  }

  return one_step;
}

} // namespace Dune::Copasi

#endif // DUNE_COPASI_MODEL_MAKE_DIFFUSION_REACTION_STEP_OPERATOR_HH
