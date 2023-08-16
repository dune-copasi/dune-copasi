#ifndef DUNE_COPASI_MODEL_LOCAL_EQUATIONS_LOCAL_EQUATIONS_HH
#define DUNE_COPASI_MODEL_LOCAL_EQUATIONS_LOCAL_EQUATIONS_HH

#include <dune/copasi/common/bit_flags.hh>
#include <dune/copasi/common/exceptions.hh>
#include <dune/copasi/concepts/compartment_tree.hh>
#include <dune/copasi/model/local_equations/functor_factory.hh>
#include <dune/copasi/parser/factory.hh>

#include <dune/pdelab/concepts/basis.hh>
#include <dune/pdelab/common/container_entry.hh>

#include <dune/common/fvector.hh>
#include <dune/common/indices.hh>
#include <dune/common/overloadset.hh>
#include <dune/common/parametertree.hh>
#include <dune/common/tuplevector.hh>

#include <function2/function2.hpp>
#include <spdlog/spdlog.h>

#include <functional>
#include <utility>
#include <variant>
#include <set>

namespace Dune::Copasi {

template<std::size_t dim>
struct LocalDomain
{
  FieldVector<double, dim> position;
  FieldVector<double, dim> normal;
  double time = 0.;
  double entity_volume = 0.;
  double cell_index = 0.;
  virtual ~LocalDomain() {}
};

// this class holds a data-structure for each equation that contains functors to be evaluated.
// Additionally, it contains the values with respect these functors may be evaluated if they are
// non-linear All the functors are require to be thread-safe!
template<std::size_t dim>
class LocalEquations : public LocalDomain<dim>
{

  using Scalar = FieldVector<double, 1>;
  using Vector = FieldVector<double, dim>;
  using Tensor = FieldMatrix<double, dim, dim>;

  using CompartmentPath = TypeTree::HybridTreePath<index_constant<0>, std::size_t, std::size_t>;
  using MembranePath = TypeTree::HybridTreePath<index_constant<1>, std::size_t, std::size_t>;

  static const Concept::CompartmentScalarLocalBasisNode auto& path_to_local_basis_node(
    CompartmentPath path,
    const Concept::CompartmentLocalBasisNode auto& lbasis)
  {
    assert(path[Indices::_1] == 0);
    return lbasis.child(path[Indices::_2]);
  }

  static const Concept::CompartmentScalarLocalBasisNode auto& path_to_local_basis_node(
    CompartmentPath path,
    const Concept::MultiCompartmentLocalBasisNode auto& lbasis)
  {
    return lbasis.child(path[Indices::_1]).child(path[Indices::_2]);
  }

  static const Concept::MembraneSubEntitiesLocalBasisNode auto& path_to_local_basis_node(
    MembranePath path,
    const Concept::MembraneLocalBasisNode auto& lbasis)
  {
    assert(path[Indices::_1] == 0);
    return lbasis.child(path[Indices::_2]);
  }

  static const Concept::MembraneSubEntitiesLocalBasisNode auto& path_to_local_basis_node(
    MembranePath path,
    const Concept::MultiCompartmentLocalBasisNode auto& lbasis)
  {
    return lbasis.child(path[Indices::_1]).child(path[Indices::_2]);
  }

  enum class FactoryFalgs
  {
    Reaction = 1 << 0,
    Diffusion = 1 << 1,
    Velocity = 1 << 2,
    Outflow = 1 << 3,
    Storage = 1 << 4
  };

public:
  class CompartmentNode;
  class MembraneNode;

  template<class Signature>
  struct CompartmentPartialDerivative : public fu2::unique_function<Signature>
  {
    CompartmentPartialDerivative(fu2::unique_function<Signature>&& callable,
                                 const CompartmentNode& _wrt)
      : fu2::unique_function<Signature>{ std::move(callable) }
      , wrt{ _wrt }
    {
    }
    const CompartmentNode& wrt;
  };

  template<class Signature>
  struct MembranePartialDerivative : public fu2::unique_function<Signature>
  {
    MembranePartialDerivative(fu2::unique_function<Signature>&& callable,
                              const MembraneNode& _wrt)
      : fu2::unique_function<Signature>{ std::move(callable) }
      , wrt{ _wrt }
    {
    }
    const MembraneNode& wrt;
  };

  template<class Signature>
  struct CompartmentDifferentiableFunction : public fu2::unique_function<Signature>
  {
    std::vector<CompartmentPartialDerivative<Signature>> compartment_jacobian;
  };

  template<class Signature>
  struct MembraneDifferentiableFunction : public fu2::unique_function<Signature>
  {
    std::vector<CompartmentPartialDerivative<Signature>> compartment_jacobian;
    std::vector<MembranePartialDerivative<Signature>> membrane_jacobian;
  };

  using CompartmentScalarFunction =
    CompartmentDifferentiableFunction<Scalar() const DUNE_COPASI_FUNCTOR_NOEXCEPT>;
  using MembraneScalarFunction =
    MembraneDifferentiableFunction<Scalar() const DUNE_COPASI_FUNCTOR_NOEXCEPT>;

  using CompartmentVectorFunction =
    CompartmentDifferentiableFunction<Vector() const DUNE_COPASI_FUNCTOR_NOEXCEPT>;
  using MembraneVectorFunction =
    MembraneDifferentiableFunction<Vector() const DUNE_COPASI_FUNCTOR_NOEXCEPT>;

  struct CompartmentDiffusionApply
    : public CompartmentDifferentiableFunction<Vector(Vector) const DUNE_COPASI_FUNCTOR_NOEXCEPT>
  {
    CompartmentDiffusionApply(
      fu2::unique_function<Vector(Vector) const DUNE_COPASI_FUNCTOR_NOEXCEPT>&& callable,
      const CompartmentNode& _wrt)
      : CompartmentDifferentiableFunction<Vector(
          Vector) const DUNE_COPASI_FUNCTOR_NOEXCEPT>{ std::move(callable) }
      , wrt{ _wrt }
    {
    }
    const CompartmentNode& wrt;
  };

  struct MembraneDiffusionApply
    : public MembraneDifferentiableFunction<Vector(Vector) const DUNE_COPASI_FUNCTOR_NOEXCEPT>
  {
    MembraneDiffusionApply(
      fu2::unique_function<Vector(Vector) const DUNE_COPASI_FUNCTOR_NOEXCEPT>&& callable,
      const MembraneNode& _wrt)
      : MembraneDifferentiableFunction<Vector(Vector)
                                         const DUNE_COPASI_FUNCTOR_NOEXCEPT>{ std::move(callable) }
      , wrt{ _wrt }
    {
    }
    const MembraneNode& wrt;
  };

  struct CompartmentNode
  {
    Scalar& value;
    Vector& gradient;

    CompartmentPath path;
    std::string name;

    CompartmentScalarFunction reaction;
    CompartmentScalarFunction storage;
    CompartmentVectorFunction velocity;

    std::vector<CompartmentDiffusionApply> cross_diffusion;
    std::vector<MembraneScalarFunction> outflow;

    const Concept::CompartmentScalarLocalBasisNode auto& to_local_basis_node(
      const PDELab::Concept::LocalBasis auto& lbasis) const
    {
      return path_to_local_basis_node(path, lbasis.tree());
    }

    void debug() const
    {
      std::cout << fmt::format("Name: {}\n", name);
      std::cout << "\tPath: " << path << std::endl;
      std::cout << fmt::format("\tValue: {}\n", value[0]);
      std::cout << "\tGradient: " << gradient << std::endl;
      if (reaction) {
        std::cout << fmt::format("\tReaction: {}\n", reaction()[0]);
        for (const auto& jac : reaction.compartment_jacobian)
          std::cout << fmt::format("\t\tJacobian wrt '{}': {}\n", jac.wrt.name, jac()[0]);
      }
      if (storage)
        std::cout << fmt::format("\tStorage: {}\n", storage()[0]);
      for (const auto& diffusion : cross_diffusion) {
        std::cout << fmt::format("\tCross Diffusion wrt '{}': ", diffusion.wrt.name)
                  << diffusion(diffusion.wrt.gradient) << std::endl;
        for (const auto& jac : diffusion.compartment_jacobian)
          std::cout << fmt::format("\t\tJacobian wrt '{}': ", jac.wrt.name) << jac(jac.wrt.gradient)
                    << std::endl;
      }
    }
  };

  struct MembraneNode
  {
    Scalar& value;
    Vector& gradient;

    MembranePath path;
    std::string name;
    bool is_linear = false;

    MembraneScalarFunction reaction;
    MembraneScalarFunction storage;
    MembraneVectorFunction velocity;

    std::vector<MembraneDiffusionApply> cross_diffusion;
    std::vector<MembraneScalarFunction> outflow;

    const Concept::MembraneSubEntitiesLocalBasisNode auto& to_local_basis_node(
      const PDELab::Concept::LocalBasis auto& lbasis) const
    {
      return path_to_local_basis_node(path, lbasis.tree());
    }

    void debug() const
    {
      std::cout << fmt::format("Name: {}\n", name);
      std::cout << "\tPath: " << path << std::endl;
      std::cout << fmt::format("\tValue: {}\n", value[0]);
      std::cout << "\tGradient: " << gradient << std::endl;
      if (reaction) {
        std::cout << fmt::format("\tReaction: {}\n", reaction()[0]);
        for (const auto& jac : reaction.compartment_jacobian)
          std::cout << fmt::format("\t\tJacobian wrt '{}': {}\n", jac.wrt.name, jac()[0]);
      }
      if (storage)
        std::cout << fmt::format("\tStorage: {}\n", storage()[0]);
      for (const auto& diffusion : cross_diffusion) {
        std::cout << fmt::format("\tCross Diffusion wrt '{}': ", diffusion.wrt.name)
                  << diffusion(diffusion.wrt.gradient) << std::endl;
        for (const auto& jac : diffusion.compartment_jacobian)
          std::cout << fmt::format("\t\tJacobian wrt '{}': ", jac.wrt.name) << jac(jac.wrt.gradient)
                    << std::endl;
      }
    }
  };

  [[nodiscard]] static std::unique_ptr<LocalEquations> make(
    const PDELab::Concept::LocalBasis auto& lbasis,
    const ParameterTree& config = {},
    FunctorFactory<dim> const * functor_factory = nullptr,
    BitFlags<FactoryFalgs> opts = BitFlags<FactoryFalgs>::no_flags())
  {
    auto local_values = std::unique_ptr<LocalEquations>(new LocalEquations{});

    // calculate how many nodes there are
    std::size_t compartment_count = 0;
    std::size_t membrane_count = 0;
    PDELab::forEachLeafNode(
      lbasis.tree(),
      overload([&](const Concept::CompartmentScalarLocalBasisNode auto&) { ++compartment_count; },
               [&](const Concept::MembraneScalarLocalBasisNode auto&) { ++membrane_count; }));
    local_values->_values.reserve(compartment_count + membrane_count);
    local_values->_gradients.reserve(compartment_count + membrane_count);

    local_values->initialize_nodes(
      lbasis.tree(), [&](auto path) { return lbasis.globalBasis().subSpace(path).name(); });

    // check that names are not repeated
    std::set<std::string_view> names;
    PDELab::forEach(local_values->_nodes, [&](auto& compartments_fncs) {
      for (const auto& compartment_fncs : compartments_fncs)
        for (const auto& component_fncs : compartment_fncs)
          if (auto [it, inserted] = names.insert(component_fncs.name); not inserted)
            throw format_exception(
              InvalidStateException{}, "\tVariable with name '{}' is repeated", *it);
    });

    if (opts.any()) {
      if (not functor_factory) {
        throw format_exception(InvalidStateException{}, "Equations cannot be configured without a functor factory");
      }
      local_values->configure(config, *functor_factory, opts);
    }
    return local_values;
  }

  static std::unique_ptr<LocalEquations> make_stiffness(
    const PDELab::Concept::LocalBasis auto& lbasis,
    const ParameterTree& config,
    const FunctorFactory<dim>& functor_factory)
  {
    return make(lbasis,
                config,
                &functor_factory,
                FactoryFalgs::Reaction | FactoryFalgs::Diffusion | FactoryFalgs::Velocity |
                  FactoryFalgs::Outflow);
  }

  static std::unique_ptr<LocalEquations> make_mass(const PDELab::Concept::LocalBasis auto& lbasis,
                                                   const ParameterTree& config,
                                                   const FunctorFactory<dim>& functor_factory)
  {
    return make(lbasis, config, &functor_factory, FactoryFalgs::Storage);
  }

  template<class Tree>
    requires Concept::CompartmentScalarLocalBasisNode<Tree> ||
             Concept::MembraneScalarLocalBasisNode<Tree>
  const auto& get_equation(const Tree& tree) const
  {
    return PDELab::containerEntry(_nodes, make_path(tree));
  }

  template<class Tree>
    requires Concept::CompartmentScalarLocalBasisNode<Tree> ||
             Concept::MembraneScalarLocalBasisNode<Tree>
  auto& get_value(const Tree& tree)
  {
    return PDELab::containerEntry(_nodes, make_path(tree)).value;
  }

  template<class Tree>
    requires Concept::CompartmentScalarLocalBasisNode<Tree> ||
             Concept::MembraneScalarLocalBasisNode<Tree>
  auto& get_gradient(const Tree& tree)
  {
    return PDELab::containerEntry(_nodes, make_path(tree)).gradient;
  }

  void clear()
  {
    std::fill(begin(_values), end(_values), 0.);
    std::fill(begin(_gradients), end(_gradients), 0.);
  }

  void debug() const
  {
    PDELab::forEach(_nodes, [&](auto& compartments) {
      for (auto& compartment_fncs : compartments)
        for (auto& func : compartment_fncs)
          func.debug();
    });
  }

  const auto& nodes() const { return _nodes; }

private:
  template<class Tree>
    requires Concept::CompartmentScalarLocalBasisNode<Tree> ||
             Concept::MembraneScalarLocalBasisNode<Tree>
  auto make_path(const Tree& tree) const
  {
    // the tree may be one of many kind of trees

    // membrane/compartment is distinsuihed because the first index is a template parameter
    // and because they have a different concept requirements
    auto compartment_type = compartment_index(tree);

    // remove compartment/membrane index (always comes from a tuple and is compile time constant)
    auto mi = [&]() {
      auto front_index = front(tree.orderingViewPath());
      // in some cases, the tree has both comartment types in the tree, then the first index of its
      // path is an integral constant
      if constexpr (IsCompileTimeConstant<decltype(front_index)>{}) {
        assert(front_index == compartment_type);
        return pop_front(tree.orderingViewPath());
      } else {
        // otherwise, the tree has only one type of compartments
        return tree.orderingViewPath();
      }
    }();
    assert(mi.size() <= 2);
    // all the compartment trees have their field id as their last path in the ordering view path
    // (different from tree path if membrane)
    auto field_id = back(mi);
    // in case of different compartments, the first index (after stripping its compartment type) is
    // the compartment id
    auto compartment_id = (mi.size() == 2) ? front(mi) : 0;
    return TypeTree::treePath(compartment_type, compartment_id, field_id);
  }

  LocalEquations() = default;

  LocalEquations(const LocalEquations&) = delete;
  LocalEquations(LocalEquations&&) = delete;

  LocalEquations& operator=(const LocalEquations&) = delete;
  LocalEquations& operator=(LocalEquations&&) = delete;

  // (bulk|membrane, compartment, component)
  TupleVector<std::vector<std::vector<CompartmentNode>>, std::vector<std::vector<MembraneNode>>>
    _nodes;
  std::array<std::vector<std::string>, 2> _compartment_names;

  std::vector<Scalar> _values;
  std::vector<Vector> _gradients;

  auto compartment_index(const Concept::CompartmentLocalBasisTree auto&) const
  {
    return Indices::_0;
  }

  auto compartment_index(const Concept::MembraneLocalBasisTree auto&) const { return Indices::_1; }

  template<class Tree>
    requires Concept::CompartmentLocalBasisNode<Tree> || Concept::MembraneLocalBasisNode<Tree>
  void initialize_nodes(const Tree& tree, auto fname, std::size_t c = 0)
  {
    auto ct = compartment_index(tree);
    _nodes[ct].resize(c + 1);
    _nodes[ct][c].clear();
    _compartment_names[ct].resize(c + 1);
    for (std::size_t i = 0; i != tree.degree(); ++i) {
      _nodes[ct][c].emplace_back(_values.emplace_back(),
                                 _gradients.emplace_back(),
                                 TypeTree::treePath(ct, c, i),
                                 fname(tree.child(i).orderingViewPath()));
      _compartment_names[ct][c] = fname(pop_back(tree.child(i).orderingViewPath()));
    }
  }

  template<class Tree>
    requires Concept::MultiCompartmentLocalBasisNode<Tree> ||
             Concept::MultiMembraneLocalBasisNode<Tree>
  void initialize_nodes(const Tree& tree, auto fname)
  {
    for (std::size_t c = 0; c != tree.degree(); ++c)
      initialize_nodes(tree.child(c), fname, c);
  }

  void initialize_nodes(const Concept::MultiCompartmentMembraneLocalBasisNode auto& tree,
                        auto fname)
  {
    using namespace Indices;
    initialize_nodes(tree.child(_0), fname);
    initialize_nodes(tree.child(_1), fname);
  }

  void initialize_nodes(const Concept::CompartmentMembraneLocalBasisNode auto& tree, auto fname)
  {
    using namespace Indices;
    initialize_nodes(tree.child(_0), fname, 0);
    initialize_nodes(tree.child(_1), fname, 0);
  }

  void configure(const ParameterTree& config,
                 const FunctorFactory<dim>& functor_factory,
                 BitFlags<FactoryFalgs> opts)
  {

    using ScalarTag = index_constant<0>;
    using VectorTag = index_constant<1>;
    using TensorApplyTag = index_constant<2>;

    auto make_functor = overload(
      [&](std::string_view prefix, const ParameterTree& config, bool membrane_expression, ScalarTag)
        -> fu2::unique_function<Scalar() const DUNE_COPASI_FUNCTOR_NOEXCEPT> {
        return functor_factory.make_scalar(
          prefix, config, std::as_const(*this), membrane_expression);
      },
      [&](std::string_view prefix, const ParameterTree& config, bool membrane_expression, VectorTag)
        -> fu2::unique_function<Vector() const DUNE_COPASI_FUNCTOR_NOEXCEPT> {
        return functor_factory.make_vector(
          prefix, config, std::as_const(*this), membrane_expression);
      },
      [&](std::string_view prefix,
          const ParameterTree& config,
          bool membrane_expression,
          TensorApplyTag)
        -> fu2::unique_function<Vector(Vector) const DUNE_COPASI_FUNCTOR_NOEXCEPT> {
        return functor_factory.make_tensor_apply(
          prefix, config, std::as_const(*this), membrane_expression);
      });

    auto set_differentiable_function = overload(
      [&]<class Signature>(CompartmentDifferentiableFunction<Signature>& function,
                           std::string_view prefix,
                           const ParameterTree& function_config,
                           auto range_tag) {
        function = CompartmentDifferentiableFunction<Signature>{ make_functor(
          prefix, function_config, false, range_tag) };
        std::set<std::string_view> debug_set;
        const auto& jac_config = function_config.sub("jacobian");
        for (const auto& compartment_fncs_jac : _nodes[Indices::_0])
          for (const auto& component_fncs_jac : compartment_fncs_jac)
            if (jac_config.hasSub(component_fncs_jac.name))
              if (auto jac =
                    make_functor(fmt::format("{}.jacobian.{}", prefix, component_fncs_jac.name),
                                 jac_config.sub(component_fncs_jac.name),
                                 false,
                                 range_tag)) {
                function.compartment_jacobian.emplace_back(std::move(jac), component_fncs_jac);
                debug_set.insert(component_fncs_jac.name);
              }
        if (jac_config.getSubKeys().size() != debug_set.size())
          spdlog::warn("Some sub-sections in \"jacobian\" section are being ignored");
      },
      [&]<class Signature>(MembraneDifferentiableFunction<Signature>& function,
                           std::string_view prefix,
                           const ParameterTree& function_config,
                           auto range_tag) {
        function = MembraneDifferentiableFunction<Signature>{ make_functor(
          prefix, function_config, true, range_tag) };
        const auto& jac_config = function_config.sub("jacobian");
        for (const auto& compartment_fncs_jac : _nodes[Indices::_0])
          for (const auto& component_fncs_jac : compartment_fncs_jac)
            if (jac_config.hasSub(component_fncs_jac.name))
              if (auto jac =
                    make_functor(fmt::format("{}.jacobian.{}", prefix, component_fncs_jac.name),
                                 jac_config.sub(component_fncs_jac.name),
                                 true,
                                 range_tag))
                function.compartment_jacobian.emplace_back(std::move(jac), component_fncs_jac);
        for (const auto& membrane_fncs_jac : _nodes[Indices::_1])
          for (const auto& component_fncs_jac : membrane_fncs_jac)
            if (jac_config.hasSub(component_fncs_jac.name))
              if (auto jac =
                    make_functor(fmt::format("{}.jacobian.{}", prefix, component_fncs_jac.name),
                                 jac_config.sub(component_fncs_jac.name),
                                 true,
                                 range_tag))
                function.membrane_jacobian.emplace_back(std::move(jac), component_fncs_jac);
      });

    PDELab::forEach(
      _nodes, [&]<class Node>(std::vector<std::vector<Node>>& compartments_fncs, auto l) {
        for (auto& compartment_fncs : compartments_fncs) {
          for (Node& component_fncs : compartment_fncs) {

            const auto& component_config = config.sub(component_fncs.name, true);

            if (opts.test(FactoryFalgs::Diffusion) and component_config.hasSub("cross_diffusion")) {
              const auto& cross_diffusion_config = component_config.sub("cross_diffusion");
              std::set<std::string_view> debug_set;
              for (const auto& compartment_fncs_cross_diff : compartments_fncs)
                for (const auto& component_fncs_cross_diff : compartment_fncs_cross_diff)
                  if (cross_diffusion_config.hasSub(component_fncs_cross_diff.name)) {
                    const auto& diffusion_config =
                      cross_diffusion_config.sub(component_fncs_cross_diff.name);
                    debug_set.insert(component_fncs_cross_diff.name);
                    auto& cross_diffusion = component_fncs.cross_diffusion.emplace_back(
                      nullptr, component_fncs_cross_diff);
                    set_differentiable_function(
                      cross_diffusion,
                      fmt::format("{}.cross_diffusion", component_fncs.name),
                      diffusion_config,
                      TensorApplyTag{});
                    if (not cross_diffusion)
                      component_fncs.cross_diffusion.pop_back();
                  }
              if (cross_diffusion_config.getSubKeys().size() != debug_set.size())
                spdlog::warn(
                  "Some sub-sections in \"{}.cross_diffusion\" section are being ignored",
                  component_fncs.name);
            }

            if (opts.test(FactoryFalgs::Reaction) and component_config.hasSub("reaction"))
              set_differentiable_function(component_fncs.reaction,
                                          fmt::format("{}.reaction", component_fncs.name),
                                          component_config.sub("reaction"),
                                          ScalarTag{});

            if (opts.test(FactoryFalgs::Velocity) and component_config.hasSub("velocity"))
              set_differentiable_function(component_fncs.velocity,
                                          fmt::format("{}.velocity", component_fncs.name),
                                          component_config.sub("velocity"),
                                          VectorTag{});

            if (opts.test(FactoryFalgs::Storage) and component_config.hasSub("storage"))
              set_differentiable_function(component_fncs.storage,
                                          fmt::format("{}.storage", component_fncs.name),
                                          component_config.sub("storage"),
                                          ScalarTag{});

            if (opts.test(FactoryFalgs::Outflow) and component_config.hasSub("outflow")) {
              if (l == 1)
                throw format_exception(NotImplemented{},
                                       "Outflow for membranes is not implemented");
              const auto& boundaries_config = component_config.sub("outflow");
              for (std::size_t i = 0; i != _compartment_names[l].size(); ++i) {
                if (boundaries_config.hasSub(_compartment_names[l][i])) {
                  component_fncs.outflow.resize(_compartment_names[l].size());
                  set_differentiable_function(
                    component_fncs.outflow[i],
                    fmt::format("{}.storage.{}", component_fncs.name, _compartment_names[l][i]),
                    boundaries_config.sub(_compartment_names[l][i]),
                    ScalarTag{});
                }
              }
            }
          }
        }
      });
  }
};

} // namespace Dune::Copasi

#endif // DUNE_COPASI_MODEL_LOCAL_EQUATIONS_LOCAL_EQUATIONS_HH
