#ifndef DUNE_COPASI_MODEL_DIFFUSION_REACTION_HH
#define DUNE_COPASI_MODEL_DIFFUSION_REACTION_HH

#include <dune/copasi/common/enum.hh>
#include <dune/copasi/concepts/grid.hh>

#include <dune/copasi/local_operator/diffusion_reaction/continuous_galerkin.hh>
#include <dune/copasi/local_operator/diffusion_reaction/finite_volume.hh>
#include <dune/copasi/local_operator/diffusion_reaction/multidomain.hh>
#include <dune/copasi/local_operator/variadic.hh>
#include <dune/copasi/model/state.hh>

#include <dune/pdelab/backend/istl.hh>
#include <dune/pdelab/constraints/conforming.hh>
#include <dune/pdelab/constraints/p0.hh>
#include <dune/pdelab/function/discretegridviewfunction.hh>
#include <dune/pdelab/gridfunctionspace/dynamicpowergridfunctionspace.hh>
#include <dune/pdelab/gridfunctionspace/gridfunctionspace.hh>
#include <dune/pdelab/gridfunctionspace/subspace.hh>
#include <dune/pdelab/gridfunctionspace/vtk.hh>
#include <dune/pdelab/gridoperator/gridoperator.hh>
#include <dune/pdelab/gridoperator/onestep.hh>
#include <dune/pdelab/finiteelementmap/p0fem.hh>
#include <dune/pdelab/finiteelementmap/pkfem.hh>

#include <dune/grid/io/file/vtk/vtksequencewriter.hh>
#include <dune/grid/io/file/vtk/vtkwriter.hh>

#include <dune/common/indices.hh>
#include <dune/common/parametertree.hh>
#include <dune/common/parametertreeparser.hh>

#include <memory>
#include <sstream>
#include <string>
#include <type_traits>
#include <vector>

namespace Dune::Copasi {

// a base class to create multi-domain multi-component models
class PkMultiDomainMultiComponentStateFactory
{
  using RF = double;

public:
  template<class GridView, class FEMOrder>
  static auto make_finite_element_map(const GridView& grid_view,
                                      FEMOrder fem_order)
  {
    static_assert(fem_order != 0,
                  "Conforming galerkin implementation does not allow 0-th "
                  "order polynomials");
    using FEMk = PDELab::PkLocalFiniteElementMap<GridView, RF, RF, fem_order>;
    return std::make_unique<FEMk>(grid_view);
  }

  template<class EntitySet, class FEMOrder>
  static auto make_component_grid_function_space(
    const std::string& component_name,
    const EntitySet& entity_set,
    FEMOrder fem_order)
  {
    auto fem = make_finite_element_map(entity_set.gridView(), fem_order);
    auto first_entity_it = entity_set.template begin<0>();
    auto is_fv = first_entity_it->type().isCube();

    using FEM = std::decay_t<decltype(*fem)>;

    //! Constraints local operator
    using CON = PDELab::ConformingDirichletConstraints;

    //! Leaf vector backend
    using LVBE = PDELab::ISTL::BackendOptions<>;

    //! Leaf grid function space
    using LGFS = PDELab::UnorderedGridFunctionSpace<EntitySet, FEM, CON, LVBE>;

    auto comp_gfs = std::make_shared<LGFS>(entity_set, std::move(fem));
    comp_gfs->name(component_name);

    using OutputType = PDELab::GridFunctionOutputParameters::Output;
    auto output_type = is_fv ? OutputType::cellData : OutputType::vertexData;
    comp_gfs->setDataSetType(output_type);
    return comp_gfs;
  }

  template<class EntitySet, class FEMOrder, class Ordering>
  static auto make_multicomponent_grid_function_space(
    const std::string& domain_name,
    const std::vector<std::string>& component_names,
    const EntitySet& entity_set,
    FEMOrder fem_order,
    Ordering ordering)
  {
    auto get_gfs = [&](const std::string& name) {
      return make_component_grid_function_space( name, entity_set, fem_order);
    };

    using LGFS = std::decay_t<decltype(*get_gfs(""))>;

    //! Vector backend
    using VBE = PDELab::ISTL::BackendOptions<>;

    //! Grid function space
    using GFS =
      PDELab::UnorderedDynamicPowerGridFunctionSpace<LGFS, VBE, Ordering>;

    auto comp_gfs_vec = typename GFS::NodeStorage{};

    for (const auto& name : component_names)
      comp_gfs_vec.push_back(get_gfs(name));

    auto gfs = std::make_shared<GFS>(comp_gfs_vec);
    gfs->name(domain_name);

    if (gfs->degree() == 0)
      DUNE_THROW(InvalidStateException,
                 "Grid function space is not correctly setup");
    return gfs;
  }

  template<class MultiDomainGrid,
           class FEMOrder,
           class SDOrderingTag>
  static auto make_multidomain_multicomponent_grid_function_space(
    std::vector<std::string> domain_names,
    const std::vector<std::vector<std::string>>& component_names,
    const MultiDomainGrid& md_grid,
    FEMOrder fem_order,
    SDOrderingTag sd_ordering)
  {
    assert(domain_names.size() == component_names.size());

    auto get_gfs = [&](const std::size_t& sub_domain) {
      auto sd_grid_view = md_grid.subDomain(sub_domain).leafGridView();
      using SubDomainGridView =
        typename MultiDomainGrid::SubDomainGrid::LeafGridView;
      using EntitySet = PDELab::NonOverlappingEntitySet<SubDomainGridView>;

      return make_multicomponent_grid_function_space(
        domain_names[sub_domain],
        component_names[sub_domain],
        EntitySet{ sd_grid_view },
        fem_order,
        sd_ordering);
    };

    using SDGFS = std::decay_t<decltype(*get_gfs(0))>;
    using VBE = PDELab::ISTL::BackendOptions<>;
    using MDOrderingTag = PDELab::LexicographicOrderingTag;
    using GFS =
      PDELab::UnorderedDynamicPowerGridFunctionSpace<SDGFS, VBE, MDOrderingTag>;

    auto gfs_vec = typename GFS::NodeStorage{};

    if (md_grid.maxAssignedSubDomainIndex() >= domain_names.size())
      DUNE_THROW(InvalidStateException,
                 "Domain names size do not match multidomain size");
    gfs_vec.resize(domain_names.size());

    for (std::size_t sub_domain = 0; sub_domain < domain_names.size();
         ++sub_domain)
      gfs_vec[sub_domain] = get_gfs(sub_domain);

    return std::make_shared<GFS>(gfs_vec, VBE{}, MDOrderingTag{});
  }

  template<class Grid,
           class FEMOrder = index_constant<1>,
           class Ordering = PDELab::EntityBlockedOrderingTag>
  static auto make_singledomain_multicomponent_state(
    std::shared_ptr<Grid> grid,
    const std::string& domain_name,
    const std::vector<std::string>& component_names,
    FEMOrder fem_order = {},
    Ordering ordering = {})
  {
    if (component_names.size() == 0)
      DUNE_THROW(InvalidStateException, "Empty compartments are not supported");
    using GridView = typename Grid::LeafGridView;
    using EntitySet = PDELab::OverlappingEntitySet<GridView>;
    EntitySet entity_set{ grid->leafGridView() };

    auto ugfs = make_multicomponent_grid_function_space(domain_name,
                                                        component_names,
                                                        entity_set,
                                                        fem_order,
                                                        ordering);
    using UGFS = std::decay_t<decltype(*ugfs)>;
    using GFS = Dune::PDELab::OrderedGridFunctionSpace<UGFS, EntitySet>;
    auto gfs = std::make_shared<GFS>(*ugfs, entity_set);
    using X = PDELab::Backend::Vector<GFS, RF>;
    auto x = std::make_shared<X>(std::shared_ptr<const GFS>{ gfs }, 0.);
    return PDELabState{ grid, gfs, x };
  }

private:
  template<class GFS, class TreePath>
  static auto subSpace(GFS& gfs, TreePath tree_path)
  {
    if constexpr (TreePath::size() > 0) {
      const auto& child_gfs = TypeTree::child(gfs, tree_path);
      auto enity_set = PDELab::Impl::first_leaf(child_gfs).entitySet();
      return PDELab::makeSubSpace(gfs, tree_path, enity_set);
    } else {
      return gfs;
    }
  }

public:
  template<class TreePath>
  static auto make_vtk_writer(TreePath tree_path)
  {
    auto log_writer = Logging::Logging::componentLogger({}, "writer");
    log_writer.debug("Setup VTK writer"_fmt);

    std::vector<double> timesteps;
    return [=](const auto& state, const fs::path& path, bool append) mutable {
      // create directory if necessary
      auto path_entry = fs::directory_entry{ path };
      if (not path_entry.exists()) {
        log_writer.info("Creating output directory '{}'"_fmt,
                        path_entry.path().string());
        std::error_code ec{ 0, std::generic_category() };
        fs::create_directories(path_entry.path(), ec);
        if (ec)
          DUNE_THROW(IOError,
                     "\n Category: " << ec.category().name() << '\n'
                                     << "Value: " << ec.value() << '\n'
                                     << "Message: " << ec.message() << '\n');
      }
      // setup writer again with old timesteps if necessary
      auto vtk_gfs = subSpace(std::as_const(state).space(), tree_path);
      using VTKGFS = std::decay_t<decltype(vtk_gfs)>;
      using VTKGridView = typename VTKGFS::Traits::GridView;

      std::string name = path.filename().string();
      if (not append) {
        timesteps.clear();
        log_writer.detail("Creating a time sequence file: '{}.pvd'"_fmt, name);
      } else {
        log_writer.trace("Overriding time sequence file: '{}.pvd'"_fmt, name);
      }

      auto writer = std::make_shared<VTKWriter<VTKGridView>>(
        vtk_gfs.gridView(), Dune::VTK::conforming);
      auto sequential_writer = VTKSequenceWriter<VTKGridView>{
        writer, name, path.string(), path.string()
      };
      Dune::PDELab::addSolutionToVTKWriter(
        sequential_writer, vtk_gfs, state.coefficients());
      sequential_writer.setTimeSteps(timesteps);

      log_writer.detail("Writing solution for {:.2f}s time stamp"_fmt,
                        state.time);
      log_writer.trace(
        "Writing vtu file: '{0}/{0}-{1:0>5}.vtu'"_fmt, name, timesteps.size());
      sequential_writer.write(state.time, Dune::VTK::base64);
      sequential_writer.vtkWriter()->clear();
      timesteps = sequential_writer.getTimeSteps();
    };
  }

  static auto make_vtk_writer()
  {
    return make_vtk_writer(TypeTree::treePath());
  }

  template<class State, class GridFunctionTree, class TreePath>
  static void interpolate(State& state,
                          const GridFunctionTree& grid_function_tree,
                          const TreePath& tree_path)
  {
    auto sub_space = subSpace(std::as_const(state).space(), tree_path);
    using namespace TypeTree::Experimental;
    const std::size_t sub_space_depth = Info::depth<decltype(sub_space)>();
    const std::size_t gf_tree_depth = Info::depth<GridFunctionTree>();
    // sub-space and grid function should have the same tree structure
    static_assert(sub_space_depth == gf_tree_depth);
    assert(sub_space.degree() == grid_function_tree.degree());
    // TODO assert grid views are the same!
    PDELab::interpolate(grid_function_tree, sub_space, state.coefficients());
  }

  template<class State, class GridFunctionTree>
  static void interpolate(State& state,
                          const GridFunctionTree& grid_function_tree)
  {
    interpolate(state, grid_function_tree, TypeTree::treePath());
  }

  // this builds recursively a grid function tree starting with root node on
  // tree_path
  template<class DataHandler, class TreePath>
  static auto as_grid_function(const std::shared_ptr<DataHandler>& data_handler,
                               const TreePath& tree_path)
  {
    const auto& child_lfs =
      TypeTree::child(data_handler->localFunctionSpace(), tree_path);
    using ChildLFS = std::decay_t<decltype(child_lfs)>;
    if constexpr (ChildLFS::isLeaf) {
      using GridView =
        typename ChildLFS::Traits::GridFunctionSpace::Traits::GridView;
      using DiscreteLeafFunction =
        PDELab::vtk::DGFTreeLeafFunction<ChildLFS, DataHandler, GridView>;
      return std::make_shared<DiscreteLeafFunction>(child_lfs, data_handler);
    } else if constexpr (ChildLFS::isPower) {
      using SubGridFunctionPtr = std::decay_t<decltype(as_grid_function(
        data_handler, push_back(tree_path, std::size_t{})))>;
      std::vector<SubGridFunctionPtr> sub_grid_functions;
      for (std::size_t i = 0; i < child_lfs.degree(); ++i) {
        // here we call recursively this function to fill the dynamic node.
        // recursion ends with leaf nodes.
        sub_grid_functions.push_back(
          as_grid_function(data_handler, push_back(tree_path, i)));
      }
      using SubGridFunction = typename SubGridFunctionPtr::element_type;
      // TODO assert lfs is dynamic power
      using GridFunctionTree =
        PDELab::DynamicPowerGridFunction<SubGridFunction>;
      return std::make_shared<GridFunctionTree>(sub_grid_functions);
    } else {
      static_assert(AlwaysFalse<DataHandler>{},
                    "Composite nodes are not implemented!");
    }
  }

  // returns a grid function tree of the same shape as the grid function space
  // rooted ath the tree_path. in case of multidomain states, the grid view will
  // be from the sub domain in case of non-empty paths.
  template<class State, class TreePath>
  static auto as_grid_function(const State& state, const TreePath& tree_path)
  {
    auto sub_gfs = subSpace(state.space(), tree_path);
    using GFS = std::decay_t<decltype(sub_gfs)>;
    using X = typename State::Coefficients;
    using DataHandler = PDELab::vtk::DGFTreeCommonData<GFS, X>;

    const auto& x = state.coefficients_storage();
    auto data_handler =
      std::make_shared<DataHandler>(std::make_shared<GFS>(sub_gfs), x);

    // This returns a shared pointer, but it is only needed for recursive
    // construction. Thus, we just return a copy of the grid function tree.
    return *as_grid_function(data_handler, TypeTree::treePath());
  }

  // The grid function is only usable with the assembly grid view. In other
  // words, In case of multidomain state, the resulting grid function will only
  // be usable by the multi-domain grid and not by the sub-domain grid nor the
  // host grid.
  template<class State>
  static auto as_grid_function(const State& state)
  {
    return as_grid_function(state, TypeTree::treePath());
  }
};

/**
 * @brief      Class for (multidomain) diffusion-reaction models.
 */
class DiffusionReactionModel : public PkMultiDomainMultiComponentStateFactory
{
  using RF = double;

  using Base = PkMultiDomainMultiComponentStateFactory;

public:
  DiffusionReactionModel(const std::string& str_config)
  {
    std::istringstream ss{ str_config };
    ParameterTreeParser::readINITree(ss, _config);
  }

  DiffusionReactionModel(const ParameterTree& config)
    : _config(config)
  {
    const auto& domain_names = _config.getSubKeys();
    for (const auto& domain_name : domain_names) {
      const auto& domain_config = _config.sub(domain_name, true);
      const auto& diff_vars = domain_config.sub("diffusion").getValueKeys();
      const auto& react_vars = domain_config.sub("reaction").getValueKeys();
      if (react_vars != diff_vars)
        DUNE_THROW(InvalidStateException,
                   "'reaction' and 'diffusion' sub-sections in '"
                     << domain_name << "' section should have the same keys");
    }
  }

  template<class Grid,
           class FEMOrder = index_constant<1>,
           class Ordering = PDELab::EntityBlockedOrderingTag>
  auto make_singledomain_multicomponent_state(std::shared_ptr<Grid> grid,
                                              const std::string& domain_name,
                                              FEMOrder fem_order = {},
                                              Ordering ordering = {}) const
  {
    const auto& vars =
      _config.sub(domain_name, true).sub("diffusion").getValueKeys();
    return Base::make_singledomain_multicomponent_state(
      grid, domain_name, vars, fem_order, ordering);
  }

  template<class Grid,
           class FEMOrder = index_constant<1>,
           class SDOrderingTag = PDELab::EntityBlockedOrderingTag>
  auto make_multidomain_multicomponent_state(
    std::shared_ptr<Grid> grid,
    FEMOrder fem_order = {},
    SDOrderingTag sd_ordering = {}) const
  {
    //! TODO: assert multidomain grid
    std::vector<std::vector<std::string>> component_names;
    auto domain_names = _config.sub("compartments",true).getValueKeys();
    sort(begin(domain_names), end(domain_names));
    for (const auto& domain_name : domain_names) {
      const auto& domain_config = _config.sub(domain_name);
      component_names.push_back(domain_config.sub("diffusion").getValueKeys());
      sort(begin(component_names), end(component_names));
    }

    auto ugfs =
      Base::make_multidomain_multicomponent_grid_function_space(domain_names,
                                                                component_names,
                                                                *grid,
                                                                fem_order,
                                                                sd_ordering);

    using GridView = typename Grid::LeafGridView;
    using EntitySet = PDELab::OverlappingEntitySet<GridView>;
    EntitySet entity_set{ grid->leafGridView() };

    using UGFS = std::decay_t<decltype(*ugfs)>;
    using GFS = Dune::PDELab::OrderedGridFunctionSpace<UGFS, EntitySet>;
    auto gfs = std::make_shared<GFS>(*ugfs, entity_set);

    using X = PDELab::Backend::Vector<GFS, RF>;
    auto x = std::make_shared<X>(std::shared_ptr<const GFS>{ gfs }, 0.);
    return PDELabState{ grid, gfs, x };
  }

public:
  struct FV_CG_Switch //! TODO: fix PDELab to allow to pass entities in
                      //! patterns, then use fem mapper
  {
    template<class T>
    std::size_t operator()(const T& fe_v) const
    {
      return fe_v.type().isCube() ? 0 : 1;
    }
    template<class T0, class T1>
    std::size_t operator()(const T0& fe_u, const T1& fe_v) const
    {
      return fe_v.type().isCube() ? 0 : 1;
    }
  };

  template<class MultiComponentGFS>
  auto make_spatial_multicomponent_local_operator(
    const MultiComponentGFS& mc_gfs) const
  {
    const std::size_t depth =
      TypeTree::Experimental::Info::depth<MultiComponentGFS>();
    static_assert(depth == 2, "MultiComponentGFS has incorrect depth");

    using ComponentGFS =
      typename MultiComponentGFS::ChildType; // (all leaf nodes have the same
                                             // type)
    using GridView = typename ComponentGFS::Traits::GridView;
    using FEM = typename ComponentGFS::Traits::FiniteElementMap;
    using BasisTraits =
      typename FEM::Traits::FiniteElement::Traits::LocalBasisType::Traits;

    GridView grid_view =
      mc_gfs.child(0).gridView(); // all components are on the same (sub-)domain
    std::string compartment = mc_gfs.name();
    const auto& sub_config = _config.sub(compartment);

    // multi-component conforming galerkin local operator
    using MCCGLOP =
      LocalOperatorDiffusionReactionCG<GridView,
                                       BasisTraits,
                                       JacobianMethod::Numerical>;
    return std::make_unique<MCCGLOP>(grid_view, sub_config);
  }

  template<class GridFunctionSpace>
  auto make_spatial_local_operator(const GridFunctionSpace& gfs) const
  {
    using Grid = typename GridFunctionSpace::Traits::GridView::Grid;
    const std::size_t depth =
      TypeTree::Experimental::Info::depth<GridFunctionSpace>();
    if constexpr (depth == 2) {
      return make_spatial_multicomponent_local_operator(gfs);
    } else {
      static_assert(depth == 3,
                    "Depth of grid function space shall be either 2 or 3");
      using MultiComponentLocalOperator =
        std::decay_t<decltype(*make_spatial_multicomponent_local_operator(
          gfs.child(0)))>;
      using LOP =
        LocalOperatorMultiDomainDiffusionReaction<Grid,
                                                  MultiComponentLocalOperator,
                                                  JacobianMethod::Numerical>;
      return std::make_unique<LOP>(gfs.gridView().grid(), _config);
    }
  }

  template<class GridFunctionSpace>
  auto make_temporal_local_operator(const GridFunctionSpace& gfs) const
  {

    using Grid = typename GridFunctionSpace::Traits::GridView::Grid;
    using ComponentGFS =
      PDELab::Impl::FirstLeaf<GridFunctionSpace>; // (all leaf nodes have the
                                                  // same type)
    using GridView = typename ComponentGFS::Traits::GridView;
    using FEM = typename ComponentGFS::Traits::FiniteElementMap;
    using BasisTraits =
      typename FEM::Traits::FiniteElement::Traits::LocalBasisType::Traits;
    using MultiComponentTemporalLocalOperator =
      TemporalLocalOperatorDiffusionReactionCG<GridView,
                                               BasisTraits,
                                               JacobianMethod::Numerical>;

    const std::size_t depth =
      TypeTree::Experimental::Info::depth<GridFunctionSpace>();

    if constexpr (depth == 2) {
      GridView grid_view =
        gfs.child(0).gridView(); // all components are on the same (sub-)domain
      const auto& sub_config = _config.sub(gfs.name());
      return std::make_unique<MultiComponentTemporalLocalOperator>(grid_view,
                                                                   sub_config);
    } else {
      static_assert(depth == 3,
                    "Depth of grid function space shall be either 2 or 3");
      using TemporalLocalOperator =
        TemporalLocalOperatorMultiDomainDiffusionReaction<
          Grid,
          MultiComponentTemporalLocalOperator,
          JacobianMethod::Numerical>;
      return std::make_unique<TemporalLocalOperator>(gfs.gridView().grid(),
                                                     _config);
    }
  }

private:
  ParameterTree _config;
};

} // namespace Dune::Copasi

#endif // DUNE_COPASI_MODEL_DIFFUSION_REACTION_HH
