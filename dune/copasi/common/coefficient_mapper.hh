#ifndef DUNE_COPASI_GLOBAL_COEFFICIENT_MAPPER_HH
#define DUNE_COPASI_GLOBAL_COEFFICIENT_MAPPER_HH

#include <dune/copasi/concepts/grid.hh>

#include <dune/pdelab/gridfunctionspace/lfsindexcache.hh>
#include <dune/pdelab/gridfunctionspace/localfunctionspace.hh>

#include <dune/common/parametertree.hh>

#include <map>
#include <memory>
#include <vector>

namespace Dune::Copasi {

/**
 * @brief      Default coefficinet mapper
 * @details    Does nothing but evaluate the second an the third arguments of
 *             the parenthesis call into the first one. Useful to keep code
 *             compatible with other coefficient mappers
 */
struct DefaultCoefficientMapper
{
  /**
   * @brief      Constructor. Arguments in this class do nothimg.
   *
   * @param[in]  operator_config  The operator configuration file
   * @param[in]  this_op          The operator for local coefficients
   */
  DefaultCoefficientMapper(const ParameterTree& operator_config,
                           std::size_t this_op)
  {}

  /**
   * @brief      Evaluates the coefficents for a give component and a degree of
   *             freedom. Either from a local vector x_view_local, or by other
   *             stored coefficients vector.
   *
   * @param[in]  x_view_local  The view on the local coefficients vector
   * @param[in]  comp          The component
   * @param[in]  dof           The degree of freedom
   *
   * @tparam     XViewLocal    The type fo the local coefficient view
   *
   * @return     The coefficient value for the requested component and the
   *             degree of freedom
   */
  template<class XViewLocal>
  auto operator()(const XViewLocal& x_view_local,
                  const std::size_t& comp,
                  const std::size_t& dof) const
  {
    return x_view_local(comp, dof);
  }

  /**
   * @brief      Bind internal coefficients to an entity
   *
   * @param[in]  <unnamed>  The entity
   *
   * @tparam     E          The entity
   */
  template<class E>
  void bind(const E&)
  {}

  /**
   * @brief      Unbind this coefficients to an entity
   */
  void unbind() {}

  /**
   * @brief      Updates the internal coefficients with a exteral source.
   *
   * @param[in]  <unnamed>  Map to the model states
   *
   * @tparam     T          Map to the model states
   */
  template<class T>
  void update(const T&)
  {}
};

/**
 * @brief      Model coefficinet mapper base
 * @details    Keeps a reference to external coefficients so that they are
 *             accessible within other local operators.
 *
 * @tparam     ModelState  The model state type
 */
template<class ModelState>
struct ModelCoefficientMapperBase
{
  const std::size_t _this_op;

  using GFS = typename ModelState::GridFunctionSpace;
  using LFS = PDELab::LocalFunctionSpace<GFS>;
  using LFSCache = PDELab::LFSIndexCache<LFS>;
  using X = typename ModelState::Coefficients;
  using XView = typename X::template ConstLocalView<LFSCache>;
  using SolutionElement = typename X::ElementType;
  using SolutionVector =
    Dune::PDELab::LocalVector<SolutionElement, Dune::PDELab::TrialSpaceTag>;

  std::map<std::size_t, std::shared_ptr<LFS>> _lfs;
  std::map<std::size_t, std::shared_ptr<LFSCache>> _lfs_cache;
  std::map<std::size_t, std::shared_ptr<XView>> _x_view;
  std::map<std::size_t, std::shared_ptr<SolutionVector>> _x;

  /**
   * @brief      Constructs the model coefficient mapper base
   * @details    The operator in the argument refers to the index to which the
   *             coefficient mapper has to use local coefficients (passed as
   *             argument during evaluation) instead of the external
   *             coefficients. The index refers to the index to which the model
   *             states are mapped from.
   *
   * @param[in]  this_op  The operator for local coefficients
   */
  ModelCoefficientMapperBase(std::size_t this_op)
    : _this_op(this_op)
  {}

  /**
   * @brief      Updates the given model states.
   *
   * @param[in]  states  The model states
   */
  void update(const std::map<std::size_t, ModelState>& states)
  {
    _x.clear();
    _x_view.clear();
    _lfs_cache.clear();
    _lfs.clear();

    for (const auto& [op, state] : states) {
      _lfs[op] = std::make_shared<LFS>(state.grid_function_space);
      if (_this_op == op)
        continue;
      _lfs_cache[op] = std::make_shared<LFSCache>(*_lfs[op]);
      _x_view[op] = std::make_shared<XView>(*(state.coefficients));
      auto max_local_size = state.grid_function_space->maxLocalSize();
      _x[op] = std::make_shared<SolutionVector>(max_local_size);
    }
  }

  /**
   * @brief      Bind internal coefficients to an entity
   * @details    This binds the entity to a internal local function space and
   *             local coefficients to the internal coefficients
   *
   * @param[in]  entity  The entity
   *
   * @tparam     E       The entity
   */
  template<class E>
  void bind(const E& entity)
  {
    for (const auto& [op, xx] : _lfs) {
      if (_this_op == op)
        continue;
      _lfs[op]->bind(entity);
      _lfs_cache[op]->update();
      _x[op]->resize(_lfs_cache[op]->size());
      _x_view[op]->bind(*_lfs_cache[op]);
      _x_view[op]->read(*_x[op]);
    }
  }

  /**
   * @brief      Unbind this coefficients to an entity
   */
  void unbind()
  {
    for (auto& [op, x_view] : _x_view)
      x_view.unbind();
  }
};

/**
 * @brief      Constructs the model coefficient mapper
 * @details    This coefficient mapper is designed to one domain local
 *             operators. That is, when the PDELab local function tree only
 *             has depth equal to 1
 *
 * @tparam     ModelState  The model state type
 */
template<class ModelState>
struct ModelCoefficientMapper : public ModelCoefficientMapperBase<ModelState>
{
  using Base = ModelCoefficientMapperBase<ModelState>;
  using Base::_lfs;
  using Base::_this_op;
  using Base::_x;

  // component -> (operator,child)
  std::unordered_map<std::size_t, std::array<std::size_t, 2>> _comp_mapper;

  /**
   * @brief      Constructor of the coefficient mapper
   *
   * @param[in]  this_op  The operator for local coefficients
   */
  ModelCoefficientMapper(std::size_t this_op = 0)
    : Base(this_op)
  {}

  /**
   * @brief      Constructor of the coefficient mapper
   *
   * @param[in]  map_operator  A vector mapping components to operators
   * @param[in]  this_op       The operator for local coefficients
   */
  ModelCoefficientMapper(const std::vector<std::size_t>& map_operator,
                         std::size_t this_op)
    : Base(this_op)
  {
    set_mapper(map_operator);
  }

  /**
   * @brief      Constructor of the coefficient mapper
   *
   * @param[in]  operator_config  The operator configuration
   * @param[in]  this_op          The operator for local coefficients
   */
  ModelCoefficientMapper(const ParameterTree& operator_config,
                         std::size_t this_op)
    : Base(this_op)
  {
    // map component -> operator
    std::vector<std::size_t> map_op;
    auto comp_names = operator_config.getValueKeys();
    std::sort(comp_names.begin(), comp_names.end());

    for (std::size_t j = 0; j < comp_names.size(); j++) {
      std::size_t k = operator_config.template get<std::size_t>(comp_names[j]);
      map_op.push_back(k);
    }

    set_mapper(map_op);
  }

  /**
   * @brief      Sets a map from component to operator and child in the PDELab
   *             tree
   *
   * @param[in]  map_op  A vector mapping components to operators
   */
  void set_mapper(const std::vector<std::size_t>& map_op)
  {
    std::size_t max_op = *std::max_element(map_op.begin(), map_op.end());
    std::vector<std::size_t> comp_child(max_op + 1, 0);

    for (std::size_t i = 0; i < map_op.size(); i++)
      _comp_mapper[i] =
        std::array<std::size_t, 2>{ map_op[i], comp_child[map_op[i]]++ };
  }

  /**
   * @brief      Evaluates the coefficents for a give component and a degree of
   *             freedom. Either from a local vector x_view_local, or by other
   *             stored coefficients vector.
   *
   * @param[in]  x_view_local  The view on the local coefficients vector
   * @param[in]  comp          The component
   * @param[in]  dof           The degree of freedom
   *
   * @tparam     XViewLocal    The type fo the local coefficient view
   *
   * @return     The coefficient value for the requested component and the
   *             degree of freedom
   */
  template<class XViewLocal>
  auto operator()(const XViewLocal& x_view_local,
                  const std::size_t& comp,
                  const std::size_t& dof) const
  {
    const auto& map = _comp_mapper.find(comp)->second;
    if (map[0] == _this_op)
      return x_view_local(map[1], dof);
    else {
      const auto& x_view_global = *_x.at(map[0]);
      const auto& lfs = _lfs.at(map[0])->child(map[1]);
      return x_view_global(lfs, dof);
    }
  }
};

/**
 * @brief      Constructs the coefficient mapper for multidomain models
 * @details    This coefficient mapper is designed to multidomain domain local
 *             operators. That is, when the PDELab local function tree has
 *             depth equal to 2
 *
 * @tparam     ModelState  The model state type
 */
template<class ModelState>
struct MultiDomainModelCoefficientMapper
  : public ModelCoefficientMapper<ModelState>
{
  using Base = ModelCoefficientMapper<ModelState>;
  using Base::_comp_mapper;
  using Base::_lfs;
  using Base::_this_op;
  using Base::_x;
  using Base::Base;
  std::size_t _domain;

  /**
   * @brief      Evaluates the coefficents for a give component and a degree of
   *             freedom. Either from a local vector x_view_local, or by other
   *             stored coefficients vector.
   *
   * @param[in]  x_view_local  The view on the local coefficients vector
   * @param[in]  comp          The component
   * @param[in]  dof           The degree of freedom
   *
   * @tparam     XViewLocal    The type fo the local coefficient view
   *
   * @return     The coefficient value for the requested component and the
   *             degree of freedom
   */
  template<class XViewLocal>
  auto operator()(const XViewLocal& x_view_local,
                  const std::size_t& comp,
                  const std::size_t& dof) const
  {
    const auto& map = _comp_mapper.find(comp)->second;
    if (map[0] == _this_op)
      return x_view_local(map[1], dof);
    else {
      const auto& x_view_global = *_x.at(map[0]);
      const auto& lfs = _lfs.at(map[0])->child(_domain).child(map[1]);
      assert(lfs.size() > 0);
      return x_view_global(lfs, dof);
    }
  }

  /**
   * @brief      Bind internal coefficients to an entity
   * @details    This procedure sets the first child in the PDELab tree by
   *             finding the domain the entity belongs to, plus further bindings
   *             for the other coefficients
   *
   * @param[in]  entity  The entity
   *
   * @tparam     E       The entity
   */
  template<class E>
  void bind(const E& entity)
  {
    Base::bind(entity);
    const auto gv = _lfs[_this_op]->gridFunctionSpace().gridView();
    using Grid = std::decay_t<decltype(gv.grid())>;
    static_assert(Concept::isMultiDomainGrid<Grid>(),
                  "MultiDomainModelCoefficientMapper should only be used with "
                  "non-ovelapping multidomain grids");
    auto domain_set = gv.indexSet().subDomains(entity);
    assert(domain_set.size() == 1);
    _domain = *(domain_set.begin());
  }
};

} // namespace Dune::Copasi

#endif // DUNE_COPASI_GLOBAL_COEFFICIENT_MAPPER_HH