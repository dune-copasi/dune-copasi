#ifndef DUNE_COPASI_GLOBAL_COEFFICIENT_MAPPER_HH
#define DUNE_COPASI_GLOBAL_COEFFICIENT_MAPPER_HH

#include <dune/pdelab/gridfunctionspace/lfsindexcache.hh>
#include <dune/pdelab/gridfunctionspace/localfunctionspace.hh>

#include <unordered_map>
#include <vector>
#include <memory>

namespace Dune::Copasi {

struct DefaultCoefficientMapper
{
  template<class XViewLocal>
  auto operator()(const XViewLocal& x_view_local,
                  const std::size_t& comp,
                  const std::size_t& dof) const
  {
    return x_view_local(comp, dof);
  }

  template<class E>
  void bind(const E&)
  {}

  void unbind() {}

  template<class T>
  void update(const T&)
  {}
};

template<class ModelState>
struct ModelCoefficientMapper
{
  // component -> (operator,child)
  std::unordered_map<std::size_t, std::array<std::size_t, 2>> _mapper;
  const std::size_t _this_op;

  using GFS = typename ModelState::GridFunctionSpace;
  using LFS = PDELab::LocalFunctionSpace<GFS>;
  using LFSCache = PDELab::LFSIndexCache<LFS>;
  using X = typename ModelState::Coefficients;
  using XView = typename X::template ConstLocalView<LFSCache>;
  using SolutionElement = typename X::ElementType;
  using SolutionVector =
    Dune::PDELab::LocalVector<SolutionElement, Dune::PDELab::TrialSpaceTag>;

  std::vector<std::shared_ptr<LFS>> _lfs;
  std::vector<std::shared_ptr<LFSCache>> _lfs_cache;
  std::vector<std::shared_ptr<XView>> _x_view;
  std::vector<std::shared_ptr<SolutionVector>> _x;

  ModelCoefficientMapper(const std::vector<std::size_t>& map_operator,
                         std::size_t this_op)
    : _this_op(this_op)
  {
    std::vector<std::size_t> comp_child(map_operator.size(), 0);

    for (std::size_t i = 0; i < map_operator.size(); i++)
      _mapper[i] = std::array<std::size_t, 2>{ map_operator[i],
                                               comp_child[map_operator[i]]++ };
  }

  void update(const std::vector<ModelState>& states)
  {
    _x.clear();
    _x_view.clear();
    _lfs_cache.clear();
    _lfs.clear();

    _lfs.resize(states.size());
    _lfs_cache.resize(states.size());
    _x_view.resize(states.size());
    _x.resize(states.size());

    for (std::size_t i = 0; i < states.size(); i++) {
      if (_this_op == i)
        continue;
      _lfs[i] = std::make_shared<LFS>(states[i].grid_function_space);
      _lfs_cache[i] = std::make_shared<LFSCache>(*_lfs[i]);
      _x_view[i] = std::make_shared<XView>(*(states[i].coefficients));
      auto max_local_size = states[i].grid_function_space->maxLocalSize();
      _x[i] = std::make_shared<SolutionVector>(max_local_size);
    }
  }

  template<class XViewLocal>
  auto operator()(const XViewLocal& x_view_local,
                  const std::size_t& comp,
                  const std::size_t& dof) const
  {
    const auto& map = _mapper.find(comp)->second;
    if (map[0] == _this_op)
      return x_view_local(map[1], dof);
    else
      return (*_x[map[0]])(_lfs[map[0]]->child(map[1]), dof);
  }

  template<class E>
  void bind(const E& entity)
  {
    for (std::size_t i = 0; i < _lfs.size(); i++) {
      if (_this_op == i)
        continue;
      _lfs[i]->bind(entity);
      _lfs_cache[i]->update();
      _x[i]->resize(_lfs_cache[i]->size());
      _x_view[i]->bind(*_lfs_cache[i]);
      _x_view[i]->read(*_x[i]);
    }
  }

  void unbind()
  {
    for (std::size_t i = 0; i < _lfs.size(); i++)
      _x_view[i].unbind();
  }
};

} // namespace Dune::Copasi

#endif // DUNE_COPASI_GLOBAL_COEFFICIENT_MAPPER_HH