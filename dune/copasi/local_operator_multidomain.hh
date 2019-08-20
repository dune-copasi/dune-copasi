#ifndef DUNE_COPASI_LOCAL_OPERATOR_MULTIDOMAIN_DIFFUSION_REACTION_HH
#define DUNE_COPASI_LOCAL_OPERATOR_MULTIDOMAIN_DIFFUSION_REACTION_HH

#include <dune/copasi/concepts/grid.hh>
#include <dune/copasi/local_operator.hh>

#include <dune/pdelab/localoperator/numericaljacobian.hh>

#include <algorithm>

namespace Dune::Copasi {

template<class Grid, class LocalFiniteElement>
class LocalOperatorMultiDomainDiffusionReaction
  : public Dune::PDELab::LocalOperatorDefaultFlags
  , public Dune::PDELab::InstationaryLocalOperatorDefaultMethods<double>
  , public Dune::PDELab::NumericalJacobianSkeleton<
      LocalOperatorMultiDomainDiffusionReaction<Grid, LocalFiniteElement>>
{
  static_assert(Concept::isMultiDomainGrid<Grid>());

  //! local basis
  using LocalBasis = typename LocalFiniteElement::Traits::LocalBasisType;

  //! domain field
  using DF = typename LocalBasis::Traits::DomainFieldType;

  //! coordinates type
  using Domain = typename LocalBasis::Traits::DomainType;

  //! range field
  using RF = typename LocalBasis::Traits::RangeFieldType;

  //! range type (for the local finite element)
  using Range = typename LocalBasis::Traits::RangeType;

  //! jacobian tpye
  using Jacobian = typename LocalBasis::Traits::JacobianType;

  //! world dimension
  static constexpr int dim = LocalBasis::Traits::dimDomain;

  static constexpr std::size_t unused_domain = ~std::size_t(0);

  using GridView = typename Grid::SubDomainGrid::LeafGridView;
  using BaseLOP = LocalOperatorDiffusionReaction<GridView, LocalFiniteElement>;

  //! number of basis per component
  const std::size_t _basis_size;

  const std::size_t _size;

  std::vector<std::shared_ptr<BaseLOP>> _local_operator;

  std::vector<std::vector<std::string>> _component_name;

  std::map<std::array<std::size_t, 3>, std::size_t> _component_offset;

public:
  //! visit skeleton from the two sides
  static constexpr bool doSkeletonTwoSided = true;

  //! pattern assembly flags
  static constexpr bool doPatternVolume = BaseLOP::doPatternVolume;

  //! pattern assembly flags
  static constexpr bool doPatternSkeleton = true;

  //! residual assembly flags
  static constexpr bool doAlphaVolume = BaseLOP::doAlphaVolume;

  //! residual assembly flags
  static constexpr bool doAlphaSkeleton = true;

  LocalOperatorMultiDomainDiffusionReaction(
    std::shared_ptr<const Grid> grid,
    const ParameterTree& config,
    const LocalFiniteElement& finite_element)
    : _basis_size(finite_element.localBasis().size())
    , _size(config.sub("compartements").getValueKeys().size())
    , _local_operator(_size)
    , _component_name(_size)
  {
    const auto& compartements = config.sub("compartements").getValueKeys();
    // warning if size > _size
    for (std::size_t i = 0; i < _size; ++i) {
      const std::string compartement = compartements[i];

      int sub_domain_id =
        config.sub("compartements").template get<int>(compartement);
      GridView sub_grid_view = grid->subDomain(sub_domain_id).leafGridView();

      auto sub_config = config.sub(compartements[i]);
      auto lp =
        std::make_shared<BaseLOP>(sub_grid_view, sub_config, finite_element);
      _local_operator[i] = lp;

      _component_name[i] = sub_config.sub("reaction").getValueKeys();
      std::sort(_component_name[i].begin(), _component_name[i].end());
    }

    for (std::size_t i = 0; i < _size; ++i) {
      for (std::size_t j = 0; j < _component_name[i].size(); j++) {
        for (std::size_t k = 0; k < _size; ++k) {
          for (std::size_t l = 0; l < _component_name[k].size(); l++) {
            if (_component_name[i][j] == _component_name[k][l]) {
              std::array<std::size_t, 3> inside_comp{ i, k, j };
              _component_offset.insert(std::make_pair(inside_comp, l));
            }
          }
        }
      }
    }
  }

  // define sparsity pattern of operator representation
  template<typename LFSU, typename LFSV, typename LocalPattern>
  void pattern_volume(const LFSU& lfsu,
                      const LFSV& lfsv,
                      LocalPattern& pattern) const
  {
    for (std::size_t i = 0; i < lfsu.degree(); ++i)
      _local_operator[i]->pattern_volume(lfsu.child(i), lfsv.child(i), pattern);
  }

  template<typename LFSU, typename LFSV, typename LocalPattern>
  void pattern_skeleton(const LFSU& lfsu_i,
                        const LFSV& lfsv_i,
                        const LFSU& lfsu_o,
                        const LFSV& lfsv_o,
                        LocalPattern& pattern_io,
                        LocalPattern& pattern_oi) const
  {
    std::size_t domain_i(unused_domain), domain_o(unused_domain);
    for (std::size_t k = 0; k < _size; k++) {
      if (lfsu_i.child(k).size() > 0)
        domain_i = k;
      if (lfsu_o.child(k).size() > 0)
        domain_o = k;
    }

    if (domain_i == domain_o)
      return;
    if ((domain_i == unused_domain) or (domain_o == unused_domain))
      return;

    assert(lfsu_i.degree() == _size);
    assert(lfsu_i.child(domain_i).degree() == lfsv_i.child(domain_i).degree());
    assert(lfsu_o.child(domain_o).degree() == lfsv_o.child(domain_o).degree());

    auto comp_size_i = lfsu_i.child(domain_i).degree();
    auto comp_size_o = lfsu_o.child(domain_o).degree();

    for (std::size_t comp_i = 0; comp_i < comp_size_i; ++comp_i) {
      std::array<std::size_t, 3> inside_comp{ domain_i, domain_o, comp_i };
      auto it = _component_offset.find(inside_comp);
      if (it != _component_offset.end()) {
        auto comp_o = it->second;
        auto lfsv_ci = lfsv_i.child(domain_i).child(comp_i);
        auto lfsu_co = lfsu_o.child(domain_o).child(comp_o);
        for (std::size_t dof_i = 0; dof_i < lfsv_ci.size(); dof_i++)
          for (std::size_t dof_o = 0; dof_o < lfsu_co.size(); dof_o++) {
            pattern_io.addLink(lfsv_ci, dof_i, lfsu_co, dof_o);
          }
      }
    }

    for (std::size_t comp_o = 0; comp_o < comp_size_o; ++comp_o) {
      std::array<std::size_t, 3> outside_comp{ domain_o, domain_i, comp_o };
      auto it = _component_offset.find(outside_comp);
      if (it != _component_offset.end()) {
        auto comp_i = it->second;
        auto lfsv_co = lfsv_o.child(domain_o).child(comp_o);
        auto lfsu_ci = lfsu_i.child(domain_i).child(comp_i);
        for (std::size_t dof_o = 0; dof_o < lfsv_co.size(); dof_o++)
          for (std::size_t dof_i = 0; dof_i < lfsu_ci.size(); dof_i++)
            pattern_oi.addLink(lfsv_co, dof_o, lfsu_ci, dof_i);
      }
    }
  }

  template<typename EG, typename LFSU, typename X, typename LFSV, typename R>
  void jacobian_apply_volume(const EG& eg,
                             const LFSU& lfsu,
                             const X& x,
                             const X& z,
                             const LFSV& lfsv,
                             R& r) const
  {
    for (std::size_t i = 0; i < _size; ++i) {
      if (lfsu.child(i).size() > 0) {
        auto sub_lfsu = lfsu.child(i);
        auto sub_lfsv = lfsv.child(i);
        if (sub_lfsu.size() == 0)
          continue;
        _local_operator[i]->jacobian_apply_volume(
          eg, sub_lfsu, x, z, sub_lfsv, r);
      }
    }
  }

  //! jacobian contribution of volume term
  template<typename EG, typename LFSU, typename X, typename LFSV, typename M>
  void jacobian_volume(const EG& eg,
                       const LFSU& lfsu,
                       const X& x,
                       const LFSV& lfsv,
                       M& mat) const
  {
    for (std::size_t i = 0; i < _size; ++i) {
      if (lfsu.child(i).size() > 0) {
        auto sub_lfsu = lfsu.child(i);
        auto sub_lfsv = lfsv.child(i);
        if (sub_lfsu.size() == 0)
          continue;
        _local_operator[i]->jacobian_volume(eg, sub_lfsu, x, sub_lfsv, mat);
      }
    }
  }

  template<typename EG, typename LFSU, typename X, typename LFSV, typename R>
  void alpha_volume(const EG& eg,
                    const LFSU& lfsu,
                    const X& x,
                    const LFSV& lfsv,
                    R& r) const
  {
    // auto subdomain = eg.subDomain();
    for (std::size_t i = 0; i < _size; ++i) {
      if (lfsu.child(i).size() > 0) {
        auto sub_lfsu = lfsu.child(i);
        auto sub_lfsv = lfsv.child(i);
        if (sub_lfsu.size() == 0)
          continue;
        _local_operator[i]->alpha_volume(eg, sub_lfsu, x, sub_lfsv, r);
      }
    }
  }

  //! apply local jacobian of the volume term -> linear variant
  template<typename EG, typename LFSU, typename X, typename LFSV, typename R>
  void jacobian_apply_volume(const EG& eg,
                             const LFSU& lfsu,
                             const X& x,
                             const LFSV& lfsv,
                             R& r) const
  {
    for (std::size_t i = 0; i < _size; ++i) {
      if (lfsu.child(i).size() > 0) {
        auto sub_lfsu = lfsu.child(i);
        auto sub_lfsv = lfsv.child(i);
        if (sub_lfsu.size() == 0)
          continue;
        _local_operator[i]->jacobian_apply_volume(eg, sub_lfsu, x, sub_lfsv, r);
      }
    }
  }

  // skeleton integral depending on test and ansatz functions
  template<typename IG, typename LFSU, typename X, typename LFSV, typename R>
  void alpha_skeleton(const IG& ig,
                      const LFSU& lfsu_i,
                      const X& x_i,
                      const LFSV& lfsv_i,
                      const LFSU& lfsu_o,
                      const X& x_o,
                      const LFSV& lfsv_o,
                      R& r_i,
                      R& r_o) const
  {
    std::size_t domain_i(unused_domain), domain_o(unused_domain);
    for (std::size_t k = 0; k < _size; k++) {
      if (lfsu_i.child(k).size() > 0)
        domain_i = k;
      if (lfsu_o.child(k).size() > 0)
        domain_o = k;
    }

    if (domain_i == domain_o)
      return;
    if ((domain_i == unused_domain) or (domain_o == unused_domain))
      return;

    assert(lfsu_i.degree() == _size);
    assert(lfsu_i.child(domain_i).size() == lfsv_i.child(domain_i).size());
    assert(lfsu_o.child(domain_o).size() == lfsv_o.child(domain_o).size());

    const auto& entity_f = ig.intersection();

    auto geo_f = entity_f.geometry();
    auto geo_in_i = entity_f.geometryInInside();
    auto geo_in_o = entity_f.geometryInOutside();

    const auto& lfsu_di = lfsu_i.child(domain_i);
    const auto& lfsu_do = lfsu_o.child(domain_o);

    std::size_t components_i = lfsu_di.degree();
    std::size_t components_o = lfsu_do.degree();

    auto lfsu_basis_i = lfsu_di.child(0).finiteElement().localBasis();
    auto lfsu_basis_o = lfsu_do.child(0).finiteElement().localBasis();

    auto x_coeff_i = [&](const std::size_t& component, const std::size_t& dof) {
      return x_i(lfsu_di.child(component), dof);
    };
    auto x_coeff_o = [&](const std::size_t& component, const std::size_t& dof) {
      return x_o(lfsu_do.child(component), dof);
    };

    auto accumulate_i = [&](const std::size_t& component,
                            const std::size_t& dof,
                            const auto& value) {
      r_i.accumulate(lfsu_di.child(component), dof, value);
    };

    auto accumulate_o = [&](const std::size_t& component,
                            const std::size_t& dof,
                            const auto& value) {
      r_o.accumulate(lfsu_do.child(component), dof, value);
    };

    typename IG::Entity::Geometry::JacobianInverseTransposed jac;

    for (const auto& it : quadratureRule(geo_f, 3)) {
      const auto& position_f = it.position();
      // position of quadrature point in local coordinates of elements
      const auto position_i = geo_in_i.global(position_f);
      const auto position_o = geo_in_o.global(position_f);

      std::vector<Range> phiu_i(lfsu_i.size());
      std::vector<Range> phiu_o(lfsu_o.size());

      // evaluate basis functions
      lfsu_basis_i.evaluateFunction(position_i, phiu_i);
      lfsu_basis_o.evaluateFunction(position_o, phiu_o);

      // evaluate concentrations at quadrature point
      DynamicVector<RF> u_i(components_i), u_o(components_o);

      for (std::size_t k = 0; k < components_i; k++) // loop over components
        for (std::size_t j = 0; j < lfsu_di.child(k).size(); j++)
          u_i[k] += x_coeff_i(k, j) * phiu_i[j];

      for (std::size_t k = 0; k < components_o; k++) // loop over components
        for (std::size_t j = 0; j < lfsu_do.child(k).size(); j++)
          u_o[k] += x_coeff_o(k, j) * phiu_o[j];

      // integration factor
      auto factor = it.weight() * geo_f.integrationElement(position_f);

      for (std::size_t k = 0; k < components_i; k++) // loop over components
      {
        std::array<std::size_t, 3> inside_comp{ domain_i, domain_o, k };
        auto it = _component_offset.find(inside_comp);
        if (it != _component_offset.end()) {
          for (std::size_t j = 0; j < lfsu_di.child(k).size(); j++)
            accumulate_i(k, j, factor * (u_i[k] - u_o[it->second]) * phiu_i[j]);
        }
      }

      for (std::size_t k = 0; k < components_o; k++) // loop over components
      {
        std::array<std::size_t, 3> outside_comp{ domain_o, domain_i, k };
        auto it = _component_offset.find(outside_comp);
        if (it != _component_offset.end()) {
          for (std::size_t j = 0; j < lfsu_do.child(k).size(); j++)
            accumulate_o(
              k, j, -factor * (u_i[it->second] - u_o[k]) * phiu_o[j]);
        }
      }
    }
  }
};

template<class Grid, class LocalFiniteElement>
class TemporalLocalOperatorMultiDomainDiffusionReaction
  : public Dune::PDELab::LocalOperatorDefaultFlags
  , public Dune::PDELab::InstationaryLocalOperatorDefaultMethods<double>
{
  static_assert(Concept::isMultiDomainGrid<Grid>());

  using GridView = typename Grid::SubDomainGrid::LeafGridView;
  using BaseLOP =
    TemporalLocalOperatorDiffusionReaction<GridView, LocalFiniteElement>;

  std::size_t _size;

  std::vector<std::shared_ptr<BaseLOP>> _local_operator;

public:
  //! pattern assembly flags
  static constexpr bool doPatternVolume = true;

  //! residual assembly flags
  static constexpr bool doAlphaVolume = true;

  TemporalLocalOperatorMultiDomainDiffusionReaction(
    std::shared_ptr<const Grid> grid,
    const ParameterTree& config,
    const LocalFiniteElement& finite_element)
    : _size(config.sub("compartements").getValueKeys().size())
    , _local_operator(_size)
  {
    const auto& compartements = config.sub("compartements").getValueKeys();

    for (std::size_t i = 0; i < _size; ++i) {
      const std::string compartement = compartements[i];

      int sub_domain_id =
        config.sub("compartements").template get<int>(compartement);
      GridView sub_grid_view = grid->subDomain(sub_domain_id).leafGridView();

      auto sub_config = config.sub(compartements[i]);
      _local_operator[i] =
        std::make_shared<BaseLOP>(sub_grid_view, sub_config, finite_element);
    }
  }

  // define sparsity pattern of operator representation
  template<typename LFSU, typename LFSV, typename LocalPattern>
  void pattern_volume(const LFSU& lfsu,
                      const LFSV& lfsv,
                      LocalPattern& pattern) const
  {
    for (std::size_t i = 0; i < _size; ++i)
      _local_operator[i]->pattern_volume(lfsu.child(i), lfsv.child(i), pattern);
  }

  template<typename EG, typename LFSU, typename X, typename LFSV, typename R>
  void alpha_volume(const EG& eg,
                    const LFSU& lfsu,
                    const X& x,
                    const LFSV& lfsv,
                    R& r) const
  {
    for (std::size_t i = 0; i < _size; ++i) {
      if (lfsu.child(i).size() > 0) {
        auto sub_lfsu = lfsu.child(i);
        auto sub_lfsv = lfsv.child(i);
        if (sub_lfsu.size() == 0)
          continue;
        _local_operator[i]->alpha_volume(eg, sub_lfsu, x, sub_lfsv, r);
      }
    }
  }

  template<typename EG, typename LFSU, typename X, typename LFSV, typename Mat>
  void jacobian_volume(const EG& eg,
                       const LFSU& lfsu,
                       const X& x,
                       const LFSV& lfsv,
                       Mat& mat) const
  {
    for (std::size_t i = 0; i < _size; ++i) {
      if (lfsu.child(i).size() > 0) {
        auto sub_lfsu = lfsu.child(i);
        auto sub_lfsv = lfsv.child(i);
        if (sub_lfsu.size() == 0)
          continue;
        _local_operator[i]->jacobian_volume(eg, sub_lfsu, x, sub_lfsv, mat);
      }
    }
  }

  template<typename EG, typename LFSU, typename X, typename LFSV, typename R>
  void jacobian_apply_volume(const EG& eg,
                             const LFSU& lfsu,
                             const X& x,
                             const X& z,
                             const LFSV& lfsv,
                             R& r) const
  {
    for (std::size_t i = 0; i < _size; ++i) {
      if (lfsu.child(i).size() > 0) {
        auto sub_lfsu = lfsu.child(i);
        auto sub_lfsv = lfsv.child(i);
        if (sub_lfsu.size() == 0)
          continue;
        _local_operator[i]->jacobian_apply_volume(
          eg, sub_lfsu, x, z, sub_lfsv, r);
      }
    }
  }

  template<typename EG, typename LFSU, typename X, typename LFSV, typename R>
  void jacobian_apply_volume(const EG& eg,
                             const LFSU& lfsu,
                             const X& x,
                             const LFSV& lfsv,
                             R& r) const
  {
    for (std::size_t i = 0; i < _size; ++i) {
      if (lfsu.child(i).size() > 0) {
        auto sub_lfsu = lfsu.child(i);
        auto sub_lfsv = lfsv.child(i);
        if (sub_lfsu.size() == 0)
          continue;
        _local_operator[i]->jacobian_apply_volume(eg, sub_lfsu, x, sub_lfsv, r);
      }
    }
  }
};

} // namespace Dune::Copasi

#endif // DUNE_COPASI_LOCAL_OPERATOR_MULTIDOMAIN_DIFFUSION_REACTION_HH