#ifndef DUNE_COPASI_LOCAL_OPERATOR_MULTIDOMAIN_DIFFUSION_REACTION_HH
#define DUNE_COPASI_LOCAL_OPERATOR_MULTIDOMAIN_DIFFUSION_REACTION_HH

#include <dune/copasi/concepts/grid.hh>
#include <dune/copasi/enum.hh>
#include <dune/copasi/local_operator.hh>

#include <dune/pdelab/localoperator/numericaljacobian.hh>
#include <dune/pdelab/localoperator/numericaljacobianapply.hh>
#include <dune/pdelab/localoperator/numericalnonlinearjacobianapply.hh>

#include <algorithm>
#include <map>
#include <vector>

namespace Dune::Copasi {

template<class Grid,
         class LocalFiniteElement,
         class CM = DefaultCoefficientMapper,
         JacobianMethod JM = JacobianMethod::Analytical>
class LocalOperatorMultiDomainDiffusionReaction
  : public Dune::PDELab::LocalOperatorDefaultFlags
  , public Dune::PDELab::InstationaryLocalOperatorDefaultMethods<double>
  , public Dune::PDELab::NumericalJacobianSkeleton<
      LocalOperatorMultiDomainDiffusionReaction<Grid,
                                                LocalFiniteElement,
                                                CM,
                                                JM>>
  , public Dune::PDELab::NumericalJacobianApplySkeleton<
      LocalOperatorMultiDomainDiffusionReaction<Grid,
                                                LocalFiniteElement,
                                                CM,
                                                JM>>
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

  //! jacobian tpye
  using IndexSet = typename Grid::LeafGridView::IndexSet;

  //! world dimension
  static constexpr int dim = LocalBasis::Traits::dimDomain;

  static constexpr std::size_t unused_domain = ~std::size_t(0);

  using GridView = typename Grid::SubDomainGrid::LeafGridView;
  using BaseLOP =
    LocalOperatorDiffusionReaction<GridView, LocalFiniteElement, CM, JM>;

  const IndexSet& _index_set;

  const LocalBasis _local_basis;

  const std::size_t _size;

  mutable std::vector<std::shared_ptr<BaseLOP>> _local_operator;

  std::vector<std::vector<std::string>> _component_name;

  // map (domain_i,domain_o,component_i) -> component_o
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
    const LocalFiniteElement& finite_element,
    std::size_t id_operator)
    : _index_set(grid->leafGridView().indexSet())
    , _local_basis(finite_element.localBasis())
    , _size(config.sub("compartments").getValueKeys().size())
    , _local_operator(_size)
    , _component_name(_size)
  {
    // std::cout << config.sub("compartments") << std::endl;
    const auto& compartments = config.sub("compartments").getValueKeys();

    // create local operators for each compartment
    for (std::size_t i = 0; i < _size; ++i) {
      const std::string compartement = compartments[i];

      int sub_domain_id =
        config.sub("compartments").template get<int>(compartement);
      GridView sub_grid_view = grid->subDomain(sub_domain_id).leafGridView();

      const auto& sub_config = config.sub(compartments[i]);
      _component_name[i] = sub_config.sub("reaction").getValueKeys();
      std::sort(_component_name[i].begin(), _component_name[i].end());

      auto lp = std::make_shared<BaseLOP>(
        sub_grid_view, sub_config, finite_element, id_operator);
      _local_operator[i] = lp;
    }

    // create mapping between all inside and outside components
    for (std::size_t domain_i = 0; domain_i < _size; ++domain_i) {
      for (std::size_t comp_i = 0; comp_i < _component_name[domain_i].size();
           comp_i++) {
        for (std::size_t domain_o = 0; domain_o < _size; ++domain_o) {
          // notice that _component_offset is agnostic on which operator is
          // given component
          for (std::size_t comp_o = 0;
               comp_o < _component_name[domain_o].size();
               comp_o++) {
            if (_component_name[domain_i][comp_i] ==
                _component_name[domain_o][comp_o]) {
              std::array<std::size_t, 3> key_i{ domain_i, domain_o, comp_i };
              _component_offset.insert(std::make_pair(key_i, comp_o));
            }
          }
        }
      }
    }
  }

  template<class States>
  void update(const States& states)
  {
    for (std::size_t i = 0; i < _local_operator.size(); ++i)
      _local_operator[i]->update(states);
  }

  // define sparsity pattern of operator representation
  template<typename LFSU, typename LFSV, typename LocalPattern>
  void pattern_volume(const LFSU& lfsu,
                      const LFSV& lfsv,
                      LocalPattern& pattern) const
  {
    for (std::size_t i = 0; i < lfsu.degree(); ++i) {
      _local_operator[i]->pattern_volume(lfsu.child(i), lfsv.child(i), pattern);
    }
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

    auto lfs_size_i = lfsu_i.child(domain_i).degree();

    for (std::size_t i = 0; i < lfs_size_i; ++i) {
      const std::size_t comp_i = _local_operator[domain_i]->_lfs_components[i];
      std::array<std::size_t, 3> inside_comp{ domain_i, domain_o, comp_i };
      auto it = _component_offset.find(inside_comp);
      if (it != _component_offset.end()) {
        auto comp_o = it->second;
        const auto& lfs_comp_o = _local_operator[domain_o]->_lfs_components;
        auto o_it = std::find(lfs_comp_o.begin(), lfs_comp_o.end(), comp_o);
        if (o_it != lfs_comp_o.end()) {
          auto lfsv_ci = lfsv_i.child(domain_i).child(i);
          std::size_t o = std::distance(lfs_comp_o.begin(), o_it);
          auto lfsu_co = lfsu_o.child(domain_o).child(o);
          for (std::size_t dof_i = 0; dof_i < lfsv_ci.size(); dof_i++) {
            for (std::size_t dof_o = 0; dof_o < lfsu_co.size(); dof_o++) {
              pattern_io.addLink(lfsv_ci, dof_i, lfsu_co, dof_o);
            }
          }
        }
      }
    }

    auto lfs_size_o = lfsu_o.child(domain_o).degree();
    for (std::size_t o = 0; o < lfs_size_o; ++o) {
      const std::size_t comp_o = _local_operator[domain_o]->_lfs_components[o];
      std::array<std::size_t, 3> outside_comp{ domain_o, domain_i, comp_o };
      auto it = _component_offset.find(outside_comp);
      if (it != _component_offset.end()) {
        auto comp_i = it->second;
        const auto& lfs_comp_i = _local_operator[domain_i]->_lfs_components;
        auto i_it = std::find(lfs_comp_i.begin(), lfs_comp_i.end(), comp_i);
        if (i_it != lfs_comp_i.end()) {
          auto lfsv_co = lfsv_o.child(domain_o).child(o);
          std::size_t i = std::distance(lfs_comp_i.begin(), i_it);
          auto lfsu_ci = lfsu_i.child(domain_i).child(i);
          for (std::size_t dof_o = 0; dof_o < lfsv_co.size(); dof_o++) {
            for (std::size_t dof_i = 0; dof_i < lfsu_ci.size(); dof_i++) {
              pattern_oi.addLink(lfsv_co, dof_o, lfsu_ci, dof_i);
            }
          }
        }
      }
    }
  }

  void setTime(double t)
  {
    Dune::PDELab::InstationaryLocalOperatorDefaultMethods<double>::setTime(t);
    for (auto lp : _local_operator)
      lp->setTime(t);
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
        const auto& sub_lfsu = lfsu.child(i);
        const auto& sub_lfsv = lfsv.child(i);
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
        const auto& sub_lfsu = lfsu.child(i);
        const auto& sub_lfsv = lfsv.child(i);
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
    // const auto& subdomain = eg.subDomain();
    for (std::size_t i = 0; i < _size; ++i) {
      if (lfsu.child(i).size() > 0) {
        const auto& sub_lfsu = lfsu.child(i);
        const auto& sub_lfsv = lfsv.child(i);
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
        const auto& sub_lfsu = lfsu.child(i);
        const auto& sub_lfsv = lfsv.child(i);
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

    const auto& entity_f = ig.intersection();
    const auto& entity_i = ig.inside();
    const auto& entity_o = ig.outside();

    assert(lfsu_i.degree() == _size);

    auto domain_set_i = _index_set.subDomains(entity_i);
    auto domain_set_o = _index_set.subDomains(entity_o);

    assert(domain_set_i.size() == 1);
    assert(domain_set_o.size() == 1);

    std::size_t domain_i = *(domain_set_i.begin());
    std::size_t domain_o = *(domain_set_o.begin());

    if (domain_i == domain_o)
      return;

    assert(lfsu_i.child(domain_i).size() == lfsv_i.child(domain_i).size());
    assert(lfsu_o.child(domain_o).size() == lfsv_o.child(domain_o).size());

    auto geo_f = entity_f.geometry();
    auto geo_in_i = entity_f.geometryInInside();
    auto geo_in_o = entity_f.geometryInOutside();

    const auto& lfsu_di = lfsu_i.child(domain_i);
    const auto& lfsu_do = lfsu_o.child(domain_o);

    std::size_t components_i = _component_name[domain_i].size();
    std::size_t components_o = _component_name[domain_o].size();

    auto x_coeff_local_i = [&](const std::size_t& component,
                               const std::size_t& dof) {
      return x_i(lfsu_di.child(component), dof);
    };
    auto x_coeff_local_o = [&](const std::size_t& component,
                               const std::size_t& dof) {
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

    auto& coefficient_mapper_i = _local_operator[domain_i]->_coefficient_mapper;
    coefficient_mapper_i.bind(entity_i);

    auto& coefficient_mapper_o = _local_operator[domain_o]->_coefficient_mapper;
    coefficient_mapper_o.bind(entity_o);

    for (const auto& it : quadratureRule(geo_f, 3)) {
      const auto& position_f = it.position();
      // position of quadrature point in local coordinates of elements
      const auto position_i = geo_in_i.global(position_f);
      const auto position_o = geo_in_o.global(position_f);

      std::vector<Range> phiu_i(lfsu_i.size());
      std::vector<Range> phiu_o(lfsu_o.size());

      // evaluate basis functions
      _local_basis.evaluateFunction(position_i, phiu_i);
      _local_basis.evaluateFunction(position_o, phiu_o);

      // evaluate concentrations at quadrature point
      DynamicVector<RF> u_i(components_i), u_o(components_o);

      for (std::size_t k = 0; k < components_i; k++) // loop over components
        for (std::size_t j = 0; j < _local_basis.size(); j++)
          u_i[k] += coefficient_mapper_i(x_coeff_local_i, k, j) * phiu_i[j];

      for (std::size_t k = 0; k < components_o; k++) // loop over components
        for (std::size_t j = 0; j < _local_basis.size(); j++)
          u_o[k] += coefficient_mapper_o(x_coeff_local_o, k, j) * phiu_o[j];

      // integration factor
      auto factor = it.weight() * geo_f.integrationElement(position_f);

      for (std::size_t i = 0; i < lfsu_di.degree(); i++) // loop over components
      {
        if (lfsu_di.child(i).size() == 0)
          continue; // child space has nothing to compute
        const std::size_t comp_i =
          _local_operator[domain_i]->_lfs_components[i];
        std::array<std::size_t, 3> inside_comp{ domain_i, domain_o, comp_i };
        auto it = _component_offset.find(inside_comp);
        if (it != _component_offset.end()) {
          for (std::size_t j = 0; j < lfsu_di.child(i).size(); j++)
            accumulate_i(
              i, j, factor * (u_i[comp_i] - u_o[it->second]) * phiu_i[j]);
        }
      }

      for (std::size_t o = 0; o < lfsu_do.degree(); o++) // loop over components
      {
        if (lfsu_do.child(o).size() == 0)
          continue; // child space has nothing to compute
        const std::size_t comp_o =
          _local_operator[domain_o]->_lfs_components[o];
        std::array<std::size_t, 3> outside_comp{ domain_o, domain_i, comp_o };
        auto it = _component_offset.find(outside_comp);
        if (it != _component_offset.end()) {
          for (std::size_t j = 0; j < lfsu_do.child(o).size(); j++)
            accumulate_o(
              o, j, -factor * (u_i[it->second] - u_o[comp_o]) * phiu_o[j]);
        }
      }
    }
  }

  template<typename IG, typename LFSU, typename X, typename LFSV, typename J>
  void jacobian_skeleton(const IG& ig,
                         const LFSU& lfsu_i,
                         const X& x_i,
                         const LFSV& lfsv_i,
                         const LFSU& lfsu_o,
                         const X& x_o,
                         const LFSV& lfsv_o,
                         J& mat_ii,
                         J& mat_io,
                         J& mat_oi,
                         J& mat_oo) const
  {
    if constexpr (JM == JacobianMethod::Numerical) {
      PDELab::NumericalJacobianSkeleton<
        LocalOperatorMultiDomainDiffusionReaction>::jacobian_skeleton(ig,
                                                                      lfsu_i,
                                                                      x_i,
                                                                      lfsv_i,
                                                                      lfsu_o,
                                                                      x_o,
                                                                      lfsv_o,
                                                                      mat_ii,
                                                                      mat_io,
                                                                      mat_oi,
                                                                      mat_oo);
      return;
    }

    const auto& entity_f = ig.intersection();
    const auto& entity_i = ig.inside();
    const auto& entity_o = ig.outside();

    assert(lfsu_i.degree() == _size);

    auto domain_set_i = _index_set.subDomains(entity_i);
    auto domain_set_o = _index_set.subDomains(entity_o);

    assert(domain_set_i.size() == 1);
    assert(domain_set_o.size() == 1);

    std::size_t domain_i = *(domain_set_i.begin());
    std::size_t domain_o = *(domain_set_o.begin());

    if (domain_i == domain_o)
      return;

    assert(lfsu_i.child(domain_i).size() == lfsv_i.child(domain_i).size());
    assert(lfsu_o.child(domain_o).size() == lfsv_o.child(domain_o).size());

    auto geo_f = entity_f.geometry();
    auto geo_in_i = entity_f.geometryInInside();
    auto geo_in_o = entity_f.geometryInOutside();

    const auto& lfsu_di = lfsu_i.child(domain_i);
    const auto& lfsu_do = lfsu_o.child(domain_o);

    auto accumulate_ii = [&](const std::size_t& component_i,
                             const std::size_t& dof_i,
                             const std::size_t& component_j,
                             const std::size_t& dof_j,
                             const auto& value) {
      mat_ii.accumulate(lfsu_di.child(component_i),
                        dof_i,
                        lfsu_di.child(component_j),
                        dof_j,
                        value);
    };

    auto accumulate_io = [&](const std::size_t& component_i,
                             const std::size_t& dof_i,
                             const std::size_t& component_j,
                             const std::size_t& dof_j,
                             const auto& value) {
      mat_io.accumulate(lfsu_di.child(component_i),
                        dof_i,
                        lfsu_do.child(component_j),
                        dof_j,
                        value);
    };

    auto accumulate_oi = [&](const std::size_t& component_i,
                             const std::size_t& dof_i,
                             const std::size_t& component_j,
                             const std::size_t& dof_j,
                             const auto& value) {
      mat_oi.accumulate(lfsu_do.child(component_i),
                        dof_i,
                        lfsu_di.child(component_j),
                        dof_j,
                        value);
    };

    auto accumulate_oo = [&](const std::size_t& component_i,
                             const std::size_t& dof_i,
                             const std::size_t& component_j,
                             const std::size_t& dof_j,
                             const auto& value) {
      mat_oo.accumulate(lfsu_do.child(component_i),
                        dof_i,
                        lfsu_do.child(component_j),
                        dof_j,
                        value);
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
      _local_basis.evaluateFunction(position_i, phiu_i);
      _local_basis.evaluateFunction(position_o, phiu_o);
      // integration factor
      auto factor = it.weight() * geo_f.integrationElement(position_f);

      for (std::size_t i = 0; i < lfsu_di.degree(); i++) // loop over components
      {
        // has child space something to compute?
        if (lfsu_di.child(i).size() == 0)
          continue;
        const std::size_t comp_i =
          _local_operator[domain_i]->_lfs_components[i];
        std::array<std::size_t, 3> inside_comp{ domain_i, domain_o, comp_i };
        auto it = _component_offset.find(inside_comp);

        if (it != _component_offset.end()) {
          // compute inside jacobian
          for (std::size_t j = 0; j < lfsu_di.child(i).size(); j++)
            for (std::size_t k = 0; k < lfsu_di.child(i).size(); k++)
              accumulate_ii(
                comp_i, j, comp_i, k, factor * phiu_i[j] * phiu_i[k]);

          // find index of outside component in lfs
          const std::size_t comp_o = it->second;
          const auto& lfs_comp_o = _local_operator[domain_o]->_lfs_components;
          auto o_it = std::find(lfs_comp_o.begin(), lfs_comp_o.end(), comp_o);
          if (o_it == lfs_comp_o.end())
            continue;
          std::size_t o = std::distance(lfs_comp_o.begin(), o_it);

          // compute inside-outside jacobian
          for (std::size_t j = 0; j < lfsu_di.child(i).size(); j++)
            for (std::size_t k = 0; k < lfsu_do.child(o).size(); k++)
              accumulate_io(
                comp_i, j, comp_o, k, -factor * phiu_i[j] * phiu_o[k]);
        }
      }

      for (std::size_t o = 0; o < lfsu_do.degree(); o++) // loop over components
      {
        // has child space something to compute?
        if (lfsu_do.child(o).size() == 0)
          continue;
        const std::size_t comp_o =
          _local_operator[domain_o]->_lfs_components[o];
        std::array<std::size_t, 3> outside_comp{ domain_o, domain_i, comp_o };
        auto it = _component_offset.find(outside_comp);

        if (it != _component_offset.end()) {
          // compute outside jacobian
          for (std::size_t j = 0; j < lfsu_do.child(o).size(); j++)
            for (std::size_t k = 0; k < lfsu_do.child(o).size(); k++)
              accumulate_oo(
                comp_o, j, comp_o, k, factor * phiu_o[j] * phiu_o[k]);

          // find index of outside component in lfs
          const std::size_t comp_i = it->second;
          const auto& lfs_comp_i = _local_operator[domain_i]->_lfs_components;
          auto i_it = std::find(lfs_comp_i.begin(), lfs_comp_i.end(), comp_i);
          if (i_it == lfs_comp_i.end())
            continue;
          std::size_t i = std::distance(lfs_comp_i.begin(), i_it);

          // compute inside-outside jacobian
          for (std::size_t j = 0; j < lfsu_do.child(o).size(); j++)
            for (std::size_t k = 0; k < lfsu_di.child(i).size(); k++)
              accumulate_oi(
                comp_o, j, comp_i, k, -factor * phiu_o[j] * phiu_i[k]);
        }
      }
    }
  }

  template<typename IG, typename LFSU, typename X, typename LFSV, typename Y>
  void jacobian_apply_skeleton(const IG& ig,
                               const LFSU& lfsu_i,
                               const X& x_i,
                               const LFSV& lfsv_i,
                               const LFSU& lfsu_o,
                               const X& x_o,
                               const LFSV& lfsv_o,
                               Y& y_i,
                               Y& y_o) const
  {
    if constexpr (JM == JacobianMethod::Numerical) {
      PDELab::NumericalJacobianApplySkeleton<
        LocalOperatorMultiDomainDiffusionReaction>::
        jacobian_apply_skeleton(
          ig, lfsu_i, x_i, lfsv_i, lfsu_o, x_o, lfsv_o, y_i, y_o);
      return;
    }
    DUNE_THROW(NotImplemented, "Analytic jacobian is not implemented");
  }

  template<typename IG, typename LFSU, typename X, typename LFSV, typename Y>
  void jacobian_apply_skeleton(const IG& ig,
                               const LFSU& lfsu_i,
                               const X& x_i,
                               const X& z_i,
                               const LFSV& lfsv_i,
                               const LFSU& lfsu_o,
                               const X& x_o,
                               const X& z_o,
                               const LFSV& lfsv_o,
                               Y& y_i,
                               Y& y_o) const
  {
    if constexpr (JM == JacobianMethod::Numerical) {
      PDELab::NumericalJacobianApplySkeleton<
        LocalOperatorMultiDomainDiffusionReaction>::
        jacobian_apply_skeleton(
          ig, lfsu_i, x_i, z_i, lfsv_i, lfsu_o, x_o, z_o, lfsv_o, y_i, y_o);
      return;
    }
    DUNE_THROW(NotImplemented, "Analytic jacobian is not implemented");
  }
};

template<class Grid,
         class LocalFiniteElement,
         JacobianMethod JM = JacobianMethod::Analytical>
class TemporalLocalOperatorMultiDomainDiffusionReaction
  : public Dune::PDELab::LocalOperatorDefaultFlags
  , public Dune::PDELab::InstationaryLocalOperatorDefaultMethods<double>
{
  static_assert(Concept::isMultiDomainGrid<Grid>());

  using GridView = typename Grid::SubDomainGrid::LeafGridView;
  using BaseLOP =
    TemporalLocalOperatorDiffusionReaction<GridView, LocalFiniteElement, JM>;

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
    const LocalFiniteElement& finite_element,
    std::size_t id_operator)
    : _size(config.sub("compartments").getValueKeys().size())
    , _local_operator(_size)
  {
    const auto& compartments = config.sub("compartments").getValueKeys();

    for (std::size_t i = 0; i < _size; ++i) {
      const std::string compartement = compartments[i];

      int sub_domain_id =
        config.sub("compartments").template get<int>(compartement);
      GridView sub_grid_view = grid->subDomain(sub_domain_id).leafGridView();

      const auto& sub_config = config.sub(compartments[i]);
      _local_operator[i] = std::make_shared<BaseLOP>(
        sub_grid_view, sub_config, finite_element, id_operator);
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
        const auto& sub_lfsu = lfsu.child(i);
        const auto& sub_lfsv = lfsv.child(i);
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
        const auto& sub_lfsu = lfsu.child(i);
        const auto& sub_lfsv = lfsv.child(i);
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
        const auto& sub_lfsu = lfsu.child(i);
        const auto& sub_lfsv = lfsv.child(i);
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
        const auto& sub_lfsu = lfsu.child(i);
        const auto& sub_lfsv = lfsv.child(i);
        if (sub_lfsu.size() == 0)
          continue;
        _local_operator[i]->jacobian_apply_volume(eg, sub_lfsu, x, sub_lfsv, r);
      }
    }
  }
};

} // namespace Dune::Copasi

#endif // DUNE_COPASI_LOCAL_OPERATOR_MULTIDOMAIN_DIFFUSION_REACTION_HH