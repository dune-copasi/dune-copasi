#ifndef DUNE_COPASI_LOCAL_OPERATOR_MULTIDOMAIN_DIFFUSION_REACTION_HH
#define DUNE_COPASI_LOCAL_OPERATOR_MULTIDOMAIN_DIFFUSION_REACTION_HH

#include <dune/copasi/common/enum.hh>
#include <dune/copasi/concepts/grid.hh>

#include <dune/pdelab/localoperator/flags.hh>
#include <dune/pdelab/localoperator/idefault.hh>
#include <dune/pdelab/localoperator/numericaljacobian.hh>
#include <dune/pdelab/localoperator/numericaljacobianapply.hh>
#include <dune/pdelab/localoperator/numericalnonlinearjacobianapply.hh>

#include <dune/common/parametertree.hh>

#include <algorithm>
#include <map>
#include <vector>

namespace Dune::Copasi {

/**
 * @brief      This class describes a PDELab local operator for multi domain
 *             diffusion reaction systems.
 * @details    This class describre the operatrions for local integrals required
 *             for diffusion reaction system. The operator is only valid for
 *             entities contained in the grid. The local finite element is used
 *             for caching shape function evaluations. The coefficient mapper is
 *             the interface to fetch date from either local or external
 *             coefficient vectors. And the jacobian method switches between
 *             numerical and analytical jacobians. This local operator creates
 *             internally an individual local operator for every subdomain in
 *             the grid
 *
 * @tparam     Grid                The grid
 * @tparam     LocalFiniteElement  Local Finite Element
 * @tparam     CM                  Coefficient Mapper
 * @tparam     JM                  Jacobian Method
 */
template<class Grid,
         class SubLocalOperator,
         JacobianMethod JM = JacobianMethod::Analytical>
class LocalOperatorMultiDomainDiffusionReaction
  : public Dune::PDELab::LocalOperatorDefaultFlags
  , public Dune::PDELab::InstationaryLocalOperatorDefaultMethods<double>
  , public Dune::PDELab::NumericalJacobianSkeleton<
      LocalOperatorMultiDomainDiffusionReaction<Grid,
                                                SubLocalOperator,
                                                JM>>
  , public Dune::PDELab::NumericalJacobianApplySkeleton<
      LocalOperatorMultiDomainDiffusionReaction<Grid,
                                                SubLocalOperator,
                                                JM>>
{
  static_assert(Concept::isMultiDomainGrid<Grid>());

  //! jacobian tpye
  using IndexSet = typename Grid::LeafGridView::IndexSet;

  static constexpr std::size_t unused_domain = ~std::size_t(0);

  using GridView = typename Grid::SubDomainGrid::LeafGridView;
  using SubLOP = SubLocalOperator;

  const IndexSet& _index_set;

  const std::size_t _size;

  mutable std::vector<std::shared_ptr<SubLOP>> _local_operator;

  std::vector<std::vector<std::string>> _component_name;

  // map (domain_i,domain_o,component_i) -> component_o
  std::map<std::array<std::size_t, 3>, std::size_t> _component_offset;

public:

  static_assert(not SubLOP::doSkeletonTwoSided);

  //! pattern assembly flags
  static constexpr bool doPatternVolume = SubLOP::doPatternVolume;

  //! pattern assembly flags
  static constexpr bool doPatternSkeleton = true;

  //! residual assembly flags
  static constexpr bool doAlphaVolume = SubLOP::doAlphaVolume;

  //! residual assembly flags
  static constexpr bool doAlphaSkeleton = true;

  /**
   * @brief      Constructs a new instance.
   *
   * @param[in]  grid            The grid
   * @param[in]  config          The configuration
   * @param[in]  finite_element  The local finite element
   * @param[in]  id_operator     The index of this operator
   */
  LocalOperatorMultiDomainDiffusionReaction(
    std::shared_ptr<const Grid> grid,
    const ParameterTree& config,
    std::size_t id_operator)
    : _index_set(grid->leafGridView().indexSet())
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

      auto lp = std::make_shared<SubLOP>(
        sub_grid_view, sub_config, id_operator);
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

  /**
   * @brief      Updates the coefficient mapper with given states.
   *
   * @param[in]  states  A map from operator index to states
   *
   * @tparam     States  Map from index to states
   */
  template<class States>
  void update(const States& states)
  {
    for (std::size_t i = 0; i < _local_operator.size(); ++i)
      _local_operator[i]->update(states);
  }

  /**
   * @copydoc LocalOperatorDiffusionReactionCG::pattern_volume
   */
  template<typename LFSU, typename LFSV, typename LocalPattern>
  void pattern_volume(const LFSU& lfsu,
                      const LFSV& lfsv,
                      LocalPattern& pattern) const
  {
    for (std::size_t i = 0; i < lfsu.degree(); ++i) {
      _local_operator[i]->pattern_volume(lfsu.child(i), lfsv.child(i), pattern);
    }
  }

  /**
   * @brief      Pattern sckeleton
   * @details    This method links degrees of freedom between trial and test
   *             spaces at entities intersection taking into account the
   *             structure of the reaction term
   *
   * @param[in]  lfsu_i        The inside trial local function space
   * @param[in]  lfsv_i        The inside test local function space
   * @param[in]  lfsu_o        The outside trial local function space
   * @param[in]  lfsv_o        The outside test local function space
   * @param      pattern_io    The inside-outside local pattern
   * @param      pattern_oi    The outside-inside local pattern
   *
   * @tparam     LFSU          The trial local function space
   * @tparam     LFSV          The test local function space
   * @tparam     LocalPattern  The local pattern
   */
  template<typename LFSU, typename LFSV, typename LocalPattern>
  void interface_pattern_skeleton(std::size_t domain_i, std::size_t domain_o,
                        const LFSU& lfsu_di,
                        const LFSV& lfsv_di,
                        const LFSU& lfsu_do,
                        const LFSV& lfsv_do,
                        LocalPattern& pattern_io,
                        LocalPattern& pattern_oi) const
  {
    auto lfs_size_i = lfsu_di.degree();
    for (std::size_t i = 0; i < lfs_size_i; ++i) {
      const std::size_t comp_i = _local_operator[domain_i]->lfs_components(lfsu_di,lfsv_di)[i];
      std::array<std::size_t, 3> inside_comp{ domain_i, domain_o, comp_i };
      auto it = _component_offset.find(inside_comp);
      if (it != _component_offset.end()) {
        auto comp_o = it->second;
        const auto& lfs_comp_o = _local_operator[domain_o]->lfs_components(lfsu_do,lfsv_do);
        auto o_it = std::find(lfs_comp_o.begin(), lfs_comp_o.end(), comp_o);
        if (o_it != lfs_comp_o.end()) {
          auto lfsv_ci = lfsv_di.child(i);
          std::size_t o = std::distance(lfs_comp_o.begin(), o_it);
          auto lfsu_co = lfsu_do.child(o);
          for (std::size_t dof_i = 0; dof_i < lfsv_ci.size(); dof_i++) {
            for (std::size_t dof_o = 0; dof_o < lfsu_co.size(); dof_o++) {
              pattern_io.addLink(lfsv_ci, dof_i, lfsu_co, dof_o);
            }
          }
        }
      }
    }

    auto lfs_size_o = lfsu_do.degree();
    for (std::size_t o = 0; o < lfs_size_o; ++o) {
      const std::size_t comp_o = _local_operator[domain_o]->lfs_components(lfsu_do,lfsv_do)[o];
      std::array<std::size_t, 3> outside_comp{ domain_o, domain_i, comp_o };
      auto it = _component_offset.find(outside_comp);
      if (it != _component_offset.end()) {
        auto comp_i = it->second;
        const auto& lfs_comp_i = _local_operator[domain_i]->lfs_components(lfsu_do,lfsv_do);
        auto i_it = std::find(lfs_comp_i.begin(), lfs_comp_i.end(), comp_i);
        if (i_it != lfs_comp_i.end()) {
          auto lfsv_co = lfsv_do.child(o);
          std::size_t i = std::distance(lfs_comp_i.begin(), i_it);
          auto lfsu_ci = lfsu_di.child(i);
          for (std::size_t dof_o = 0; dof_o < lfsv_co.size(); dof_o++) {
            for (std::size_t dof_i = 0; dof_i < lfsu_ci.size(); dof_i++) {
              pattern_oi.addLink(lfsv_co, dof_o, lfsu_ci, dof_i);
            }
          }
        }
      }
    }
  }

  /**
   * @brief      Pattern sckeleton
   * @details    This method links degrees of freedom between trial and test
   *             spaces at entities intersection taking into account the
   *             structure of the reaction term
   *
   * @param[in]  lfsu_i        The inside trial local function space
   * @param[in]  lfsv_i        The inside test local function space
   * @param[in]  lfsu_o        The outside trial local function space
   * @param[in]  lfsv_o        The outside test local function space
   * @param      pattern_io    The inside-outside local pattern
   * @param      pattern_oi    The outside-inside local pattern
   *
   * @tparam     LFSU          The trial local function space
   * @tparam     LFSV          The test local function space
   * @tparam     LocalPattern  The local pattern
   */
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

    if ((domain_i == unused_domain) or (domain_o == unused_domain))
      return;

    const auto& lfsu_di = lfsu_i.child(domain_i);
    const auto& lfsv_di = lfsv_i.child(domain_i);
    const auto& lfsu_do = lfsu_o.child(domain_o);
    const auto& lfsv_do = lfsv_o.child(domain_o);

    assert(lfsu_i.degree() == _size);
    assert(lfsu_di.degree() == lfsv_di.degree());
    assert(lfsu_do.degree() == lfsv_do.degree());

    if (domain_i != domain_o)
      interface_pattern_skeleton(domain_i,domain_o,lfsu_di,lfsv_di,lfsu_do,lfsv_do,pattern_io,pattern_oi);
    else if constexpr (SubLOP::doPatternSkeleton)
      _local_operator[domain_i]->pattern_skeleton(lfsu_di,lfsv_di,lfsu_do,lfsv_do,pattern_io,pattern_oi);
  }

  /**
   * @brief      Sets the time.
   *
   * @param[in]  t     The new time
   */
  void setTime(double t)
  {
    Dune::PDELab::InstationaryLocalOperatorDefaultMethods<double>::setTime(t);
    for (auto lp : _local_operator)
      lp->setTime(t);
  }

  /**
   * @copydoc LocalOperatorDiffusionReactionCG::jacobian_apply_volume
   * @details    This particular operator does a jacobian apply volume for the
   *             LocalOperatorDiffusionReactionCG corresponding to incoming
   * entity
   */
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
        _local_operator[i]->jacobian_apply_volume(
          eg, sub_lfsu, x, z, sub_lfsv, r);
      }
    }
  }

  /**
   * @copydoc LocalOperatorDiffusionReactionCG::jacobian_apply_volume
   * @details    This particular operator does a jacobian apply volume for the
   *             LocalOperatorDiffusionReactionCG corresponding to incoming
   * entity
   */
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
        _local_operator[i]->jacobian_apply_volume(eg, sub_lfsu, x, sub_lfsv, r);
      }
    }
  }

  /**
   * @copydoc LocalOperatorDiffusionReactionCG::jacobian_volume
   * @details    This particular operator does a jacobian volume for the
   *             LocalOperatorDiffusionReactionCG corresponding to incoming
   * entity
   */
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
        _local_operator[i]->jacobian_volume(eg, sub_lfsu, x, sub_lfsv, mat);
      }
    }
  }

  /**
   * @copydoc LocalOperatorDiffusionReactionCG::alpha_volume
   * @details    This particular operator does a alpha volume for the
   *             LocalOperatorDiffusionReactionCG corresponding to incoming
   * entity
   */
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
        _local_operator[i]->alpha_volume(eg, sub_lfsu, x, sub_lfsv, r);
      }
    }
  }

  /**
   * @brief      The skeleton integral
   * @details    This integral is only performed at the interface between
   *             different domains. Currently it has the form of
   *             dichlet-dirichlet boundary condition between domains
   *
   * @param[in]  ig      The intersection
   * @param[in]  lfsu_i  The inside trial local function space
   * @param[in]  x_i     The inside local coefficient vector
   * @param[in]  lfsv_i  The inside test local function space
   * @param[in]  lfsu_o  The outside trial local function space
   * @param[in]  x_o     The outside local coefficient vector
   * @param[in]  lfsv_o  The outside test local function space
   * @param      r_i     The inside residual vector
   * @param      r_o     The outside residual vector
   *
   * @tparam     IG      The indersection
   * @tparam     LFSU    The trial local function space
   * @tparam     X       The local coefficient vector
   * @tparam     LFSV    The test local function space
   * @tparam     R       The residual vector
   */
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
    const auto& entity_i = ig.inside();
    const auto& entity_o = ig.outside();

    assert(_size >= lfsu_i.degree());
    assert(_size >= lfsu_o.degree());
    assert(lfsu_i.degree() == lfsv_i.degree());
    assert(lfsu_o.degree() == lfsv_o.degree());

    auto domain_set_i = _index_set.subDomains(entity_i);
    auto domain_set_o = _index_set.subDomains(entity_o);

    assert(domain_set_i.size() == 1);
    assert(domain_set_o.size() == 1);

    std::size_t domain_i = *(domain_set_i.begin());
    std::size_t domain_o = *(domain_set_o.begin());

    const auto& lfsu_di = lfsu_i.child(domain_i);
    const auto& lfsv_di = lfsv_i.child(domain_i);
    const auto& lfsu_do = lfsu_o.child(domain_o);
    const auto& lfsv_do = lfsv_o.child(domain_o);

    assert(lfsu_di.size() == lfsv_di.size());
    assert(lfsu_do.size() == lfsv_do.size());

    if (domain_i != domain_o)
      interface_alpha_skeleton(domain_i,domain_o,ig,lfsu_di,x_i,lfsv_di,lfsu_do,x_o,lfsv_do,r_i,r_o);
    else if constexpr (SubLocalOperator::doAlphaSkeleton)
      _local_operator[domain_i]->alpha_skeleton(ig,lfsu_di,x_i,lfsv_di,lfsu_do,x_o,lfsv_do,r_i,r_o);
  }

  /**
   * @brief      The skeleton integral
   * @details    This integral is only performed at the interface between
   *             different domains. Currently it has the form of
   *             dichlet-dirichlet boundary condition between domains
   *
   * @param[in]  ig      The intersection
   * @param[in]  lfsu_i  The inside trial local function space
   * @param[in]  x_i     The inside local coefficient vector
   * @param[in]  lfsv_i  The inside test local function space
   * @param[in]  lfsu_o  The outside trial local function space
   * @param[in]  x_o     The outside local coefficient vector
   * @param[in]  lfsv_o  The outside test local function space
   * @param      r_i     The inside residual vector
   * @param      r_o     The outside residual vector
   *
   * @tparam     IG      The indersection
   * @tparam     LFSU    The trial local function space
   * @tparam     X       The local coefficient vector
   * @tparam     LFSV    The test local function space
   * @tparam     R       The residual vector
   */
  template<typename IG, typename LFSU, typename X, typename LFSV, typename R>
  void interface_alpha_skeleton(std::size_t domain_i, std::size_t domain_o,
                      const IG& ig,
                      const LFSU& lfsu_di,
                      const X& x_i,
                      const LFSV& lfsv_di,
                      const LFSU& lfsu_do,
                      const X& x_o,
                      const LFSV& lfsv_do,
                      R& r_i,
                      R& r_o) const
  {
    assert(lfsu_di.size() == lfsv_di.size());
    assert(lfsu_do.size() == lfsv_do.size());

    const auto& entity_f = ig.intersection();
    const auto& entity_i = ig.inside();
    const auto& entity_o = ig.outside();

    auto geo_f = entity_f.geometry();
    auto geo_in_i = entity_f.geometryInInside();
    auto geo_in_o = entity_f.geometryInOutside();

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

    const auto& coeff_mapper_i = _local_operator[domain_i]->coefficient_mapper_inside(entity_i);
    const auto& coeff_mapper_o = _local_operator[domain_o]->coefficient_mapper_outside(entity_o);

    const auto& local_basis_i = lfsu_di.child(0).finiteElement().localBasis();
    const auto& local_basis_o = lfsu_do.child(0).finiteElement().localBasis();

    using LocalBasis = std::decay_t<decltype(local_basis_i)>;
    using Range = typename LocalBasis::Traits::RangeType;
    using RF = typename LocalBasis::Traits::RangeFieldType;

    std::vector<Range> phiu_i(lfsu_di.size());
    std::vector<Range> phiu_o(lfsu_do.size());

    std::vector<RF> u_i(components_i), u_o(components_o);

    for (const auto& it : quadratureRule(geo_f, 3)) {
      const auto& position_f = it.position();
      // position of quadrature point in local coordinates of elements
      const auto position_i = geo_in_i.global(position_f);
      const auto position_o = geo_in_o.global(position_f);

      // evaluate basis functions
      phiu_i.clear(), phiu_o.clear();
      local_basis_i.evaluateFunction(position_i, phiu_i);
      local_basis_o.evaluateFunction(position_o, phiu_o);

      // evaluate concentrations at quadrature point

      u_i.clear(), u_o.clear();
      for (std::size_t k = 0; k < components_i; k++) // loop over components
        for (std::size_t j = 0; j < local_basis_i.size(); j++)
          u_i[k] += coeff_mapper_i(x_coeff_local_i, k, j) * phiu_i[j];

      for (std::size_t k = 0; k < components_o; k++) // loop over components
        for (std::size_t j = 0; j < local_basis_o.size(); j++)
          u_o[k] += coeff_mapper_o(x_coeff_local_o, k, j) * phiu_o[j];

      // integration factor
      auto factor = it.weight() * geo_f.integrationElement(position_f);

      for (std::size_t i = 0; i < lfsu_di.degree(); i++) // loop over components
      {
        if (lfsu_di.child(i).size() == 0)
          continue; // child space has nothing to compute
        const std::size_t comp_i =
          _local_operator[domain_i]->lfs_components(lfsu_di,lfsv_di)[i];
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
          _local_operator[domain_o]->lfs_components(lfsu_do,lfsv_do)[o];
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

  /**
   * @brief      The jacobian skeleton integral
   * @copydetails alpha_skeleton
   *
   * @param[in]  ig      The intersection
   * @param[in]  lfsu_i  The inside trial local function space
   * @param[in]  x_i     The inside local coefficient vector
   * @param[in]  lfsv_i  The inside test local function space
   * @param[in]  lfsu_o  The outside trial local function space
   * @param[in]  x_o     The outside local coefficient vector
   * @param[in]  lfsv_o  The outside test local function space
   * @param      mat_ii  The local inside-inside matrix
   * @param      mat_io  The local inside-outside matrix
   * @param      mat_oi  The local outside-inside matrix
   * @param      mat_oo  The local outside-outside matrix
   *
   * @tparam     IG      The indersection
   * @tparam     LFSU    The trial local function space
   * @tparam     X       The local coefficient vector
   * @tparam     LFSV    The test local function space
   * @tparam     J       The local jacobian matrix
   */
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
    assert(_size >= lfsu_i.degree());
    assert(_size >= lfsu_o.degree());
    assert(lfsu_i.degree() == lfsv_i.degree());
    assert(lfsu_o.degree() == lfsv_o.degree());

    const auto& entity_i = ig.inside();
    const auto& entity_o = ig.outside();

    auto domain_set_i = _index_set.subDomains(entity_i);
    auto domain_set_o = _index_set.subDomains(entity_o);

    assert(domain_set_i.size() == 1);
    assert(domain_set_o.size() == 1);

    std::size_t domain_i = *(domain_set_i.begin());
    std::size_t domain_o = *(domain_set_o.begin());

    const auto& lfsu_di = lfsu_i.child(domain_i);
    const auto& lfsv_di = lfsv_i.child(domain_i);
    const auto& lfsu_do = lfsu_o.child(domain_o);
    const auto& lfsv_do = lfsv_o.child(domain_o);

    assert(lfsu_di.size() == lfsv_di.size());
    assert(lfsu_do.size() == lfsv_do.size());

    if (domain_i == domain_o)
    {
      if constexpr (SubLocalOperator::doAlphaSkeleton)
      {
        _local_operator[domain_i]->jacobian_skeleton(ig,lfsu_di,x_i,lfsv_di,lfsu_do,x_o,lfsv_do,mat_ii,mat_io,mat_oi,mat_oo);
      }
    }
    else if constexpr (JM == JacobianMethod::Numerical)
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
    else
      interface_jacobian_skeleton(domain_i,domain_o,ig,lfsu_di,x_i,lfsv_di,lfsu_do,x_o,lfsv_do,mat_ii,mat_io,mat_oi,mat_oo);
  }

  /**
   * @brief      The jacobian skeleton integral
   * @copydetails alpha_skeleton
   *
   * @param[in]  ig      The intersection
   * @param[in]  lfsu_i  The inside trial local function space
   * @param[in]  x_i     The inside local coefficient vector
   * @param[in]  lfsv_i  The inside test local function space
   * @param[in]  lfsu_o  The outside trial local function space
   * @param[in]  x_o     The outside local coefficient vector
   * @param[in]  lfsv_o  The outside test local function space
   * @param      mat_ii  The local inside-inside matrix
   * @param      mat_io  The local inside-outside matrix
   * @param      mat_oi  The local outside-inside matrix
   * @param      mat_oo  The local outside-outside matrix
   *
   * @tparam     IG      The indersection
   * @tparam     LFSU    The trial local function space
   * @tparam     X       The local coefficient vector
   * @tparam     LFSV    The test local function space
   * @tparam     J       The local jacobian matrix
   */
  template<typename IG, typename LFSU, typename X, typename LFSV, typename J>
  void interface_jacobian_skeleton(std::size_t domain_i, std::size_t domain_o,
                         const IG& ig,
                         const LFSU& lfsu_di,
                         const X& x_i,
                         const LFSV& lfsv_di,
                         const LFSU& lfsu_do,
                         const X& x_o,
                         const LFSV& lfsv_do,
                         J& mat_ii,
                         J& mat_io,
                         J& mat_oi,
                         J& mat_oo) const
  {
    assert(lfsu_di.size() == lfsv_di.size());
    assert(lfsu_do.size() == lfsv_do.size());

    const auto& entity_f = ig.intersection();

    auto geo_f = entity_f.geometry();
    auto geo_in_i = entity_f.geometryInInside();
    auto geo_in_o = entity_f.geometryInOutside();

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

    const auto& local_basis_i = lfsu_di.child(0).finiteElement().localBasis();
    const auto& local_basis_o = lfsu_do.child(0).finiteElement().localBasis();

    using LocalBasis = std::decay_t<decltype(local_basis_i)>;
    using Range = typename LocalBasis::Traits::RangeType;

    std::vector<Range> phiu_i(lfsu_di.size());
    std::vector<Range> phiu_o(lfsu_do.size());

    for (const auto& it : quadratureRule(geo_f, 3)) {
      const auto& position_f = it.position();
      // position of quadrature point in local coordinates of elements
      const auto position_i = geo_in_i.global(position_f);
      const auto position_o = geo_in_o.global(position_f);

      // evaluate basis functions
      phiu_i.clear(), phiu_o.clear();
      local_basis_i.evaluateFunction(position_i, phiu_i);
      local_basis_o.evaluateFunction(position_o, phiu_o);
      // integration factor
      auto factor = it.weight() * geo_f.integrationElement(position_f);

      for (std::size_t i = 0; i < lfsu_di.degree(); i++) // loop over components
      {
        // has child space something to compute?
        if (lfsu_di.child(i).size() == 0)
          continue;
        const std::size_t comp_i =
          _local_operator[domain_i]->lfs_components(lfsu_di,lfsv_di)[i];
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
          const auto& lfs_comp_o = _local_operator[domain_o]->lfs_components(lfsu_do,lfsv_do);
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
          _local_operator[domain_o]->lfs_components(lfsu_do,lfsv_do)[o];
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
          const auto& lfs_comp_i = _local_operator[domain_i]->lfs_components(lfsu_di,lfsv_di);
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
};

/**
 * @brief      This class describes a PDELab temporal local operator for multi
 *             domain diffusion reaction systems.
 * @details    This class describre the operatrions for local integrals required
 *             for diffusion reaction system. The operator is only valid for
 *             entities contained in the grid. The local finite element is used
 *             for caching shape function evaluations.And the jacobian method
 *             switches between numerical and analytical jacobians. This local
 *             operator creates internally an individual local operator for
 *             every subdomain in the grid
 * @todo       Add numerical jacobian methods
 *
 * @tparam     Grid                The grid
 * @tparam     LocalFiniteElement  The local finite element
 * @tparam     JM                  The jacobian method
 */
template<class Grid,
         class SubLocalOperator,
         JacobianMethod JM = JacobianMethod::Analytical>
class TemporalLocalOperatorMultiDomainDiffusionReaction
  : public Dune::PDELab::LocalOperatorDefaultFlags
  , public Dune::PDELab::InstationaryLocalOperatorDefaultMethods<double>
{
  static_assert(Concept::isMultiDomainGrid<Grid>());

  using GridView = typename Grid::SubDomainGrid::LeafGridView;
  using SubLOP = SubLocalOperator;

  std::size_t _size;

  std::vector<std::shared_ptr<SubLOP>> _local_operator;

public:
  //! pattern assembly flags
  static constexpr bool doPatternVolume = true;

  //! residual assembly flags
  static constexpr bool doAlphaVolume = true;

  /**
   * @brief      Constructs a new instance.
   *
   * @param[in]  grid            The grid
   * @param[in]  config          The configuration
   * @param[in]  finite_element  The local finite element
   * @param[in]  id_operator     The index of this operator
   */
  TemporalLocalOperatorMultiDomainDiffusionReaction(
    std::shared_ptr<const Grid> grid,
    const ParameterTree& config,
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
      _local_operator[i] = std::make_shared<SubLOP>(
        sub_grid_view, sub_config, id_operator);
    }
  }

  /**
   * @copydoc TemporalLocalOperatorDiffusionReactionCG::pattern_volume
   */
  template<typename LFSU, typename LFSV, typename LocalPattern>
  void pattern_volume(const LFSU& lfsu,
                      const LFSV& lfsv,
                      LocalPattern& pattern) const
  {
    for (std::size_t i = 0; i < _size; ++i)
      _local_operator[i]->pattern_volume(lfsu.child(i), lfsv.child(i), pattern);
  }

  /**
   * @copydoc LocalOperatorMultiDomainDiffusionReaction::alpha_volume
   */
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
        _local_operator[i]->alpha_volume(eg, sub_lfsu, x, sub_lfsv, r);
      }
    }
  }

  /**
   * @copydoc LocalOperatorMultiDomainDiffusionReaction::jacobian_volume
   */
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
        _local_operator[i]->jacobian_volume(eg, sub_lfsu, x, sub_lfsv, mat);
      }
    }
  }

  /**
   * @copydoc LocalOperatorMultiDomainDiffusionReaction::jacobian_apply_volume
   */
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
        _local_operator[i]->jacobian_apply_volume(
          eg, sub_lfsu, x, z, sub_lfsv, r);
      }
    }
  }

  /**
   * @copydoc LocalOperatorMultiDomainDiffusionReaction::jacobian_apply_volume
   */
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
        _local_operator[i]->jacobian_apply_volume(eg, sub_lfsu, x, sub_lfsv, r);
      }
    }
  }
};

} // namespace Dune::Copasi

#endif // DUNE_COPASI_LOCAL_OPERATOR_MULTIDOMAIN_DIFFUSION_REACTION_HH