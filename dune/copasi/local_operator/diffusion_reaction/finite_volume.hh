#ifndef DUNE_COPASI_LOCAL_OPERATOR_DIFFUSION_REACTION_FV_HH
#define DUNE_COPASI_LOCAL_OPERATOR_DIFFUSION_REACTION_FV_HH

#include <dune/copasi/common/coefficient_mapper.hh>
#include <dune/copasi/common/enum.hh>
#include <dune/copasi/common/pdelab_expression_adapter.hh>
#include <dune/copasi/local_operator/diffusion_reaction/base.hh>

#include <dune/pdelab/localoperator/numericaljacobian.hh>
#include <dune/pdelab/localoperator/numericaljacobianapply.hh>
#include <dune/pdelab/localoperator/pattern.hh>
#include <dune/pdelab/localoperator/flags.hh>
#include <dune/pdelab/localoperator/idefault.hh>

namespace Dune::Copasi {

// todo add the correct patterns
template<class GV,
         class LBT,
         class CM = DefaultCoefficientMapper,
         JacobianMethod JM = JacobianMethod::Analytical>
class LocalOperatorDiffusionReactionFV
  : public Dune::PDELab::LocalOperatorDefaultFlags
  , public Dune::PDELab::InstationaryLocalOperatorDefaultMethods<double>
  , protected LocalOperatorDiffusionReactionBase<GV,LBT,CM>
  , public PDELab::NumericalJacobianVolume<
      LocalOperatorDiffusionReactionFV<GV, LBT, CM, JM>>
  , public PDELab::NumericalJacobianApplyVolume<
      LocalOperatorDiffusionReactionFV<GV, LBT, CM, JM>>
  , public PDELab::NumericalJacobianSkeleton<
      LocalOperatorDiffusionReactionFV<GV, LBT, CM, JM>>
  , public PDELab::NumericalJacobianApplySkeleton<
      LocalOperatorDiffusionReactionFV<GV, LBT, CM, JM>>
  , public PDELab::FullSkeletonPattern
{
  using GridView = GV;
  using LOPBase = LocalOperatorDiffusionReactionBase<GV,LBT,CM>;

  using RF = typename LOPBase::RF;
  using DF = typename LOPBase::DF;
  using LOPBase::_component_pattern;
  using LOPBase::_components;
  using LOPBase::dim;
  using LOPBase::_diffusion_gf;
  using LOPBase::_reaction_gf;
  using LOPBase::_jacobian_gf;

public:

  //! pattern assembly flags
  static constexpr bool doPatternVolume = true;

  //! residual assembly flags
  static constexpr bool doAlphaVolume = true;

  //! pattern assembly flags
  static constexpr bool doPatternSkeleton = true;

  //! residual assembly flags
  static constexpr bool doAlphaSkeleton = true;

  using LOPBase::lfs_components;

  using LOPBase::coefficient_mapper_inside;
  using LOPBase::coefficient_mapper_outside;
  using LOPBase::update;

  LocalOperatorDiffusionReactionFV(GridView grid_view,
                                   const ParameterTree& config,
                                   std::size_t id_operator)
    : LOPBase(config,id_operator,grid_view)
  {}



  void setTime(double t)
  {
    Dune::PDELab::InstationaryLocalOperatorDefaultMethods<double>::setTime(t);
    LOPBase::setTime(t);
  }

  /**
   * @brief      Pattern volume
   * @details    This method links degrees of freedom between trial and test
   *             spaces taking into account the structure of the reaction term
   *
   * @param[in]  lfsu          The trial local function space
   * @param[in]  lfsv          The test local function space
   * @param      pattern       The local pattern
   *
   * @tparam     LFSU          The trial local function space
   * @tparam     LFSV          The test local function space
   * @tparam     LocalPattern  The local pattern
   */
  template<typename LFSU, typename LFSV, typename LocalPattern>
  void pattern_volume(const LFSU& lfsu,
                      const LFSV& lfsv,
                      LocalPattern& pattern) const
  {
    // assert(this->cached());
    auto do_link = [&](std::size_t comp_i, std::size_t comp_j) {
      auto it = _component_pattern.find(std::make_pair(comp_i, comp_j));
      return (it != _component_pattern.end());
    };
    for (std::size_t i = 0; i < lfsv.degree(); ++i)
      // if (this->cached(lfsu.child(i).finiteElement()) and this->cached(lfsv.child(i).finiteElement()))
        for (std::size_t j = 0; j < lfsu.degree(); ++j)
          if (do_link(i, j))
            for (std::size_t k = 0; k < lfsv.child(i).size(); ++k)
              for (std::size_t l = 0; l < lfsu.child(j).size(); ++l)
                pattern.addLink(lfsv.child(i), k, lfsu.child(j), l);
  }

  /**
   * @brief      The jacobian volume integral for matrix free operations
   * @details    This only switches between the actual implementation (in
   *             _jacobian_apply_volume) and the numerical jacobian
   *
   * @param[in]  eg    The entity
   * @param[in]  lfsu  The trial local function space
   * @param[in]  x     The local coefficient vector
   * @param[in]  z     The local position in the trial space to which to apply
   *                   the Jacobian.
   * @param[in]  lfsv  The test local function space
   * @param      r     The resulting vector
   *
   * @tparam     EG    The entity
   * @tparam     LFSU  The trial local function space
   * @tparam     X     The local coefficient vector type
   * @tparam     LFSV  The test local function space
   * @tparam     R     The resulting vector
   */
  template<typename EG, typename LFSU, typename X, typename LFSV, typename R>
  void jacobian_apply_volume(const EG& eg,
                             const LFSU& lfsu,
                             const X& x,
                             const X& z,
                             const LFSV& lfsv,
                             R& r) const
  {
    if constexpr (JM == JacobianMethod::Numerical) {
      PDELab::NumericalJacobianApplyVolume<LocalOperatorDiffusionReactionFV>::
        jacobian_apply_volume(eg, lfsu, x, z, lfsv, r);
    } else {
      _jacobian_apply_volume(eg, lfsu, x, z, lfsv, r);
    }
  }

  /**
   * @brief      The jacobian volume integral for matrix free operations (linear
   *             variant)
   * @details    This only switches between the actual implementation (in
   *             _jacobian_apply_volume) and the numerical jacobian
   *
   * @param[in]  eg    The entity
   * @param[in]  lfsu  The trial local function space
   * @param[in]  x     The local coefficient vector
   * @param[in]  lfsv  The test local function space
   * @param      r     The resulting vector
   *
   * @tparam     EG    The entity
   * @tparam     LFSU  The trial local function space
   * @tparam     X     The local coefficient vector type
   * @tparam     LFSV  The test local function space
   * @tparam     R     The resulting vector
   */
  template<typename EG, typename LFSU, typename X, typename LFSV, typename R>
  void jacobian_apply_volume(const EG& eg,
                             const LFSU& lfsu,
                             const X& x,
                             const LFSV& lfsv,
                             R& r) const
  {
    if constexpr (JM == JacobianMethod::Numerical) {
      PDELab::NumericalJacobianApplyVolume<LocalOperatorDiffusionReactionFV>::
        jacobian_apply_volume(eg, lfsu, x, lfsv, r);
    } else {
      _jacobian_apply_volume(eg, lfsu, x, x, lfsv, r);
    }
  }

  /**
   * @brief      The volume integral
   *
   * @param[in]  eg    The entity
   * @param[in]  lfsu  The trial local function space
   * @param[in]  x     The local coefficient vector
   * @param[in]  lfsv  The test local function space
   * @param      r     The local residual vector
   *
   * @tparam     EG    The entity
   * @tparam     LFSU  The trial local function space
   * @tparam     X     The local coefficient vector type
   * @tparam     LFSV  The test local function space
   * @tparam     R     The local residual vector
   */
  template<typename EG, typename LFSU, typename X, typename LFSV, typename R>
  void alpha_volume(const EG& eg,
                    const LFSU& lfsu,
                    const X& x,
                    const LFSV& lfsv,
                    R& r) const
  {
    _jacobian_apply_volume(eg, lfsu, x, x, lfsv, r);
  }

  /**
   * @brief      The jacobian volume integral for matrix free operations
   *
   * @param[in]  eg    The entity
   * @param[in]  lfsu  The trial local function space
   * @param[in]  x     The local coefficient vector
   * @param[in]  z     The local position in the trial space to which to apply
   *                   the Jacobian.
   * @param[in]  lfsv  The test local function space
   * @param      r     The resulting vector
   *
   * @tparam     EG    The entity
   * @tparam     LFSU  The trial local function space
   * @tparam     X     The local coefficient vector type
   * @tparam     LFSV  The test local function space
   * @tparam     R     The resulting vector
   */
  template<typename EG, typename LFSU, typename X, typename LFSV, typename R>
  void _jacobian_apply_volume(const EG& eg,
                              const LFSU& lfsu,
                              const X& x,
                              const X& z,
                              const LFSV& lfsv,
                              R& r) const
  {
    // assume we receive a power local finite element!
    auto x_coeff_local = [&](const std::size_t& component,
                             const std::size_t& dof) {
      return x(lfsu.child(component), dof);
    };

    auto accumulate = [&](const std::size_t& component,
                          const std::size_t& dof,
                          const auto& value) {
      r.accumulate(lfsu.child(component), dof, value);
    };

    // get entity
    const auto entity = eg.entity();

    // get geometry
    const auto geo = eg.geometry();

    const auto& coeff_mapper = this->coefficient_mapper_inside(entity);
    const auto& lfs_components = this->lfs_components(lfsu,lfsv);

    DynamicVector<RF> u(_components);
    DynamicVector<RF> reaction(lfs_components.size());

    std::fill(u.begin(), u.end(), 0.);
    std::fill(reaction.begin(), reaction.end(), 0.);

    // get center in local coordinates
    const auto ref_el = referenceElement(geo);
    const auto position = ref_el.position(0, 0);

    // evaluate concentrations at quadrature point
    for (std::size_t k = 0; k < _components; k++)
      u[k] += coeff_mapper(x_coeff_local, k, 0);

    // get reaction term
    for (std::size_t k = 0; k < lfs_components.size(); k++) {
      _reaction_gf[k]->update(u);
      _reaction_gf[k]->evaluate(entity, position, reaction[k]);
    }

    // contribution for each component
    for (std::size_t k = 0; k < lfs_components.size(); k++) {
      accumulate(k, 0, -reaction[k] * geo.volume());
    }
  }

  /**
   * @brief      The jacobian volume integral
   *
   * @param[in]  eg    The entity
   * @param[in]  lfsu  The trial local function space
   * @param[in]  x     The local coefficient vector
   * @param[in]  lfsv  The test local function space
   * @param      mat   The local matrix
   *
   * @tparam     EG    The entity
   * @tparam     LFSU  The trial local function space
   * @tparam     X     The local coefficient vector type
   * @tparam     LFSV  The test local function space
   * @tparam     M     The local matrix
   */
  template<typename EG, typename LFSU, typename X, typename LFSV, typename M>
  void jacobian_volume(const EG& eg,
                       const LFSU& lfsu,
                       const X& x,
                       const LFSV& lfsv,
                       M& mat) const
  {
    if constexpr (JM == JacobianMethod::Numerical) {
      PDELab::NumericalJacobianVolume<LocalOperatorDiffusionReactionFV>::
        jacobian_volume(eg, lfsu, x, lfsv, mat);
      return;
    }
    // assume we receive a power local finite element!
    auto x_coeff_local = [&](const std::size_t& component,
                             const std::size_t& dof) {
      return x(lfsu.child(component), dof);
    };

    auto accumulate = [&](const std::size_t& component_i,
                          const std::size_t& dof_i,
                          const std::size_t& component_j,
                          const std::size_t& dof_j,
                          const auto& value) {
      mat.accumulate(
        lfsv.child(component_i), dof_i, lfsu.child(component_j), dof_j, value);
    };

    // get entity
    const auto entity = eg.entity();

    // get geometry
    const auto geo = eg.geometry();

    const auto& coeff_mapper = this->coefficient_mapper_inside(entity);
    const auto& lfs_components = this->lfs_components(lfsu,lfsv);

    DynamicVector<RF> u(_components);
    DynamicVector<RF> jacobian(lfs_components.size() * lfs_components.size());

    // get center in local coordinates
    const auto ref_el = referenceElement(geo);
    const auto position = ref_el.position(0, 0);

    std::fill(u.begin(), u.end(), 0.);
    std::fill(jacobian.begin(), jacobian.end(), 0.);

    // evaluate concentrations at quadrature point
    for (std::size_t k = 0; k < _components; k++)
      u[k] += coeff_mapper(x_coeff_local, k, 0);

    // evaluate reaction term
    for (std::size_t k = 0; k < lfs_components.size(); k++) {
      for (std::size_t l = 0; l < lfs_components.size(); l++) {
        const auto j = lfs_components.size() * k + l;
        _jacobian_gf[j]->update(u);
        _jacobian_gf[j]->evaluate(entity, position, jacobian[j]);
      }
    }

    auto do_link = [&](std::size_t comp_i, std::size_t comp_j) {
      auto it = _component_pattern.find(std::make_pair(comp_i, comp_j));
      return (it != _component_pattern.end());
    };

    for (std::size_t k = 0; k < lfs_components.size(); k++) {
      for (std::size_t l = 0; l < lfs_components.size(); l++) {
        if (not do_link(k, l))
          continue;
        const auto j = lfs_components.size() * k + l;
        accumulate(k, 0, l, 0, -jacobian[j] * geo.volume());
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
    auto x_coeff_local_i = [&](const std::size_t& component,
                               const std::size_t& dof) {
      return x_i(lfsu_i.child(component), dof);
    };
    auto x_coeff_local_o = [&](const std::size_t& component,
                               const std::size_t& dof) {
      return x_o(lfsu_o.child(component), dof);
    };

    auto accumulate_i = [&](const std::size_t& component,
                            const std::size_t& dof,
                            const auto& value) {
      r_i.accumulate(lfsu_i.child(component), dof, value);
    };

    auto accumulate_o = [&](const std::size_t& component,
                            const std::size_t& dof,
                            const auto& value) {
      r_o.accumulate(lfsu_o.child(component), dof, value);
    };

    // get cell entities from both sides of the intersection
    const auto& entity_i = ig.inside();
    const auto& entity_o = ig.outside();
    const auto& entity_f = ig.intersection();

    const auto& lfs_components_i = this->lfs_components(lfsu_i,lfsv_i);
    [[maybe_unused]] const auto& lfs_components_o = this->lfs_components(lfsu_o,lfsv_o);

    assert(&lfs_components_i == &lfs_components_o); // this lop only works for same lfs componets in the skeleton
    assert(lfsu_i.degree() == lfs_components_i.size());
    assert(lfsu_o.degree() == lfs_components_o.size());

    auto geo_i = entity_i.geometry();
    auto geo_o = entity_o.geometry();
    auto geo_f = entity_f.geometry();
    auto geo_in_i = entity_f.geometryInInside();
    auto geo_in_o = entity_f.geometryInOutside();

    const auto& coeff_mapper_i = this->coefficient_mapper_inside(entity_i);
    const auto& coeff_mapper_o = this->coefficient_mapper_outside(entity_o);

    // cell centers in codim 1 reference element (facet)
    auto ref_el_i = referenceElement(geo_in_i);
    auto ref_el_o = referenceElement(geo_in_o);
    auto position_f_i = ref_el_i.position(0, 0);
    auto position_f_o = ref_el_o.position(0, 0);

    // cell centers in codim 0 reference element (cube)
    auto position_i = geo_in_i.global(position_f_i);
    auto position_o = geo_in_o.global(position_f_o);

    // cell centers in global coordinates reference element (cube)
    auto position_g_i = geo_i.center();
    auto position_g_o = geo_o.center();
    // distance between the two cell centers
    position_g_o -= position_g_i;
    auto distance = position_g_o.two_norm();

    // get diffusion coefficient
    for (std::size_t k = 0; k < lfs_components_i.size(); k++) {
      RF diffusion_i(0.), diffusion_o(0.);
      _diffusion_gf[k]->evaluate(entity_i, position_i, diffusion_i);
      _diffusion_gf[k]->evaluate(entity_o, position_o, diffusion_o);
      RF diffusion =
        2.0 / (1.0 / (diffusion_i + 1E-30) + 1.0 / (diffusion_o + 1E-30));
      RF gradu = coeff_mapper_o(x_coeff_local_o, k, 0) -
                 coeff_mapper_i(x_coeff_local_i, k, 0);
      gradu /= distance;
      // contribution to residual on inside element, other residual is computed
      // by symmetric call
      accumulate_i(k, 0, - diffusion * gradu * geo_f.volume());
      accumulate_o(k, 0,   diffusion * gradu * geo_f.volume());
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
    if constexpr (JM == JacobianMethod::Numerical) {
      PDELab::NumericalJacobianSkeleton<
        LocalOperatorDiffusionReactionFV>::jacobian_skeleton(ig,
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

    auto accumulate_ii = [&](const std::size_t& component_i,
                             const std::size_t& dof_i,
                             const std::size_t& component_j,
                             const std::size_t& dof_j,
                             const auto& value) {
      mat_ii.accumulate(lfsu_i.child(component_i),
                        dof_i,
                        lfsu_i.child(component_j),
                        dof_j,
                        value);
    };

    auto accumulate_io = [&](const std::size_t& component_i,
                             const std::size_t& dof_i,
                             const std::size_t& component_j,
                             const std::size_t& dof_j,
                             const auto& value) {
      mat_io.accumulate(lfsu_i.child(component_i),
                        dof_i,
                        lfsu_o.child(component_j),
                        dof_j,
                        value);
    };

    auto accumulate_oi = [&](const std::size_t& component_i,
                             const std::size_t& dof_i,
                             const std::size_t& component_j,
                             const std::size_t& dof_j,
                             const auto& value) {
      mat_oi.accumulate(lfsu_o.child(component_i),
                        dof_i,
                        lfsu_i.child(component_j),
                        dof_j,
                        value);
    };

    auto accumulate_oo = [&](const std::size_t& component_i,
                             const std::size_t& dof_i,
                             const std::size_t& component_j,
                             const std::size_t& dof_j,
                             const auto& value) {
      mat_oo.accumulate(lfsu_o.child(component_i),
                        dof_i,
                        lfsu_o.child(component_j),
                        dof_j,
                        value);
    };
    // get cell entities from both sides of the intersection
    const auto& entity_i = ig.inside();
    const auto& entity_o = ig.outside();
    const auto& entity_f = ig.intersection();

    const auto& lfs_components_i = this->lfs_components(lfsu_i,lfsv_i);
    [[maybe_unused]] const auto& lfs_components_o = this->lfs_components(lfsu_o,lfsv_o);

    assert(&lfs_components_i == &lfs_components_o); // this lop only works for same lfs componets in the skeleton
    assert(lfsu_i.degree() == lfs_components_i.size());
    assert(lfsu_o.degree() == lfs_components_o.size());

    auto geo_i = entity_i.geometry();
    auto geo_o = entity_o.geometry();
    auto geo_f = entity_f.geometry();
    auto geo_in_i = entity_f.geometryInInside();
    auto geo_in_o = entity_f.geometryInOutside();

    // cell centers in codim 1 reference element (facet)
    auto ref_el_i = referenceElement(geo_in_i);
    auto ref_el_o = referenceElement(geo_in_o);
    auto position_f_i = ref_el_i.position(0, 0);
    auto position_f_o = ref_el_o.position(0, 0);

    // cell centers in codim 0 reference element (cube)
    auto position_i = geo_in_i.global(position_f_i);
    auto position_o = geo_in_o.global(position_f_o);

    // cell centers in global coordinates reference element (cube)
    auto position_g_i = geo_i.center();
    auto position_g_o = geo_o.center();
    // distance between the two cell centers
    position_g_o -= position_g_i;
    auto distance = position_g_o.two_norm();

    // get diffusion coefficient
    for (std::size_t k = 0; k < lfs_components_i.size(); k++) {
      RF diffusion_i(0.), diffusion_o(0.);
      _diffusion_gf[k]->evaluate(entity_i, position_i, diffusion_i);
      _diffusion_gf[k]->evaluate(entity_o, position_o, diffusion_o);
      RF diffusion =
        2.0 / (1.0 / (diffusion_i + 1E-30) + 1.0 / (diffusion_o + 1E-30));

      // contribution to residual on inside element, other residual is computed
      // by symmetric call
      accumulate_ii(k, 0, k, 0,   diffusion * geo_f.volume() / distance);
      accumulate_io(k, 0, k, 0, - diffusion * geo_f.volume() / distance);
      accumulate_oo(k, 0, k, 0,   diffusion * geo_f.volume() / distance);
      accumulate_oi(k, 0, k, 0, - diffusion * geo_f.volume() / distance);
    }
  }
};

} // namespace Dune::Copasi

#endif // DUNE_COPASI_LOCAL_OPERATOR_DIFFUSION_REACTION_FV_HH