#ifndef DUNE_COPASI_LOCAL_OPERATOR_DIFFUSION_REACTION_HH
#define DUNE_COPASI_LOCAL_OPERATOR_DIFFUSION_REACTION_HH

#include <dune/copasi/common/coefficient_mapper.hh>
#include <dune/copasi/common/enum.hh>
#include <dune/copasi/common/pdelab_expression_adapter.hh>

namespace Dune::Copasi {

template<class GV,
         class CM = DefaultCoefficientMapper,
         JacobianMethod JM = JacobianMethod::Numerical>
class LocalOperatorDiffusionReactionFV
  : public Dune::PDELab::LocalOperatorDefaultFlags
  , public Dune::PDELab::InstationaryLocalOperatorDefaultMethods<double>
  , public PDELab::NumericalJacobianVolume<
      LocalOperatorDiffusionReactionFV<GV, LFE, CM, JM>>
  , public PDELab::NumericalJacobianApplyVolume<
      LocalOperatorDiffusionReactionFV<GV, LFE, CM, JM>>
{
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
      return;
    } else {
      _jacobian_apply_volume(eg, lfsu, x, x, lfsv, r);
    }
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

    _coefficient_mapper.bind(entity);

    DynamicVector<RF> u(_components);
    DynamicVector<RF> reaction(_lfs_components.size());

    std::fill(u.begin(), u.end(), 0.);
    std::fill(reaction.begin(), reaction.end(), 0.);

    // get center in local coordinates
    const auto ref_el = referenceElement(geo);
    const auto position = ref_el.position(0, 0);

    // get diffusion coefficient
    for (std::size_t k = 0; k < _lfs_components.size(); k++)
      _diffusion_gf[k]->evaluate(entity, position, diffusion[k]);

    // evaluate concentrations at quadrature point
    for (std::size_t k = 0; k < _components; k++)
      u[k] += _coefficient_mapper(x_coeff_local, k, 0);

    // get reaction term
    for (std::size_t k = 0; k < _lfs_components.size(); k++) {
      _reaction_gf[k]->update(u);
      _reaction_gf[k]->evaluate(entity, position, reaction[k]);
    }

    // contribution for each component
    for (std::size_t k = 0; k < _lfs_components.size(); k++) {
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

    assert(lfsu_i.degree() == 0);

    // get entity
    const auto entity = eg.entity();

    // get geometry
    const auto geo = eg.geometry();

    _coefficient_mapper.bind(entity);

    DynamicVector<RF> u(_components);
    DynamicVector<RF> jacobian(_lfs_components.size() * _lfs_components.size());

    // get center in local coordinates
    const auto ref_el = referenceElement(geo);
    const auto position = ref_el.position(0, 0);

    std::fill(u.begin(), u.end(), 0.);
    std::fill(jacobian.begin(), jacobian.end(), 0.);

    // evaluate concentrations at quadrature point
    for (std::size_t k = 0; k < _components; k++)
      u[k] += _coefficient_mapper(x_coeff_local, k, 0);

    // evaluate reaction term
    for (std::size_t k = 0; k < _lfs_components.size(); k++) {
      for (std::size_t l = 0; l < _lfs_components.size(); l++) {
        const auto j = _lfs_components.size() * k + l;
        _jacobian_gf[j]->update(u);
        _jacobian_gf[j]->evaluate(entity, position, jacobian[j]);
      }
    }

    auto do_link = [&](std::size_t comp_i, std::size_t comp_j) {
      auto it = _component_pattern.find(std::make_pair(comp_i, comp_j));
      return (it != _component_pattern.end());
    };

    for (std::size_t k = 0; k < _lfs_components.size(); k++) {
      for (std::size_t l = 0; l < _lfs_components.size(); l++) {
        if (not do_link(k, l))
          continue;
        const auto j = _lfs_components.size() * k + l;
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
    const auto& entity_f = ig.intersection();
    const auto& entity_i = ig.inside();
    const auto& entity_o = ig.outside();

    assert(lfsu_i.degree() == 0);

    auto geo_f = entity_f.geometry();
    auto geo_in_i = entity_f.geometryInInside();
    auto geo_in_o = entity_f.geometryInOutside();

    // get center in local coordinates
    const auto ref_el_f = referenceElement(geo);
    const auto position_f = ref_el.position(0, 0);

    // cell centers in references elements
    auto ref_el_i = referenceElement(geo_in_i);
    auto ref_el_o = referenceElement(geo_in_o);
    auto position_i = ref_el_i.position(0, 0);
    auto position_o = ref_el_o.position(0, 0);

    // cell centers in global coordinates
    auto position_g_i = geo_in_i.global(position_i);
    auto position_g_o = geo_in_o.global(position_o);

    // distance between the two cell centers
    position_g_i -= position_g_o;
    auto distance = position_g_i.two_norm();

    // get diffusion coefficient
    for (std::size_t k = 0; k < _lfs_components.size(); k++) {
      RF diffusion_i(0.), diffusion_o(0.);
      _diffusion_gf[k]->evaluate(entity_i, position_i, diffusion_i);
      _diffusion_gf[k]->evaluate(entity_o, position_o, diffusion_o);
      RF diffusion =
        2.0 / (1.0 / (diffusion_i + 1E-30) + 1.0 / (diffusion_o + 1E-30));
      RF gradu = coefficient_mapper_i(x_coeff_local_i, k, 0) -
                 coefficient_mapper_o(x_coeff_local_o, k, 0);
      gradu /= distance;
      // contribution to residual on inside element, other residual is computed
      // by symmetric call
      accumulate_i(k, 0, diffusion * gradu * geo_f.volume());
      accumulate_o(k, 0, diffusion * gradu * geo_f.volume());
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

    const auto& entity_f = ig.intersection();
    const auto& entity_i = ig.inside();
    const auto& entity_o = ig.outside();

    assert(lfsu_i.degree() == 0);

    auto geo_f = entity_f.geometry();
    auto geo_in_i = entity_f.geometryInInside();
    auto geo_in_o = entity_f.geometryInOutside();

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

    auto geo_f = entity_f.geometry();
    auto geo_in_i = entity_f.geometryInInside();
    auto geo_in_o = entity_f.geometryInOutside();

    // get center in local coordinates
    const auto ref_el_f = referenceElement(geo);
    const auto position_f = ref_el.position(0, 0);

    // cell centers in references elements
    auto ref_el_i = referenceElement(geo_in_i);
    auto ref_el_o = referenceElement(geo_in_o);
    auto position_i = ref_el_i.position(0, 0);
    auto position_o = ref_el_o.position(0, 0);

    // cell centers in global coordinates
    auto position_g_i = geo_in_i.global(position_i);
    auto position_g_o = geo_in_o.global(position_o);

    // distance between the two cell centers
    position_g_i -= position_g_o;
    auto distance = position_g_i.two_norm();

    // get diffusion coefficient
    for (std::size_t k = 0; k < _lfs_components.size(); k++) {
      RF diffusion_i(0.), diffusion_o(0.);
      _diffusion_gf[k]->evaluate(entity_i, position_i, diffusion_i);
      _diffusion_gf[k]->evaluate(entity_o, position_o, diffusion_o);
      RF diffusion =
        2.0 / (1.0 / (diffusion_i + 1E-30) + 1.0 / (diffusion_o + 1E-30));

      // contribution to residual on inside element, other residual is computed
      // by symmetric call
      accumulate_ii(k, 0, k, 0, diffusion * geo_f.volume());
      accumulate_io(k, 0, k, 0, -diffusion * geo_f.volume());
      accumulate_oo(k, 0, k, 0, diffusion * geo_f.volume());
      accumulate_oi(k, 0, k, 0, -diffusion * geo_f.volume());
    }
  }
};

} // namespace Dune::Copasi

#endif // DUNE_COPASI_LOCAL_OPERATOR_DIFFUSION_REACTION_HH