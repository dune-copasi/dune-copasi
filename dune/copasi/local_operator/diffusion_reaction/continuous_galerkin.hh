#ifndef DUNE_COPASI_LOCAL_OPERATOR_DIFFUSION_REACTION_CG_HH
#define DUNE_COPASI_LOCAL_OPERATOR_DIFFUSION_REACTION_CG_HH


#include <dune/copasi/finite_element/local_basis_cache.hh>
#include <dune/copasi/common/enum.hh>
#include <dune/copasi/local_operator/diffusion_reaction/base.hh>

#include <dune/assembler/common/tree_traversal.hh>

#include <dune/pdelab/localoperator/pattern.hh>
#include <dune/pdelab/localoperator/flags.hh>
#include <dune/pdelab/localoperator/idefault.hh>
#include <dune/pdelab/common/quadraturerules.hh>

#include <dune/geometry/referenceelements.hh>
#include <dune/geometry/type.hh>

#include <dune/common/fvector.hh>

namespace Dune::Copasi {

/**
 * @brief      This class describes a PDELab local operator for diffusion
 *             reaction systems.
 * @details    This class describre the operatrions for local integrals required
 *             for diffusion reaction system. The operator is only valid for
 *             entities contained in the entity set. The local finite element is
 *             used for caching shape function evaluations.
 *
 * @tparam     ES    Entity Set
 * @tparam     LBT   Local basis traits
 */
template<class ES, class LBT>
class LocalOperatorDiffusionReactionCG
  : public Dune::PDELab::LocalOperatorDefaultFlags
  , public Dune::PDELab::InstationaryLocalOperatorDefaultMethods<double>
  , protected LocalOperatorDiffusionReactionBase<ES,LBT>
{
  using EntitySet = ES;

  using LOPBase = LocalOperatorDiffusionReactionBase<ES,LBT>;

  using RF = typename LOPBase::RF;
  using DF = typename LOPBase::DF;
  using LOPBase::_component_pattern;
  using LOPBase::_components;
  using LOPBase::dim;
  using LOPBase::_diffusion_gf;
  using LOPBase::_reaction_gf;
  using LOPBase::_jacobian_gf;

  mutable LocalBasisCache<LBT> _fe_cache;
  mutable LocalBasisCache<LBT> _test_cache;

  mutable DynamicVector<RF> _u;
  mutable DynamicVector<RF> _diffusion;
  mutable DynamicVector<RF> _reaction;
  mutable DynamicVector<RF> _jacobian;
  mutable DynamicVector<FieldVector<RF, dim>> _gradphi;

public:

  //! selective assembly flags
  static constexpr bool doSkipEntity = false;
  static constexpr bool doSkipIntersection = false;

  //! pattern assembly flags
  static constexpr bool doPatternVolume = true;
  static constexpr bool doPatternSkeleton = false;

  //! residual assembly flags
  static constexpr bool doAlphaVolume = true;
  static constexpr bool doAlphaSkeleton = false;

  /**
   * @brief      Constructs a new instance.
   *
   * @todo       Make integration order variable depending on user requirements
   *             and polynomail order of the local finite element
   *
   * @param[in]  entity_set      The entity set where this local operator is valid
   * @param[in]  config          The configuration tree
   * @param[in]  finite_element  The local finite element
   */
  LocalOperatorDiffusionReactionCG(EntitySet entity_set,
                                   const ParameterTree& config)
    : LOPBase(config,entity_set)
    , _fe_cache()
  {
    _u.resize(_components);
    _diffusion.resize(_components);
    _reaction.resize(_components);
    _jacobian.resize(_components*_components);
  }

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
  template<class Entity, typename LFSU, typename LFSV, typename LocalPattern>
  void pattern_volume(const Entity& entity,
                      const LFSU& lfsu,
                      const LFSV& lfsv,
                      LocalPattern& pattern) const
  {
    for (std::size_t k = 0; k != lfsv.degree(); ++k)
      for (std::size_t l : _component_pattern[k])
        for (std::size_t i = 0; i != lfsu.child(k).size(); ++i)
          for (std::size_t j = 0; j != lfsv.child(l).size(); ++j)
            pattern.addLink(lfsv.child(k), i, lfsu.child(l), j);
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
  void jacobian_apply_volume(const EG& entity,
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
    auto z_coeff_local = [&](const std::size_t& component,
                             const std::size_t& dof) {
      return z(lfsu.child(component), dof);
    };

    auto accumulate = [&](const std::size_t& component,
                          const std::size_t& dof,
                          const auto& value) {
      r.accumulate(lfsu.child(component), dof, value);
    };

    // get geometry
    const auto& geo = entity.geometry();

    // TODO: check geometry type
    const auto& trial_finite_element = lfsu.child(0).finiteElement();

    _fe_cache.bind(trial_finite_element);
    _gradphi.resize(trial_finite_element.size());

    const auto& rule = _fe_cache.rule();

    assert(geo.affine());
    const FieldMatrix<DF, dim, dim>& jac = geo.jacobianInverseTransposed( rule[0].position());

    // loop over quadrature points
    for (std::size_t q = 0; q != rule.size(); ++q) {
      const auto& position = rule[q].position();

      std::fill(_u.begin(), _u.end(), 0.);
      std::fill(_diffusion.begin(), _diffusion.end(), 0.);
      std::fill(_reaction.begin(), _reaction.end(), 0.);
      std::fill(_gradphi.begin(), _gradphi.end(), 0.);

      const auto& phi = _fe_cache.evaluateFunction(q);
      const auto& jacphi = _fe_cache.evaluateJacobian(q);

      for (std::size_t i=0; i != _gradphi.size(); ++i)
        jac.mv(jacphi[i][0], _gradphi[i]);

      // galerkin!
      const auto& psi = phi;
      const auto& gradpsi = _gradphi;

      // get diffusion coefficient
      for (std::size_t k = 0; k != _components; ++k)
        _diffusion_gf[k]->evaluate(entity, position, _diffusion[k]);

      // evaluate concentrations at quadrature point
      for (std::size_t comp = 0; comp != _components; ++comp)
        for (std::size_t dof = 0; dof != phi.size(); ++dof)
          _u[comp] += x_coeff_local(comp, dof) * phi[dof];

      RF factor = rule[q].weight() * geo.integrationElement(position);

      // contribution for each component
      for (std::size_t k = 0; k < _components; k++) {
        // get reaction term
        _reaction_gf[k]->update(_u);
        _reaction_gf[k]->evaluate(entity, position, _reaction[k]);
        // compute gradient u_h
        FieldVector<RF, dim> graduh(.0);
        for (std::size_t d = 0; d != dim; ++d)
          for (std::size_t j = 0; j != _gradphi.size(); ++j)
            graduh[d] += _gradphi[j][d] * z_coeff_local(k, j);

        // scalar products
        for (std::size_t i = 0; i != psi.size(); ++i) // test func. loop
        {
          auto rhs = -_reaction[k] * psi[i];
          for (std::size_t d = 0; d != dim; ++d) // rows of grad
            rhs += _diffusion[k] * gradpsi[i][d] * graduh[d];
          accumulate(k, i, rhs * factor);
        }
      }
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
  void jacobian_volume(const EG& entity,
                       const LFSU& lfsu,
                       const X& x,
                       const LFSV& lfsv,
                       M& mat) const
  {
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

    // get geometry
    const auto& geo = entity.geometry();

    // TODO: check geometry type
    const auto& trial_finite_element = lfsu.child(0).finiteElement();

    _fe_cache.bind(trial_finite_element);
    _gradphi.resize(trial_finite_element.size());

    const auto& rule = _fe_cache.rule();

    assert(geo.affine());
    const FieldMatrix<DF, dim, dim>& jac = geo.jacobianInverseTransposed( rule[0].position());

    // loop over quadrature points
    for (std::size_t q = 0; q != rule.size(); ++q) {
      const auto& position = rule[q].position();

      const auto& phi = _fe_cache.evaluateFunction(q);
      const auto& jacphi = _fe_cache.evaluateJacobian(q);

      std::fill(_u.begin(), _u.end(), 0.);
      std::fill(_diffusion.begin(), _diffusion.end(), 0.);
      std::fill(_jacobian.begin(), _jacobian.end(), 0.);
      std::fill(_gradphi.begin(), _gradphi.end(), 0.);

      // get diffusion coefficient
      for (std::size_t k = 0; k < lfsv.degree(); k++)
        _diffusion_gf[k]->evaluate(entity, position, _diffusion[k]);

      // evaluate concentrations at quadrature point
      for (std::size_t comp = 0; comp < lfsu.degree(); comp++)
        for (std::size_t dof = 0; dof < lfsu.child(comp).size(); dof++) //  ansatz func. loop
          _u[comp] += x_coeff_local(comp, dof) * phi[dof];

      RF factor = rule[q].weight() * geo.integrationElement(position);

      for (std::size_t i=0; i != _gradphi.size(); ++i)
        jac.mv(jacphi[i][0], _gradphi[i]);

      const auto& gradpsi = _gradphi;
      const auto& psi = phi;

      for (std::size_t k = 0; k != lfsv.degree(); ++k) {
        // diffusion part
        for (std::size_t m = 0; m != lfsv.child(k).size(); ++m) {
          std::size_t l = k;
          for (std::size_t n = 0; n != lfsu.child(l).size(); ++n) {
            double ldiff = 0.;
            for (std::size_t d = 0; d < dim; d++)
              ldiff += _diffusion[k] * _gradphi[m][d] * gradpsi[n][d];
            accumulate(k, m, l, n, ldiff * factor);
          }
        }

        // reaction part
        for (std::size_t l : _component_pattern[k]) {
          const auto j = _components * k + l;
          // evaluate reaction term
          _jacobian_gf[j]->update(_u);
          _jacobian_gf[j]->evaluate(entity, position, _jacobian[j]);
          for (std::size_t m = 0; m != lfsv.child(k).size(); ++m) {
            for (std::size_t n = 0; n != lfsu.child(l).size(); ++n) {
              accumulate(k, m, l, n, -_jacobian[j] * phi[m] * psi[n] * factor);
            }
          }
        }
      }
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
    jacobian_apply_volume(eg, lfsu, x, x, lfsv, r);
  }
};

/**
 * @brief      This class describes a PDELab local operator for temporal part
 *             for diffusion reaction systems.
 * @details    This class describre the operatrions for local integrals required
 *             for diffusion reaction system. The operator is only valid for
 *             entities contained in the entity set. The local finite element is
 *             used for caching shape function evaluations.
 *
 * @tparam     ES    Entity Set
 * @tparam     LBT   Local Finite Element
 */
template<class ES, class LBT>
class TemporalLocalOperatorDiffusionReactionCG
  : public Dune::PDELab::LocalOperatorDefaultFlags
  , public Dune::PDELab::InstationaryLocalOperatorDefaultMethods<double>
  , protected LocalOperatorDiffusionReactionBase<ES,LBT>
{
  //! entity set
  using EntitySet = ES;

  using LOPBase = LocalOperatorDiffusionReactionBase<ES,LBT>;

  using RF = typename LOPBase::RF;
  using LOPBase::_components;
  using LOPBase::dim;


public:
  //! pattern assembly flags
  static constexpr bool doPatternVolume = true;

  //! residual assembly flags
  static constexpr bool doAlphaVolume = true;

  template<class Entity, typename LFSU, typename LFSV, typename LocalPattern>
  void pattern_volume(const Entity& entity,
                      const LFSU& lfsu,
                      const LFSV& lfsv,
                      LocalPattern& pattern) const
  {
    Assembler::forEachLeafNode(lfsv, [&](auto& lfsv_leaf, auto){
      Assembler::forEachLeafNode(lfsu, [&](auto& lfsu_leaf, auto){
        for (size_t i=0; i != lfsv_leaf.size(); ++i)
          for (size_t j=0; j != lfsu_leaf.size(); ++j)
            pattern.addLink(lfsv_leaf,i,lfsu_leaf,j);
      });
    });
  }


  /**
   * @brief      Constructs a new instance.
   *
   * @todo       Make integration order variable depending on user requirements
   *             and polynomail order of the local finite element
   *
   * @param[in]  entity_set       The entity set where this local operator is
   *                             valid
   * @param[in]  config          The configuration tree
   * @param[in]  finite_element  The local finite element
   */
  TemporalLocalOperatorDiffusionReactionCG(
    EntitySet entity_set,
    const ParameterTree& config)
    : LOPBase(config,entity_set)
  {}

  void setTime(double t)
  {
    Dune::PDELab::InstationaryLocalOperatorDefaultMethods<double>::setTime(t);
    LOPBase::setTime(t);
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
  void alpha_volume(const EG& entity,
                    const LFSU& lfsu,
                    const X& x,
                    const LFSV& lfsv,
                    R& r) const
  {
    auto x_coeff_local = [&](const std::size_t& component,
                             const std::size_t& dof) {
      return x(lfsu.child(component), dof);
    };
    auto accumulate = [&](const std::size_t& component,
                          const std::size_t& dof,
                          const auto& value) {
      r.accumulate(lfsu.child(component), dof, value);
    };

    const auto& geo = entity.geometry();

    // TODO: check geometry type
    const auto& trial_finite_element = lfsu.child(0).finiteElement();

    _fe_cache.bind(trial_finite_element);

    // TODO: check geometry type
    const auto& rule = _fe_cache.rule();

    // loop over quadrature points
    for (std::size_t q = 0; q != rule.size(); ++q) {
      const auto& position = rule[q].position();

      // get Jacobian and determinant
      RF factor = rule[q].weight() * geo.integrationElement(position);

      const auto& phi = _fe_cache.evaluateFunction(q);

      // contribution for each component
      for (std::size_t k = 0; k != lfsv.degree(); ++k) // loop over components
      {
        // compute value of component
        double u = 0.0;
        for (std::size_t j = 0; j < lfsu.child(k).size(); j++) // ansatz func. loop
          u += x_coeff_local(k, j) * phi[j];

        for (std::size_t i = 0; i < lfsv.child(k).size(); i++) // test func. loop
          accumulate(k, i, u * phi[i] * factor);
      }
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
    auto accumulate = [&](const std::size_t& component_i,
                          const std::size_t& dof_i,
                          const std::size_t& component_j,
                          const std::size_t& dof_j,
                          const auto& value) {
      mat.accumulate(
        lfsv.child(component_i), dof_i, lfsu.child(component_j), dof_j, value);
    };

    // get geometry
    const auto& geo = eg.geometry();

    // TODO: check geometry type
    const auto& trial_finite_element = lfsu.child(0).finiteElement();

    _fe_cache.bind(trial_finite_element);

    // TODO: check geometry type
    const auto& rule = _fe_cache.rule();

    // loop over quadrature points
    for (std::size_t q = 0; q != rule.size(); ++q) {
      const auto& position = rule[q].position();

      // get Jacobian and determinant
      RF factor = rule[q].weight() * geo.integrationElement(position);

      const auto& phi = _fe_cache.evaluateFunction(q);

      // integrate mass matrix
      for (std::size_t k = 0; k < _components; k++) // loop over components
      {
        for (std::size_t i = 0; i != lfsv.child(k).size(); ++i)
          for (std::size_t j = 0; j != lfsu.child(k).size(); ++j)
            accumulate(k, i, k, j, phi[i] * phi[j] * factor);
      }
    }
  }

  /**
   * @brief      The jacobian volume integral for matrix free operations
   * @details    This only switches between the actual implementation (in
   *             alpha_volume) and the numerical jacobian
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
    alpha_volume(eg, lfsu, z, lfsv, r);
  }

private:
   mutable LocalBasisCache<LBT> _fe_cache;
};

} // namespace Dune::Copasi

#endif // DUNE_COPASI_LOCAL_OPERATOR_DIFFUSION_REACTION_CG_HH
