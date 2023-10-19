#ifndef DUNE_COPASI_LOCAL_OPERATOR_DIFFUSION_REACTION_CG_HH
#define DUNE_COPASI_LOCAL_OPERATOR_DIFFUSION_REACTION_CG_HH


#include <dune/copasi/finite_element/local_basis_cache.hh>
#include <dune/copasi/common/enum.hh>
#include <dune/copasi/local_operator/diffusion_reaction/base.hh>

#include <dune/assembler/common/tree_traversal.hh>

#include <dune/pdelab/localoperator/flags.hh>
#include <dune/pdelab/localoperator/idefault.hh>

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
 * @tparam     Space    Space
 * @tparam     LBT   Local basis traits
 */
template<Assembler::Concept::DiscreteFunctionSpace Space, class LBT>
class LocalOperatorDiffusionReactionCG
  : public Dune::PDELab::LocalOperatorDefaultFlags
  , public Dune::PDELab::InstationaryLocalOperatorDefaultMethods<double>
  , protected LocalOperatorDiffusionReactionBase<Space,LBT>
{
  using LOPBase = LocalOperatorDiffusionReactionBase<Space,LBT>;

  using RF = typename LOPBase::RF;
  using DF = typename LOPBase::DF;
  using LOPBase::_component_pattern;
  using LOPBase::dim;
  using LOPBase::_u;
  using LOPBase::_diffusion_gf;
  using LOPBase::_reaction_gf;
  using LOPBase::_jacobian_gf;

  mutable LocalBasisCache<LBT> _fe_cache;

  mutable std::vector<RF> _diffusion;
  mutable std::vector<RF> _reaction;
  mutable std::vector<RF> _jacobian;
  mutable std::vector<FieldVector<RF, dim>> _gradphi;

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
   * @param[in]  space
   * @param[in]  config          The configuration tree
   * @param[in]  finite_element  The local finite element
   */
  LocalOperatorDiffusionReactionCG(const Space& space,
                                   const ParameterTree& config)
    : LOPBase(space, config)
    , _fe_cache()
  {
    _diffusion.resize(_diffusion_gf.size());
    _reaction.resize(_diffusion_gf.size());
    _jacobian.resize(_diffusion_gf.size()*_diffusion_gf.size());
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

  template<class R, class X>
  struct PseudoJacobian {
    void accumulate(const auto& ltest, auto test_dof, const auto& ltrial, auto trial_dof, auto value) {
      _r.accumulate(ltest, test_dof, _z(ltrial, trial_dof) * value);
    }

    R& _r;
    const X& _z;
  };


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
                              R& r)
  {
    PseudoJacobian<R,X> mat{r,z};
    jacobian_volume(entity, lfsu, x, lfsv, mat);
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
                       M& mat)
  {
    if (lfsu.degree() == 0)
      return;

    // get geometry
    const auto& geo = entity.geometry();

    // TODO: check geometry type
    const auto& trial_finite_element = lfsu.child(0).finiteElement();

    _fe_cache.bind(trial_finite_element);
    _gradphi.resize(trial_finite_element.size());

    const auto& rule = _fe_cache.rule();

    assert(geo.affine());
    const FieldMatrix<DF, dim, dim>& jac = geo.jacobianInverse( rule[0].position() );

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
      for (std::size_t k = 0; k != lfsv.degree(); ++k)
        _diffusion_gf[k]->evaluate(entity, position, _diffusion[k]);

      // evaluate concentrations at quadrature point
      for (std::size_t comp = 0; comp != lfsu.degree(); ++comp)
        for (std::size_t dof = 0; dof != lfsu.child(comp).size(); ++dof)
          _u[comp] += x(lfsu.child(comp), dof) * phi[dof];

      for (std::size_t dof=0; dof != _gradphi.size(); ++dof)
        _gradphi[dof] = (jacphi[dof] * jac)[0]; // assert (1 x dim)

      RF factor = rule[q].weight() * geo.integrationElement(position);
      const auto& gradpsi = _gradphi;
      const auto& psi = phi;

      for (std::size_t comp_i = 0; comp_i != lfsv.degree(); ++comp_i) {
        // diffusion part
        for (std::size_t dof_i = 0; dof_i != lfsv.child(comp_i).size(); ++dof_i) {
          std::size_t comp_j = comp_i;
          for (std::size_t dof_j = 0; dof_j != lfsu.child(comp_j).size(); ++dof_j) {
            mat.accumulate(lfsv.child(comp_i), dof_i, lfsu.child(comp_j), dof_j, _diffusion[comp_i] * dot(_gradphi[dof_i], gradpsi[dof_j]) * factor);
          }
        }

        // reaction part
        for (std::size_t comp_j : _component_pattern[comp_i]) {
          // evaluate reaction term
          const auto jac_index = lfsv.degree() * comp_i + comp_j;
          _jacobian_gf[jac_index]->evaluate(entity, position, _jacobian[jac_index]);
          for (std::size_t dof_i = 0; dof_i != lfsv.child(comp_i).size(); ++dof_i)
            for (std::size_t dof_j = 0; dof_j != lfsu.child(comp_j).size(); ++dof_j)
              mat.accumulate(lfsv.child(comp_i), dof_i, lfsu.child(comp_j), dof_j, -_jacobian[jac_index] * phi[dof_i] * psi[dof_j] * factor);
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
  void alpha_volume(const EG& entity,
                    const LFSU& lfsu,
                    const X& x,
                    const LFSV& lfsv,
                    R& r)
  {
    if (lfsu.degree() == 0)
      return;

    // get geometry
    const auto& geo = entity.geometry();

    // TODO: check geometry type
    const auto& trial_finite_element = lfsu.child(0).finiteElement();

    _fe_cache.bind(trial_finite_element);
    _gradphi.resize(trial_finite_element.size());

    const auto& rule = _fe_cache.rule();

    assert(geo.affine());
    const FieldMatrix<DF, dim, dim>& jac = geo.jacobianInverse( rule[0].position() );

    // loop over quadrature points
    for (std::size_t q = 0; q != rule.size(); ++q) {
      const auto& position = rule[q].position();

      std::fill(_u.begin(), _u.end(), 0.);
      std::fill(_diffusion.begin(), _diffusion.end(), 0.);
      std::fill(_reaction.begin(), _reaction.end(), 0.);
      std::fill(_gradphi.begin(), _gradphi.end(), 0.);

      const auto& phi = _fe_cache.evaluateFunction(q);
      const auto& jacphi = _fe_cache.evaluateJacobian(q);

      for (std::size_t dof=0; dof != _gradphi.size(); ++dof)
        _gradphi[dof] = (jacphi[dof] * jac)[0]; // assert (1 x dim)

      // galerkin!
      const auto& psi = phi;
      const auto& gradpsi = _gradphi;

      // evaluate concentrations at quadrature point
      for (std::size_t comp = 0; comp != lfsu.degree(); ++comp)
        for (std::size_t dof = 0; dof != phi.size(); ++dof)
          _u[comp] += x(lfsu.child(comp), dof) * phi[dof];

      auto factor = rule[q].weight() * geo.integrationElement(position);

      // contribution for each component
      for (std::size_t comp = 0; comp != lfsv.degree(); ++comp) {
        _reaction_gf[comp]->evaluate(entity, position, _reaction[comp]);
        _diffusion_gf[comp]->evaluate(entity, position, _diffusion[comp]);
        // compute gradient u_h
        FieldVector<RF, dim> graduh(.0);
        for (std::size_t dof = 0; dof != _gradphi.size(); ++dof)
          graduh += _gradphi[dof] * x(lfsu.child(comp), dof);

        // scalar products
        for (std::size_t dof = 0; dof != psi.size(); ++dof) {
          r.accumulate(lfsv.child(comp), dof, (_diffusion[comp] * dot(gradpsi[dof], graduh) - _reaction[comp] * psi[dof]) * factor);
        }
      }
    }
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
template<Assembler::Concept::DiscreteFunctionSpace Space, class LBT>
class TemporalLocalOperatorDiffusionReactionCG
  : public Dune::PDELab::LocalOperatorDefaultFlags
  , public Dune::PDELab::InstationaryLocalOperatorDefaultMethods<double>
{
  using LocalBasisTraits = LBT;
  using RF = typename LocalBasisTraits::RangeFieldType;
  static constexpr int dim = LocalBasisTraits::dimDomain;


public:
  //! pattern assembly flags
  static constexpr bool doPatternVolume = true;

  //! residual assembly flags
  static constexpr bool doAlphaVolume = true;

 /**
   * @brief      Constructs a new instance.
   *
   * @todo       Make integration order variable depending on user requirements
   *             and polynomail order of the local finite element
   *
   * @param[in]  space
   *                             valid
   * @param[in]  config          The configuration tree
   * @param[in]  finite_element  The local finite element
   */
  TemporalLocalOperatorDiffusionReactionCG(
    const Space& space,
    const ParameterTree& config)
  {}

  template<class Entity, typename LFSU, typename LFSV, typename LocalPattern>
  void pattern_volume(const Entity& entity,
                      const LFSU& lfsu,
                      const LFSV& lfsv,
                      LocalPattern& pattern) const
  {
    assert(lfsv.degree() == lfsu.degree());
    for (std::size_t k = 0; k != lfsv.degree(); ++k)
      for (std::size_t i = 0; i != lfsu.child(k).size(); ++i)
        for (std::size_t j = 0; j != lfsv.child(k).size(); ++j)
          pattern.addLink(lfsv.child(k), i, lfsu.child(k), j);
  }

  void setTime(double t)
  {
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
                    R& r)
  {
    if (lfsu.degree() == 0)
      return;

    const auto& geo = entity.geometry();
    const auto& trial_finite_element = lfsu.child(0).finiteElement();
    _fe_cache.bind(trial_finite_element);
    const auto& rule = _fe_cache.rule();

    // loop over quadrature points
    for (std::size_t q = 0; q != rule.size(); ++q) {
      const auto& position = rule[q].position();
      auto factor = rule[q].weight() * geo.integrationElement(position);
      const auto& phi = _fe_cache.evaluateFunction(q);

      // contribution for each component
      for (std::size_t comp = 0; comp != lfsv.degree(); ++comp){
        double u = 0.0;
        for (std::size_t dof = 0; dof != lfsu.child(comp).size(); ++dof)
          u += x(lfsu.child(comp), dof) * phi[dof];

        for (std::size_t dof = 0; dof != lfsv.child(comp).size(); ++dof)
          r.accumulate(lfsv.child(comp), dof, u * phi[dof] * factor);
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
                       M& mat)
  {
    if (lfsu.degree() == 0)
      return;

    const auto& geo = eg.geometry();
    const auto& trial_finite_element = lfsu.child(0).finiteElement();
    _fe_cache.bind(trial_finite_element);
    const auto& rule = _fe_cache.rule();

    // loop over quadrature points
    for (std::size_t q = 0; q != rule.size(); ++q) {
      const auto& position = rule[q].position();

      auto factor = rule[q].weight() * geo.integrationElement(position);
      const auto& phi = _fe_cache.evaluateFunction(q);

      // integrate mass matrix
      for (std::size_t comp = 0; comp != lfsv.degree(); ++comp) {
        for (std::size_t dof_i = 0; dof_i != lfsv.child(comp).size(); ++dof_i)
          for (std::size_t dof_j = 0; dof_j != lfsu.child(comp).size(); ++dof_j)
            mat.accumulate(lfsv.child(comp), dof_i, lfsu.child(comp), dof_j, phi[dof_i] * phi[dof_j] * factor);
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
                             R& r)
  {
    alpha_volume(eg, lfsu, z, lfsv, r);
  }

private:
   mutable LocalBasisCache<LBT> _fe_cache;
};

} // namespace Dune::Copasi

#endif // DUNE_COPASI_LOCAL_OPERATOR_DIFFUSION_REACTION_CG_HH
