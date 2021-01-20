#ifndef DUNE_COPASI_LOCAL_OPERATOR_DIFFUSION_REACTION_CG_HH
#define DUNE_COPASI_LOCAL_OPERATOR_DIFFUSION_REACTION_CG_HH

#include <dune/copasi/finite_element/local_basis_cache.hh>
#include <dune/copasi/common/enum.hh>
#include <dune/copasi/local_operator/diffusion_reaction/base.hh>

#include <dune/pdelab/localoperator/numericaljacobian.hh>
#include <dune/pdelab/localoperator/numericaljacobianapply.hh>
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
 *             entities contained in the grid view. The local finite element is
 *             used for caching shape function evaluations. And the jacobian
 *             method switches between numerical and analytical jacobians.
 *
 * @tparam     GV    Grid View
 * @tparam     LBT   Local basis traits
 * @tparam     JM    Jacobian Method
 */
template<class GV,
         class LBT,
         JacobianMethod JM = JacobianMethod::Analytical>
class LocalOperatorDiffusionReactionCG
  : public Dune::PDELab::LocalOperatorDefaultFlags
  , public Dune::PDELab::InstationaryLocalOperatorDefaultMethods<double>
  , protected LocalOperatorDiffusionReactionBase<GV,LBT>
  , public PDELab::NumericalJacobianVolume<
      LocalOperatorDiffusionReactionCG<GV, LBT, JM>>
  , public PDELab::NumericalJacobianApplyVolume<
      LocalOperatorDiffusionReactionCG<GV, LBT, JM>>
{
  using GridView = GV;

  using LOPBase = LocalOperatorDiffusionReactionBase<GV,LBT>;

  using RF = typename LOPBase::RF;
  using DF = typename LOPBase::DF;
  using LOPBase::_component_pattern;
  using LOPBase::_components;
  using LOPBase::dim;
  using LOPBase::_diffusion_gf;
  using LOPBase::_reaction_gf;
  using LOPBase::_jacobian_gf;

  mutable LocalBasisCache<LBT> _trial_cache;
  mutable LocalBasisCache<LBT> _test_cache;
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
   * @param[in]  grid_view       The grid view where this local operator is
   *                             valid
   * @param[in]  config          The configuration tree
   * @param[in]  finite_element  The local finite element
   */
  LocalOperatorDiffusionReactionCG(GridView grid_view,
                                   const ParameterTree& config)
    : LOPBase(config,grid_view)
    , _trial_cache()
    , _test_cache(_trial_cache)
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
    auto do_link = [&](std::size_t comp_i, std::size_t comp_j) {
      auto it = _component_pattern.find(std::make_pair(comp_i, comp_j));
      return (it != _component_pattern.end());
    };
    for (std::size_t i = 0; i < lfsv.degree(); ++i)
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
      PDELab::NumericalJacobianApplyVolume<LocalOperatorDiffusionReactionCG>::
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
      PDELab::NumericalJacobianApplyVolume<LocalOperatorDiffusionReactionCG>::
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
    auto z_coeff_local = [&](const std::size_t& component,
                             const std::size_t& dof) {
      return z(lfsu.child(component), dof);
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

    // TODO: check geometry type
    const auto& rule = quadratureRule(geo,3);
    const auto& trial_finite_element = lfsu.child(0).finiteElement();
    const auto& test_finite_element = lfsv.child(0).finiteElement();
    const auto& trial_basis = trial_finite_element.localBasis();
    const auto& test_basis = test_finite_element.localBasis();

    _trial_cache.bind(trial_finite_element);
    _test_cache.bind(test_finite_element);

    DynamicVector<RF> u(_components);
    DynamicVector<RF> diffusion(_components);
    DynamicVector<RF> reaction(_components);
    DynamicVector<FieldVector<RF, dim>> gradphi(trial_basis.size());
    DynamicVector<FieldVector<RF, dim>> gradpsi(test_basis.size());

    // loop over quadrature points
    for (const auto& point : rule) {
      const auto& position = point.position();

      std::fill(u.begin(), u.end(), 0.);
      std::fill(diffusion.begin(), diffusion.end(), 0.);
      std::fill(reaction.begin(), reaction.end(), 0.);
      std::fill(gradphi.begin(), gradphi.end(), 0.);
      std::fill(gradpsi.begin(), gradpsi.end(), 0.);

      const auto& phi = _trial_cache.evaluateFunction(position);
      const auto& jacphi = _trial_cache.evaluateJacobian(position);

      const auto& psi = _trial_cache.evaluateFunction(position);
      const auto& jacpsi = _trial_cache.evaluateJacobian(position);

      FieldMatrix<DF, dim, dim> jac = geo.jacobianInverseTransposed(position);

      for (std::size_t i=0; i<gradphi.size(); i++)
        jac.mv(jacphi[i][0],gradphi[i]);

      for (std::size_t i=0; i<gradpsi.size(); i++)
        jac.mv(jacpsi[i][0],gradpsi[i]);

      // get diffusion coefficient
      for (std::size_t k = 0; k < _components; k++)
        _diffusion_gf[k]->evaluate(entity, position, diffusion[k]);

      // evaluate concentrations at quadrature point
      for (std::size_t comp = 0; comp < _components; comp++)
        for (std::size_t dof = 0; dof < phi.size(); dof++)
          u[comp] += x_coeff_local(comp, dof) * phi[dof];

      RF factor = point.weight() * geo.integrationElement(position);

      // contribution for each component
      for (std::size_t k = 0; k < _components; k++) {
        // get reaction term
        _reaction_gf[k]->update(u);
        _reaction_gf[k]->evaluate(entity, position, reaction[k]);
        // compute gradient u_h
        FieldVector<RF, dim> graduh(.0);
        for (std::size_t d = 0; d < dim; d++)
          for (std::size_t j = 0; j < gradphi.size(); j++)
            graduh[d] += gradphi[j][d] * z_coeff_local(k, j);

        // scalar products
        for (std::size_t i = 0; i < psi.size(); i++) // test func. loop
        {
          typename R::value_type rhs = -reaction[k] * psi[i];
          for (std::size_t d = 0; d < dim; d++) // rows of grad
            rhs += diffusion[k] * gradpsi[i][d] * graduh[d];
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
  void jacobian_volume(const EG& eg,
                       const LFSU& lfsu,
                       const X& x,
                       const LFSV& lfsv,
                       M& mat) const
  {
    if constexpr (JM == JacobianMethod::Numerical) {
      PDELab::NumericalJacobianVolume<LocalOperatorDiffusionReactionCG>::
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

    const auto& rule = QuadratureRules<RF, dim>::rule(geo.type(),3);

    const auto& trial_basis = lfsu.child(0).finiteElement().localBasis();
    const auto basis_size = trial_basis.size();
    std::vector<typename LBT::RangeType> phi(basis_size);
    std::vector<typename LBT::JacobianType> jac(basis_size);

    DynamicVector<RF> u(_components);
    DynamicVector<RF> diffusion(_components);
    DynamicVector<RF> jacobian(_components * _components);
    DynamicVector<FieldVector<RF, dim>> grad(basis_size);

    // loop over quadrature points
    for (const auto& point : rule) {
      const auto& position = point.position();

      trial_basis.evaluateFunction(position, phi);
      trial_basis.evaluateJacobian(position, jac);

      std::fill(u.begin(), u.end(), 0.);
      std::fill(diffusion.begin(), diffusion.end(), 0.);
      std::fill(jacobian.begin(), jacobian.end(), 0.);
      std::fill(grad.begin(), grad.end(), FieldVector<RF, dim>(0.));

      // get diffusion coefficient
      for (std::size_t k = 0; k < _components; k++)
        _diffusion_gf[k]->evaluate(entity, position, diffusion[k]);

      // evaluate concentrations at quadrature point
      for (std::size_t comp = 0; comp < _components; comp++)
        for (std::size_t dof = 0; dof < basis_size; dof++) //  ansatz func. loop
          u[comp] += x_coeff_local(comp, dof) * phi[dof];

      // get jacobian and determinant
      FieldMatrix<DF, dim, dim> S = geo.jacobianInverseTransposed(position);
      RF factor = point.weight() * geo.integrationElement(position);

      // compute gradients of basis functions in transformed element
      // (independent of component)
      for (std::size_t i = 0; i < dim; i++)             // rows of S
        for (std::size_t k = 0; k < dim; k++)           // columns of S
          for (std::size_t j = 0; j < basis_size; j++) // columns of _gradhat
            grad[j][i] += S[i][k] * jac[j][0][k];

      auto do_link = [&](std::size_t comp_i, std::size_t comp_j) {
        auto it = _component_pattern.find(std::make_pair(comp_i, comp_j));
        return (it != _component_pattern.end());
      };

      // compute grad^T * grad
      for (std::size_t k = 0; k < _components; k++) {
        for (std::size_t l = 0; l < _components; l++) {
          if (not do_link(k, l))
            continue;
          const auto j = _components * k + l;
          // evaluate reaction term
          _jacobian_gf[j]->update(u);
          _jacobian_gf[j]->evaluate(entity, position, jacobian[j]);
          for (std::size_t m = 0; m < basis_size; m++) {
            for (std::size_t n = 0; n < basis_size; n++) {
              typename M::value_type ljac =
                -jacobian[j] * phi[m] * phi[n];
              if (l == k)
                for (std::size_t d = 0; d < dim; d++)
                  ljac += diffusion[k] * grad[m][d] * grad[n][d];
              accumulate(k, m, l, n, ljac * factor);
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
    _jacobian_apply_volume(eg, lfsu, x, x, lfsv, r);
  }
};

/**
 * @brief      This class describes a PDELab local operator for temporal part
 *             for diffusion reaction systems.
 * @details    This class describre the operatrions for local integrals required
 *             for diffusion reaction system. The operator is only valid for
 *             entities contained in the grid view. The local finite element is
 *             used for caching shape function evaluations. And the jacobian
 *             method switches between numerical and analytical jacobians.
 *
 * @tparam     GV    Grid View
 * @tparam     LBT   Local Finite Element
 * @tparam     JM    Jacobian Method
 */
template<class GV, class LBT, JacobianMethod JM = JacobianMethod::Analytical>
class TemporalLocalOperatorDiffusionReactionCG
  : public Dune::PDELab::LocalOperatorDefaultFlags
  , public Dune::PDELab::FullVolumePattern
  , public Dune::PDELab::InstationaryLocalOperatorDefaultMethods<double>
  , protected LocalOperatorDiffusionReactionBase<GV,LBT>
  , public Dune::PDELab::NumericalJacobianVolume<
      TemporalLocalOperatorDiffusionReactionCG<GV, LBT, JM>>
  , public Dune::PDELab::NumericalJacobianApplyVolume<
      TemporalLocalOperatorDiffusionReactionCG<GV, LBT, JM>>
{
  //! grid view
  using GridView = GV;

  using LOPBase = LocalOperatorDiffusionReactionBase<GV,LBT>;

  using RF = typename LOPBase::RF;
  using LOPBase::_components;
  using LOPBase::dim;


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
   * @param[in]  grid_view       The grid view where this local operator is
   *                             valid
   * @param[in]  config          The configuration tree
   * @param[in]  finite_element  The local finite element
   */
  TemporalLocalOperatorDiffusionReactionCG(
    GridView grid_view,
    const ParameterTree& config)
    : LOPBase(config,grid_view)
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
  void alpha_volume(const EG& eg,
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

    // get geometry
    const auto geo = eg.geometry();

    // TODO: check geometry type
    const auto& rule = QuadratureRules<RF, dim>::rule(geo.type(),3);

    // loop over quadrature points
    for (const auto& point : rule) {
      const auto& position = point.position();

      // get Jacobian and determinant
      RF factor = point.weight() * geo.integrationElement(position);

      // contribution for each component
      for (std::size_t k = 0; k < _components;
           k++) // loop over components
      {
        const auto& trial_basis = lfsu.child(k).finiteElement().localBasis();
        const auto basis_size = trial_basis.size();
        std::vector<typename LBT::RangeType> phi(basis_size);
          trial_basis.evaluateFunction(position, phi);

        // compute value of component
        double u = 0.0;
        for (std::size_t j = 0; j < basis_size; j++) // ansatz func. loop
          u += x_coeff_local(k, j) * phi[j];

        for (std::size_t i = 0; i < basis_size; i++) // test func. loop
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
    if constexpr (JM == JacobianMethod::Numerical) {
      PDELab::NumericalJacobianVolume<
        TemporalLocalOperatorDiffusionReactionCG>::jacobian_volume(eg,
                                                                   lfsu,
                                                                   x,
                                                                   lfsv,
                                                                   mat);
      return;
    }

    auto accumulate = [&](const std::size_t& component_i,
                          const std::size_t& dof_i,
                          const std::size_t& component_j,
                          const std::size_t& dof_j,
                          const auto& value) {
      mat.accumulate(
        lfsv.child(component_i), dof_i, lfsu.child(component_j), dof_j, value);
    };

    // get geometry
    const auto geo = eg.geometry();

    // TODO: check geometry type
    const auto& rule = QuadratureRules<RF, dim>::rule(geo.type(),3);

    // loop over quadrature points
    for (const auto& point : rule) {
      const auto& position = point.position();


      // get Jacobian and determinant
      RF factor = point.weight() * geo.integrationElement(position);

      // integrate mass matrix
      for (std::size_t k = 0; k < _components;
           k++) // loop over components
      {
        const auto& trial_basis = lfsu.child(k).finiteElement().localBasis();
        const auto basis_size = trial_basis.size();
        std::vector<typename LBT::RangeType> phi(basis_size);
          trial_basis.evaluateFunction(position, phi);
        for (std::size_t i = 0; i < basis_size; i++)
          for (std::size_t j = 0; j < basis_size; j++)
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
    if constexpr (JM == JacobianMethod::Numerical) {
      PDELab::NumericalJacobianApplyVolume<
        TemporalLocalOperatorDiffusionReactionCG>::jacobian_apply_volume(eg,
                                                                         lfsu,
                                                                         x,
                                                                         z,
                                                                         lfsv,
                                                                         r);
      return;
    }
    alpha_volume(eg, lfsu, z, lfsv, r);
  }

  /**
   * @brief      The jacobian volume integral for matrix free operations (linear
   *             variant)
   * @details    This only switches between the actual implementation (in
   *             alpha_volume) and the numerical jacobian
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
      PDELab::NumericalJacobianApplyVolume<
        TemporalLocalOperatorDiffusionReactionCG>::jacobian_apply_volume(eg,
                                                                         lfsu,
                                                                         x,
                                                                         lfsv,
                                                                         r);
      return;
    }
    alpha_volume(eg, lfsu, x, lfsv, r);
  }
};

} // namespace Dune::Copasi

#endif // DUNE_COPASI_LOCAL_OPERATOR_DIFFUSION_REACTION_CG_HH
