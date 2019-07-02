#ifndef DUNE_COPASI_LOCAL_OPERATOR_DIFFUSION_REACTION_HH
#define DUNE_COPASI_LOCAL_OPERATOR_DIFFUSION_REACTION_HH

// #include<dune/pdelab/common/geometrywrapper.hh>
// #include<dune/pdelab/common/quadraturerules.hh>
// #include<dune/pdelab/localoperator/defaultimp.hh>
// #include<dune/pdelab/localoperator/pattern.hh>
// #include<dune/pdelab/localoperator/flags.hh>
// #include<dune/pdelab/localoperator/idefault.hh>
// #include<dune/pdelab/finiteelement/localbasiscache.hh>

// #include<dune/common/exceptions.hh>

// #include<dune/geometry/referenceelements.hh>
// #include<dune/geometry/type.hh>
// #include<dune/common/fvector.hh>
// #include<dune/common/power.hh>

#include <dune/copasi/pdelab_expression_adapter.hh>


namespace Dune::Copasi::exp {

template<class GridView, class LocalFiniteElement>
class LocalOperatorDiffusionReaction :
  public Dune::PDELab::FullVolumePattern,
  public Dune::PDELab::LocalOperatorDefaultFlags,
  public Dune::PDELab::InstationaryLocalOperatorDefaultMethods<double>
{
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

  //! range dimension
  static constexpr int dim_range = LocalBasis::Traits::dimRange;

  // this operator only support scalar ranges
  static_assert(dim_range == 1);

  //! number of basis per component
  const std::size_t _basis_size;
  //! number of total components
  const std::size_t _components;
  //! reference to a quadrature rule
  const QuadratureRule<RF,dim>& _rule;

  //! basis functions at quadrature points
  std::vector<std::vector<Range>>     _phihat;
  //! basis function gradients at quadrature points
  std::vector<std::vector<Jacobian>>  _gradhat;

  ExpressionToGridFunctionAdapter<GridView,RF> _diffusion_gf;
  mutable ExpressionToGridFunctionAdapter<GridView,RF> _reaction_gf;
  mutable ExpressionToGridFunctionAdapter<GridView,RF> _jacobian_gf;

  Logging::Logger  _logger;

public:
  //! pattern assembly flags
  static constexpr bool doPatternVolume = true;

  //! residual assembly flags
  static constexpr bool doAlphaVolume = true;

  LocalOperatorDiffusionReaction(GridView& grid_view, 
                                 const ParameterTree& config, 
                                 const LocalFiniteElement& finite_element)
    : _basis_size(finite_element.size())
    , _components(config.sub("reaction").getValueKeys().size())
    , _rule(QuadratureRules<RF,dim>::rule(finite_element.type(),3)) // TODO: make order variable
    , _diffusion_gf(grid_view,config.sub("diffusion"))
    , _reaction_gf(grid_view,config.sub("reaction"),config.sub("reaction"))
    , _jacobian_gf(grid_view,config.sub("jacobian"),config.sub("reaction"))
    , _logger(Logging::Logging::componentLogger(config,"default"))
  {
    assert(_components == config.sub("diffusion").getValueKeys().size());
    assert(std::pow(_components,2) == config.sub("jacobian").getValueKeys().size());

    _logger.trace("cache finite element evaluations on reference element"_fmt);
    std::vector<Range> phi(_basis_size);
    std::vector<Jacobian> jac(_basis_size);
    for (const auto& point : _rule) 
    {
      const auto& position = point.position();
      _logger.trace("position: {}"_fmt,position);

      const auto& local_basis = finite_element.localBasis();
      local_basis.evaluateFunction(position,phi);
      local_basis.evaluateJacobian(position,jac);

      for (int i = 0; i < _basis_size; ++i)
      {
        _logger.trace(" value[{}]: {}"_fmt,i,phi[i]);
        _logger.trace(" jacobian[{}]: {}"_fmt,i,jac[i]);
      }
      _phihat.push_back(phi);  phi.clear();
      _gradhat.push_back(jac); jac.clear();
    }

    _logger.debug("LocalOperatorDiffusionReaction constructed"_fmt);
  }

  template<typename EG, typename LFSU, typename X,
           typename LFSV, typename R>
  void jacobian_apply_volume (const EG& eg, const LFSU& lfsu,
                              const X& x, const X& z, const LFSV& lfsv,
                              R& r) const
  {
    // assume we receive a power local finite element!
    auto x_coeff = [&](const std::size_t& component, const std::size_t& dof)
    {
      return x(lfsu,component*_basis_size+dof);
    };
    auto z_coeff = [&](const std::size_t& component, const std::size_t& dof)
    {
      return z(lfsu,component*_basis_size+dof);
    };

    auto accumulate = [&](const std::size_t& component, const std::size_t& dof, const auto& value)
    {
      r.accumulate(lfsu,component*_basis_size+dof,value);
    };

    // get entity
    const auto entity = eg.entity();

    // get geometry
    const auto geo = eg.geometry();

    // loop over quadrature points
    for (std::size_t q = 0; q < _rule.size(); q++)
    {
      const auto& position = _rule[q].position();
      
      // get diffusion coefficient
      DynamicVector<RF> diffusion;
      _diffusion_gf.evaluate(entity,position,diffusion);

      // get jacobian and determinant
      FieldMatrix<DF,dim,dim> S = geo.jacobianInverseTransposed(position);
      RF factor = _rule[q].weight() * geo.integrationElement(position);

      // evaluate concentrations at quadrature point
      DynamicVector<RF> u(_components);
      for (std::size_t k=0; k<_components; k++) // loop over components
        for (std::size_t j=0; j<_basis_size; j++) // loop over ansatz functions
          u[k] += x_coeff(k,j)*_phihat[q][j];

      // get diffusion coefficient
      DynamicVector<RF> reaction;
      _reaction_gf.bind(entity,u);
      _reaction_gf.evaluate(entity,position,reaction);

      // compute gradients of basis functions in transformed element (independent of component)
      DynamicVector<FieldVector<RF,dim>> grad(_basis_size);
      for (int i=0; i<dim; i++) // rows of S
        for (int k=0; k<dim; k++) // columns of S
          for (std::size_t j=0; j<_basis_size; j++) // columns of _gradhat
            grad[j][i] += S[i][k] * _gradhat[q][j][0][k];

      // contribution for each component
      for (std::size_t k=0; k<_components; k++) // loop over components
      {
        // compute gradient u_h
        FieldVector<RF,dim> graduh(.0);
        for (int d=0; d<dim; d++) // rows of grad
          for (int j=0; j<_basis_size; j++) // columns of grad
            graduh[d] += grad[j][d]*z_coeff(k,j);

        // scalar products
        for (int d=0; d<dim; d++) // rows of grad
          for (int i=0; i<_basis_size; i++) // loop over test functions
            accumulate(k,i,diffusion[k]*grad[i][d]*graduh[d]*factor);

        // reaction term
        for (int i=0; i<_basis_size; i++) // loop over test functions
          accumulate(k,i,reaction[k]*_phihat[q][i]*factor);
      }
    }
  }

  //! jacobian contribution of volume term
  template<typename EG, typename LFSU, typename X,
           typename LFSV, typename M>
  void jacobian_volume (const EG& eg, const LFSU& lfsu, const X& x,
                        const LFSV& lfsv, M& mat) const
  {

    // assume we receive a power local finite element!
    auto x_coeff = [&](const std::size_t& component, const std::size_t& dof)
    {
      return x(lfsu,component*_basis_size+dof);
    };

    auto accumulate = [&](const std::size_t& component_i, const std::size_t& dof_i, 
                          const std::size_t& component_j, const std::size_t& dof_j, 
                          const auto& value)
    {
      mat.accumulate(lfsv,component_i*_basis_size+dof_i,lfsu,component_j*_basis_size+dof_j,value);
    };

    // get entity
    const auto entity = eg.entity();

    // get geometry
    const auto geo = eg.geometry();

    // loop over quadrature points
    for (std::size_t q = 0; q < _rule.size(); q++)
    {
      // local stiffness matrix (independent of component)
      std::vector<std::vector<RF>> A(_basis_size);
      std::fill(A.begin(),A.end(),std::vector<RF>(_basis_size));

      const auto& position = _rule[q].position();
      
      // get diffusion coefficient
      DynamicVector<RF> diffusion;
      _diffusion_gf.evaluate(entity,position,diffusion);

      // get jacobian and determinant
      FieldMatrix<DF,dim,dim> S = geo.jacobianInverseTransposed(position);
      RF factor = _rule[q].weight() * geo.integrationElement(position);

      // evaluate concentrations at quadrature point
      DynamicVector<RF> u(_components);
      for (std::size_t k=0; k<_components; k++) // loop over components
        for (std::size_t j=0; j<_basis_size; j++) // loop over ansatz functions
          u[k] += x_coeff(k,j)*_phihat[q][j];

      // evaluate reaction term
      DynamicVector<RF> jacobian;
      _jacobian_gf.bind(entity,u);
      _jacobian_gf.evaluate(entity,position,jacobian);

      // compute gradients of basis functions in transformed element (independent of component)
      DynamicVector<FieldVector<RF,dim>> grad(_basis_size);
      for (int i=0; i<dim; i++) // rows of S
        for (int k=0; k<dim; k++) // columns of S
          for (int j=0; j<_basis_size; j++) // columns of _gradhat
            grad[j][i] += S[i][k] * _gradhat[q][j][0][k];

      // compute grad^T * grad
      for (int k=0; k<_components; k++)
        for (int i=0; i<_basis_size; i++)
          for (int d=0; d<dim; d++)
            for (int j=0; j<_basis_size; j++)
              accumulate(k,i,k,j,diffusion[k]*grad[i][d]*grad[j][d]*factor);

      int count = 0;
      for (int k=0; k<_components; k++)
        for (int l=0; l<_components; l++, count++)
          for (int i=0; i<_basis_size; i++)
            for (int j=0; j<_basis_size; j++)
              accumulate(k,i,l,j,_phihat[q][i]*jacobian[count]*_phihat[q][j]*factor);
    }
  }

  template<typename EG, typename LFSU, typename X,
           typename LFSV, typename R>
  void alpha_volume (const EG& eg, const LFSU& lfsu, const X& x,
                     const LFSV& lfsv, R& r) const
  {
    jacobian_apply_volume(eg,lfsu,x,x,lfsv,r);
  }

  //! apply local jacobian of the volume term -> linear variant
  template<typename EG, typename LFSU, typename X,
           typename LFSV, typename R>
  void jacobian_apply_volume (const EG& eg, const LFSU& lfsu,
                              const X& x, const LFSV& lfsv,
                              R& r) const
  {
    jacobian_apply_volume(eg,lfsu,x,x,lfsv,r);
  }
};


/******************************************************/
/** a local operator for the spatial part of a diffusion-reaction system
 *
 * We assume that both components use the same finite element space.
 * Beware of line number changes, they may corrupt docu!
 * Template parameters are:
 * \tparam FiniteElementMap type for finite element map
 * \tparam components       size of the system
 */
/******************************************************/
template<typename FiniteElementMap, int components>
class TemporalLocalOperatorDiffusionReaction :
  // public Dune::PDELab::NumericalJacobianVolume<TemporalLocalOperatorDiffusionReaction<FiniteElementMap,components> >,
  public Dune::PDELab::FullVolumePattern,
  public Dune::PDELab::LocalOperatorDefaultFlags,
  public Dune::PDELab::InstationaryLocalOperatorDefaultMethods<typename FiniteElementMap::Traits::FiniteElementType::Traits::LocalBasisType::Traits::RangeFieldType>
{
//   // define useful types
//   typedef typename FiniteElementMap::Traits::FiniteElementType
//      FiniteElementType;
//   typedef typename FiniteElementType::Traits::LocalBasisType
//      LocalBasisType;
//   typedef typename LocalBasisType::Traits::DomainFieldType
//      DF;
//   typedef typename LocalBasisType::Traits::DomainType
//      DomainType;
//   typedef typename LocalBasisType::Traits::RangeFieldType
//      RF;
//   typedef typename LocalBasisType::Traits::RangeType
//      RangeType;
//   typedef typename LocalBasisType::Traits::JacobianType
//      JacobianType;

//   // sizes
//   enum {dim=LocalBasisType::Traits::dimDomain}; // dimension of domain
//   enum {n=4};        // number of basis functions per components
//   enum {m1d=2};           // number of quadrature points per direction
//   enum {m = Dune::StaticPower<m1d,dim>::power }; // number of quadrature points

//   // data members
//   DomainType qp[m];       // quadrature points
//   double weight[m];       // quadrature weights
//   double _phihat[m][n];    // basis functions at quadrature points

// public:
//   // pattern assembly flags
//   enum { doPatternVolume = true };

//   // residual assembly flags
//   enum { doAlphaVolume = true };

//   //! constructor stores the speed of sound
//   // Constructor precomputes element independent data
//   TemporalLocalOperatorDiffusionReaction (const FiniteElementType& fel)
//   {
//     // get quadrature rule of order q
//     Dune::GeometryType gt = fel.type();
//     const Dune::QuadratureRule<RF,dim>&
//       rule = Dune::QuadratureRules<RF,dim>::rule(gt,3);
//     std::cout << "number of quadrature points:"
//               << " wanted " << m
//               << " got " << rule.size()
//               << std::endl;
//     if (rule.size()!=m) {
//       std::cout << "mismatch in number of quadrature points:"
//                 << " wanted " << m
//                 << " got " << rule.size()
//                 << std::endl;
//       exit(1);
//     }

//     // position and weight of the quadrature point
//     for (int k=0; k<m; k++) {
//       weight[k] = rule[k].weight();
//       qp[k] = rule[k].position();
//     }

//     // check size of the basis
//     std::cout << "number of basis functions"
//               << " wanted " << n
//               << " got " << int(fel.localBasis().size()/components)
//               << std::endl;
//     if (int(fel.localBasis().size()/components)!=n) {
//       std::cout << "mismatch in number of basis functions"
//                 << " wanted " << n
//                 << " got " << int(fel.localBasis().size()/components)
//                 << std::endl;
//       exit(1);
//     }

//     // evaluate basis functions on refelem
//     std::vector<RangeType> phi(n);
//     for (int k=0; k<m; k++) {
//       fel.localBasis().evaluateFunction(qp[k],phi);
//       for (int i=0; i<n; i++)
//         _phihat[k][i] = phi[i];
//     }
//   }

//   //! volume integral depending on test and ansatz functions
//   template<typename EG, typename LFSU, typename X,
//            typename LFSV, typename R>
//   void alpha_volume (const EG& eg, const LFSU& lfsu, const X& x,
//                      const LFSV& lfsv, R& r) const
//   {
//     // select the two components (but assume Galerkin scheme U=V)
//     // assert(LFSU::CHILDREN==components);

//     // extract coefficients
//     double xc[components][n];
//     for (int k=0; k<components; k++) // loop over components
//       {
//         for (int j=0; j<n; j++) // loop of ansatz functions
//           xc[k][j] = x(lfsu,k*n+j); // read coeffs
//       }

//     // the result
//     double a[components][n] = {{0.0}};

//     // get geometry
//     const auto geo = eg.geometry();

//     // loop over quadrature points
//     for (int q=0; q<m; q++)
//       {
//         // get Jacobian and determinant
//         RF factor = weight[q]*geo.integrationElement(qp[q]);

//         // contribution for each component
//         for (int k=0; k<components; k++) // loop over components
//           {
//             // compute value of component
//             double u=0.0;
//             for (int j=0; j<n; j++) // loop over ansatz functions
//               u += xc[k][j]*_phihat[q][j];

//             // reaction term
//             for (int i=0; i<n; i++) // loop over test functions
//               a[k][i] += u*_phihat[q][i]*factor;
//           }
//       }

//       // store result
//       for (int k=0; k<components; k++) // loop over components
//         {
//           for (int i=0; i<n; i++)
//             r.accumulate(lfsu,k*n+i,a[k][i]);
//         }
//   }

//   //! jacobian contribution of volume term
//   template<typename EG, typename LFSU, typename X,
//            typename LFSV, typename Mat>
//   void jacobian_volume (const EG& eg, const LFSU& lfsu, const X& x,
//                         const LFSV& lfsv, Mat& mat) const
//   {
//     // select the two components (but assume Galerkin scheme U=V)
//     // assert(LFSU::CHILDREN==components);

//     // local stiffness matrix (independent of component)
//     double M[n][n] = {{0.0}};

//     // get geometry
//     const auto geo = eg.geometry();

//     // loop over quadrature points
//     for (int q=0; q<m; q++)
//       {
//         // get Jacobian and determinant
//         RF factor = weight[q]*geo.integrationElement(qp[q]);

//         // integrate mass matrix
//         for (int i=0; i<n; i++)
//           for (int j=0; j<n; j++)
//             M[i][j] += _phihat[q][i]*_phihat[q][j]*factor;
//       }

//     // store in result
//     for (int k=0; k<components; k++) // loop over components
//       {
//         for (int i=0; i<n; i++)
//           for (int j=0; j<n; j++)
//             mat.accumulate(lfsu,k*n+i,lfsu,k*n+j,M[i][j]);
//       }
//   }

//   //! apply local jacobian of the volume term -> nonlinear variant
//   template<typename EG, typename LFSU, typename X,
//            typename LFSV, typename R>
//   void jacobian_apply_volume (const EG& eg, const LFSU& lfsu,
//                               const X& x, const X& z, const LFSV& lfsv,
//                               R& r) const
//   {
//     alpha_volume(eg,lfsu,z,lfsv,r);
//   }

//   //! apply local jacobian of the volume term -> linear variant
//   template<typename EG, typename LFSU, typename X,
//            typename LFSV, typename R>
//   void jacobian_apply_volume (const EG& eg, const LFSU& lfsu,
//                               const X& x, const LFSV& lfsv,
//                               R& r) const
//   {
//     alpha_volume(eg,lfsu,x,lfsv,r);
//   }
};

} // Dune::Copasi namespace

#endif // DUNE_COPASI_LOCAL_OPERATOR_DIFFUSION_REACTION_HH