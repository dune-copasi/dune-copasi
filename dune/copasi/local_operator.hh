#include<dune/pdelab/common/geometrywrapper.hh>
#include<dune/pdelab/common/quadraturerules.hh>
#include<dune/pdelab/localoperator/defaultimp.hh>
#include<dune/pdelab/localoperator/pattern.hh>
#include<dune/pdelab/localoperator/flags.hh>
#include<dune/pdelab/localoperator/idefault.hh>
#include<dune/pdelab/finiteelement/localbasiscache.hh>

#include<dune/common/exceptions.hh>
#include<dune/common/fvector.hh>
#include<dune/geometry/referenceelements.hh>
#include<dune/geometry/type.hh>

#include<dune/common/power.hh>

namespace Dune {
namespace Copasi {

/**
 * @brief      Local operator for the spatial part of a diffusion-reaction
 *             system
 * @details    We assume that both components use the same finite element space.
 *
 * @tparam     Param             Parameter class
 * @tparam     FiniteElementMap  Finite element map class
 * @tparam     order             Order of the quadrature rule
 */
template<typename Param, typename FiniteElementMap, int order=3>
class LocalOperatorDiffusionReaction :
  // public Dune::PDELab::NumericalJacobianVolume<LocalOperatorDiffusionReaction<Param,FiniteElementMap,order> >,
  public Dune::PDELab::FullVolumePattern,
  public Dune::PDELab::LocalOperatorDefaultFlags,
  public Dune::PDELab::InstationaryLocalOperatorDefaultMethods<typename FiniteElementMap::Traits::FiniteElementType::Traits::LocalBasisType::Traits::RangeFieldType>
{
  // define useful types
  typedef typename FiniteElementMap::Traits::FiniteElementType
     FiniteElementType;
  typedef typename FiniteElementType::Traits::LocalBasisType
     LocalBasisType;
  typedef typename LocalBasisType::Traits::DomainFieldType
     DF;
  typedef typename LocalBasisType::Traits::DomainType
     DomainType;
  typedef typename LocalBasisType::Traits::RangeFieldType
     RF;
  typedef typename LocalBasisType::Traits::RangeType
     RangeType;
  typedef typename LocalBasisType::Traits::JacobianType
     JacobianType;

  // sizes
  enum {dim=LocalBasisType::Traits::dimDomain}; // dimension of domain
  enum {n=4};        // number of basis functions per components
  enum {m1d=order/2+1};   // number of quadrature points per direction
  enum {m = Dune::StaticPower<m1d,dim>::power }; // number of quadrature points
  enum {components=Param::components}; // number of components in diffusion-reaction system

  // data members
  Param& param;              // parameter functions
  DomainType center;         // center of reference element
  DomainType qp[m];          // quadrature points
  double weight[m];          // quadrature weights
  double phihat[m][n];       // basis functions at quadrature points
  double gradhat[m][dim][n]; // gradient of basis functions at quadrature points
  RF time;                   // current time

public:
  // pattern assembly flags
  enum { doPatternVolume = true };

  // residual assembly flags
  enum { doAlphaVolume = true };

  // Constructor precomputes element independent data
  LocalOperatorDiffusionReaction (Param& param_, const FiniteElementType& fel)
    : param(param_)
  {
    for (int i=0; i<dim; i++)
      center[i] = 0.5;

    // get quadrature rule of order q
    Dune::GeometryType gt = fel.type();
    const Dune::QuadratureRule<RF,dim>&
      rule = Dune::QuadratureRules<RF,dim>::rule(gt,order);
    std::cout << "number of quadrature points:"
              << " wanted " << m
              << " got " << rule.size()
              << std::endl;
    if (rule.size()!=m) {
      std::cout << "mismatch in number of quadrature points:"
                << " wanted " << m
                << " got " << rule.size()
                << std::endl;
      exit(1);
    }

    // position and weight of the quadrature point
    for (int k=0; k<m; k++) {
      weight[k] = rule[k].weight();
      qp[k] = rule[k].position();
    }

    // check size of the basis
    std::cout << "number of basis functions"
              << " wanted " << n
              << " got " << int(fel.localBasis().size()/components)
              << std::endl;
    if (int(fel.localBasis().size()/components)!=n) {
      std::cout << "mismatch in number of basis functions"
                << " wanted " << n
                << " got " << int(fel.localBasis().size()/components)
                << std::endl;
      exit(1);
    }

    // evaluate basis functions on refelem
    std::vector<RangeType> phi(n);
    for (int k=0; k<m; k++) {
      fel.localBasis().evaluateFunction(qp[k],phi);
      for (int i=0; i<n; i++)
        phihat[k][i] = phi[i];
    }

    // evaluate gradients of basis functions on refelem
    std::vector<JacobianType> js(n);
    for (int k=0; k<m; k++) {
      fel.localBasis().evaluateJacobian(qp[k],js);
      for (int i=0; i<n; i++)
        for (int j=0; j<dim; j++)
          gradhat[k][j][i] = js[i][0][j];
    }
  }

  void setTime (RF t_)
  {
    time = t_;
    param.setTime(time);
  }

  //! volume integral depending on test and ansatz functions
  template<typename EG, typename LFSU, typename X,
           typename LFSV, typename R>
  void jacobian_apply_volume (const EG& eg, const LFSU& lfsu,
                              const X& x, const X& z, const LFSV& lfsv,
                              R& r) const
  {
    // select the two components (but assume Galerkin scheme U=V)
    // assert(LFSU::CHILDREN==components);

    // extract coefficients
    double xc[components][n];
    double zc[components][n];
    for (int k=0; k<components; k++) // loop over components
      {
        for (int j=0; j<n; j++) // loop of ansatz functions
        {
          xc[k][j] = x(lfsu,k*n+j);
          zc[k][j] = z(lfsu,k*n+j);  // read coeffs
        }
      }

    // the result
    double a[components][n] = {{0.0}};

    // get geometry
    const auto geo = eg.geometry();

    // get diffusion coefficient
    auto D=param.D(eg.entity(),center);

    // loop over quadrature points
    for (int q=0; q<m; q++)
      {
        // get Jacobian and determinant
        Dune::FieldMatrix<DF,dim,dim> S = geo.jacobianInverseTransposed(qp[q]);
        RF factor = weight[q]*geo.integrationElement(qp[q]);

        // evaluate concentrations at quadrature point
        double u[components] = {0.0};
        for (int k=0; k<components; k++) // loop over components
          for (int j=0; j<n; j++) // loop over ansatz functions
            u[k] += xc[k][j]*phihat[q][j];

        // evaluate reaction term
        auto f=param.f(eg.entity(),qp[q],u);

        // compute gradients of basis functions in transformed element (independent of component)
        double grad[dim][n] = {{0.0}};  // coordinate x #basisfct
        for (int i=0; i<dim; i++) // rows of S
          for (int k=0; k<dim; k++) // columns of S
            for (int j=0; j<n; j++) // columns of gradhat
              grad[i][j] += S[i][k] * gradhat[q][k][j];

        // contribution for each component
        for (int k=0; k<components; k++) // loop over components
          {
            // compute gradient u_h
            double graduh[dim] = {0.0};
            for (int d=0; d<dim; d++) // rows of grad
              for (int j=0; j<n; j++) // columns of grad
                graduh[d] += grad[d][j]*zc[k][j];

            // scalar products
            for (int d=0; d<dim; d++) // rows of grad
              for (int i=0; i<n; i++) // loop over test functions
                a[k][i] += D[k]*grad[d][i]*graduh[d]*factor;

            // reaction term
            for (int i=0; i<n; i++) // loop over test functions
              a[k][i] += f[k]*phihat[q][i]*factor;
          }
      }

      // store result
      for (int k=0; k<components; k++) // loop over components
        {
          for (int i=0; i<n; i++)
            r.accumulate(lfsu,k*n+i,a[k][i]);
        }
  }

  //! jacobian contribution of volume term
  template<typename EG, typename LFSU, typename X,
           typename LFSV, typename M>
  void jacobian_volume (const EG& eg, const LFSU& lfsu, const X& x,
                        const LFSV& lfsv, M& mat) const
  {
    // select the two components (but assume Galerkin scheme U=V)
    // assert(LFSU::CHILDREN==components);

    // extract coefficients
    double xc[components][n];
    for (int k=0; k<components; k++) // loop over components
      {
        for (int j=0; j<n; j++) // loop of ansatz functions
          xc[k][j] = x(lfsu,k*n+j); // read coeffs
      }

    // local stiffness matrix (independent of component)
    double A[n][n] = {{0.0}};

    // get geometry
    const auto geo = eg.geometry();

    // get diffusion coefficient
    auto D=param.D(eg.entity(),center);

    // loop over quadrature points
    for (int q=0; q<m; q++)
      {
        // get Jacobian and determinant
        Dune::FieldMatrix<DF,dim,dim> S = geo.jacobianInverseTransposed(qp[q]);
        RF factor = weight[q]*geo.integrationElement(qp[q]);

        // compute gradients of basis functions in transformed element at qp
        double grad[dim][n] = {{0.0}}; // coordinate x #basisfct
        for (int i=0; i<dim; i++) // rows of S
          for (int d=0; d<dim; d++) // columns of S
            for (int j=0; j<n; j++) // columns of gradhat
              grad[i][j] += S[i][d] * gradhat[q][d][j];

        // compute grad^T * grad
        for (int i=0; i<n; i++)
          for (int d=0; d<dim; d++)
            for (int j=0; j<n; j++)
              A[i][j] += grad[d][i]*grad[d][j]*factor;

        // evaluate concentrations at quadrature point
        double u[components] = {0.0};
        for (int k=0; k<components; k++) // loop over components
          for (int j=0; j<n; j++) // loop over ansatz functions
            u[k] += xc[k][j]*phihat[q][j];

        // evaluate reaction term
        auto nablaf=param.nablaf(eg.entity(),qp[q],u);
        for (int e=0; e<Param::nonzeroes; e++)
          {
            auto childrow = nablaf[e].row;
            auto childcol = nablaf[e].col;
            for (int i=0; i<n; i++)
              for (int j=0; j<n; j++)
                mat.accumulate(lfsv,childrow*n+i,lfsu,childcol*n+j,phihat[q][i]*nablaf[e].value*phihat[q][j]*factor);
          }
      }

    // store in result
    for (int k=0; k<components; k++) // loop over components
      {
        for (int i=0; i<n; i++)
          for (int j=0; j<n; j++)
            mat.accumulate(lfsu,k*n+i,lfsu,k*n+j,D[k]*A[i][j]);
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
  // define useful types
  typedef typename FiniteElementMap::Traits::FiniteElementType
     FiniteElementType;
  typedef typename FiniteElementType::Traits::LocalBasisType
     LocalBasisType;
  typedef typename LocalBasisType::Traits::DomainFieldType
     DF;
  typedef typename LocalBasisType::Traits::DomainType
     DomainType;
  typedef typename LocalBasisType::Traits::RangeFieldType
     RF;
  typedef typename LocalBasisType::Traits::RangeType
     RangeType;
  typedef typename LocalBasisType::Traits::JacobianType
     JacobianType;

  // sizes
  enum {dim=LocalBasisType::Traits::dimDomain}; // dimension of domain
  enum {n=4};        // number of basis functions per components
  enum {m1d=2};           // number of quadrature points per direction
  enum {m = Dune::StaticPower<m1d,dim>::power }; // number of quadrature points

  // data members
  DomainType qp[m];       // quadrature points
  double weight[m];       // quadrature weights
  double phihat[m][n];    // basis functions at quadrature points

public:
  // pattern assembly flags
  enum { doPatternVolume = true };

  // residual assembly flags
  enum { doAlphaVolume = true };

  //! constructor stores the speed of sound
  // Constructor precomputes element independent data
  TemporalLocalOperatorDiffusionReaction (const FiniteElementType& fel)
  {
    // get quadrature rule of order q
    Dune::GeometryType gt = fel.type();
    const Dune::QuadratureRule<RF,dim>&
      rule = Dune::QuadratureRules<RF,dim>::rule(gt,3);
    std::cout << "number of quadrature points:"
              << " wanted " << m
              << " got " << rule.size()
              << std::endl;
    if (rule.size()!=m) {
      std::cout << "mismatch in number of quadrature points:"
                << " wanted " << m
                << " got " << rule.size()
                << std::endl;
      exit(1);
    }

    // position and weight of the quadrature point
    for (int k=0; k<m; k++) {
      weight[k] = rule[k].weight();
      qp[k] = rule[k].position();
    }

    // check size of the basis
    std::cout << "number of basis functions"
              << " wanted " << n
              << " got " << int(fel.localBasis().size()/components)
              << std::endl;
    if (int(fel.localBasis().size()/components)!=n) {
      std::cout << "mismatch in number of basis functions"
                << " wanted " << n
                << " got " << int(fel.localBasis().size()/components)
                << std::endl;
      exit(1);
    }

    // evaluate basis functions on refelem
    std::vector<RangeType> phi(n);
    for (int k=0; k<m; k++) {
      fel.localBasis().evaluateFunction(qp[k],phi);
      for (int i=0; i<n; i++)
        phihat[k][i] = phi[i];
    }
  }

  //! volume integral depending on test and ansatz functions
  template<typename EG, typename LFSU, typename X,
           typename LFSV, typename R>
  void alpha_volume (const EG& eg, const LFSU& lfsu, const X& x,
                     const LFSV& lfsv, R& r) const
  {
    // select the two components (but assume Galerkin scheme U=V)
    // assert(LFSU::CHILDREN==components);

    // extract coefficients
    double xc[components][n];
    for (int k=0; k<components; k++) // loop over components
      {
        for (int j=0; j<n; j++) // loop of ansatz functions
          xc[k][j] = x(lfsu,k*n+j); // read coeffs
      }

    // the result
    double a[components][n] = {{0.0}};

    // get geometry
    const auto geo = eg.geometry();

    // loop over quadrature points
    for (int q=0; q<m; q++)
      {
        // get Jacobian and determinant
        RF factor = weight[q]*geo.integrationElement(qp[q]);

        // contribution for each component
        for (int k=0; k<components; k++) // loop over components
          {
            // compute value of component
            double u=0.0;
            for (int j=0; j<n; j++) // loop over ansatz functions
              u += xc[k][j]*phihat[q][j];

            // reaction term
            for (int i=0; i<n; i++) // loop over test functions
              a[k][i] += u*phihat[q][i]*factor;
          }
      }

      // store result
      for (int k=0; k<components; k++) // loop over components
        {
          for (int i=0; i<n; i++)
            r.accumulate(lfsu,k*n+i,a[k][i]);
        }
  }

  //! jacobian contribution of volume term
  template<typename EG, typename LFSU, typename X,
           typename LFSV, typename Mat>
  void jacobian_volume (const EG& eg, const LFSU& lfsu, const X& x,
                        const LFSV& lfsv, Mat& mat) const
  {
    // select the two components (but assume Galerkin scheme U=V)
    // assert(LFSU::CHILDREN==components);

    // local stiffness matrix (independent of component)
    double M[n][n] = {{0.0}};

    // get geometry
    const auto geo = eg.geometry();

    // loop over quadrature points
    for (int q=0; q<m; q++)
      {
        // get Jacobian and determinant
        RF factor = weight[q]*geo.integrationElement(qp[q]);

        // integrate mass matrix
        for (int i=0; i<n; i++)
          for (int j=0; j<n; j++)
            M[i][j] += phihat[q][i]*phihat[q][j]*factor;
      }

    // store in result
    for (int k=0; k<components; k++) // loop over components
      {
        for (int i=0; i<n; i++)
          for (int j=0; j<n; j++)
            mat.accumulate(lfsu,k*n+i,lfsu,k*n+j,M[i][j]);
      }
  }

  //! apply local jacobian of the volume term -> nonlinear variant
  template<typename EG, typename LFSU, typename X,
           typename LFSV, typename R>
  void jacobian_apply_volume (const EG& eg, const LFSU& lfsu,
                              const X& x, const X& z, const LFSV& lfsv,
                              R& r) const
  {
    alpha_volume(eg,lfsu,z,lfsv,r);
  }

  //! apply local jacobian of the volume term -> linear variant
  template<typename EG, typename LFSU, typename X,
           typename LFSV, typename R>
  void jacobian_apply_volume (const EG& eg, const LFSU& lfsu,
                              const X& x, const LFSV& lfsv,
                              R& r) const
  {
    alpha_volume(eg,lfsu,x,lfsv,r);
  }
};

} // Copasi namespace
} // Dune namespace
