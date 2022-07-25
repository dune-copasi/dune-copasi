#ifndef DUNE_COPASI_LOCAL_FINITE_ELEMENT_PK_HH
#define DUNE_COPASI_LOCAL_FINITE_ELEMENT_PK_HH

#include <dune/copasi/common/factory.hh>

#include<dune/localfunctions/lagrange/lagrangesimplex.hh>

#include <dune/geometry/type.hh>

namespace Dune::Copasi {

/**
 * @brief      Factory for LagrangeSimplexLocalFiniteElement instances
 * @ingroup    Factory, FiniteElement
 * @tparam     <unnamed>  Template paramenters of the PkLocalFiniteElement
 */
template<class D, class R, int d, int k>
struct Factory<Dune::LagrangeSimplexLocalFiniteElement<D,R,d,k>>
{
  /**
   * @brief      Create method
   *
   * @param      ctx   Empty @ref DataContext
   *
   * @tparam     Ctx   Universal reference to the @ref DataContext
   *
   * @return     Instance of PkLocalFiniteElement
   */
  template<class Ctx>
  static auto create(Ctx&& ctx)
  {
    return std::make_unique<Dune::LagrangeSimplexLocalFiniteElement<D,R,d,k>>();
  }
};

} // namespace Dune::Copasi

#endif // DUNE_COPASI_LOCAL_FINITE_ELEMENT_PK_HH
