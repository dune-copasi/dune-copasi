#ifndef DUNE_COPASI_PARAMETER_BASE_HH
#define DUNE_COPASI_PARAMETER_BASE_HH

namespace Dune {
namespace Copasi {
/**
 * @brief      Base class for parameters.
 * @details    This class is an interface that every parameter class must
 *             follow.
 *
 * @tparam     Imp   Parameter implementation
 */
template<class Imp>
class ParameterBase
{
  /**
   * @brief      Evaluation of the grid parameter in a grid entity with respect
   *             to a Parametrization at position x
   *
   * @param[in]  e                    The entity to evaluate on
   * @param[in]  x                    The position in entity-local coordinates
   * @param[in]  p                    The parametrization object
   *
   * @tparam     EtityType            The type of the entity
   * @tparam     DomainType           The type of the domain
   * @tparam     ParametrizationType  The type of the parametrization
   *
   * @return     The result of the evaluation of the parameter with respect to
   *             the parametrization
   */
  template<class EtityType, class DomainType, class ParametrizationType>
  inline
  auto evaluate ( const EtityType& e,
                  const DomainType& x,
                  const ParametrizationType& p) const
  {
    return asImp().evaluate(e,x,p);
  }

private:
  Imp& asImp () {return static_cast<Imp &> (*this);}
  const Imp& asImp () const {return static_cast<const Imp &>(*this);}
};

} // Copasi namespace
} // Dune namespace

#endif // DUNE_COPASI_PARAMETER_BASE_HH
