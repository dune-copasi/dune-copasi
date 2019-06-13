#ifndef DUNE_COPASI_PARAMETERERIZATION_HH
#define DUNE_COPASI_PARAMETERERIZATION_HH

#include <dune/copasi/common/util.hh>

namespace Dune {
namespace Copasi {

/**
 * @brief      Parameter container
 *
 * @tparam     Parameters  A variadic set of parameters to store.
 */
template<class... Parameters>
class Parametererization
{

  // Assign an id for each parameter
  template<class Parameter>
  struct parameter_id : public argument_id<Parameter,0,Parameters...> {};

  // Parameter container
  using ParameterContainer = std::tuple<std::shared_ptr<Parameters>...>;

public:


  /**
   * @brief      Gets the parameter
   *
   * @tparam     Parameter  Parameter type
   *
   * @return     Shared pointer to the parameter
   */
  template<class Parameter>
  std::shared_ptr<Parameter> get_parameter()
  {
    constexpr int id = parameter_id<Parameter>::value;
    return std::get<id>(parameter_container);
  }

  /**
   * @brief      Gets the parameter
   *
   * @tparam     Parameter  Parameter type
   *
   * @return     Shared pointer to the parameter
   */
  template<class Parameter>
  std::shared_ptr<const Parameter> get_parameter() const
  {
    constexpr int id = parameter_id<Parameter>::value;
    return std::get<id>(parameter_container);
  }

  /**
   * @brief      Evaluation of the grid parameter in a grid entity with respect
   *             to a parameterization p at position x
   *
   * @param[in]  e                    The entity to evaluate on
   * @param[in]  x                    The position in entity-local coordinates
   * @param[in]  p                    The parametrization object
   *
   * @tparam     Parameter            The type of the parameter
   * @tparam     EntityType           The type of the entity
   * @tparam     DomainType           The type of the domain
   * @tparam     ParametrizationType  The type of the parametrization
   *
   * @return     The result of the evaluation of the parameter with respect to
   *             the parameterization p
   */
  template<class Parameter,  class EntityType, 
           class DomainType, class ParametrizationType>
  inline
  auto evaluate ( const EntityType& e,
                  const DomainType& x,
                  const ParametrizationType& p) const
  {
    auto parameter = this->get_parameter<Parameter>();
    return parameter->evaluate(e,x,p);
  }

private:
  ParameterContainer parameter_container;
};

} // Copasi namespace
} // Dune namespace

#endif // DUNE_COPASI_PARAMETERERIZATION_HH
