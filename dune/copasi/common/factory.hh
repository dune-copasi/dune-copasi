#ifndef DUNE_COPASI_FACTORY_HH
#define DUNE_COPASI_FACTORY_HH

#include <memory>

#include <dune/common/typetraits.hh>

namespace Dune::Copasi {

/**
@defgroup Factory Factory
@{
  @brief Class instance creator

  A factory should be able to create an instance of a type out of @ref DataContext.
  This is done by the definition of an specialization of the class Factory and a static function create.
  One thing to take into account is that data context can only store unique values of certain data type. 
  This is quite restrictive is the contructor of T contains repeated or very common types (e.g. `T{int,int,double}`),
  in shuch case, is best to wrap these values in a unique struct that contains all of these.
@}
*/


/**
 * @brief      Class Factory
 *
 * @tparam     T     Type to be created
 */
template<class T>
struct Factory {
  // Factory for type T has not been instantiated
  static_assert(Dune::AlwaysFalse<T>::value, "Factory does not exist");
};

} // namespace Dune::Copasi

#endif // DUNE_COPASI_FACTORY_HH