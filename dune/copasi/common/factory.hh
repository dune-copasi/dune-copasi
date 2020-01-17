#ifndef DUNE_COPASI_FACTORY_HH
#define DUNE_COPASI_FACTORY_HH

namespace Dune::Copasi {

// class has to be specializated to be a factory of T!
template<class T>
struct Factory {
  static_assert(Dune::AlwaysFalse<T>::value, "Factory does not exist");
};

} // namespace Dune::Copasi

#endif // DUNE_COPASI_FACTORY_HH