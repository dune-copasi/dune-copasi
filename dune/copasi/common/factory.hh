#ifndef DUNE_COPASI_FACTORY_HH
#define DUNE_COPASI_FACTORY_HH

namespace Dune::Copasi {

// class has to be specializated to be a factory of T!
template<class T>
struct Factory {};

} // namespace Dune::Copasi

#endif // DUNE_COPASI_FACTORY_HH