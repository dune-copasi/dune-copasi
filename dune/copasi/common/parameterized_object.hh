#ifndef DUNE_COPASI_COMMON_PARAMETERIZED_OBJECT_HH
#define DUNE_COPASI_COMMON_PARAMETERIZED_OBJECT_HH

#include <dune/common/parameterizedobject.hh>

#include <set>
#include <string>

namespace Dune::Copasi {

//! simple wrapper that captures and exposes the keys of Dune::ParameterizedObjectFactory
template<typename Signature,
         typename Key = std::string>
class ParameterizedObjectFactoryWrapper : public Dune::ParameterizedObjectFactory<Signature, Key> {
  using Super = Dune::ParameterizedObjectFactory<Signature, Key>;
public:

  using Super::Super;

  template<class T>
  void define(Key const& key, T&& t)
  {
    Super::define(key, std::forward<T>(t));
    _keys.insert(key);
  }

  const std::set<Key>& keys() const {
    return _keys;
  }

private:
  std::set<Key> _keys;
};

} // namespace Dune::Copasi

#endif // DUNE_COPASI_COMMON_PARAMETERIZED_OBJECT_HH
