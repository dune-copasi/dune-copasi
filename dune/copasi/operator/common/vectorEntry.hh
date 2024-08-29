#ifndef DUNE_COPASI_VECTORENTRY_HH
#define DUNE_COPASI_VECTORENTRY_HH

#include <dune/copasi/operator/common/concepts.hh>
#include <dune/pdelab/basis/backend/istl.hh> 

namespace Dune::Copasi{

  namespace impl{

    template<std::size_t l>
    struct level_v
        : public std::integral_constant<std::size_t,l>
    {};

    template<class Container, std::size_t l>
    constexpr decltype(auto) _vectorEntry(Container&& container, level_v<l>, const std::vector<std::size_t>& index)
    requires (Concept::BracketIndexable<decltype(container),decltype(index.front())>)
    {

        return _vectorEntry(std::forward<Container>(container)[index[l]], level_v<l+1>(), index);
    }

    template<std::size_t l>
    constexpr decltype(auto) _vectorEntry(Dune::BlockVector<double>& container, level_v<l>, const std::vector<std::size_t>& index)
    requires (Concept::BracketIndexable<decltype(container),decltype(index.front())>)
    {
        return container[index[l]];
    }

    template<std::size_t l>
    constexpr decltype(auto) _vectorEntry(Dune::BlockVector<int>& container, level_v<l>, const std::vector<std::size_t>& index)
    requires (Concept::BracketIndexable<decltype(container),decltype(index.front())>)
    {
        return container[index[l]];
    }
  }

  /**
   * @brief Access an entry in a container
   * @details Use front-most index to return a bracket-indexed entry from the container
   *
   * @tparam Container        Container type
   * @tparam MultiIndex       Multi-index type
   * @param container         Container with the desired entry
   * @param rowindex          Multi-index to reach the target row
   * @param colindex          Multi-index to reach the target col
   * @return decltype(auto)   Entry at multi-index
   */
  template<class Container>
  constexpr decltype(auto) vectorEntry(Container&& container, const std::vector<std::size_t>& index)
  requires (Concept::BracketIndexable<decltype(container),decltype(index.front())>)
  {
    return impl::_vectorEntry(std::forward<Container>(container), impl::level_v<0>() , index);
  }
}


#endif // DUNE_COPASI_VECTORENTRY_HH
