#ifndef DUNE_COPASI_MATRIXENTRY_HH
#define DUNE_COPASI_MATRIXENTRY_HH

#include <dune/copasi/operator/common/concepts.hh>

namespace Dune::Copasi{

  namespace impl{

    template<std::size_t l>
    struct level
        : public std::integral_constant<std::size_t,l>
    {};

    template<class Container, std::size_t l>
    constexpr decltype(auto) _matrixEntry(Container&& container, level<l>, const std::vector<std::size_t>& rowindex, const std::vector<std::size_t>& colindex)
    requires (Concept::DoubleBracketIndexable<decltype(container),decltype(rowindex.front())>)
    {
        return _matrixEntry(std::forward<Container>(container)[rowindex[l]][colindex[l]], level<l+1>(), rowindex, colindex);
    }

    template<std::size_t l>
    constexpr decltype(auto) _matrixEntry(Dune::BCRSMatrix<double>& container, level<l>, const std::vector<std::size_t>& rowindex, const std::vector<std::size_t>& colindex)
    requires (Concept::DoubleBracketIndexable<decltype(container),decltype(rowindex.front())>)
    {
        return container.entry(rowindex[l],colindex[l]);
    }

    template<class Container, std::size_t l>
    constexpr decltype(auto) _matrixContainer(Container&& container, level<l>, const std::vector<std::size_t>& rowindex, const std::vector<std::size_t>& colindex)
    requires (Concept::DoubleBracketIndexable<decltype(container),decltype(rowindex.front())>)
    {
        return _matrixContainer(std::forward<Container>(container)[rowindex[l]][colindex[l]], level<l+1>(), rowindex, colindex);
    }

    template<std::size_t l>
    constexpr decltype(auto) _matrixContainer(Dune::BCRSMatrix<double>& container, level<l>, const std::vector<std::size_t>& rowindex, const std::vector<std::size_t>& colindex)
    requires (Concept::DoubleBracketIndexable<decltype(container),decltype(rowindex.front())>)
    {
        return container;
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
  constexpr decltype(auto) matrixEntry(Container&& container, const std::vector<std::size_t>& rowindex, const std::vector<std::size_t>& colindex)
  requires (Concept::DoubleBracketIndexable<decltype(container),decltype(rowindex.front())>)
  {
    return impl::_matrixEntry(std::forward<Container>(container), impl::level<0>() ,rowindex, colindex);
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
  constexpr decltype(auto) matrixContainer(Container&& container, const std::vector<std::size_t>& rowindex, const std::vector<std::size_t>& colindex)
  requires (Concept::DoubleBracketIndexable<decltype(container),decltype(rowindex.front())>)
  {
    return impl::_matrixContainer(std::forward<Container>(container), impl::level<0>() ,rowindex, colindex);
  }


}

#endif // DUNE_COPASI_MATRIXENTRY_HH
