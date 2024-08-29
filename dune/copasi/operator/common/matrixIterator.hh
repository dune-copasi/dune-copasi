#ifndef DUNE_COPASI_OPERATOR_COMMON_MATRIXITERATOR_HH
#define DUNE_COPASI_OPERATOR_COMMON_MATRIXITERATOR_HH

#include <dune/istl/bcrsmatrix.hh>

#include <cassert>
#include <tuple>

namespace Dune::Copasi{
  namespace impl{

    template<typename T, bool is_const, typename... Iterators>
    struct _extract_matrix_iterators;

    template<typename T, typename... Iterators>
    struct _extract_matrix_iterators<T,true,Iterators...>
        : public _extract_matrix_iterators<typename T::block_type,
                                    true,
                                    Iterators..., typename T::row_type::const_iterator
                                    >
    {};

    template<typename T, typename... Iterators>
    struct _extract_matrix_iterators<T,false,Iterators...>
        : public _extract_matrix_iterators<typename T::block_type,
                                    false,
                                    Iterators..., typename T::row_type::iterator
                                    >
    {};

    template<typename... Iterators>
    struct _extract_matrix_iterators< Dune::BCRSMatrix<double>,true, Iterators...>
    {
        typedef std::tuple<Iterators...,typename Dune::BCRSMatrix<double>::row_type::const_iterator> type;
    };

    template<typename... Iterators>
    struct _extract_matrix_iterators<Dune::BCRSMatrix<double>,false,Iterators...>
    {
        typedef std::tuple<Iterators...,typename Dune::BCRSMatrix<double>::row_type::iterator> type;
    };

    template<typename M>
    struct extract_matrix_iterators
        : public _extract_matrix_iterators<M,false>
    {};

    template<typename M>
    struct extract_matrix_iterators<const M>
        : public _extract_matrix_iterators<M,true>
    {};
    

    /**
     * @brief implementation: matrix iterator base type
     **/
    template<typename M>
    struct matrix_iterator_base
    {
      // iterator types
      using iterator_category = std::forward_iterator_tag;
      using value_type = double; //M::field_type;
      using difference_type = std::ptrdiff_t;
      using pointer = double*; //M::field_type*;
      using reference = double&; //M::field_type&;

      using matrix = M;
      using matrix_reference = M&;
      static const bool is_const = false;
    };

    /**
     * @brief implementation: const matrix iterator base type
     **/
    template<typename M>
    struct matrix_iterator_base<const M>
    {
      // iterator types
      using iterator_category = std::forward_iterator_tag;
      using value_type = double; //M::field_type;
      using difference_type = std::ptrdiff_t;
      using pointer = const double*; //M::field_type*;
      using reference = const double&; //M::field_type&;

      using matrix = M;
      using matrix_reference = const M&;
      static const bool is_const = true;
    };

  }

  template<typename M>
  class matrix_row_iterator 
        : public impl::matrix_iterator_base<M>
  {

    using Base = impl::matrix_iterator_base<M>;
    using matrix = typename Base::matrix;
    using matrix_reference = typename Base::matrix_reference;
    using Iterators = typename impl::extract_matrix_iterators<M>::type;
    static const bool is_const = Base::is_const;
    static const std::size_t depth = std::tuple_size<Iterators>::value;

    template<typename>
    friend class matrix_row_iterator;

  public: 
    /**
     * @brief Basic constructor taking a reference to the matrix and the index of the row (as std::array<std::size_t, N>)
     */
    matrix_row_iterator(M& m, std::array<std::size_t, depth> row_index) : _matrix_reference(m) , _row_index(row_index)
    {
      if (!_at_end)
        if (!start_row())
          _at_end = true;
    }

    /**
     * @brief Basic constructor taking a reference to the matrix and the index of the row (as std::vector<std::size_t>)
     */
    matrix_row_iterator(M& m, std::vector<std::size_t> row_index) : _matrix_reference(m)
    {
      std::move(row_index.begin(), row_index.begin() + depth, _row_index.begin());
      if (!_at_end)
        if (!start_row())
          _at_end = true;
    }

    /**
     * @brief returns bool indicating if referenced matrix is constant 
     * @return is constant bool
     */
    static bool constant(){
      return is_const;
    }

    /**
     * @brief returns bool indicating if iterator is at end
     * @return is at_end bool
     */
    bool at_end(){
      return _at_end;
    }

    /**
     * @brief arrow operator for pointer 
     * @return the current value by pointer
     */
    typename Base::pointer operator->() const
    {
      assert(!_at_end);
      return _current;
    }

    /**
     * @brief dereference the current value 
     * @return the current value by reference
     */
    typename Base::reference operator*() const
    {
      assert(!_at_end);
      return *_current;
    }

    /**
     * @brief increments the iterator 
     * @return the current iterator
     */
    matrix_row_iterator& operator++()
    {
      increment();
      return *this;
    }

  /*******************************************
   *  private section with member function
   *******************************************/
  private:
    template<std::size_t l>
    struct level
      : public std::integral_constant<std::size_t,l>
    {};

    /**
     * @brief Increment the iterator to the next element.
     * 
     * Asserts that the iterator is not at the end and attempts to advance to the next element.
     */    
    void increment()
    {
      assert(!_at_end);
      if (!advance(level<0>()))
        _at_end = true;
    }

    /**
     * @brief Advance the iterator at a given level.
     * 
     * @tparam l The level in the hierarchy.
     * @param level The level object.
     * @return true if the iterator was successfully advanced (at this level).
     * @return false if the iterator could not be advanced (at this level).
     */
    template<std::size_t l>
    bool advance(level<l>)
    {
      using iterator = typename std::tuple_element<l,Iterators>::type;
      iterator& it = std::get<l>(_iterators);
      iterator& end = std::get<l>(_end);

      if (advance(level<l+1>()))
        return true;

      ++it;

      while (it != end)
        {
          if (start(level<l+1>(),*it))
            return true;

          ++it;
        }

      return false;
    }

    /**
     * @brief Specialization of advance function for the deepest level.
     * 
     * @param level The level object.
     * @return true if the iterator was successfully advanced.
     * @return false if the iterator could not be advanced.
     */
    bool advance(level<depth-1>)
    {
      using iterator = typename std::tuple_element<depth-1,Iterators>::type;
      iterator& it = std::get<depth-1>(_iterators);
      const iterator& end = std::get<depth-1>(_end);

      ++it;

      if (it == end)
        return false;

      _current = &(*it);

      return true;
    }

    /**
     * @brief Start the iterator at a given level with a block.
     * 
     * @tparam l The level in the hierarchy.
     * @tparam Block The type of the block.
     * @param level The level object.
     * @param block The block to search for a valid starting iterator.
     * @return true if a valid starting iterator was found at this level.
     * @return false if no valid starting iterator was found at this level (-> end of level).
     */
    template<std::size_t l, typename Block>
    bool start(level<l>, Block& block)
    {
      using iterator =  typename std::tuple_element<l,Iterators>::type;
      iterator& it = std::get<l>(_iterators);
      iterator& end = std::get<l>(_end);

      it = block[_row_index[l]].begin();
      end = block[_row_index[l]].end();

      while (it != end)
        {
          if (start(level<l+1>(),*it))
            return true;

          ++it;
        }

      return false;
    }

    /**
     * @brief Start the iterator at the deepest level with a block.
     * 
     * @tparam Block The type of the block.
     * @param level The level object.
     * @param block The block to start the iterator with.
     * @return true if the begin iterator was successfully defined.
     * @return false if the begin iterator could not be defined (-> end of level).
     */
    template<typename Block>
    bool start(level<depth-1>, Block& block)
    {
      using iterator = std::tuple_element<depth-1,Iterators>::type;
      iterator& it = std::get<depth-1>(_iterators);
      iterator& end = std::get<depth-1>(_end);

      it = block[_row_index[depth-1]].begin();
      end = block[_row_index[depth-1]].end();

      if (it == end)
        return false;

      _current = &(*it);

      return true;
    }

    /**
     * @brief Start the iterator to run over the row of the matrix.
     * 
     * @return true if the begin iterator was successfully defined.
     * @return false if the begin iterator could not be defined (-> at end).
     */
    bool start_row(){
      return start(level<0>{}, _matrix_reference);
    }

  /*******************************************
   *  private section with member variables
   *******************************************/
  private: 
    // row index
    std::array<std::size_t, depth> _row_index;

    // iterator member variables
    bool _at_end = false;
    typename Base::pointer _current;
    matrix_reference _matrix_reference;
    Iterators _iterators;
    Iterators _end;
  };
}

#endif