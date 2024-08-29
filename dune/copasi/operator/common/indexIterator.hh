#ifndef DUNE_COPASI_OPERATOR_COMMON_INDEXITERATOR_HH
#define DUNE_COPASI_OPERATOR_COMMON_INDEXITERATOR_HH

#include <dune/copasi/operator/common/concepts.hh>
#include <dune/pdelab/common/multiindex.hh>

namespace Dune::Copasi{


  namespace impl{

    template<typename V>
    struct index_iterator_base
    {
      // iterator types
      using iterator_category = std::forward_iterator_tag;
      using value_type = V::field_type;
      using difference_type = std::ptrdiff_t;
      using pointer = V::field_type*;
      using reference = V::field_type&;

      // specific types
      typedef V vector;
      typedef V& vector_reference;
      static const bool is_const = false;
    };

    template<typename V>
    struct index_iterator_base<const V>
    {
      // iterator types
      using iterator_category = std::forward_iterator_tag;
      using value_type = V::field_type;
      using difference_type = std::ptrdiff_t;
      using pointer = const V::field_type*;
      using reference = const V::field_type&;

      // specific types
      typedef V vector;
      typedef const V& vector_reference;
      static const bool is_const = true;
    };

  }

  template<typename V>
  class index_iterator : public impl::index_iterator_base<V>
  {
    typedef impl::index_iterator_base<V> BaseT;
    typedef typename BaseT::vector vector;
    typedef typename BaseT::vector_reference vector_reference;
    static const bool is_const = BaseT::is_const;

    using MultiIndex = typename Dune::PDELab::MultiIndex<std::size_t, 10>;

  public:

    index_iterator(vector_reference vector)
      : _vector(vector)
    {
      initialize();
    }

    std::vector<std::size_t> index(){
      return _index;
    }

    MultiIndex multiIndex(){
      MultiIndex mi(_index.begin(), _index.end());
      return mi;
    }

    index_iterator& operator++()
    {
      increment();
      return *this;
    }

    typename BaseT::pointer operator->() const
    {
      assert(!_at_end);
      return &get_current();
    }

    typename BaseT::reference operator*() const
    {
      assert(!_at_end);
      return get_current();
    }

    /**
     * @brief returns bool indicating if iterator is at end
     * @return is at_end bool
     */
    bool at_end(){
      return _at_end;
    }

    void print(){
      std::cout << "---- index_iterator ---- \n";
      std::cout << "index: ";
      if( _index.empty() )
        std::cout << "** empty ** -> reached end? " << _at_end << std::endl;
      else{
        for(std::size_t l = 0; l < _index.size() ; ++l){
          if( l < _index.size()-1)
            std::cout <<  _index[l] << ", ";
          else
            std::cout << _index[l] << std::endl;
        }
        std::cout << "value: " << get_current() << std::endl;
      }
      std::cout << std::endl;
    }
      
  private:

    void initialize(){

      // clear the index vectors
      _index.clear();
 
      // check if not empty and flag if so
      if(not start(0, _vector))
        _at_end = true;
    }

    typename BaseT::reference get_current() const {
      return get_current(0, _vector);
    }

    template<typename Block>
    typename BaseT::reference get_current(std::size_t l, Block& block) const{
      return block;
    }

    template<typename Block>
    typename BaseT::reference get_current(std::size_t l, Block& block) const
    requires Dune::Copasi::Concept::BracketIndexable<Block, std::size_t> 
    {
      return get_current(l+1, block[_index[l]] );
    }


    // Default fallback if Block is no longer indexable -> basic entry
    template<typename Block>
    bool start(std::size_t l, Block& block)
    { 
      return true;
    }

    // If indexable then add a zero and parse the front of the block
    template<typename Block>
    bool start(std::size_t l, Block& block)
    requires Dune::Copasi::Concept::BracketIndexable<Block, std::size_t> 
    {
      // add the zero index 
      _index.emplace_back(0);
      while (_index[l] != block.N())
      {
        if (start(l+1,block[_index[l]]))
          return true;

        ++_index[l];
      }
      // remove last index
      _index.pop_back();
      return false;
    }

    void increment()
    {
      assert(!_at_end); // the end of the vector is reached.
      if (!advance(0, _vector))
        _at_end = true;
    }

    // if the container Block is not indexable then move
    // to the next element in the container (we do not
    // check here if it exists ! )
    template<typename Block>
    bool advance(std::size_t l, Block& block){
      ++_index[l-1];
      return true;
    }

    template<typename Block>
    bool advance(std::size_t l, Block& block)
    requires Concept::BracketIndexable<Block, std::size_t> 
    {
      // Check if the container is non-empty
      if(block.N() == 0){
        _index.pop_back();
        return false;
      }

      if (advance(l+1, block[ _index[l] ])){
        if( _index[l] < block.N() ){
          return true;
        } else {
          _index.pop_back();
          return false;
        }
      }

      ++_index[l];

      while ( _index[l] != block.N())
        {
          if (start(l+1,block[ _index[l] ]))
            return true;

          ++_index[l];
        }
      _index.pop_back();
      return false;
    }

    // members variables 
    std::vector<std::size_t> _index;
    vector_reference _vector;
    bool _at_end = false;
  };
  
}


#endif // DUNE_ECGI_INDEXITERATROR_HH
