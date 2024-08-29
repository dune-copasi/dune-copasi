#ifndef DUNE_COPASI_OPERATOR_COMMON_CONCEPTS_HH
#define DUNE_COPASI_OPERATOR_COMMON_CONCEPTS_HH

namespace Dune::Copasi{


  namespace Concept{

    template<class C, class Index>
    concept BracketIndexable = requires(C&& c, Index index)
    {
      c[index];
    };

    template<class C, class Index>
    concept DoubleBracketIndexable = requires(C&& c, Index index1, Index index2)
    {
      c[index1][index2];
    };

  }

}

#endif