#ifndef DUNE_COPASI_UTIL_META_HH
#define DUNE_COPASI_UTIL_META_HH

#include <dune/pdelab/function/callableadapter.hh>
#include <dune/pdelab/common/function.hh>

#include <utility>

namespace Dune::Copasi {

// TODO: Try this with concepts

template<class GV, class F, class = void>
struct is_pdelab_callable : std::false_type {};

template<class GV, class F>
struct is_pdelab_callable<GV,F,
                      std::enable_if_t<
                        AlwaysTrue <
                          decltype(std::declval<PDELab::makeGridFunctionFromCallable(std::declval<GV>(),
                                                                                     std::declval<F>()  )>() )
                        >::value >>
      : std::true_type {};

template<class GF, class = void>
struct is_grid_function : std::false_type {};


template<class GF>
struct is_grid_function<GF,
            Dune::PDELab::GridFunctionInterface<typename GF::Traits,GF>>
      : std::true_type {};

} // Dune::Copasi namespace

#endif // DUNE_COPASI_UTIL_META_HH