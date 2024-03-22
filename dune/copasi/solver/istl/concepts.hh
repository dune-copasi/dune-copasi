#ifndef DUNE_COPASI_SOLVER_ISTL_CONCEPTS_HH
#define DUNE_COPASI_SOLVER_ISTL_CONCEPTS_HH

#include <dune/istl/solver.hh>
#include <dune/istl/preconditioner.hh>


namespace Dune::Copasi::ISTL::Concept {

template<class O>
concept LinearOperator = std::derived_from<O,Dune::LinearOperator<typename O::domain_type, typename O::range_type>>;

template<class O>
concept AssembledLinearOperator = std::derived_from<O,Dune::AssembledLinearOperator<typename O::matrix_type, typename O::domain_type, typename O::range_type>>;

} // Dune::Copasi::ISTL::Concept

#endif // DUNE_COPASI_SOLVER_ISTL_CONCEPTS_HH
