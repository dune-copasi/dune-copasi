#ifndef DUNE_COPASI_MODEL_LOCAL_EQUATIONS_FUNCTOR_FACTORY_HH
#define DUNE_COPASI_MODEL_LOCAL_EQUATIONS_FUNCTOR_FACTORY_HH

#include <dune-copasi-config.hh>

#include <dune/common/fmatrix.hh>
#include <dune/common/fvector.hh>
#include <dune/common/parametertree.hh>

#include <function2/function2.hpp>

#include <functional>
#include <string>

namespace Dune::Copasi {

template<std::size_t dim>
class LocalDomain;

// the owner of the resulting function must hold the local equations for the lifetime of created
// functors
template<std::size_t dim>
class FunctorFactory
{
public:
  using Scalar = FieldVector<double, 1>;
  using Vector = FieldVector<double, dim>;
  using Tensor = FieldMatrix<double, dim, dim>;

  using ScalarFunctor = fu2::unique_function<Scalar() const noexcept>;
  using VectorFunctor = fu2::unique_function<Vector() const noexcept>;
  using TensorApplyFunctor =
    fu2::unique_function<Vector(Vector) const noexcept>;

  FunctorFactory() = default;
  FunctorFactory(const FunctorFactory&) = delete;
  FunctorFactory(FunctorFactory&&) = delete;

  virtual ~FunctorFactory() = default;

  [[nodiscard]] virtual ScalarFunctor make_scalar(std::string_view,
                                                  const ParameterTree&,
                                                  const LocalDomain<dim>&,
                                                  int /*codim*/ = 0) const = 0;

  [[nodiscard]] virtual VectorFunctor make_vector(std::string_view,
                                                  const ParameterTree&,
                                                  const LocalDomain<dim>&,
                                                  int /*codim*/ = 0) const = 0;

  [[nodiscard]] virtual TensorApplyFunctor make_tensor_apply(std::string_view,
                                                             const ParameterTree&,
                                                             const LocalDomain<dim>&,
                                                             int /*codim*/ = 0) const = 0;
};

} // namespace Dune::Copasi

#endif // DUNE_COPASI_MODEL_LOCAL_EQUATIONS_FUNCTOR_FACTORY_HH
