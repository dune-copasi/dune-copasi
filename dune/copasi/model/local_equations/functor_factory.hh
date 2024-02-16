#ifndef DUNE_COPASI_MODEL_LOCAL_EQUATIONS_FUNCTOR_FACTORY_HH
#define DUNE_COPASI_MODEL_LOCAL_EQUATIONS_FUNCTOR_FACTORY_HH

#include <dune-copasi-config.hh>
#include <dune/copasi/concepts/grid.hh>
#include <dune/copasi/parser/context.hh>
#include <dune/copasi/parser/grid_context.hh>

#include <dune/common/fmatrix.hh>
#include <dune/common/fvector.hh>
#include <dune/common/parametertree.hh>

#include <function2/function2.hpp>

#include <functional>
#include <string>

namespace Dune::Copasi {

template<std::size_t dim>
struct LocalDomain;

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
  using TensorApplyFunctor = fu2::unique_function<Vector(Vector) const noexcept>;

  FunctorFactory() = default;
  FunctorFactory(const FunctorFactory&) = delete;
  FunctorFactory(FunctorFactory&&) = delete;

  FunctorFactory& operator=(const FunctorFactory&) = delete;
  FunctorFactory& operator=(FunctorFactory&&) = delete;

  virtual ~FunctorFactory() = default;

  // ---------------------------------------------------------------------------
  // Defines the interface for making functors that can be evaluated
  // ---------------------------------------------------------------------------
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


  virtual std::shared_ptr<const ParserContext> parser_context() const = 0;

  // ---------------------------------------------------------------------------
  // Defines the interface for updating grid entity specific values
  // ---------------------------------------------------------------------------
  virtual double get_gmsh_id( std::any entity) const = 0;

  virtual void update_grid_data(std::unordered_map<std::string, double>& cell_data, std::any entity) const = 0;

  virtual const std::unordered_map<std::string, std::function<double*(std::size_t)>>& get_cell_functor() const = 0;

};

} // namespace Dune::Copasi

#endif // DUNE_COPASI_MODEL_LOCAL_EQUATIONS_FUNCTOR_FACTORY_HH
