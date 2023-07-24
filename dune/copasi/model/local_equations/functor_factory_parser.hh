#ifndef DUNE_COPASI_MODEL_LOCAL_EQUATIONS_FUNCTOR_FACTORY_PARSER_HH
#define DUNE_COPASI_MODEL_LOCAL_EQUATIONS_FUNCTOR_FACTORY_PARSER_HH

#include <dune/copasi/model/local_equations/functor_factory.hh>
#include <dune/copasi/parser/context.hh>

#include <memory>

namespace Dune::Copasi {

template<std::size_t dim>
class FunctorFactoryParser final : public FunctorFactory<dim>
{
public:
  using Scalar = FieldVector<double, 1>;
  using Vector = FieldVector<double, dim>;
  using Tensor = FieldMatrix<double, dim, dim>;

  using ScalarFunctor = typename FunctorFactory<dim>::ScalarFunctor;
  using VectorFunctor = typename FunctorFactory<dim>::VectorFunctor;
  using TensorApplyFunctor = typename FunctorFactory<dim>::TensorApplyFunctor;

  explicit FunctorFactoryParser(ParserType parser_type = default_parser,
                                std::shared_ptr<const ParserContext> parser_context = nullptr)
    : FunctorFactory<dim>()
    , _parser_type{ parser_type }
    , _parser_context{ std::move(parser_context) }
  {
  }

  ~FunctorFactoryParser() override = default;

  [[nodiscard]] ScalarFunctor make_scalar(std::string_view /*prefix*/,
                                          const ParameterTree& /*config*/,
                                          const LocalDomain<dim>& /*local_domain*/,
                                          bool /*is_membrane_expression*/) const noexcept override;

  [[nodiscard]] VectorFunctor make_vector(std::string_view /*prefix*/,
                                          const ParameterTree& /*config*/,
                                          const LocalDomain<dim>& /*local_domain*/,
                                          bool /*is_membrane_expression*/) const noexcept override;

  [[nodiscard]] TensorApplyFunctor make_tensor_apply(
    std::string_view /*prefix*/,
    const ParameterTree& /*config*/,
    const LocalDomain<dim>& /*local_domain*/,
    bool /*is_membrane_expression*/) const noexcept override;


private:
  [[nodiscard]] ScalarFunctor parse_scalar_expression(
    const ParameterTree& /*config*/,
    const LocalDomain<dim>& /*local_values*/,
    bool /*is_membrane_expression*/) const noexcept;

  ParserType _parser_type;
  std::shared_ptr<const ParserContext> _parser_context;
};

} // namespace Dune::Copasi

#ifndef DUNE_COPASI_PRECOMPILED_MODE
#include <dune/copasi/model/local_equations/functor_factory_parser.impl.hh>
#endif

#endif // DUNE_COPASI_MODEL_LOCAL_EQUATIONS_FUNCTOR_FACTORY_PARSER_HH
