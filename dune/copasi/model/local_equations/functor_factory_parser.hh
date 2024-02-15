#ifndef DUNE_COPASI_MODEL_LOCAL_EQUATIONS_FUNCTOR_FACTORY_PARSER_HH
#define DUNE_COPASI_MODEL_LOCAL_EQUATIONS_FUNCTOR_FACTORY_PARSER_HH

#include <dune/copasi/model/local_equations/functor_factory.hh>


#include <memory>
#include <fstream>
#include <iostream>

namespace Dune::Copasi {

template<Dune::Concept::Grid Grid>
class FunctorFactoryParser final : public FunctorFactory<Grid::dimensionworld>
{
public:

  static constexpr int dim = Grid::dimensionworld;

  using Scalar = FieldVector<double, 1>;
  using Vector = FieldVector<double, dim>;
  using Tensor = FieldMatrix<double, dim, dim>;

  using ScalarFunctor = typename FunctorFactory<dim>::ScalarFunctor;
  using VectorFunctor = typename FunctorFactory<dim>::VectorFunctor;
  using TensorApplyFunctor = typename FunctorFactory<dim>::TensorApplyFunctor;

  explicit FunctorFactoryParser(ParserType parser_type = default_parser,
                                std::shared_ptr<const ParserContext> parser_context = nullptr,
                                std::shared_ptr<const ParserGridContext<Grid>> parser_grid_context = nullptr)
    : FunctorFactory<dim>()
    , _parser_type{ parser_type }
    , _parser_context{ std::move(parser_context)}
    , _parser_grid_context{ std::move(parser_grid_context) }
  {
  }

  FunctorFactoryParser(const FunctorFactoryParser&) = delete;
  FunctorFactoryParser(FunctorFactoryParser&&) = delete;

  FunctorFactoryParser& operator=(const FunctorFactoryParser&) = delete;
  FunctorFactoryParser& operator=(FunctorFactoryParser&&) = delete;

  ~FunctorFactoryParser() override = default;

  [[nodiscard]] ScalarFunctor make_scalar(std::string_view /*prefix*/,
                                          const ParameterTree& /*config*/,
                                          const LocalDomain<dim>& /*local_domain*/,
                                          int /*codim*/ = 0) const override;

  [[nodiscard]] VectorFunctor make_vector(std::string_view /*prefix*/,
                                          const ParameterTree& /*config*/,
                                          const LocalDomain<dim>& /*local_domain*/,
                                          int /*codim*/ = 0) const override;

  [[nodiscard]] TensorApplyFunctor make_tensor_apply( std::string_view /*prefix*/,
                                                      const ParameterTree& /*config*/,
                                                      const LocalDomain<dim>& /*local_domain*/,
                                                      int /*codim*/ = 0) const override;

  std::shared_ptr<const ParserContext> parser_context() const { return _parser_context; };

  double get_gmsh_id( std::any entity) const{
    return _parser_grid_context->get_gmsh_id(entity);
  }

  void update_grid_data(std::unordered_map<std::string, double>& cell_data, std::any entity) const {
    _parser_grid_context->update_grid_data( cell_data, entity);
  }

  const std::unordered_map<std::string, std::function<double*(std::size_t)>>& get_cell_functor() const {
    return _parser_grid_context->get_cell_functor();
  }

private:
  [[nodiscard]] ScalarFunctor parse_scalar_expression(const ParameterTree& /*config*/,
                                                      const LocalDomain<dim>& /*local_values*/,
                                                      int /*codim*/) const;

  ParserType _parser_type;
  std::shared_ptr<const ParserContext> _parser_context;
  std::shared_ptr<const ParserGridContext<Grid>> _parser_grid_context;
};

} // namespace Dune::Copasi

#ifndef DUNE_COPASI_PRECOMPILED_MODE
#include <dune/copasi/model/local_equations/functor_factory_parser.impl.hh>
#endif

#endif // DUNE_COPASI_MODEL_LOCAL_EQUATIONS_FUNCTOR_FACTORY_PARSER_HH
