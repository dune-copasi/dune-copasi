#ifndef DUNE_COPASI_LOCAL_FINITE_ELEMENT_CACHE_HH
#define DUNE_COPASI_LOCAL_FINITE_ELEMENT_CACHE_HH

#include <dune-copasi-config.h>

#include <dune/geometry/quadraturerules.hh>

#include <span>
#include <typeindex>

namespace Dune::Copasi {

/**
 * @brief      This class describes a local basis cache.
 * @details    This class cache results of the local basis for different local basis with same
 traits.
 * @ingroup    FiniteElement
 *
 * @tparam     LocalBasisTraits  Local basis traits
 */
template<class LocalBasisTraits>
class LocalBasisCache
{
  using Domain = typename LocalBasisTraits::DomainType;
  using DomainField = typename LocalBasisTraits::DomainFieldType;

  using Key = std::type_index;
  using QuadratureRule = Dune::QuadratureRule<DomainField, LocalBasisTraits::dimDomain>;

  using Range = typename LocalBasisTraits::RangeType;
  using RangeField = typename LocalBasisTraits::RangeFieldType;
  using Jacobian = typename LocalBasisTraits::JacobianType;

public:

  /**
   * @brief      Constructs a new instance with empty values.
   */
  LocalBasisCache()
    : _bound_basis{ typeid(void) }
  {
    unbind();
  }

  /**
   * @brief      Binds a finite element to this cache
   * @details    The finite element gets evaluated if necessary
   *
   * @param[in]  finite_element  The finite element
   *
   * @tparam     FiniteElement   Type of the finite element
   */
  template<class FiniteElement>
  void bind(const FiniteElement& finite_element)
  {
    // assume that the type uniquely identifies the finite element
    const std::type_index& basis_id = typeid(finite_element.localBasis());

    if (_bound_basis == basis_id) {
      return;
    }

    _fe_size = finite_element.size();
    _bound_basis = basis_id;

    std::size_t const order = 4;
    _bound_quadrature = &QuadratureRules<DomainField, LocalBasisTraits::dimDomain>::rule(
      finite_element.type(), order);

    auto rit = _range.find(_bound_basis);
    if (rit != _range.end()) {
      _bound_range = rit->second.data();
    } else {
      _bound_range = cache_evaluate(finite_element.localBasis());
    }

    auto jit = _jacobian.find(_bound_basis);
    if (jit != _jacobian.end()) {
      _bound_jacobian = jit->second.data();
    } else {
      _bound_jacobian = cache_jacobian(finite_element.localBasis());
    }
  }

  /**
   * @brief      Unbind finite element from this cache
   */
  void unbind()
  {
    _bound_basis = typeid(void);
    _bound_quadrature = nullptr;
    _bound_range = nullptr;
    _bound_jacobian = nullptr;
  }

  /**
   * @brief      Quadrature rule used to cache the finite element
   *
   * The index of the positions in the quadrature rule are the position_id
   * used to fetch the cached values.
   *
   * @return     Quadrature rule used to cache the finite element
   */
  const QuadratureRule& rule() const {
    assert(_bound_quadrature);
    return *_bound_quadrature;
  }

  /**
   * @brief      Evaluate function with the bound finite element
   *
   * @param[in]  position_id  The position id to be evaluated
   *
   * @return     The evaluation for each of the local basis at position_id
   */
  std::span<const Range> evaluateFunction(std::size_t position_id) const
  {
    assert(_bound_range);
    assert(position_id < _bound_quadrature->size());
    return {_bound_range + _fe_size * position_id, _fe_size};
  }

  /**
   * @brief      Evaluate jacobian with the bound finite element
   *
   * @param[in]  position_id  The position id to be evaluated
   *
   * @return     The jacobian for each of the local basis at position_id
   */
  std::span<const Jacobian> evaluateJacobian(std::size_t position_id) const
  {
    assert(_bound_jacobian);
    assert(position_id < _bound_quadrature->size());
    return {_bound_jacobian + _fe_size * position_id, _fe_size};
  }

private:

  template<class Basis>
  Range* cache_evaluate(const Basis& basis)
  {
    std::vector<Range> quad_range;
    std::vector<Range> range;
    for (const auto& quad : *_bound_quadrature) {
      basis.evaluateFunction(quad.position(), quad_range);
      range.insert(end(range), begin(quad_range), end(quad_range));
    }
    assert(range.size() == _fe_size * _bound_quadrature->size());
    const auto insert_pair = _range.emplace(_bound_basis, std::move(range));
    assert(insert_pair.second);
    return insert_pair.first->second.data();
  }

  template<class Basis>
  Jacobian* cache_jacobian(const Basis& basis)
  {
    std::vector<Jacobian> quad_jacobian;
    std::vector<Jacobian> jacobian;
    for (const auto& quad : *_bound_quadrature) {
      basis.evaluateJacobian(quad.position(), quad_jacobian);
      jacobian.insert(end(jacobian),begin(quad_jacobian), end(quad_jacobian));
    }
    assert(jacobian.size() == _fe_size * _bound_quadrature->size());
    const auto insert_pair = _jacobian.emplace(_bound_basis, std::move(jacobian));
    assert(insert_pair.second);
    return insert_pair.first->second.data();
  }

  std::type_index _bound_basis;
  std::size_t _fe_size{};
  QuadratureRule const * _bound_quadrature;
  Range const * _bound_range;
  Jacobian const * _bound_jacobian;
  std::unordered_map<Key,std::vector<Range>> _range;
  std::unordered_map<Key,std::vector<Jacobian>> _jacobian;
};

} // namespace Dune::Copasi

#endif // DUNE_COPASI_LOCAL_FINITE_ELEMENT_CACHE_HH
