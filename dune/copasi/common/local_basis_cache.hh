#ifndef DUNE_COPASI_LOCAL_FINITE_ELEMENT_CACHE_HH
#define DUNE_COPASI_LOCAL_FINITE_ELEMENT_CACHE_HH

#include <dune/localfunctions/common/virtualinterface.hh>
#include <dune/localfunctions/common/virtualwrappers.hh>

#include <dune/common/hash.hh>

#include <functional>
#include <unordered_map>
#include <vector>
#include <memory>
#include <typeindex>

namespace Dune::Copasi {

/**
 * @brief      This class describes a quadrature point cache key.
 *
 * @tparam     Domain  Domain type to cache
 */
template<class Domain>
class QuadraturePointCacheKey
{
public:

  /**
   * @brief      Constructs a new instance.
   *
   * @param[in]  position  The position
   * @param[in]  basis_id  The basis identifier
   */
  QuadraturePointCacheKey(const Domain& position, const std::type_index& basis_id)
    : _position(position)
    , _basis_id(basis_id)
  {}

  /**
   * @brief      Constructs a new instance.
   *
   * @param[in]  position    The position
   * @param[in]  basis       The local basis
   *
   * @tparam     LocalBasis  Local basis type
   */
  template<class LocalBasis>
  QuadraturePointCacheKey(const Domain& position, const LocalBasis& basis)
    : QuadraturePointCacheKey(position,typeid(basis))
  {}

  /**
   * @brief      Returns the hash code for this object.
   *
   * @return     The hash code value for this object
   */
  std::size_t hash_code() const noexcept
  {
    std::size_t seed = 46'255'207; // some prime number
    Dune::hash_combine(seed,_basis_id);
    for (auto&& i : _position)
      Dune::hash_combine(seed,i);
    return seed;
  }

  /**
   * @brief      Equality operator.
   *
   * @param[in]  lhs   The left hand side
   * @param[in]  rhs   The right hand side
   *
   * @return     The result of the equality
   */
  friend inline bool
  operator==(const QuadraturePointCacheKey& lhs, const QuadraturePointCacheKey& rhs)
  {
    return (lhs._basis_id == rhs._basis_id) and (lhs._position == rhs._position);
  }

  /**
   * @brief      Inequality operator.
   *
   * @param[in]  lhs   The left hand side
   * @param[in]  rhs   The right hand side
   *
   * @return     The result of the inequality
   */
  friend inline bool
  operator!=(const QuadraturePointCacheKey& lhs, const QuadraturePointCacheKey& rhs)
  {
    return not (lhs == rhs);
  }

private:
  const Domain _position;
  const std::type_index _basis_id;
};

} // namespace Dune::Copasi

namespace std {

  /**
   * @brief      Standard overload of the hash function
   *
   * @tparam     Domain  Domain type
   */
  template<typename Domain>
  struct hash<Dune::Copasi::QuadraturePointCacheKey<Domain>>
  {
    std::size_t operator()(const Dune::Copasi::QuadraturePointCacheKey<Domain>& key) const
    {
      return key.hash_code();
    }
  };
} // namespace std

namespace Dune::Copasi {

/**
 * @brief      This class describes a local basis cache.
 * @details    This class cache results of the local basis for different local basis with same traits.
 * @warning    First tests show that this is relatively slower that the local basis cache provided by PDELab.
 * @todo       Improve preformance!
 * @ingroup    FiniteElement
 *
 * @tparam     LocalBasisTraits  Local basis traits
 */
template<class LocalBasisTraits>
class LocalBasisCache
{
  using Domain = typename LocalBasisTraits::DomainType;

  using BasisKey = std::type_index;
  using PointKey = QuadraturePointCacheKey<Domain>;

  using Range = typename LocalBasisTraits::RangeType;
  using RangeField = typename LocalBasisTraits::RangeFieldType;
  using Jacobian = typename LocalBasisTraits::JacobianType;

  using FiniteElementInterface = Dune::LocalFiniteElementVirtualInterface<LocalBasisTraits>;

public:
  /**
   * @brief      Constructs a new instance.
   */
  LocalBasisCache()
    : _bound_basis(typeid(void))
    , _range(std::make_shared<std::unordered_map<PointKey,std::vector<Range>>>())
    , _jacobian(std::make_shared<std::unordered_map<PointKey,std::vector<Jacobian>>>())
  {}

  /**
   * @brief      Copy constructor
   * @details    This makes caches to share underlying data. This is useful to bind to different 
   *             finite elements but still using the same cache
   * 
   * @param[in]  other  Local basis cache
   */
  LocalBasisCache(const LocalBasisCache& other)
    : _bound_basis(typeid(void))
  {
    // copy pointers to actual cache data
    _range = other._range;
    _jacobian = other._jacobian;
  }

  /**
   * @brief      Binds a finite element to this cache
   * @details    The finite element gets virtualized to be evaluated if necessary
   * 
   * @param[in]  finite_element  The finite element
   *
   * @tparam     FiniteElement   Type of the finite element
   */
  template<class FiniteElement>
  void bind(const FiniteElement& finite_element)
  {
    const std::type_index basis_id = typeid(finite_element.localBasis());

    if (_bound_basis == basis_id)
      return;

    _bound_basis = basis_id;

    using BasisTraits = typename FiniteElement::Traits::LocalBasisType::Traits;
    static_assert(std::is_same_v<LocalBasisTraits,BasisTraits>);
    using FiniteElementWrapper = Dune::LocalFiniteElementVirtualImp<FiniteElement>;

    _finite_element = std::make_unique<FiniteElementWrapper>(finite_element);
  }

  /**
   * @brief      Unbind finite element from this cache
   */
  void unbind()
  {
    _bound_basis = typeid(void);
    _finite_element.release();
  }

  /**
   * @brief      Evaluate function with the bound finite element
   *
   * @param[in]  position  The position to be evaluated
   *
   * @return     The evaluation for each of the local basis
   */
  inline const std::vector<Range>& evaluateFunction(const Domain& position) const
  {
    assert(_bound_basis != typeid(void));
    return evaluateFunction(position,_bound_basis);
  }

  /**
   * @brief      Evaluate jacobian with the bound finite element
   *
   * @param[in]  position  The position to be evaluated
   *
   * @return     The jacobian evaluation for each of the local basis
   */
  inline const std::vector<Jacobian>& evaluateJacobian(const Domain& position) const
  {
    assert(_bound_basis != typeid(void));
    return evaluateJacobian(position,_bound_basis);
  }

private:

  /**
   * @brief      Cache the evaluate method
   *
   * @param[in]  position  The position
   * @param[in]  basis     The local basis
   * @param[in]  basis_id  The local basis identifier
   *
   * @tparam     Basis     Local basis type
   *
   * @return     A reference to the cached evaluations
   */
  template<class Basis>
  const std::vector<Range>& cache_evaluate(const Domain& position, const Basis& basis, const std::type_index& basis_id) const
  {
    std::vector<Range> range;
    basis.evaluateFunction(position,range);
    auto&& range_pair = std::make_pair(PointKey{position,basis_id},std::move(range));
    const auto insert_pair = _range->insert(range_pair);
    assert(insert_pair.second);
    return insert_pair.first->second;
  }

  /**
   * @brief      Cache the jacobian method
   *
   * @param[in]  position  The position
   * @param[in]  basis     The local basis
   * @param[in]  basis_id  The local basis identifier
   *
   * @tparam     Basis     Local basis type
   *
   * @return     A reference to the cached jacobian evaluations
   */
  template<class Basis>
  const std::vector<Jacobian>& cache_jacobian(const Domain& position, const Basis& basis, const std::type_index& basis_id) const
  {
    std::vector<Jacobian> jacobian;
    basis.evaluateJacobian(position,jacobian);
    auto&& jacobian_pair = std::make_pair(PointKey{position,basis_id},std::move(jacobian));
    const auto insert_pair = _jacobian->insert(jacobian_pair);
    assert(insert_pair.second);
    return insert_pair.first->second;
  }

  /**
   * @brief      Cache the evaluate method
   * @details    If basis_id is already cached returs the stored value, 
   *             otherwise evaluates and caches it and returns the new value
   *
   * @param[in]  position  The position
   * @param[in]  basis_id  The local basis identifier
   *
   * @return     A reference to the cached evaluations
   */
  const std::vector<Range>& evaluateFunction(const Domain& position, const std::type_index& basis_id) const
  {
    auto it = _range->find(PointKey{position,basis_id});
    if (it != _range->end())
      return it->second;
    else
      return cache_evaluate(position,_finite_element->localBasis(),basis_id);
  }
  
  /**
   * @brief      Cache the jacobian method
   * @details    If basis_id is already cached returs the stored value, 
  *              otherwise evaluates and caches it and returns the new value
   *
   * @param[in]  position  The position
   * @param[in]  basis_id  The local basis identifier
   *
   * @return     A reference to the cached jacobian
   */
  const std::vector<Jacobian>& evaluateJacobian(const Domain& position, const std::type_index& basis_id) const
  {
    auto it = _jacobian->find(PointKey{position,basis_id});
    if (it != _jacobian->end())
      return it->second;
    else
      return cache_jacobian(position,_finite_element->localBasis(),basis_id);
  }

  std::type_index _bound_basis;
  std::unique_ptr<FiniteElementInterface> _finite_element;
  std::shared_ptr<std::unordered_map<PointKey,std::vector<Range>>> _range;
  std::shared_ptr<std::unordered_map<PointKey,std::vector<Jacobian>>> _jacobian;
};

} // namespace Dune::Copasi

#endif // DUNE_COPASI_LOCAL_FINITE_ELEMENT_CACHE_HH