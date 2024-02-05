#ifndef DUNE_COPASI_LOCAL_FINITE_ELEMENT_CACHE_HH
#define DUNE_COPASI_LOCAL_FINITE_ELEMENT_CACHE_HH

#include <dune-copasi-config.hh>

#include <dune/copasi/common/exceptions.hh>

#include <dune/geometry/quadraturerules.hh>
#include <dune/geometry/referenceelements.hh>

#include <dune/common/exceptions.hh>
#include <dune/common/hash.hh>

#include <span>
#include <typeindex>

namespace Dune::Copasi::Impl {

//! Key for the local basis cache
template<typename T>
struct LocalBasisCahceKey
{
  LocalBasisCahceKey(std::type_index fem_id = typeid(void),
                     void const* quad_id = nullptr,
                     T proj_id = {})
    : _fem_id{ fem_id }
    , _quad_id{ quad_id }
    , _proj_id{ proj_id }
  {
  }
  std::type_index _fem_id;
  void const* _quad_id;
  T _proj_id;

  friend inline std::size_t hash_value(const LocalBasisCahceKey& key)
  {
    std::size_t hash = 233;
    hash_combine(hash, key._fem_id);
    hash_combine(hash, key._quad_id);
    hash_range(hash, key._proj_id.begin(), key._proj_id.end());
    return hash;
  }

  bool operator==(const LocalBasisCahceKey&) const = default;
};

} // namespace Dune::Copasi::Impl

DUNE_DEFINE_HASH(DUNE_HASH_TEMPLATE_ARGS(typename T),
                 DUNE_HASH_TYPE(Dune::Copasi::Impl::LocalBasisCahceKey<T>))

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
  static constexpr int dimDomain = LocalBasisTraits::dimDomain;

  using Range = typename LocalBasisTraits::RangeType;
  using RangeField = typename LocalBasisTraits::RangeFieldType;
  using Jacobian = typename LocalBasisTraits::JacobianType;

public:
  using Key = Impl::LocalBasisCahceKey<Domain>;

  /**
   * @brief      Constructs a new instance with empty values.
   */
  LocalBasisCache() {}

  /**
   * @brief Key to match previews applications of this cache
   * @details The key is composed with the std::type_index of the finite element,
   * the address of the quadrature rule, and the projected center of the quadrature
   * geometry type. This means that:
   * - Finite elements where shape functions change dynamically cannot be uniquely represented.
   * - Quadrature rules on non-conforming intersections cannot be uniquely represented.
   *
   * @tparam dim            Dimension of the quadrature rule to apply
   * @param finite_element  Finite element to apply cache to
   * @param quad_rule       Quadrature rule to evaluate the finite element
   * @param quad_proj       Quadrature positions projection into the finite element reference element
   * @return Key            Key with the parameter ids
   */
  template<int dim>
  Key key(const auto& finite_element,
          const Dune::QuadratureRule<DomainField, dim>& quad_rule,
          auto&& quad_proj) const noexcept
  {
    const auto& fem_id = typeid(finite_element.localBasis());
    const auto quad_id = &quad_rule;
    const auto proj_id =
      quad_proj(referenceElement<DomainField, dim>(quad_rule.type()).position(0, 0));
    return { fem_id, quad_id, proj_id };
  }

  /**
   * @brief Binds a finite element to this cache
   */
  void bind(const auto& finite_element,
            const Dune::QuadratureRule<DomainField, dimDomain>& quad_rule,
            bool force = false)
  {
    bind(finite_element, quad_rule, std::identity{}, force);
  }

  /**
   * @brief Binds a finite element to this cache
   * @details The values of `evaluateFunction` and `evaluateJacobian` of the
   * local basis in the finite element will be reused wrt a previous usage of
   * this cache if the key has already been used. See key function for details.
   * The force parameter will force recomputation of the values. This is useful
   * if the key cannot uniquely the local basis evaluation. In either case the values
   * resulting from  `evaluateFunction` and `evaluateJacobian` will be stored within
   * this cache. The quadrature projection is useful when the quadrature rule belongs
   * to a sub-entity or an intersection of the finite element reference element.
   *
   * @tparam dim            Dimension of the quadrature rule to apply
   * @param finite_element  Finite element to apply cache to
   * @param quad_rule       Quadrature rule to evaluate the finite element
   * @param quad_proj       Quadrature positions projection into the finite element reference element
   * @param force           Whether to force a recomputation of the results
   */
  template<int dim, typename ProjQuad>
    requires std::is_invocable_r_v<Domain, ProjQuad, FieldVector<DomainField, dim>>
  void bind(const auto& finite_element,
            const Dune::QuadratureRule<DomainField, dim>& quad_rule,
            ProjQuad&& quad_proj,
            bool force = false)
  {
    if (force) {
      _key = Key{};
    } else {
      auto new_key = key(finite_element, quad_rule, std::forward<ProjQuad>(quad_proj));
      if (std::exchange(_key, new_key) == new_key) {
        return;
      }
    }

    _fe_size = finite_element.size();
    _quad_size = quad_rule.size();

    auto rit = (_key == Key{}) ? _range.end() : _range.find(_key);
    if (rit != _range.end()) {
      _bound_range = rit->second.data();
    } else {
      _bound_range = cache_evaluate(finite_element.localBasis(), quad_rule, quad_proj, _key);
    }

    auto jit = (_key == Key{}) ? _jacobian.end() : _jacobian.find(_key);
    if (jit != _jacobian.end()) {
      _bound_jacobian = jit->second.data();
    } else {
      _bound_jacobian = cache_jacobian(finite_element.localBasis(), quad_rule, quad_proj, _key);
    }
  }

  /**
   * @brief      Unbind finite element from this cache
   */
  void unbind() noexcept
  {
    assert(std::exchange(_key, Key{}) != Key{});
    assert(std::exchange(_bound_range, nullptr));
    assert(std::exchange(_bound_jacobian, nullptr));
  }

  /**
   * @brief      Evaluate function with the bound finite element and quadrature rule
   *
   * @param[in]  position_id  The position id wrt quadrature rule
   *
   * @return     The evaluation for each of the local basis at position_id
   */
  std::span<const Range> evaluateFunction(std::size_t position_id) const
  {
    assert(_bound_range);
    assert(position_id < _quad_size);
    return { _bound_range + _fe_size * position_id, _fe_size };
  }

  /**
   * @brief      Evaluate jacobian with the bound finite element and quadrature rule
   *
   * @param[in]  position_id  The position id wrt quadrature rule
   *
   * @return     The jacobian for each of the local basis at position_id
   */
  std::span<const Jacobian> evaluateJacobian(std::size_t position_id) const
  {
    assert(_bound_jacobian);
    assert(position_id < _quad_size);
    return { _bound_jacobian + _fe_size * position_id, _fe_size };
  }

private:
  Range* cache_evaluate(const auto& basis, const auto& quad_rule, const auto& projection, Key key)
  {
    auto [it, inserted] = _range.emplace(key, std::vector<Range>{});
    assert(inserted or key == Key{});
    auto& range = it->second;
    range.clear();
    for (const auto& quad : quad_rule) {
      basis.evaluateFunction(projection(quad.position()), _quad_range);
      range.insert(end(range), begin(_quad_range), end(_quad_range));
    }
    assert(range.size() == _fe_size * _quad_size);
    return range.data();
  }

  Jacobian* cache_jacobian(const auto& basis,
                           const auto& quad_rule,
                           const auto& projection,
                           Key key)
  {
    auto [it, inserted] = _jacobian.emplace(key, std::vector<Jacobian>{});
    assert(inserted or key == Key{});
    auto& jacobian = it->second;
    jacobian.clear();
    for (const auto& quad : quad_rule) {
      basis.evaluateJacobian(projection(quad.position()), _quad_jacobian);
      jacobian.insert(end(jacobian), begin(_quad_jacobian), end(_quad_jacobian));
    }
    assert(jacobian.size() == _fe_size * _quad_size);
    return jacobian.data();
  }

  Key _key = {};
  std::size_t _fe_size = 0, _quad_size = 0;
  Range const* _bound_range = nullptr;
  Jacobian const* _bound_jacobian = nullptr;
  std::unordered_map<Key, std::vector<Range>> _range;
  std::unordered_map<Key, std::vector<Jacobian>> _jacobian;
  std::vector<Range> _quad_range;
  std::vector<Jacobian> _quad_jacobian;
};

} // namespace Dune::Copasi

#endif // DUNE_COPASI_LOCAL_FINITE_ELEMENT_CACHE_HH
