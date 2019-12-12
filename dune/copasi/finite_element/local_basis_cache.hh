#ifndef DUNE_COPASI_LOCAL_FINITE_ELEMENT_CACHE_HH
#define DUNE_COPASI_LOCAL_FINITE_ELEMENT_CACHE_HH

#include <functional>
#include <unordered_map>
#include <vector>
#include <typeindex>

#include <dune/common/hash.hh>

namespace Dune::Copasi {

template<class QuadratureRule>
class QuadratureRuleCacheKey
{
public:

  QuadratureRuleCacheKey(const QuadratureRule& rule, const std::type_index& basis_id)
    : _rule(rule)
    , _basis_id(basis_id)
  {}

  template<class LocalBasis>
  QuadratureRuleCacheKey(const QuadratureRule& rule, const LocalBasis& basis)
    : _rule(rule)
    , _basis_id(typeid(basis))
  {}

  std::size_t hash_code() const noexcept
  {
    std::size_t seed = 1'226'688'347; // some prime number
    Dune::hash_combine(seed,_rule);
    Dune::hash_combine(seed,_basis_id);
    return seed;
  }

  friend inline bool
  operator==(const QuadratureRuleCacheKey& lhs, const QuadratureRuleCacheKey& rhs)
  {
    return (lhs._basis_id == rhs._basis_id) and (lhs._rule == rhs._rule);
  }

  friend inline bool
  operator!=(const QuadratureRuleCacheKey& lhs, const QuadratureRuleCacheKey& rhs)
  {
    return not (lhs == rhs);
  }

private:
  const QuadratureRule _rule;
  const std::type_index _basis_id;
};

template<class Domain>
class QuadraturePointCacheKey
{
public:

  QuadraturePointCacheKey(const Domain& domain, const std::type_index& basis_id)
    : _domain(domain)
    , _basis_id(basis_id)
  {}

  template<class LocalBasis>
  QuadraturePointCacheKey(const Domain& domain, const LocalBasis& basis)
    : QuadraturePointCacheKey(domain,typeid(basis))
  {}

  std::size_t hash_code() const noexcept
  {
    std::size_t seed = 46'255'207; // some prime number
    Dune::hash_combine(seed,_basis_id);
    for (auto&& i : _domain)
      Dune::hash_combine(seed,i);
    return seed;
  }

  friend inline bool
  operator==(const QuadraturePointCacheKey& lhs, const QuadraturePointCacheKey& rhs)
  {
    return (lhs._basis_id == rhs._basis_id) and (lhs._domain == rhs._domain);
  }

  friend inline bool
  operator!=(const QuadraturePointCacheKey& lhs, const QuadraturePointCacheKey& rhs)
  {
    return not (lhs == rhs);
  }

private:
  const Domain _domain;
  const std::type_index _basis_id;
};

} // namespace Dune::Copasi

namespace std {

  template<typename Rule>
  struct hash<Dune::Copasi::QuadratureRuleCacheKey<Rule>>
  {
    std::size_t operator()(const Dune::Copasi::QuadratureRuleCacheKey<Rule>& key) const
    {
      return key.hash_code();
    }
  };

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

// @warning This caches assumes that the quadrature rule (from dune-geometry) is using the singleton
// pattern and its address is always the same, and that the local basis are the same iff its
// (maybe polymorphic) type is the same.
template<class LocalBasisTraits>
class LocalBasisCache
{
  using DomainField = typename LocalBasisTraits::DomainFieldType;
  static constexpr std::size_t dimDomain = LocalBasisTraits::dimDomain;
  using Domain = typename LocalBasisTraits::DomainType;

  using Rule = Dune::QuadratureRule<DomainField,dimDomain>;
  using RuleWrapper = Dune::PDELab::QuadratureRuleWrapper<Rule>;

  using BasisKey = std::type_index;
  using RuleKey = QuadratureRuleCacheKey<RuleWrapper>;
  using PointKey = QuadraturePointCacheKey<Domain>;

  using Range = typename LocalBasisTraits::RangeType;
  using RangeField = typename LocalBasisTraits::RangeFieldType;
  using Jacobian = typename LocalBasisTraits::JacobianType;

public:
  LocalBasisCache(bool do_function = true, bool do_jacobian = true)
    : _do_function(do_function)
    , _do_jacobian(do_jacobian)
    , _bound_basis(typeid(void))
  {}

  template<class Basis>
  void bind(const Rule& rule, const Basis& basis)
  {
    bind(RuleWrapper(rule),basis);
  }

  template<class Basis>
  void bind(const RuleWrapper& rule, const Basis& basis)
  {
    _bound_basis = typeid(basis);
    RuleKey rule_key(rule,_bound_basis);

    // is this pair (rule,basis) already cached?
    if (_cached_rules.find(rule_key) != _cached_rules.end())
      return;

    for(auto&& point : rule)
    {
      const auto& position = point.position();
      if (_do_function)
        cache_evaluate(position,basis,_bound_basis);
      if (_do_jacobian)
        cache_jacobian(position,basis,_bound_basis);
    }
    _cached_rules.insert(std::move(rule_key));
  }

  void unbind()
  {
    _bound_basis = typeid(void);
  }

  inline const std::vector<Range>& evaluateFunction(const Domain& domain) const
  {
    assert(_bound_basis != typeid(void));
    return evaluateFunction(PointKey{domain,_bound_basis});
  }

  inline const std::vector<Jacobian>& evaluateJacobian(const Domain& domain) const
  {
    assert(_bound_basis != typeid(void));
    return evaluateJacobian(PointKey{domain,_bound_basis});
  }

private:

  template<class Basis>
  void cache_evaluate(const Domain& position, const Basis& basis, const std::type_index& basis_id)
  {
    std::vector<Range> range;
    basis.evaluateFunction(position,range);
    auto&& range_pair = std::make_pair(PointKey{position,basis_id},std::move(range));
    auto insert_pair = _range.insert(range_pair);
    assert(insert_pair.second);
  }

  template<class Basis>
  void cache_jacobian(const Domain& position, const Basis& basis, const std::type_index& basis_id)
  {
    std::vector<Jacobian> jacobian;
    basis.evaluateJacobian(position,jacobian);
    auto&& jacobian_pair = std::make_pair(PointKey{position,basis_id},std::move(jacobian));
    auto insert_pair = _jacobian.insert(jacobian_pair);
    assert(insert_pair.second);
  }

  const std::vector<Range>& evaluateFunction(const PointKey& point_key) const
  {
    assert(_do_function);
    auto it = _range.find(point_key);
    assert(it != _range.end());
    return it->second;
  }

  const std::vector<Jacobian>& evaluateJacobian(const PointKey& point_key) const
  {
    assert(_do_jacobian);
    auto it = _jacobian.find(point_key);
    assert(it != _jacobian.end());
    return it->second;
  }

  const bool _do_function;
  const bool _do_jacobian;
  std::type_index _bound_basis;
  std::unordered_set<RuleKey> _cached_rules;
  std::unordered_map<PointKey,std::vector<Range>> _range;
  std::unordered_map<PointKey,std::vector<Jacobian>> _jacobian;
};

} // namespace Dune::Copasi

#endif // DUNE_COPASI_LOCAL_FINITE_ELEMENT_CACHE_HH