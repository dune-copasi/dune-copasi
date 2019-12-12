#ifndef DUNE_COPASI_DYNAMIC_POWER_LOCAL_FINITE_ELEMENT_MAP_HH
#define DUNE_COPASI_DYNAMIC_POWER_LOCAL_FINITE_ELEMENT_MAP_HH

#include <dune/copasi/common/factory.hh>
#include <dune/copasi/finite_element/dynamic_power.hh>

#include <dune/pdelab/finiteelementmap/finiteelementmap.hh>

namespace Dune::Copasi {

/**
 * @brief      This class describes a dynamic power local finite element map.
 * @details    This class wrapps a usual PDELab finite element map into a
 *             dynamic power finite element map.
 *
 * @tparam     FiniteElementMap  The original finite element map to wrap
 */
template<class FiniteElementMap>
class DynamicPowerLocalFiniteElementMap
  : public PDELab::LocalFiniteElementMapInterface<
      PDELab::LocalFiniteElementMapTraits<DynamicPowerLocalFiniteElement<
        typename FiniteElementMap::Traits::FiniteElement>>,
      DynamicPowerLocalFiniteElementMap<FiniteElementMap>>
{
  using BaseFiniteElement = typename FiniteElementMap::Traits::FiniteElement;
  using FiniteElement = DynamicPowerLocalFiniteElement<BaseFiniteElement>;

public:
  static const int dimension = FiniteElementMap::dimension;

  DynamicPowerLocalFiniteElementMap(std::unique_ptr<FiniteElementMap>&& fem,
                                    std::size_t power_size = 1)
    : _power_size(power_size)
    , _fem(fem)
  {}

  /**
   * @brief      Constructs a new instance.
   *
   * @param[in]  fem                  The finite element map to wrap
   * @param[in]  base_finite_element  The base finite element of the finite
   *                                  element map
   * @param[in]  power_size           The power size for the sub domain part
   *                                  (this is always 0 outside of the
   *                                  subdomain)
   */
  DynamicPowerLocalFiniteElementMap(const FiniteElementMap& fem,
                                    std::size_t power_size = 1)
    : _power_size(power_size)
  {
    if constexpr (std::is_polymorphic_v<FiniteElementMap>)
      _fem = std::make_unique<FiniteElementMap>(fem.clone());
    else
      _fem = std::make_unique<FiniteElementMap>(fem);
  }

  DynamicPowerLocalFiniteElementMap(const DynamicPowerLocalFiniteElementMap& other)
    : _power_size(other._power_size)
  {
    if constexpr (std::is_polymorphic_v<FiniteElementMap>)
      _fem = std::make_unique<FiniteElementMap>(other._fem->clone());
    else
      _fem = std::make_unique<FiniteElementMap>(*other._fem);
  }

  /**
   * @brief      Searches for the finite element for entity e.
   * @warning    The return value is valid until the next entity
   *             tries to find another finite element. Hence, not
   *             suitable for concurrency
   * @todo       Cache more than one finite element
   *
   * @param[in]  e           The entity
   *
   * @tparam     EntityType  The entity
   *
   * @return     A dynamic power local finite element
   *             @DynamicPowerLocalFiniteElement
   */
  template<class EntityType>
  const FiniteElement& find(const EntityType& e) const
  {
    // get base finite element
    const auto& base_fe = _fem->find(e);

    // wrap base finite element into a dynamic power finite element and cache it
    if (_fe_cache.find(&base_fe) == _fe_cache.end())
      _fe_cache[&base_fe] = std::make_unique<const FiniteElement>(base_fe, _power_size);

    return *_fe_cache[&base_fe];
  }

  /**
   * @brief      Returns true if this finite element map has a fixed size
   *
   * @return     Always the underlaying fixed size method
   */
  bool fixedSize() const { return _fem->fixedSize(); }

  /**
   * @brief      Size for a given geometry type
   *
   * @param[in]  gt    The geometry type
   *
   * @return     The size for the given geometry type
   */
  std::size_t size(GeometryType gt) const
  {
    return _power_size * _fem->size(gt);
  }

  /**
   * @brief      Describes wheter a codim has associated degrees of freedom
   *
   * @param[in]  codim  The codim
   *
   * @return     True if codim has dregees of freedom
   */
  bool hasDOFs(int codim) const
  {
    return (_power_size != 0) && _fem->hasDOFs(codim);
  }

  /**
   * @brief      Maximum local size for all finite elements
   *
   * @return     The maximim size
   */
  std::size_t maxLocalSize() const { return _power_size * _fem->maxLocalSize(); }

private:
  std::size_t _power_size;
  std::unique_ptr<const FiniteElementMap> _fem;
  mutable BaseFiniteElement const * _base_fe_cache;
  mutable std::map<BaseFiniteElement const *,std::unique_ptr<const FiniteElement>> _fe_cache;
};

template<class BaseLocalFiniteElementMap>
struct Factory<DynamicPowerLocalFiniteElementMap<BaseLocalFiniteElementMap>>
{
  template<class Context>
  static auto create(const Context& ctx)
  {
    using FEM = DynamicPowerLocalFiniteElementMap<BaseLocalFiniteElementMap>;
    // todo add power size
    return std::make_unique<FEM>(Factory<BaseLocalFiniteElementMap>::create(ctx));
  }
};

} // namespace Dune::Copasi

#endif // DUNE_COPASI_DYNAMIC_POWER_LOCAL_FINITE_ELEMENT_MAP_HH