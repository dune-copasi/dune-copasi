#ifndef DUNE_COPASI_DYNAMIC_POWER_LOCAL_FINITE_ELEMENT_MAP_HH
#define DUNE_COPASI_DYNAMIC_POWER_LOCAL_FINITE_ELEMENT_MAP_HH

#include <dune/copasi/common/factory.hh>
#include <dune/copasi/finite_element/dynamic_power.hh>

#include <dune/pdelab/finiteelementmap/finiteelementmap.hh>

namespace Dune::Copasi {

/**
 * @brief      This class describes a dynamic power local finite element map.
 * @details    This class wraps a usual PDELab finite element map into a
 *             dynamic power finite element map.
 * @ingroup    FiniteElementMap
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
   * @return     Always the underlying fixed size method
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
   * @brief      Describes whether a codim has associated degrees of freedom
   *
   * @param[in]  codim  The codim
   *
   * @return     True if codim has degrees of freedom
   */
  template<class = void>
  bool hasDOFs(int codim) const
  {
    return (_power_size != 0) && _fem->hasDOFs(codim);
  }

  /**
   * @brief      Maximum local size for all finite elements
   *
   * @return     The maximum size
   */
  std::size_t maxLocalSize() const { return _power_size * _fem->maxLocalSize(); }

private:
  std::size_t _power_size;
  std::unique_ptr<const FiniteElementMap> _fem;
  mutable std::map<BaseFiniteElement const *,std::unique_ptr<const FiniteElement>> _fe_cache;
};

/**
 * @brief      Factory for DynamicPowerLocalFiniteElementMap instances
 * @ingroup    Factory, FiniteElementMap
 * @tparam     <unnamed>  Template parameters of the DynamicPowerLocalFiniteElementMap
 */
template<class BaseLocalFiniteElementMap>
struct Factory<DynamicPowerLocalFiniteElementMap<BaseLocalFiniteElementMap>>
{
  /**
   * @brief      Create method
   * @todo       Add power size
   *
   * @param      ctx   @ref DataContext containing sufficient data to construct a BaseLocalFiniteElementMap
   *
   * @tparam     Ctx   Universal reference to the @ref DataContext
   *
   * @return     Instance of DynamicPowerLocalFiniteElementMap
   */
  template<class Ctx>
  static auto create(Ctx&& ctx)
  {
    using FEM = DynamicPowerLocalFiniteElementMap<BaseLocalFiniteElementMap>;
    return std::make_unique<FEM>(Factory<BaseLocalFiniteElementMap>::create(std::forward<Ctx>(ctx)));
  }
};

} // namespace Dune::Copasi

#endif // DUNE_COPASI_DYNAMIC_POWER_LOCAL_FINITE_ELEMENT_MAP_HH
