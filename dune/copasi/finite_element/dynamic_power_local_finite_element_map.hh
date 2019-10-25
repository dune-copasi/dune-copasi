#ifndef DUNE_COPASI_DYNAMIC_POWER_LOCAL_FINITE_ELEMENT_MAP_HH
#define DUNE_COPASI_DYNAMIC_POWER_LOCAL_FINITE_ELEMENT_MAP_HH

#include <dune/copasi/finite_element/dynamic_power_local_finite_element.hh>

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
  DynamicPowerLocalFiniteElementMap(FiniteElementMap fem,
                                   std::size_t power_size = 1)
    : _power_size(power_size)
    , _fem(fem)
    , _fe_cache(NULL)
  {}

  /**
   * @brief      Destroys the object.
   */
  ~DynamicPowerLocalFiniteElementMap() { delete _fe_cache; }

  /**
   * @brief      Searches for the finite element for entity e.
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
    auto base_fe = _fem.find(e);

    // cache the last used base finite element
    if (_base_fe_cache != &base_fe) {
      _base_fe_cache = &base_fe;
      if (_fe_cache != NULL)
        delete _fe_cache;
      _fe_cache = new FiniteElement(base_fe, _power_size);
    }
    return *_fe_cache;
  }

  /**
   * @brief      Returns true if this finite element map has a fixed size
   *
   * @return     Always the underlaying fixed size method
   */
  bool fixedSize() const { return _fem.fixedSize(); }

  /**
   * @brief      Size for a given geometry type
   *
   * @param[in]  gt    The geometry type
   *
   * @return     The size for the given geometry type
   */
  std::size_t size(GeometryType gt) const
  {
    return _power_size * _fem.size(gt);
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
    return (_power_size != 0) && _fem.hasDOFs(codim);
  }

  /**
   * @brief      Maximum local size for all finite elements
   *
   * @return     The maximim size
   */
  std::size_t maxLocalSize() const { return _power_size * _fem.maxLocalSize(); }

private:
  std::size_t _power_size;
  FiniteElementMap _fem;
  mutable BaseFiniteElement* _base_fe_cache;
  mutable FiniteElement* _fe_cache;
};

} // namespace Dune::Copasi

#endif // DUNE_COPASI_DYNAMIC_POWER_LOCAL_FINITE_ELEMENT_MAP_HH
