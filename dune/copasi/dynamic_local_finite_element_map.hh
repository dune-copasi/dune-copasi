#ifndef DUNE_COPASI_DYNAMIC_LOCAL_FINITE_ELEMENT_MAP_HH
#define DUNE_COPASI_DYNAMIC_LOCAL_FINITE_ELEMENT_MAP_HH

#include <dune/pdelab/finiteelementmap/finiteelementmap.hh>

#include <dune/copasi/dynamic_local_finite_element.hh>

namespace Dune::Copasi {

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

  DynamicPowerLocalFiniteElementMap(FiniteElementMap fem,
                                    std::size_t power_size)
    : _power_size(power_size)
    , _fem(fem)
  {}

  DynamicPowerLocalFiniteElementMap(std::size_t power_size)
    : DynamicPowerLocalFiniteElementMap(FiniteElementMap{}, power_size)
  {}

  ~DynamicPowerLocalFiniteElementMap()
  {
    delete _fe_cache;
  }

  template<class EntityType>
  const FiniteElement& find(const EntityType& e) const
  {
    auto base_fe = _fem.find(e);
    // cache the last used base finite element
    if (_base_fe_cache != &base_fe)
    {
      if (_fe_cache) delete _fe_cache;
      _fe_cache = new FiniteElement(base_fe,_power_size);
    }
    return *_fe_cache;
  }

  bool fixedSize() const { return _fem.fixedSize(); }

  std::size_t size(GeometryType gt) const
  {
    return _power_size * _fem.size(gt);
  }

  bool hasDOFs(int codim) const { return _fem.hasDOFs(codim); }

  std::size_t maxLocalSize() const { return _power_size * _fem.maxLocalSize(); }

private:
  std::size_t _power_size;
  FiniteElementMap _fem;
  mutable BaseFiniteElement* _base_fe_cache;
  mutable FiniteElement* _fe_cache;
};

} // namespace Dune::Copasi

#endif // DUNE_COPASI_DYNAMIC_LOCAL_FINITE_ELEMENT_MAP_HH
