#ifndef DUNE_COPASI_DYNAMIC_LOCAL_FINITE_ELEMENT_MAP_HH
#define DUNE_COPASI_DYNAMIC_LOCAL_FINITE_ELEMENT_MAP_HH

#include <dune/pdelab/finiteelementmap/finiteelementmap.hh>

#include <dune/copasi/dynamic_local_finite_element.hh>

namespace Dune::Copasi {

template<class FiniteElementMap>
class DynamicPowerLocalFiniteElementMap 
  : public PDELab::LocalFiniteElementMapInterface<
          PDELab::LocalFiniteElementMapTraits<
            DynamicPowerLocalFiniteElement<typename FiniteElementMap::Traits::FiniteElement>
          >, 
          DynamicPowerLocalFiniteElementMap<FiniteElementMap>>
{
  // static_assert(std::is_base_of_v<FiniteElementMap,PDELab::SimpleLocalFiniteElementMap<FiniteElementMap>>);

  using BaseFiniteElement = typename FiniteElementMap::Traits::FiniteElement;
  using FiniteElement = DynamicPowerLocalFiniteElement<BaseFiniteElement>;

public:

  static const int dimension = FiniteElementMap::dimension;

  DynamicPowerLocalFiniteElementMap(FiniteElementMap fem, std::size_t power_size)
    : _power_size(power_size)
    , _fem(fem)
    , _fe(_fem.find(0),power_size)
  {}

  DynamicPowerLocalFiniteElementMap(std::size_t power_size)
    : DynamicPowerLocalFiniteElementMap(FiniteElementMap{},power_size)
  {}

  template<class EntityType>
  const FiniteElement& find (const EntityType& e) const
  {
    return _fe;
  }

  bool fixedSize() const 
  {
    return _fem.fixedSize();
  }

  std::size_t size(GeometryType gt) const
  {
    return _power_size * _fem.size(gt);
  }

  bool hasDOFs(int codim) const
  {
    return _fem.hasDOFs(codim);
  }

  std::size_t maxLocalSize() const
  {
    return _power_size * _fem.maxLocalSize();
  }

private:
  std::size_t       _power_size;
  FiniteElementMap  _fem;
  FiniteElement     _fe;
};

} // Dune::Copasi namespace

#endif // DUNE_COPASI_DYNAMIC_LOCAL_FINITE_ELEMENT_MAP_HH
