#ifndef DUNE_COPASI_VIRTUAL_LOCAL_FINITE_ELEMENT_HH
#define DUNE_COPASI_VIRTUAL_LOCAL_FINITE_ELEMENT_HH

#include <dune/copasi/common/factory.hh>

#include <dune/pdelab/finiteelementmap/finiteelementmap.hh>

#include <dune/localfunctions/common/virtualinterface.hh>
#include <dune/localfunctions/common/virtualwrappers.hh>

#include <array>
#include <memory>
#include <type_traits>

namespace Dune::Copasi {

// LBT Local basis traits
// Entity grid entity type
template<class LBT, class Entity>
struct VirtualLocalFiniteElementMapTraits
  : public PDELab::LocalFiniteElementMapTraits<LocalFiniteElementVirtualInterface<LBT>>
{
  using EntityType = Entity;
};

// LBT Local basis traits
// Entity grid entity type
template<class LBT, class Entity>
struct VirtualLocalFiniteElementMapInterface
  : public PDELab::LocalFiniteElementMapInterface<VirtualLocalFiniteElementMapTraits<LBT,Entity>,
              VirtualLocalFiniteElementMapInterface<LBT,Entity>>
{
  using Traits = VirtualLocalFiniteElementMapTraits<LBT,Entity>;

  virtual ~VirtualLocalFiniteElementMapInterface() {}

  virtual std::unique_ptr<VirtualLocalFiniteElementMapInterface<LBT,Entity>> clone() const = 0;

  virtual const typename Traits::FiniteElementType&
  find (const typename Traits::EntityType& e) const = 0;

  virtual bool fixedSize() const = 0;

  virtual std::size_t size(GeometryType gt) const = 0;

  virtual bool hasDOFs(int codim) const = 0;

  virtual std::size_t maxLocalSize() const = 0;
};

template<class FEM, class GV>
class VirtualLocalFiniteElementMapWrapper
  : public VirtualLocalFiniteElementMapInterface<
      typename FEM::Traits::FiniteElementType::Traits::LocalBasisType::Traits, typename GV::template Codim<0>::Entity>
{
  static_assert(not std::is_polymorphic_v<FEM>);

  using BaseFE = typename FEM::Traits::FiniteElementType;
  using VirtualFE = LocalFiniteElementVirtualImp<BaseFE>;
public:

  using Interface = VirtualLocalFiniteElementMapInterface<typename BaseFE::Traits::LocalBasisType::Traits, typename GV::template Codim<0>::Entity>;
  using Traits = typename Interface::Traits;
  static constexpr int dimension = FEM::dimension;

  VirtualLocalFiniteElementMapWrapper(std::unique_ptr<const FEM>&& fem)
    : _fem(std::move(fem))
  {}

  VirtualLocalFiniteElementMapWrapper(const FEM& fem)
    : VirtualLocalFiniteElementMapWrapper(std::make_unique<FEM>(fem))
  {}


  virtual ~VirtualLocalFiniteElementMapWrapper()
  {}

public:
  virtual std::unique_ptr<Interface> clone() const override
  {
    return std::make_unique<VirtualLocalFiniteElementMapWrapper<FEM,GV>>(*_fem);
  }

  virtual const typename Traits::FiniteElementType&
  find (const typename Traits::EntityType& e) const override
  {
    const auto& base_fe = _fem->find(e);

    // cache the last used base finite elements
    if (_virtual_fe.find(&base_fe) == _virtual_fe.end())
      _virtual_fe[&base_fe] = std::make_unique<const VirtualFE>(base_fe);
    return *_virtual_fe[&base_fe];
  }

  virtual bool fixedSize() const override
  {
    return _fem->fixedSize();
  };

  virtual std::size_t size(GeometryType gt) const override
  {
    return _fem->size(gt);
  }

  virtual bool hasDOFs(int codim) const override
  {
    return _fem->hasDOFs(codim);
  }

  virtual std::size_t maxLocalSize() const override
  {
    return _fem->maxLocalSize();
  }

protected:
  const std::unique_ptr<const FEM> _fem;
  mutable BaseFE const * _base_fe_cache;
    mutable std::map<BaseFE const *,std::unique_ptr<const VirtualFE>> _virtual_fe;
};

template<class FEM, class GV>
struct Factory<VirtualLocalFiniteElementMapWrapper<FEM,GV>>
{
  template<class Ctx>
  static auto create(Ctx&& ctx)
  {
    auto fem = Factory<FEM>::create(std::forward<Ctx>(ctx));
    return std::make_unique<VirtualLocalFiniteElementMapWrapper<FEM,GV>>(*fem);
  }
};

} // namespace Dune::Copasi

#endif // DUNE_COPASI_VIRTUAL_LOCAL_FINITE_ELEMENT_HH