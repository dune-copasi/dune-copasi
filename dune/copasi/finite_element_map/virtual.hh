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

/**
 * @brief      Virtual local finite element traits
 * @details    Injects the entity type to the LocalFiniteElementMapTraits
 * @ingroup    FiniteElementMap
 *
 * @tparam     LBT     Local Basis Traits
 * @tparam     Entity  Entity Type
 */
template<class LBT, class Entity>
struct VirtualLocalFiniteElementMapTraits
  : public PDELab::LocalFiniteElementMapTraits<LocalFiniteElementVirtualInterface<LBT>>
{
  using EntityType = Entity;
};

/**
 * @brief      Virtual local finite element interface
 * @details    Defines a base class for virtual finite element maps
 *
 * @tparam     LBT     Local Basis Traits
 * @tparam     Entity  Entity Type
 */
template<class LBT, class Entity>
struct VirtualLocalFiniteElementMapInterface
  : public PDELab::LocalFiniteElementMapInterface<VirtualLocalFiniteElementMapTraits<LBT,Entity>,
              VirtualLocalFiniteElementMapInterface<LBT,Entity>>
{
  //! Traits for the local finite element
  using Traits = VirtualLocalFiniteElementMapTraits<LBT,Entity>;

  //! Destructor
  virtual ~VirtualLocalFiniteElementMapInterface() = default;

  /**
   * @brief      Creates a new instance of the object with same properties than original.
   *
   * @return     Copy of this object.
   */
  virtual std::unique_ptr<VirtualLocalFiniteElementMapInterface<LBT,Entity>> clone() const = 0;

  /**
   * @brief      Searches for the first match.
   *
   * @param[in]  e     Entity
   *
   * @return     A virtual finite element
   */
  virtual const typename Traits::FiniteElementType&
  find (const typename Traits::EntityType& e) const = 0;

  /**
   * @brief      Is fixed size
   *
   * @return     true if size does not depend on the entity
   */
  virtual bool fixedSize() const = 0;

  /**
   * @brief      Size for a give geometry type
   *
   * @param[in]  gt    The geometry type
   *
   * @return     The size for geometry type gt
   */
  virtual std::size_t size(GeometryType gt) const = 0;

  /**
   * @brief      Determines if codim has degrees of freedom
   *
   * @param[in]  codim  The codim
   *
   * @return     true if codim has degrees of freedom
   */
  virtual bool hasDOFs(int codim) const = 0;

  /**
   * @brief      Max local size for this finite element
   *
   * @return     Max local size for this finite element
   */
  virtual std::size_t maxLocalSize() const = 0;
};

/**
 * @brief      This class describes a virtual local finite element map wrapper.
 *
 * @tparam     FEM     The wrapped finite element map
 * @tparam     Entity  Entity type
 */
template<class FEM, class Entity>
class VirtualLocalFiniteElementMapWrapper
  : public VirtualLocalFiniteElementMapInterface<
      typename FEM::Traits::FiniteElementType::Traits::LocalBasisType::Traits, Entity>
{
  // ensure to do not wrap polymorphic types
  static_assert(not std::is_polymorphic_v<FEM>);

  using BaseFE = typename FEM::Traits::FiniteElementType;
  using VirtualFE = LocalFiniteElementVirtualImp<BaseFE>;
public:

  // Export interface, traits, and dimension
  using Interface = VirtualLocalFiniteElementMapInterface<typename BaseFE::Traits::LocalBasisType::Traits, Entity>;
  using Traits = typename Interface::Traits;
  static constexpr int dimension = FEM::dimension;

  /**
   * @brief      Constructs a new instance.
   *
   * @param      fem   A unique pointer to the wrapped finite element
   */
  VirtualLocalFiniteElementMapWrapper(std::unique_ptr<const FEM>&& fem)
    : _fem(std::move(fem))
  {}

  /**
   * @brief      Constructs a new instance.
   *
   * @param[in]  fem   The fem
   */
  VirtualLocalFiniteElementMapWrapper(const FEM& fem)
    : VirtualLocalFiniteElementMapWrapper(std::make_unique<FEM>(fem))
  {}

  //! Default destructor
  ~VirtualLocalFiniteElementMapWrapper() override = default;

  /**
   * @brief      Creates a new instance of the object with same properties than original.
   *
   * @return     Copy of this object.
   */
  std::unique_ptr<Interface> clone() const override
  {
    return std::make_unique<VirtualLocalFiniteElementMapWrapper<FEM,Entity>>(*_fem);
  }

    /**
   * @brief      Searches for the first match.
   *
   * @param[in]  e     Entity
   *
   * @return     A virtual finite element
   */
  const typename Traits::FiniteElementType&
  find (const typename Traits::EntityType& e) const override
  {
    const auto& base_fe = _fem->find(e);

    // cache the last used base finite element with the address of the wrapped result of find
    // this assumes that if the address of the base fem is unchanged, then the resulting finite element is also unchanged
    if (_virtual_fe.find(&base_fe) == _virtual_fe.end())
      _virtual_fe[&base_fe] = std::make_unique<const VirtualFE>(base_fe);

    // return value associated to the base finite element address
    return *_virtual_fe[&base_fe];
  }

  /**
   * @brief      Is fixed size
   *
   * @return     true if size does not depend on the entity
   */
  bool fixedSize() const override
  {
    return _fem->fixedSize();
  };

  /**
   * @brief      Size for a give geometry type
   *
   * @param[in]  gt    The geometry type
   *
   * @return     The size for geometry type gt
   */
  std::size_t size(GeometryType gt) const override
  {
    return _fem->size(gt);
  }

  /**
   * @brief      Determines if codim has degrees of freedom
   *
   * @param[in]  codim  The codim
   *
   * @return     true if codim has degrees of freedom
   */
  bool hasDOFs(int codim) const override
  {
    return _fem->hasDOFs(codim);
  }


  /**
   * @brief      Max local size for this finite element
   *
   * @return     Max local size for this finite element
   */
  std::size_t maxLocalSize() const override
  {
    return _fem->maxLocalSize();
  }

protected:
  const std::unique_ptr<const FEM> _fem;
  mutable std::map<BaseFE const *,std::unique_ptr<const VirtualFE>> _virtual_fe;
};

/**
 * @brief      Factory for VirtualLocalFiniteElementMapWrapper instances
 * @ingroup    Factory, FiniteElementMap
 * @tparam     <unnamed>  Template paramenters of the VirtualLocalFiniteElementMapWrapper
 */
template<class FEM, class Entity>
struct Factory<VirtualLocalFiniteElementMapWrapper<FEM,Entity>>
{
  /**
   * @brief      Create method
   *
   * @param      ctx   @ref DataContext containing sufficient data to 
   *                   create base finite element map of the type FEM from
   *                   another factory.
   *
   * @tparam     Ctx   Universal reference to the @ref DataContext
   *
   * @return     Instance of VirtualLocalFiniteElementMapWrapper
   */
  template<class Ctx>
  static auto create(Ctx&& ctx)
  {
    auto fem = Factory<FEM>::create(std::forward<Ctx>(ctx));
    return std::make_unique<VirtualLocalFiniteElementMapWrapper<FEM,Entity>>(*fem);
  }
};

} // namespace Dune::Copasi

#endif // DUNE_COPASI_VIRTUAL_LOCAL_FINITE_ELEMENT_HH