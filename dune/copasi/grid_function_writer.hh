#ifndef DUNE_COPASI_GRID_FUNCTION_WRITER_HH
#define DUNE_COPASI_GRID_FUNCTION_WRITER_HH

#include <dune/pdelab/gridfunctionspace/gridfunctionspaceutilities.hh>

#include <dune/grid/io/file/vtk/subsamplingvtkwriter.hh>
#include <dune/grid/io/file/vtk/vtksequencewriter.hh>

#include <dune/common/exceptions.hh>
#include <dune/common/parametertree.hh>

#include <memory>
#include <string>
#include <type_traits>

namespace Dune::Copasi {

template<class GFS, typename X>
class DiscreteGridFunction
  : public PDELab::GridFunctionBase<
      PDELab::GridFunctionTraits<
        typename GFS::Traits::GridViewType,
        typename BasisInterfaceSwitch<typename FiniteElementInterfaceSwitch<
          typename GFS::Traits::FiniteElementType>::Basis>::RangeField,
        BasisInterfaceSwitch<typename FiniteElementInterfaceSwitch<
          typename GFS::Traits::FiniteElementType>::Basis>::dimRange,
        typename BasisInterfaceSwitch<typename FiniteElementInterfaceSwitch<
          typename GFS::Traits::FiniteElementType>::Basis>::Range>,
      DiscreteGridFunction<GFS, X>>
{
  typedef
    typename Dune::BasisInterfaceSwitch<typename FiniteElementInterfaceSwitch<
      typename GFS::Traits::FiniteElementType>::Basis>
      BasisSwitch;
  typedef PDELab::GridFunctionBase<
    PDELab::GridFunctionTraits<typename GFS::Traits::GridViewType,
                               typename BasisSwitch::RangeField,
                               BasisSwitch::dimRange,
                               typename BasisSwitch::Range>,
    DiscreteGridFunction<GFS, X>>
    BaseT;

public:
  typedef typename BaseT::Traits Traits;

  /** \brief Construct a DiscreteGridFunction
   *
   * \param gfs The GridFunctionsSpace
   * \param x_  The coefficients vector
   */
  DiscreteGridFunction(const GFS& gfs,
                       const X& x_,
                       std::size_t id_begin,
                       std::size_t id_end)
    : pgfs(stackobject_to_shared_ptr(gfs))
    , lfs(gfs)
    , lfs_cache(lfs)
    , x_view(x_)
    , xl(gfs.maxLocalSize())
    , yb(gfs.maxLocalSize())
    , _id_begin(id_begin)
    , _id_end(id_end)
  {}

  /** \brief Construct a DiscreteGridFunction
   *
   * \param gfs shared pointer to the GridFunctionsSpace
   * \param x_  shared pointer to the coefficients vector
   */
  DiscreteGridFunction(std::shared_ptr<const GFS> gfs,
                       std::shared_ptr<const X> x_,
                       std::size_t id_begin,
                       std::size_t id_end)
    : pgfs(gfs)
    , lfs(*gfs)
    , lfs_cache(lfs)
    , x_view(*x_)
    , xl(gfs->maxLocalSize())
    , yb(gfs->maxLocalSize())
    , px(x_) // FIXME: The LocalView should handle a shared_ptr correctly!
    , _id_begin(id_begin)
    , _id_end(id_end)
  {}

  // Evaluate
  inline void evaluate(const typename Traits::ElementType& e,
                       const typename Traits::DomainType& x,
                       typename Traits::RangeType& y) const
  {
    typedef FiniteElementInterfaceSwitch<
      typename Dune::PDELab::LocalFunctionSpace<GFS>::Traits::FiniteElementType>
      FESwitch;
    lfs.bind(e);
    lfs_cache.update();
    x_view.bind(lfs_cache);
    x_view.read(xl);
    x_view.unbind();
    FESwitch::basis(lfs.finiteElement()).evaluateFunction(x, yb);
    y = 0;
    for (unsigned int i = _id_begin; i < _id_end; i++)
      y.axpy(xl[i], yb[i]);
  }

  //! get a reference to the GridView
  inline const typename Traits::GridViewType& getGridView() const
  {
    return pgfs->gridView();
  }

private:
  typedef PDELab::LocalFunctionSpace<GFS> LFS;
  typedef PDELab::LFSIndexCache<LFS> LFSCache;
  typedef typename X::template ConstLocalView<LFSCache> XView;

  std::shared_ptr<GFS const> pgfs;
  mutable LFS lfs;
  mutable LFSCache lfs_cache;
  mutable XView x_view;
  mutable std::vector<typename Traits::RangeFieldType> xl;
  mutable std::vector<typename Traits::RangeType> yb;
  std::shared_ptr<const X>
    px; // FIXME: dummy pointer to make sure we take ownership of X
  const std::size_t _id_begin, _id_end;
};

/**
 * @brief      Converts PDELab grid function to match a dune-function.
 * @note       This adapter differs from the one provided in PDELab in that this
 *             one has the current interface of dune-grid to write vtk files.
 *             The main difference using this adapter is that it allows to write
 *             more than 3 components.
 *
 * @tparam     GF    PDELab grid function type.
 */
template<class GF>
class VTKGridFunctionAdapter
{
  using GridFunction = typename std::decay<GF>::type;
  using GridView = typename GF::Traits::GridViewType;
  using Entity = typename GridView::template Codim<0>::Entity;
  using Range = typename GF::Traits::RangeType;
  using Domain = typename GF::Traits::DomainType;

public:
  /**
   * @brief      Constructs the object
   *
   * @param[in]  gf    pointer to a constant grid function
   */
  VTKGridFunctionAdapter(std::shared_ptr<const GF> gf)
    : _gf(gf)
    , _entity(nullptr)
  {}

  /**
   * @brief      Binds an entity to the grid function
   *
   * @param[in]  e     Entity
   */
  void bind(const Entity& e) { _entity = &e; }

  /**
   * @brief      Unbinds the previously binded entity
   */
  void unbind() { _entity = nullptr; }

  /**
   * @brief      Evaluates the grid function at a given coordinate
   *
   * @param[in]  x     Local coordinate with respect to the last binded entity
   */
  Range operator()(const Domain& x) const
  {
    Range y;
    assert(_entity);
    _gf->evaluate(*_entity, x, y);
    return y;
  }

private:
  const Entity* _entity;
  std::shared_ptr<const GF> _gf;
};

/*-------------------------------------------------------------------------*/ /**
                                                                               * @brief      Class for grid fucntion vtk sequence writer. This class works
                                                                               *             exacly the same way as a VTKSequence writer but it receives
                                                                               *             directly smart pointers to grid functions
                                                                               *
                                                                               * @tparam     GV    GridView of the dune grid
                                                                               */
template<class GV>
class GridFunctionVTKSequenceWriter : public Dune::VTKSequenceWriter<GV>
{
public:
  // export constructors of vtk sequence writer
  using Dune::VTKSequenceWriter<GV>::VTKSequenceWriter;

  /*-----------------------------------------------------------------------*/ /**
                                                                               * @brief      Adds a vertex data.
                                                                               *
                                                                               * @param[in]  gf_ptr  A pointer to a grid function of the type GF
                                                                               * @param[in]  name    Name of the variable
                                                                               *
                                                                               * @tparam     GF      Type of the grid function
                                                                               */
  template<class GF>
  void addVertexData(std::shared_ptr<const GF> gf_ptr, std::string name)
  {
    static_assert(std::is_same<typename GF::Traits::GridViewType, GV>::value,
                  "GridFunctionVTKSequenceWriter only works with only one "
                  "GridView type.");

    using Range = typename GF::Traits::RangeType;
    using RF = typename GF::Traits::RangeFieldType;

    VTKGridFunctionAdapter<GF> vtk_gf(gf_ptr);

    const int vtk_ncomps = GF::Traits::dimRange;
    const int dim = GF::Traits::dimDomain;
    auto vtk_type = (vtk_ncomps == dim) ? VTK::FieldInfo::Type::vector
                                        : VTK::FieldInfo::Type::scalar;
    auto vtk_info = VTK::FieldInfo(name, vtk_type, vtk_ncomps);

    this->vtkWriter()->addVertexData(vtk_gf, vtk_info);
  }

  /*-----------------------------------------------------------------------*/ /**
                                                                               * @brief      Adds a cell data.
                                                                               *
                                                                               * @param[in]  gf_ptr  A pointer to a grid function of the type GF
                                                                               * @param[in]  name    Name of the variable
                                                                               *
                                                                               * @tparam     GF      Type of the grid function
                                                                               */
  template<class GF>
  void addCellData(std::shared_ptr<const GF> gf_ptr, std::string name)
  {
    static_assert(std::is_same<typename GF::Traits::GridViewType, GV>::value,
                  "GridFunctionVTKSequenceWriter only works with only one "
                  "GridView type.");

    VTKGridFunctionAdapter<GF> vtk_gf(gf_ptr);

    const int vtk_ncomps = GF::Traits::dimRange;
    const int dim = GF::Traits::dimDomain;
    auto vtk_type = (vtk_ncomps == dim) ? VTK::FieldInfo::Type::vector
                                        : VTK::FieldInfo::Type::scalar;
    auto vtk_info = VTK::FieldInfo(name, vtk_type, vtk_ncomps);

    this->vtkWriter()->addCellData(vtk_gf, vtk_info);
  }

  /*-----------------------------------------------------------------------*/ /**
                                                                               * @brief      Clear all the attached VTKFunctions.
                                                                               * @details    This is necessary for each iteration if VTKFunctions are
                                                                               *             managed internaly with pointers instead of references
                                                                               */
  void clear() { this->vtkWriter()->clear(); }
};

} // namespace Dune::Copasi

#endif // DUNE_COPASI_GRID_FUNCTION_WRITER_HH
