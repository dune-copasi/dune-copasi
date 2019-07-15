#ifndef DUNE_COPASI_LOCAL_OPERATOR_MULTIDOMAIN_DIFFUSION_REACTION_HH
#define DUNE_COPASI_LOCAL_OPERATOR_MULTIDOMAIN_DIFFUSION_REACTION_HH

#include <dune/copasi/concepts.hh>
#include <dune/copasi/local_operator.hh>

#include <dune/pdelab/localoperator/numericaljacobian.hh>


#include <algorithm>

namespace Dune::Copasi {

template<class Grid, class LocalFiniteElement, int max_subdomains>
class LocalOperatorMultiDomainDiffusionReaction
  : public Dune::PDELab::FullSkeletonPattern
  , public Dune::PDELab::LocalOperatorDefaultFlags
  , public Dune::PDELab::InstationaryLocalOperatorDefaultMethods<double>
  , public Dune::PDELab::NumericalJacobianSkeleton<LocalOperatorMultiDomainDiffusionReaction<Grid,LocalFiniteElement,max_subdomains>>
{
  static_assert(Concept::isMultiDomainGrid<Grid>());

  //! local basis
  using LocalBasis = typename LocalFiniteElement::Traits::LocalBasisType;

  //! domain field
  using DF = typename LocalBasis::Traits::DomainFieldType;

  //! coordinates type
  using Domain = typename LocalBasis::Traits::DomainType;

  //! range field
  using RF = typename LocalBasis::Traits::RangeFieldType;

  //! range type (for the local finite element)
  using Range = typename LocalBasis::Traits::RangeType;

  //! jacobian tpye
  using Jacobian = typename LocalBasis::Traits::JacobianType;

  //! world dimension
  static constexpr int dim = LocalBasis::Traits::dimDomain;

  using GridView = typename Grid::SubDomainGrid::LeafGridView;
  using BaseLOP = LocalOperatorDiffusionReaction<GridView, LocalFiniteElement>;

  std::array<std::shared_ptr<BaseLOP>, max_subdomains> _local_operator;

  //! number of basis per component
  const std::size_t _basis_size;

public:
  //! pattern assembly flags
  static constexpr bool doPatternVolume = true;

  //! pattern assembly flags
  static constexpr bool doPatternSkeleton = true;

  //! visit skeleton from the two sides
  static constexpr bool doSkeletonTwoSided = true;

  //! residual assembly flags
  static constexpr bool doAlphaVolume = true;

  //! residual assembly flags
  static constexpr bool doAlphaSkeleton = true;

  LocalOperatorMultiDomainDiffusionReaction(
    std::shared_ptr<const Grid> grid,
    const ParameterTree& config,
    const LocalFiniteElement& finite_element)
    : _basis_size(finite_element.localBasis().size())
  {
    const auto& compartements = config.sub("compartements").getValueKeys();
    int size = compartements.size();
    // warning if size > max_subdomains
    for (int i = 0; i < std::min(size, max_subdomains); ++i) {
      const std::string compartement = compartements[i];
      auto& model_config = config.sub(compartement);

      int sub_domain_id =
        config.sub("compartements").template get<int>(compartement);
      GridView sub_grid_view = grid->subDomain(sub_domain_id).leafGridView();

      auto sub_config = config.sub(compartements[i]);
      auto lp =
        std::make_shared<BaseLOP>(sub_grid_view, sub_config, finite_element);
      _local_operator[i] = lp;
    }
  }

  // define sparsity pattern of operator representation
  template<typename LFSU, typename LFSV, typename LocalPattern>
  void pattern_volume(const LFSU& lfsu,
                      const LFSV& lfsv,
                      LocalPattern& pattern) const
  {
    for (int i = 0; i < max_subdomains; ++i)
      _local_operator[i]->pattern_volume(lfsu.child(i), lfsv.child(i), pattern);
  }

  template<typename LFSU, typename LFSV, typename LocalPattern>
  void pattern_skeleton(const LFSU& lfsu_i,
                        const LFSV& lfsv_i,
                        const LFSU& lfsu_o,
                        const LFSV& lfsv_o,
                        LocalPattern& pattern_io,
                        LocalPattern& pattern_oi) const
  {
    int domain_i(-99), domain_o(-99);
    for (int k = 0; k < max_subdomains; k++) {
      if (lfsu_i.child(k).size() > 0)
        domain_i = k;
      if (lfsu_o.child(k).size() > 0)
        domain_o = k;
    }

    if (domain_i == domain_o)
      return;
    if ((domain_i == -99) or (domain_o == -99))
      return;

    assert(lfsu_i.child(domain_i).size() == lfsv_i.child(domain_i).size());
    assert(lfsu_o.child(domain_o).size() == lfsv_o.child(domain_o).size());

    // TODO: Remove off diagonal on interdomain skeleton
    for (unsigned int i = 0; i < lfsv_i.child(domain_i).size(); ++i)
      for (unsigned int j = 0; j < lfsu_o.child(domain_o).size(); ++j)
        pattern_io.addLink(
          lfsv_i.child(domain_i), i, lfsu_o.child(domain_o), j);

    for (unsigned int i = 0; i < lfsv_o.child(domain_o).size(); ++i)
      for (unsigned int j = 0; j < lfsu_i.child(domain_i).size(); ++j)
        pattern_oi.addLink(
          lfsv_o.child(domain_o), i, lfsu_i.child(domain_i), j);
  }

  template<typename EG, typename LFSU, typename X, typename LFSV, typename R>
  void jacobian_apply_volume(const EG& eg,
                             const LFSU& lfsu,
                             const X& x,
                             const X& z,
                             const LFSV& lfsv,
                             R& r) const
  {
    for (int i = 0; i < max_subdomains; ++i) {
      if (lfsu.child(i).size() > 0) {
        auto sub_lfsu = lfsu.child(i);
        auto sub_lfsv = lfsv.child(i);
        if (sub_lfsu.size() == 0)
          continue;
        _local_operator[i]->jacobian_apply_volume(
          eg, sub_lfsu, x, z, sub_lfsv, r);
      }
    }
  }

  //! jacobian contribution of volume term
  template<typename EG, typename LFSU, typename X, typename LFSV, typename M>
  void jacobian_volume(const EG& eg,
                       const LFSU& lfsu,
                       const X& x,
                       const LFSV& lfsv,
                       M& mat) const
  {
    for (int i = 0; i < max_subdomains; ++i) {
      if (lfsu.child(i).size() > 0) {
        auto sub_lfsu = lfsu.child(i);
        auto sub_lfsv = lfsv.child(i);
        if (sub_lfsu.size() == 0)
          continue;
        _local_operator[i]->jacobian_volume(eg, sub_lfsu, x, sub_lfsv, mat);
      }
    }
  }

  template<typename EG, typename LFSU, typename X, typename LFSV, typename R>
  void alpha_volume(const EG& eg,
                    const LFSU& lfsu,
                    const X& x,
                    const LFSV& lfsv,
                    R& r) const
  {
    // auto subdomain = eg.subDomain();
    for (int i = 0; i < max_subdomains; ++i) {
      if (lfsu.child(i).size() > 0) {
        auto sub_lfsu = lfsu.child(i);
        auto sub_lfsv = lfsv.child(i);
        if (sub_lfsu.size() == 0)
          continue;
        _local_operator[i]->alpha_volume(eg, sub_lfsu, x, sub_lfsv, r);
      }
    }
  }

  //! apply local jacobian of the volume term -> linear variant
  template<typename EG, typename LFSU, typename X, typename LFSV, typename R>
  void jacobian_apply_volume(const EG& eg,
                             const LFSU& lfsu,
                             const X& x,
                             const LFSV& lfsv,
                             R& r) const
  {
    for (int i = 0; i < max_subdomains; ++i) {
      if (lfsu.child(i).size() > 0) {
        auto sub_lfsu = lfsu.child(i);
        auto sub_lfsv = lfsv.child(i);
        if (sub_lfsu.size() == 0)
          continue;
        _local_operator[i]->jacobian_apply_volume(eg, sub_lfsu, x, sub_lfsv, r);
      }
    }
  }

  // skeleton integral depending on test and ansatz functions
  // each face is only visited ONCE!
  template<typename IG, typename LFSU, typename X, typename LFSV, typename R>
  void alpha_skeleton(const IG& ig,
                      const LFSU& lfsu_i,
                      const X& x_i,
                      const LFSV& lfsv_i,
                      const LFSU& lfsu_o,
                      const X& x_o,
                      const LFSV& lfsv_o,
                      R& r_i,
                      R& r_o) const
  {
    const int dim = IG::Entity::dimension;
    const int order = 3;

    int domain_i(-99), domain_o(-99);
    for (int k = 0; k < max_subdomains; k++) {
      if (lfsu_i.child(k).size() > 0)
        domain_i = k;
      if (lfsu_o.child(k).size() > 0)
        domain_o = k;
    }

    if ((domain_i == domain_o))
      return;
    if ((domain_i == -99) or (domain_o == -99))
      return;

    assert(lfsu_i.child(domain_i).size() == lfsv_i.child(domain_i).size());
    assert(lfsu_o.child(domain_o).size() == lfsv_o.child(domain_o).size());

    const auto& entity_f = ig.intersection();
    const auto& entity_i = ig.inside();
    const auto& entity_o = ig.outside();

    auto geo_f = entity_f.geometry();
    auto geo_i = entity_i.geometry();
    auto geo_o = entity_o.geometry();
    auto geo_in_i = entity_f.geometryInInside();
    auto geo_in_o = entity_f.geometryInOutside();

    const auto& lfsu_basis_i = lfsu_i.child(domain_i).finiteElement().localBasis();
    const auto& lfsu_basis_o = lfsu_o.child(domain_o).finiteElement().localBasis();

    const auto& lfsv_basis_i = lfsv_i.child(domain_i).finiteElement().localBasis();
    const auto& lfsv_basis_o = lfsv_o.child(domain_o).finiteElement().localBasis();

    std::size_t components_i = lfsu_basis_i.size()/_basis_size;
    std::size_t components_o = lfsu_basis_o.size()/_basis_size;

    auto x_coeff_i = [&](const std::size_t& component, const std::size_t& dof) {
      return x_i(lfsu_i.child(domain_i), component * _basis_size + dof);
    };
    auto x_coeff_o = [&](const std::size_t& component, const std::size_t& dof) {
      return x_o(lfsu_o.child(domain_o), component * _basis_size + dof);
    };

    auto accumulate_i = [&](const std::size_t& component,
                          const std::size_t& dof,
                          const auto& value) {
      r_i.accumulate(lfsu_i.child(domain_i), component * _basis_size + dof, value);
    };

    auto accumulate_o = [&](const std::size_t& component,
                          const std::size_t& dof,
                          const auto& value) {
      r_o.accumulate(lfsu_o.child(domain_o), component * _basis_size + dof, value);
    };

    typename IG::Entity::Geometry::JacobianInverseTransposed jac;


    for (const auto& it : quadratureRule(geo_f,3))
    {
      const auto& position_f = it.position();
      // position of quadrature point in local coordinates of elements
      const auto position_i = geo_in_i.global(position_f);
      const auto position_o = geo_in_o.global(position_f);

      // face normal vector
      const auto normal_f = entity_f.unitOuterNormal(position_f);

      std::vector<Range> phiu_i(lfsu_i.size());
      std::vector<Range> phiu_o(lfsu_o.size());

      // evaluate basis functions
      lfsu_basis_i.evaluateFunction(position_i,phiu_i);
      lfsu_basis_o.evaluateFunction(position_o,phiu_o);

      // evaluate concentrations at quadrature point
      DynamicVector<RF> u_i(components_i), u_o(components_o);

      int count=0;
      for (std::size_t k = 0; k < components_i; k++) // loop over components
        for (std::size_t j = 0; j < _basis_size; j++, count++)
          u_i[k] += x_coeff_i(k, j) * phiu_i[count];
      
      count=0;
      for (std::size_t k = 0; k < components_o; k++) // loop over components
        for (std::size_t j = 0; j < _basis_size; j++, count++)
          u_o[k] += x_coeff_o(k, j) * phiu_o[count];

      // integration factor
      auto factor = it.weight() * geo_f.integrationElement(position_f);

      count=0;
      for (std::size_t k = 0; k < components_i; k++) // loop over components
        for (std::size_t j = 0; j < _basis_size; j++, count++)
            accumulate_i(k, j, factor*(u_i[k]-u_o[k]) * phiu_i[count]);

      count=0;
      for (std::size_t k = 0; k < components_o; k++) // loop over components
        for (std::size_t j = 0; j < _basis_size; j++, count++)
            accumulate_o(k, j, -factor*(u_i[k]-u_o[k]) * phiu_o[count]);
    }
  }
};

template<class Grid, class LocalFiniteElement, int max_subdomains>
class TemporalLocalOperatorMultiDomainDiffusionReaction
  : public Dune::PDELab::LocalOperatorDefaultFlags
  , public Dune::PDELab::InstationaryLocalOperatorDefaultMethods<double>
{
  static_assert(Concept::isMultiDomainGrid<Grid>());

  using GridView = typename Grid::SubDomainGrid::LeafGridView;
  using BaseLOP =
    TemporalLocalOperatorDiffusionReaction<GridView, LocalFiniteElement>;

  std::array<std::shared_ptr<BaseLOP>, max_subdomains> _local_operator;

public:
  //! pattern assembly flags
  static constexpr bool doPatternVolume = true;

  //! residual assembly flags
  static constexpr bool doAlphaVolume = true;

  TemporalLocalOperatorMultiDomainDiffusionReaction(
    std::shared_ptr<const Grid> grid,
    const ParameterTree& config,
    const LocalFiniteElement& finite_element)
  {
    const auto& compartements = config.sub("compartements").getValueKeys();
    int size = compartements.size();
    // warning if size > max_subdomains
    for (int i = 0; i < std::min(size, max_subdomains); ++i) {
      const std::string compartement = compartements[i];
      auto& model_config = config.sub(compartement);

      int sub_domain_id =
        config.sub("compartements").template get<int>(compartement);
      GridView sub_grid_view = grid->subDomain(sub_domain_id).leafGridView();

      auto sub_config = config.sub(compartements[i]);
      _local_operator[i] =
        std::make_shared<BaseLOP>(sub_grid_view, sub_config, finite_element);
    }
  }

  // define sparsity pattern of operator representation
  template<typename LFSU, typename LFSV, typename LocalPattern>
  void pattern_volume(const LFSU& lfsu,
                      const LFSV& lfsv,
                      LocalPattern& pattern) const
  {
    for (int i = 0; i < max_subdomains; ++i)
      _local_operator[i]->pattern_volume(lfsu.child(i), lfsv.child(i), pattern);
  }

  template<typename EG, typename LFSU, typename X, typename LFSV, typename R>
  void alpha_volume(const EG& eg,
                    const LFSU& lfsu,
                    const X& x,
                    const LFSV& lfsv,
                    R& r) const
  {
    for (int i = 0; i < max_subdomains; ++i) {
      if (lfsu.child(i).size() > 0) {
        auto sub_lfsu = lfsu.child(i);
        auto sub_lfsv = lfsv.child(i);
        if (sub_lfsu.size() == 0)
          continue;
        _local_operator[i]->alpha_volume(eg, sub_lfsu, x, sub_lfsv, r);
      }
    }
  }

  template<typename EG, typename LFSU, typename X, typename LFSV, typename Mat>
  void jacobian_volume(const EG& eg,
                       const LFSU& lfsu,
                       const X& x,
                       const LFSV& lfsv,
                       Mat& mat) const
  {
    for (int i = 0; i < max_subdomains; ++i) {
      if (lfsu.child(i).size() > 0) {
        auto sub_lfsu = lfsu.child(i);
        auto sub_lfsv = lfsv.child(i);
        if (sub_lfsu.size() == 0)
          continue;
        _local_operator[i]->jacobian_volume(eg, sub_lfsu, x, sub_lfsv, mat);
      }
    }
  }

  template<typename EG, typename LFSU, typename X, typename LFSV, typename R>
  void jacobian_apply_volume(const EG& eg,
                             const LFSU& lfsu,
                             const X& x,
                             const X& z,
                             const LFSV& lfsv,
                             R& r) const
  {
    for (int i = 0; i < max_subdomains; ++i) {
      if (lfsu.child(i).size() > 0) {
        auto sub_lfsu = lfsu.child(i);
        auto sub_lfsv = lfsv.child(i);
        if (sub_lfsu.size() == 0)
          continue;
        _local_operator[i]->jacobian_apply_volume(
          eg, sub_lfsu, x, z, sub_lfsv, r);
      }
    }
  }

  template<typename EG, typename LFSU, typename X, typename LFSV, typename R>
  void jacobian_apply_volume(const EG& eg,
                             const LFSU& lfsu,
                             const X& x,
                             const LFSV& lfsv,
                             R& r) const
  {
    for (int i = 0; i < max_subdomains; ++i) {
      if (lfsu.child(i).size() > 0) {
        auto sub_lfsu = lfsu.child(i);
        auto sub_lfsv = lfsv.child(i);
        if (sub_lfsu.size() == 0)
          continue;
        _local_operator[i]->jacobian_apply_volume(eg, sub_lfsu, x, sub_lfsv, r);
      }
    }
  }
};

} // namespace Dune::Copasi

#endif // DUNE_COPASI_LOCAL_OPERATOR_MULTIDOMAIN_DIFFUSION_REACTION_HH