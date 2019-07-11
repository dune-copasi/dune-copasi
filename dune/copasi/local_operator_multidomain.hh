#ifndef DUNE_COPASI_LOCAL_OPERATOR_MULTIDOMAIN_DIFFUSION_REACTION_HH
#define DUNE_COPASI_LOCAL_OPERATOR_MULTIDOMAIN_DIFFUSION_REACTION_HH

#include <dune/copasi/concepts.hh>
#include <dune/copasi/local_operator.hh>

#include <algorithm>

namespace Dune::Copasi {

template<class Grid, class LocalFiniteElement, int max_subdomains>
class LocalOperatorMultiDomainDiffusionReaction
  : public Dune::PDELab::FullVolumePattern
  , public Dune::PDELab::LocalOperatorDefaultFlags
  , public Dune::PDELab::InstationaryLocalOperatorDefaultMethods<double>
{
  static_assert(Concept::isMultiDomainGrid<Grid>());

  using GridView = typename Grid::SubDomainGrid::LeafGridView;
  using BaseLOP = LocalOperatorDiffusionReaction<GridView, LocalFiniteElement>;

  std::array<std::shared_ptr<BaseLOP>, max_subdomains> _local_operator;

public:
  //! pattern assembly flags
  static constexpr bool doPatternVolume = true;

  //! residual assembly flags
  static constexpr bool doAlphaVolume = true;

  LocalOperatorMultiDomainDiffusionReaction(
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
      auto lp =
        std::make_shared<BaseLOP>(sub_grid_view, sub_config, finite_element);
      _local_operator[i] = lp;
    }
  }

  // define sparsity pattern of operator representation
  // template<typename LFSU, typename LFSV, typename LocalPattern>
  // void pattern_volume (const LFSU& lfsu, const LFSV& lfsv,
  //                      LocalPattern& pattern) const
  // {
  //   for (int i = 0; i < max_subdomains; ++i)
  //     for (size_t j=0; j<lfsv.child(i).size(); ++j)
  //       for (size_t k=0; k<lfsu.child(i).size(); ++k)
  //         pattern.addLink(lfsv.child(i),j,lfsu.child(i),k);
  // }

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
        if (sub_lfsu.size() == 0) continue;
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
        if (sub_lfsu.size() == 0) continue;
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
        if (sub_lfsu.size() == 0) continue;
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
        if (sub_lfsu.size() == 0) continue;
        _local_operator[i]->jacobian_apply_volume(eg, sub_lfsu, x, sub_lfsv, r);
      }
    }
  }
};

template<class Grid, class LocalFiniteElement, int max_subdomains>
class TemporalLocalOperatorMultiDomainDiffusionReaction
  : public Dune::PDELab::FullVolumePattern
  , public Dune::PDELab::LocalOperatorDefaultFlags
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
  // template<typename LFSU, typename LFSV, typename LocalPattern>
  // void pattern_volume (const LFSU& lfsu, const LFSV& lfsv,
  //                      LocalPattern& pattern) const
  // {
  //   for (int i = 0; i < max_subdomains; ++i)
  //     for (size_t j=0; j<lfsv.child(i).size(); ++j)
  //       for (size_t k=0; k<lfsu.child(i).size(); ++k)
  //         pattern.addLink(lfsv.child(i),j,lfsu.child(i),k);
  // }

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
        if (sub_lfsu.size() == 0) continue;
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
        if (sub_lfsu.size() == 0) continue;
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
        if (sub_lfsu.size() == 0) continue;
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
        if (sub_lfsu.size() == 0) continue;
        _local_operator[i]->jacobian_apply_volume(eg, sub_lfsu, x, sub_lfsv, r);
      }
    }
  }
};

} // namespace Dune::Copasi

#endif // DUNE_COPASI_LOCAL_OPERATOR_MULTIDOMAIN_DIFFUSION_REACTION_HH