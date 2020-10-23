#ifndef DUNE_COPASI_LOCAL_OPERATOR_MULTIDOMAIN_DIFFUSION_REACTION_HH
#define DUNE_COPASI_LOCAL_OPERATOR_MULTIDOMAIN_DIFFUSION_REACTION_HH

#include <dune/copasi/common/enum.hh>
#include <dune/copasi/concepts/grid.hh>

#include <dune/pdelab/localoperator/flags.hh>
#include <dune/pdelab/localoperator/idefault.hh>
#include <dune/pdelab/localoperator/pattern.hh>
#include <dune/pdelab/localoperator/numericaljacobian.hh>
#include <dune/pdelab/localoperator/numericaljacobianapply.hh>
#include <dune/pdelab/localoperator/numericalnonlinearjacobianapply.hh>

#include <dune/logging.hh>

#include <dune/common/parametertree.hh>

#include <muParser.h>

#include <algorithm>
#include <map>
#include <vector>

namespace Dune::Copasi {

/**
 * @brief      This class describes a PDELab local operator for multi domain
 *             diffusion reaction systems.
 * @details    This class describre the operatrions for local integrals required
 *             for diffusion reaction system. The operator is only valid for
 *             entities contained in the grid. The local finite element is used
 *             for caching shape function evaluations. And the jacobian method
 *             switches between numerical and analytical jacobians. This local
 *             operator creates internally an individual local operator for
 *             every subdomain in the grid
 *
 * @tparam     Grid                The grid
 * @tparam     LocalFiniteElement  Local Finite Element
 * @tparam     JM                  Jacobian Method
 */
template<class Grid,
         class SubLocalOperator,
         JacobianMethod JM = JacobianMethod::Analytical>
class LocalOperatorMultiDomainDiffusionReaction
  : public Dune::PDELab::LocalOperatorDefaultFlags
  , public Dune::PDELab::InstationaryLocalOperatorDefaultMethods<double>
  , public Dune::PDELab::NumericalJacobianSkeleton<
      LocalOperatorMultiDomainDiffusionReaction<Grid,
                                                SubLocalOperator,
                                                JM>>
  , public Dune::PDELab::NumericalJacobianApplySkeleton<
      LocalOperatorMultiDomainDiffusionReaction<Grid,
                                                SubLocalOperator,
                                                JM>>
{
  static_assert(Concept::isMultiDomainGrid<Grid>());

  //! jacobian tpye
  using IndexSet = typename Grid::LeafGridView::IndexSet;

  static constexpr std::size_t unused_domain = ~std::size_t(0);

  using GridView = typename Grid::SubDomainGrid::LeafGridView;
  using SubLOP = SubLocalOperator;

  Logging::Logger _logger;

  const IndexSet& _index_set;

  // interior domain local operators
  std::vector<std::shared_ptr<SubLOP>> _lops;

  // map (domain_i,domain_o,component_i) -> component_o
  std::map<std::array<std::size_t, 3>, std::size_t> _component_offset;

  /// map: (domain_i,domain_o) -> vector of outflow_i persers
  std::map<std::array<std::size_t,2>,std::vector<mu::Parser>> _outflow_parser;
  std::map<std::array<std::size_t,2>,std::vector<mu::Parser>> _outflow_jac_parser;
  std::map<std::array<std::size_t,2>,std::vector<mu::Parser>> _outflow_cross_jac_parser;

  // value of trial spaces (for parser evaluation)
  mutable std::vector<std::vector<double>> _components;
  mutable std::vector<std::vector<double>> _components_dn;

public:

  static_assert(not SubLOP::doSkeletonTwoSided);
  static_assert(not SubLOP::doSkipEntity);

  //! selective assembly flags
  static constexpr bool doSkipEntity = false;
  static constexpr bool doSkipIntersection = true;

  //! pattern assembly flags
  static constexpr bool doPatternVolume = SubLOP::doPatternVolume;
  static constexpr bool doPatternSkeleton = true;

  //! residual assembly flags
  static constexpr bool doAlphaVolume = SubLOP::doAlphaVolume;
  static constexpr bool doAlphaSkeleton = true;

  /**
   * @brief      Output information on the parser error and throw DUNE exception
   *
   * @param[in]  e     Exception thrown by the parser
   */
  void handle_parser_error(const mu::Parser::exception_type& e) const
  {
    _logger.error("Evaluating muParser expression failed:"_fmt);
    _logger.error("  Parsed expression:   {}"_fmt, e.GetExpr());
    _logger.error("  Token:               {}"_fmt, e.GetToken());
    _logger.error("  Error position:      {}"_fmt, e.GetPos());
    _logger.error("  Error code:          {}"_fmt, int(e.GetCode()));
    _logger.error("  Error message:       {}"_fmt, e.GetMsg());
    DUNE_THROW(IOError, "Error evaluating muParser expression");
  }

  /**
   * @brief      Constructs a new instance.
   *
   * @param[in]  grid            The grid
   * @param[in]  config          The configuration
   * @param[in]  finite_element  The local finite element
   */
  LocalOperatorMultiDomainDiffusionReaction(
    std::shared_ptr<const Grid> grid,
    const ParameterTree& config)
    : _logger(Logging::Logging::componentLogger({}, "model"))
    , _index_set(grid->leafGridView().indexSet())
  {
    auto& compartment_name = config.sub("compartments",true).getValueKeys();
    _lops.resize(compartment_name.size());
    std::vector<std::vector<std::string>> component_name{_lops.size()};

    // create local operators for each compartment
    for (std::size_t i = 0; i < _lops.size(); ++i) {
      int sub_domain_id =
        config.sub("compartments",true).template get<int>(compartment_name[i]);
      GridView sub_grid_view = grid->subDomain(sub_domain_id).leafGridView();

      const auto& sub_config = config.sub(compartment_name[i],true);
      component_name[i] = sub_config.sub("reaction",true).getValueKeys();
      std::sort(component_name[i].begin(), component_name[i].end());

      auto lp = std::make_shared<SubLOP>(sub_grid_view, sub_config);
      _lops[i] = lp;
    }

    // create mapping between all inside and outside components
    _components.resize(_lops.size());
    _components_dn.resize(_lops.size());
    for (std::size_t domain_i = 0; domain_i < _lops.size(); ++domain_i) {
      _components[domain_i].assign(component_name[domain_i].size(),std::numeric_limits<double>::quiet_NaN());
      _components_dn[domain_i].assign(component_name[domain_i].size(),std::numeric_limits<double>::quiet_NaN());
      for (std::size_t comp_i = 0; comp_i < component_name[domain_i].size();
           comp_i++) {
        for (std::size_t domain_o = 0; domain_o < _lops.size(); ++domain_o) {
          if (domain_i == domain_o)
            continue;
          for (std::size_t comp_o = 0;
               comp_o < component_name[domain_o].size();
               comp_o++) {
            if (component_name[domain_i][comp_i] ==
                component_name[domain_o][comp_o]) {
              std::array<std::size_t, 3> key_i{ domain_i, domain_o, comp_i };
              _component_offset.insert(std::make_pair(key_i, comp_o));
            }
          }
        }
      }
    }

    // Helper lambda to register variables into the parser
    auto add_var = [&](auto& parser, const auto& expr, const auto& name, auto& var) {
      if (expr.find(name) != std::string::npos){
        _logger.trace(4,"Adding variable: {}"_fmt, name);
        parser.DefineVar(name, &var);
      }
    };

    // Create outflow parsers
    for (std::size_t domain_i = 0; domain_i < _lops.size(); ++domain_i) {
      for (std::size_t domain_o = 0; domain_o < _lops.size(); ++domain_o) {
        _logger.trace("Transmission condition: {} - {}"_fmt, compartment_name[domain_i], compartment_name[domain_o]);
        if (domain_i == domain_o)
          continue;
        auto& outflow_config = config.sub(compartment_name[domain_i],true).sub("outflow", true);
        if (not outflow_config.hasSub(compartment_name[domain_o]))
          continue;
        auto& outflow_config_o = outflow_config.sub(compartment_name[domain_o],true);

        // std::string outflow_jac_section = compartment_name[domain_i]
        //                                   + ".outflow."
        //                                   + compartment_name[domain_o]
        //                                   + ".jacobian";
        // auto& outflow_jac_config = config.sub(outflow_jac_section);

        auto& parser = _outflow_parser[{domain_i,domain_o}];
        auto& parser_jac = _outflow_jac_parser[{domain_i,domain_o}];
        auto& parser_cross_jac = _outflow_cross_jac_parser[{domain_i,domain_o}];
        std::size_t comp_size_i = component_name[domain_i].size();
        std::size_t comp_size_o = component_name[domain_o].size();
        parser.resize(comp_size_i);
        parser_jac.resize(comp_size_i * comp_size_i);
        parser_cross_jac.resize(comp_size_i * comp_size_o);

        // Do parser
        for (std::size_t outflow_i = 0; outflow_i < comp_size_i; ++outflow_i) {
          std::string expr = outflow_config_o.template get<std::string>(component_name[domain_i][outflow_i]);
          _logger.debug(2,"Setup expression ({}): {}"_fmt, component_name[domain_i][outflow_i], expr);

          auto& parser_i = parser.at(outflow_i);
          for (std::size_t comp_i = 0; comp_i < comp_size_i; ++comp_i) {
            add_var(parser_i,expr,component_name[domain_i][comp_i]+"_i", _components[domain_i][comp_i]);
            add_var(parser_i,expr,"d"+component_name[domain_i][comp_i]+"_i__dn", _components_dn[domain_i][comp_i]);
          }
          for (std::size_t comp_o = 0; comp_o < comp_size_o; ++comp_o) {
            add_var(parser_i,expr,component_name[domain_o][comp_o]+"_o", _components[domain_o][comp_o]);
            add_var(parser_i,expr,"d"+component_name[domain_o][comp_o]+"_o__dn", _components_dn[domain_o][comp_o]);
          }

          try {
            parser_i.SetExpr(expr);
            // try to compile expression
            parser_i.Eval();
          } catch (mu::Parser::exception_type& e) {
            handle_parser_error(e);
          }

          if (JM == JacobianMethod::Numerical)
            continue;

        // std::string outflow_jac_section = compartment_name[domain_i]
        //                                   + ".outflow."
        //                                   + compartment_name[domain_o]
        //                                   + ".jacobian";
        // auto& outflow_jac_config = config.sub(outflow_jac_section);

          // // Do self jacobian
          // for (std::size_t outflow_ii = 0; outflow_ii < comp_size_i; ++outflow_ii) {
          //   std::string jac_name = "d"+ component_name[domain_i][outflow_i] + "__d" + component_name[domain_i][outflow_ii] + "_i";
          //   std::string default_flux = "0";
          //   std::string expr = outflow_jac_config.get(jac_name,default_flux);
          //     _logger.debug(2,"Setup self jacobian expression ({}): {}"_fmt, jac_name, expr);

          //   std::size_t jac_index = outflow_i*comp_size_i+outflow_ii;
          //   auto& parser_i = parser_jac.at(jac_index);
          //   for (std::size_t comp_i = 0; comp_i < comp_size_i; ++comp_i) {
          //     add_var(parser_i,expr,component_name[domain_i][comp_i]+"_i", _components[domain_i][comp_i]);
          //     add_var(parser_i,expr,"d"+component_name[domain_i][comp_i]+"_i__dn", _components_dn[domain_i][comp_i]);
          //   }
          //   for (std::size_t comp_o = 0; comp_o < _components[domain_o].size(); ++comp_o) {
          //     if (_component_offset.find({domain_o,domain_i,comp_o}) == _component_offset.end())
          //       continue;
          //     add_var(parser_i,expr,component_name[domain_o][comp_o]+"_o", _components[domain_o][comp_o]);
          //     add_var(parser_i,expr,"d"+component_name[domain_o][comp_o]+"_o__dn", _components_dn[domain_o][comp_o]);
          //   }

          //   try {
          //     parser_i.SetExpr(expr);
          //     // try to compile expression
          //     parser_i.Eval();
          //   } catch (mu::Parser::exception_type& e) {
          //     handle_parser_error(e);
          //   }
          // }

          // // Do cross jacobian
          // for (std::size_t outflow_io = 0; outflow_io < comp_size_o; ++outflow_io) {

          //   std::string jac_name = "d"+ component_name[domain_i][outflow_i] + "__d" + component_name[domain_o][outflow_io]+"_o";
          //   bool outflux = component_name[domain_i][outflow_i] == component_name[domain_o][outflow_io];
          //   std::string default_flux = outflux ? "-1" : "0";
          //   std::string expr = outflow_jac_config.get(jac_name,default_flux);
          //   _logger.debug(2,"Setup cross jacobian expression ({}): {}"_fmt, jac_name, expr);

          //   std::size_t jac_index = outflow_i*comp_size_o+outflow_io;
          //   auto& parser_i = parser_cross_jac.at(jac_index);
          //   for (std::size_t comp_i = 0; comp_i < comp_size_i; ++comp_i) {
          //     add_var(parser_i,expr,component_name[domain_i][comp_i]+"_i", _components[domain_i][comp_i]);
          //     add_var(parser_i,expr,"d"+component_name[domain_i][comp_i]+"_i__dn", _components_dn[domain_i][comp_i]);
          //   }
          //   for (std::size_t comp_o = 0; comp_o < comp_size_o; ++comp_o) {
          //     if (_component_offset.find({domain_o,domain_i,comp_o}) == _component_offset.end())
          //       continue;
          //     add_var(parser_i,expr,component_name[domain_o][comp_o]+"_o", _components[domain_o][comp_o]);
          //     add_var(parser_i,expr,"d"+component_name[domain_o][comp_o]+"_o__dn", _components_dn[domain_o][comp_o]);
          //   }

          //   try {
          //     parser_i.SetExpr(expr);
          //     // try to compile expression
          //     parser_i.Eval();
          //   } catch (mu::Parser::exception_type& e) {
          //     handle_parser_error(e);
          //   }
          // }
        }
      }
    }

  }

  template<class EG>
  std::size_t subDomain(const EG& entity) const
  {
    auto domain_set = _index_set.subDomains(entity);
    assert(domain_set.size() == 1);
    return *(domain_set.begin());
  }

  template<class IG>
  void skip_intersection(const IG& ig, bool& skip) const
  {
    const auto& entity_i = ig.inside();
    const auto& entity_o = ig.outside();

    std::size_t domain_i = subDomain(entity_i);
    std::size_t domain_o = subDomain(entity_o);

    if (domain_i == domain_o) {
      if constexpr (SubLOP::doSkipIntersection)
        _local_operator[domain_i]->skip_intersection(ig,skip);
      else
        skip |= true;
    }
  }

  /**
   * @copydoc LocalOperatorDiffusionReactionCG::pattern_volume
   */
  template<typename LFSU, typename LFSV, typename LocalPattern>
  void pattern_volume(const LFSU& lfsu,
                      const LFSV& lfsv,
                      LocalPattern& pattern) const
  {
    for (std::size_t i = 0; i < lfsu.degree(); ++i) {
      _lops[i]->pattern_volume(lfsu.child(i), lfsv.child(i), pattern);
    }
  }

  /**
   * @brief      Pattern sckeleton
   * @details    This method links degrees of freedom between trial and test
   *             spaces at entities intersection taking into account the
   *             structure of the reaction term
   *
   * @param[in]  lfsu_i        The inside trial local function space
   * @param[in]  lfsv_i        The inside test local function space
   * @param[in]  lfsu_o        The outside trial local function space
   * @param[in]  lfsv_o        The outside test local function space
   * @param      pattern_io    The inside-outside local pattern
   * @param      pattern_oi    The outside-inside local pattern
   *
   * @tparam     LFSU          The trial local function space
   * @tparam     LFSV          The test local function space
   * @tparam     LocalPattern  The local pattern
   */
  // template<typename LFSU, typename LFSV, typename LocalPattern>
  // void interface_pattern_skeleton(std::size_t domain_i, std::size_t domain_o,
  //                       const LFSU& lfsu_di,
  //                       const LFSV& lfsv_di,
  //                       const LFSU& lfsu_do,
  //                       const LFSV& lfsv_do,
  //                       LocalPattern& pattern_io,
  //                       LocalPattern& pattern_oi) const
  // {
  //   auto lfs_size_i = lfsu_di.degree();
  //   for (std::size_t comp_i = 0; comp_i < lfs_size_i; ++comp_i) {
  //     std::array<std::size_t, 3> inside_comp{ domain_i, domain_o, comp_i };
  //     auto it = _component_offset.find(inside_comp);
  //     if (it != _component_offset.end()) {
  //       auto comp_o = it->second;
  //       auto& lfsv_ci = lfsv_di.child(comp_i);
  //       auto& lfsu_co = lfsu_do.child(comp_o);
  //       for (std::size_t dof_i = 0; dof_i < lfsv_ci.size(); dof_i++) {
  //         for (std::size_t dof_o = 0; dof_o < lfsu_co.size(); dof_o++) {
  //           pattern_io.addLink(lfsv_ci, dof_i, lfsu_co, dof_o);
  //         }
  //       }
  //     }
  //   }

  //   auto lfs_size_o = lfsu_do.degree();
  //   for (std::size_t comp_o = 0; comp_o < lfs_size_o; ++comp_o) {
  //     std::array<std::size_t, 3> outside_comp{ domain_o, domain_i, comp_o };
  //     auto it = _component_offset.find(outside_comp);
  //     if (it != _component_offset.end()) {
  //       auto comp_i = it->second;
  //       auto& lfsv_co = lfsv_do.child(comp_o);
  //       auto& lfsu_ci = lfsu_di.child(comp_i);
  //       for (std::size_t dof_o = 0; dof_o < lfsv_co.size(); dof_o++) {
  //         for (std::size_t dof_i = 0; dof_i < lfsu_ci.size(); dof_i++) {
  //           pattern_oi.addLink(lfsv_co, dof_o, lfsu_ci, dof_i);
  //         }
  //       }
  //     }
  //   }
  // }

  /**
   * @brief      Pattern sckeleton
   * @details    This method links degrees of freedom between trial and test
   *             spaces at entities intersection taking into account the
   *             structure of the reaction term
   *
   * @param[in]  lfsu_i        The inside trial local function space
   * @param[in]  lfsv_i        The inside test local function space
   * @param[in]  lfsu_o        The outside trial local function space
   * @param[in]  lfsv_o        The outside test local function space
   * @param      pattern_io    The inside-outside local pattern
   * @param      pattern_oi    The outside-inside local pattern
   *
   * @tparam     LFSU          The trial local function space
   * @tparam     LFSV          The test local function space
   * @tparam     LocalPattern  The local pattern
   */
  template<typename LFSU, typename LFSV, typename LocalPattern>
  void pattern_skeleton(const LFSU& lfsu_i,
                        const LFSV& lfsv_i,
                        const LFSU& lfsu_o,
                        const LFSV& lfsv_o,
                        LocalPattern& pattern_io,
                        LocalPattern& pattern_oi) const
  {
    std::size_t domain_i(unused_domain), domain_o(unused_domain);
    for (std::size_t k = 0; k < _lops.size(); k++) {
      if (lfsu_i.child(k).size() > 0)
        domain_i = k;
      if (lfsu_o.child(k).size() > 0)
        domain_o = k;
    }

    if ((domain_i == unused_domain) or (domain_o == unused_domain))
      return;

    const auto& lfsu_di = lfsu_i.child(domain_i);
    const auto& lfsv_di = lfsv_i.child(domain_i);
    const auto& lfsu_do = lfsu_o.child(domain_o);
    const auto& lfsv_do = lfsv_o.child(domain_o);

    assert(lfsu_i.degree() == _lops.size());
    assert(lfsu_di.degree() == lfsv_di.degree());
    assert(lfsu_do.degree() == lfsv_do.degree());

    if (domain_i != domain_o)
      Dune::PDELab::FullSkeletonPattern{}.pattern_skeleton(lfsu_di,lfsv_di,lfsu_do,lfsv_do,pattern_io,pattern_oi);
      // interface_pattern_skeleton(domain_i,domain_o,lfsu_di,lfsv_di,lfsu_do,lfsv_do,pattern_io,pattern_oi);
    else if constexpr (SubLOP::doPatternSkeleton)
      _lops[domain_i]->pattern_skeleton(lfsu_di,lfsv_di,lfsu_do,lfsv_do,pattern_io,pattern_oi);
  }

  /**
   * @brief      Sets the time.
   *
   * @param[in]  t     The new time
   */
  void setTime(double t)
  {
    Dune::PDELab::InstationaryLocalOperatorDefaultMethods<double>::setTime(t);
    for (auto lp : _lops)
      lp->setTime(t);
  }

  /**
   * @copydoc LocalOperatorDiffusionReactionCG::jacobian_apply_volume
   * @details    This particular operator does a jacobian apply volume for the
   *             LocalOperatorDiffusionReactionCG corresponding to incoming
   * entity
   */
  template<typename EG, typename LFSU, typename X, typename LFSV, typename R>
  void jacobian_apply_volume(const EG& eg,
                             const LFSU& lfsu,
                             const X& x,
                             const X& z,
                             const LFSV& lfsv,
                             R& r) const
  {
    for (std::size_t i = 0; i < _lops.size(); ++i) {
      if (lfsu.child(i).size() > 0) {
        const auto& sub_lfsu = lfsu.child(i);
        const auto& sub_lfsv = lfsv.child(i);
        _lops[i]->jacobian_apply_volume(
          eg, sub_lfsu, x, z, sub_lfsv, r);
      }
    }
  }

  /**
   * @copydoc LocalOperatorDiffusionReactionCG::jacobian_apply_volume
   * @details    This particular operator does a jacobian apply volume for the
   *             LocalOperatorDiffusionReactionCG corresponding to incoming
   * entity
   */
  template<typename EG, typename LFSU, typename X, typename LFSV, typename R>
  void jacobian_apply_volume(const EG& eg,
                             const LFSU& lfsu,
                             const X& x,
                             const LFSV& lfsv,
                             R& r) const
  {
    for (std::size_t i = 0; i < _lops.size(); ++i) {
      if (lfsu.child(i).size() > 0) {
        const auto& sub_lfsu = lfsu.child(i);
        const auto& sub_lfsv = lfsv.child(i);
        _lops[i]->jacobian_apply_volume(eg, sub_lfsu, x, sub_lfsv, r);
      }
    }
  }

  /**
   * @copydoc LocalOperatorDiffusionReactionCG::jacobian_volume
   * @details    This particular operator does a jacobian volume for the
   *             LocalOperatorDiffusionReactionCG corresponding to incoming
   * entity
   */
  template<typename EG, typename LFSU, typename X, typename LFSV, typename M>
  void jacobian_volume(const EG& eg,
                       const LFSU& lfsu,
                       const X& x,
                       const LFSV& lfsv,
                       M& mat) const
  {
    for (std::size_t i = 0; i < _lops.size(); ++i) {
      if (lfsu.child(i).size() > 0) {
        const auto& sub_lfsu = lfsu.child(i);
        const auto& sub_lfsv = lfsv.child(i);
        _lops[i]->jacobian_volume(eg, sub_lfsu, x, sub_lfsv, mat);
      }
    }
  }

  /**
   * @copydoc LocalOperatorDiffusionReactionCG::alpha_volume
   * @details    This particular operator does a alpha volume for the
   *             LocalOperatorDiffusionReactionCG corresponding to incoming
   * entity
   */
  template<typename EG, typename LFSU, typename X, typename LFSV, typename R>
  void alpha_volume(const EG& eg,
                    const LFSU& lfsu,
                    const X& x,
                    const LFSV& lfsv,
                    R& r) const
  {
    // const auto& subdomain = eg.subDomain();
    for (std::size_t i = 0; i < _lops.size(); ++i) {
      if (lfsu.child(i).size() > 0) {
        const auto& sub_lfsu = lfsu.child(i);
        const auto& sub_lfsv = lfsv.child(i);
        _lops[i]->alpha_volume(eg, sub_lfsu, x, sub_lfsv, r);
      }
    }
  }

  /**
   * @brief      The skeleton integral
   * @details    This integral is only performed at the interface between
   *             different domains. Currently it has the form of
   *             dichlet-dirichlet boundary condition between domains
   *
   * @param[in]  ig      The intersection
   * @param[in]  lfsu_i  The inside trial local function space
   * @param[in]  x_i     The inside local coefficient vector
   * @param[in]  lfsv_i  The inside test local function space
   * @param[in]  lfsu_o  The outside trial local function space
   * @param[in]  x_o     The outside local coefficient vector
   * @param[in]  lfsv_o  The outside test local function space
   * @param      r_i     The inside residual vector
   * @param      r_o     The outside residual vector
   *
   * @tparam     IG      The indersection
   * @tparam     LFSU    The trial local function space
   * @tparam     X       The local coefficient vector
   * @tparam     LFSV    The test local function space
   * @tparam     R       The residual vector
   */
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
    const auto& entity_i = ig.inside();
    const auto& entity_o = ig.outside();

    assert(_lops.size() >= lfsu_i.degree());
    assert(_lops.size() >= lfsu_o.degree());
    assert(lfsu_i.degree() == lfsv_i.degree());
    assert(lfsu_o.degree() == lfsv_o.degree());

    std::size_t domain_i = subDomain(entity_i);
    std::size_t domain_o = subDomain(entity_o);

    const auto& lfsu_di = lfsu_i.child(domain_i);
    const auto& lfsv_di = lfsv_i.child(domain_i);
    const auto& lfsu_do = lfsu_o.child(domain_o);
    const auto& lfsv_do = lfsv_o.child(domain_o);

    assert(lfsu_di.size() == lfsv_di.size());
    assert(lfsu_do.size() == lfsv_do.size());

    if (domain_i != domain_o)
      interface_alpha_skeleton(domain_i,domain_o,ig,lfsu_di,x_i,lfsv_di,lfsu_do,x_o,lfsv_do,r_i,r_o);
    else if constexpr (SubLocalOperator::doAlphaSkeleton)
      _lops[domain_i]->alpha_skeleton(ig,lfsu_di,x_i,lfsv_di,lfsu_do,x_o,lfsv_do,r_i,r_o);
  }

  /**
   * @brief      The skeleton integral
   * @details    This integral is only performed at the interface between
   *             different domains. Currently it has the form of
   *             dichlet-dirichlet boundary condition between domains
   *
   * @param[in]  ig      The intersection
   * @param[in]  lfsu_i  The inside trial local function space
   * @param[in]  x_i     The inside local coefficient vector
   * @param[in]  lfsv_i  The inside test local function space
   * @param[in]  lfsu_o  The outside trial local function space
   * @param[in]  x_o     The outside local coefficient vector
   * @param[in]  lfsv_o  The outside test local function space
   * @param      r_i     The inside residual vector
   * @param      r_o     The outside residual vector
   *
   * @tparam     IG      The indersection
   * @tparam     LFSU    The trial local function space
   * @tparam     X       The local coefficient vector
   * @tparam     LFSV    The test local function space
   * @tparam     R       The residual vector
   */
  template<typename IG, typename LFSU, typename X, typename LFSV, typename R>
  void interface_alpha_skeleton(std::size_t domain_i, std::size_t domain_o,
                      const IG& ig,
                      const LFSU& lfsu_di,
                      const X& x_i,
                      const LFSV& lfsv_di,
                      const LFSU& lfsu_do,
                      const X& x_o,
                      const LFSV& lfsv_do,
                      R& r_i,
                      R& r_o) const
  {
    assert(lfsu_di.size() == lfsv_di.size());
    assert(lfsu_do.size() == lfsv_do.size());

    const auto& entity_f = ig.intersection();
    const auto& entity_i = ig.inside();
    const auto& entity_o = ig.outside();

    auto geo_f = entity_f.geometry();
    auto geo_i = entity_i.geometry();
    auto geo_o = entity_o.geometry();
    auto geo_in_i = entity_f.geometryInInside();
    auto geo_in_o = entity_f.geometryInOutside();

    std::size_t components_i = _components[domain_i].size();
    std::size_t components_o = _components[domain_o].size();
    assert(components_i == lfsu_di.degree());
    assert(components_o == lfsu_do.degree());

    auto x_coeff_local_i = [&](const std::size_t& component,
                               const std::size_t& dof) {
      return x_i(lfsu_di.child(component), dof);
    };
    auto x_coeff_local_o = [&](const std::size_t& component,
                               const std::size_t& dof) {
      return x_o(lfsu_do.child(component), dof);
    };

    auto accumulate_i = [&](const std::size_t& component,
                            const std::size_t& dof,
                            const auto& value) {
      r_i.accumulate(lfsu_di.child(component), dof, value);
    };

    auto accumulate_o = [&](const std::size_t& component,
                            const std::size_t& dof,
                            const auto& value) {
      r_o.accumulate(lfsu_do.child(component), dof, value);
    };

    const auto& local_basis_i = lfsu_di.child(0).finiteElement().localBasis();
    const auto& local_basis_o = lfsu_do.child(0).finiteElement().localBasis();

    using LocalBasis = std::decay_t<decltype(local_basis_i)>;
    using RF = typename LocalBasis::Traits::RangeFieldType;
    using Range = typename LocalBasis::Traits::RangeType;
    using Jacobian = typename LocalBasis::Traits::JacobianType;
    constexpr int dim = Grid::dimension;

    std::vector<Range> phiu_i;
    std::vector<Range> phiu_o;

    std::vector<Jacobian> jacphi_i;
    std::vector<Jacobian> jacphi_o;

    std::vector<FieldVector<RF, dim>> gradphi_i(local_basis_i.size());
    std::vector<FieldVector<RF, dim>> gradphi_o(local_basis_o.size());

    auto normal_f = ig.centerUnitOuterNormal();

    for (const auto& it : quadratureRule(geo_f, 3)) {
      const auto& position_f = it.position();

      // position of quadrature point in local coordinates of elements
      const auto position_i = geo_in_i.global(position_f);
      const auto position_o = geo_in_o.global(position_f);

      // evaluate basis functions
      local_basis_i.evaluateFunction(position_i, phiu_i);
      local_basis_o.evaluateFunction(position_o, phiu_o);

      // evaluate concentrations at quadrature point
      std::fill(_components[domain_i].begin(),_components[domain_i].end(),0.);
      std::fill(_components[domain_o].begin(),_components[domain_o].end(),0.);
      for (std::size_t comp = 0; comp < components_i; comp++)
        for (std::size_t dof = 0; dof < local_basis_i.size(); dof++)
          _components[domain_i][comp] += x_coeff_local_i(comp, dof) * phiu_i[dof];

      for (std::size_t comp = 0; comp < components_o; comp++)
        for (std::size_t dof = 0; dof < local_basis_o.size(); dof++)
          _components[domain_o][comp] += x_coeff_local_o(comp, dof) * phiu_o[dof];

      std::fill(_components_dn[domain_i].begin(),_components_dn[domain_i].end(),0.);
      std::fill(_components_dn[domain_o].begin(),_components_dn[domain_o].end(),0.);

      std::fill(gradphi_i.begin(), gradphi_i.end(), 0.);
      std::fill(gradphi_o.begin(), gradphi_o.end(), 0.);

      if (local_basis_i.order() == 0)
      {
        RF dn = (geo_f.center() - geo_i.center()).two_norm();
        for (std::size_t comp_i = 0; comp_i < components_i; comp_i++)
        {
          auto comp_o_it = _component_offset.find({domain_i,domain_o,comp_i});
          if (comp_o_it != _component_offset.end())
          {
            std::size_t comp_o = comp_o_it->second;
            _components_dn[domain_i][comp_i] = (_components[domain_o][comp_o] - _components[domain_i][comp_i])/dn;
          }
        }
      } else {
        local_basis_i.evaluateJacobian(position_i,jacphi_i);
        auto jac_i = geo_i.jacobianInverseTransposed(position_i);
        for (std::size_t i=0; i<gradphi_i.size(); i++)
          jac_i.mv(jacphi_i[i][0],gradphi_i[i]);
        for (std::size_t comp_i = 0; comp_i < components_i; comp_i++)
          for (std::size_t j = 0; j < gradphi_i.size(); j++)
            _components_dn[domain_i][comp_i] += (normal_f * gradphi_i[j]) * x_coeff_local_i(comp_i, j);
      }

      if (local_basis_o.order() == 0)
      {
        RF dn = (geo_f.center() - geo_o.center()).two_norm();
        for (std::size_t comp_o = 0; comp_o < components_i; comp_o++)
        {
          auto comp_i_it = _component_offset.find({domain_o,domain_i,comp_o});
          if (comp_i_it != _component_offset.end())
          {
            std::size_t comp_i = comp_i_it->second;
            _components_dn[domain_o][comp_o] = (_components[domain_o][comp_o] - _components[domain_i][comp_i])/dn;
          }
        }
      } else {
        local_basis_o.evaluateJacobian(position_o,jacphi_o);
        auto jac_o = geo_o.jacobianInverseTransposed(position_o);
        for (std::size_t i=0; i<gradphi_o.size(); i++)
          jac_o.mv(jacphi_o[i][0],gradphi_o[i]);

        for (std::size_t comp_o = 0; comp_o < components_o; comp_o++)
          for (std::size_t j = 0; j < gradphi_i.size(); j++)
            _components_dn[domain_o][comp_o] += (normal_f * gradphi_o[j]) * x_coeff_local_o(comp_o, j);
      }

      // integration factor
      auto factor = it.weight() * geo_f.integrationElement(position_f);

      auto outflow_parser_io_it = _outflow_parser.find({domain_i,domain_o});
      if (outflow_parser_io_it == _outflow_parser.end())
        DUNE_THROW(IOError,"Outflux section is missing");
      const auto& parser_io = outflow_parser_io_it->second;
      for (std::size_t comp_i = 0; comp_i < components_i; comp_i++) {
        double value = parser_io[comp_i].Eval();
        for (std::size_t dof = 0; dof < lfsu_di.child(comp_i).size(); dof++)
          accumulate_i(comp_i, dof, factor * value * phiu_i[dof]);
      }

      auto outflow_parser_oi_it = _outflow_parser.find({domain_o,domain_i});
      if (outflow_parser_oi_it == _outflow_parser.end())
        DUNE_THROW(IOError,"Outflux section is missing");
      const auto& parser_oi = outflow_parser_oi_it->second;
      for (std::size_t comp_o = 0; comp_o < components_o; comp_o++) {
       double value = parser_oi[comp_o].Eval();
        for (std::size_t dof = 0; dof < lfsu_do.child(comp_o).size(); dof++)
          accumulate_o(comp_o, dof, factor * value * phiu_o[dof]);
      }
    }
  }

  /**
   * @brief      The jacobian skeleton integral
   * @copydetails alpha_skeleton
   *
   * @param[in]  ig      The intersection
   * @param[in]  lfsu_i  The inside trial local function space
   * @param[in]  x_i     The inside local coefficient vector
   * @param[in]  lfsv_i  The inside test local function space
   * @param[in]  lfsu_o  The outside trial local function space
   * @param[in]  x_o     The outside local coefficient vector
   * @param[in]  lfsv_o  The outside test local function space
   * @param      mat_ii  The local inside-inside matrix
   * @param      mat_io  The local inside-outside matrix
   * @param      mat_oi  The local outside-inside matrix
   * @param      mat_oo  The local outside-outside matrix
   *
   * @tparam     IG      The indersection
   * @tparam     LFSU    The trial local function space
   * @tparam     X       The local coefficient vector
   * @tparam     LFSV    The test local function space
   * @tparam     J       The local jacobian matrix
   */
  template<typename IG, typename LFSU, typename X, typename LFSV, typename J>
  void jacobian_skeleton(const IG& ig,
                         const LFSU& lfsu_i,
                         const X& x_i,
                         const LFSV& lfsv_i,
                         const LFSU& lfsu_o,
                         const X& x_o,
                         const LFSV& lfsv_o,
                         J& mat_ii,
                         J& mat_io,
                         J& mat_oi,
                         J& mat_oo) const
  {
    assert(_lops.size() >= lfsu_i.degree());
    assert(_lops.size() >= lfsu_o.degree());
    assert(lfsu_i.degree() == lfsv_i.degree());
    assert(lfsu_o.degree() == lfsv_o.degree());

    const auto& entity_i = ig.inside();
    const auto& entity_o = ig.outside();

    std::size_t domain_i = subDomain(entity_i);
    std::size_t domain_o = subDomain(entity_o);

    const auto& lfsu_di = lfsu_i.child(domain_i);
    const auto& lfsv_di = lfsv_i.child(domain_i);
    const auto& lfsu_do = lfsu_o.child(domain_o);
    const auto& lfsv_do = lfsv_o.child(domain_o);

    assert(lfsu_di.size() == lfsv_di.size());
    assert(lfsu_do.size() == lfsv_do.size());

    if (domain_i == domain_o)
    {
      if constexpr (SubLocalOperator::doAlphaSkeleton)
      {
        _lops[domain_i]->jacobian_skeleton(ig,lfsu_di,x_i,lfsv_di,lfsu_do,x_o,lfsv_do,mat_ii,mat_io,mat_oi,mat_oo);
      }
    }
    else if constexpr (JM == JacobianMethod::Numerical)
      PDELab::NumericalJacobianSkeleton<
        LocalOperatorMultiDomainDiffusionReaction>::jacobian_skeleton(ig,
                                                                      lfsu_i,
                                                                      x_i,
                                                                      lfsv_i,
                                                                      lfsu_o,
                                                                      x_o,
                                                                      lfsv_o,
                                                                      mat_ii,
                                                                      mat_io,
                                                                      mat_oi,
                                                                      mat_oo);
    else
      interface_jacobian_skeleton(domain_i,domain_o,ig,lfsu_di,x_i,lfsv_di,lfsu_do,x_o,lfsv_do,mat_ii,mat_io,mat_oi,mat_oo);
  }

  /**
   * @brief      The jacobian skeleton integral
   * @copydetails alpha_skeleton
   *
   * @param[in]  ig      The intersection
   * @param[in]  lfsu_i  The inside trial local function space
   * @param[in]  x_i     The inside local coefficient vector
   * @param[in]  lfsv_i  The inside test local function space
   * @param[in]  lfsu_o  The outside trial local function space
   * @param[in]  x_o     The outside local coefficient vector
   * @param[in]  lfsv_o  The outside test local function space
   * @param      mat_ii  The local inside-inside matrix
   * @param      mat_io  The local inside-outside matrix
   * @param      mat_oi  The local outside-inside matrix
   * @param      mat_oo  The local outside-outside matrix
   *
   * @tparam     IG      The indersection
   * @tparam     LFSU    The trial local function space
   * @tparam     X       The local coefficient vector
   * @tparam     LFSV    The test local function space
   * @tparam     J       The local jacobian matrix
   */
  template<typename IG, typename LFSU, typename X, typename LFSV, typename J>
  void interface_jacobian_skeleton(std::size_t domain_i, std::size_t domain_o,
                         const IG& ig,
                         const LFSU& lfsu_di,
                         const X& x_i,
                         const LFSV& lfsv_di,
                         const LFSU& lfsu_do,
                         const X& x_o,
                         const LFSV& lfsv_do,
                         J& mat_ii,
                         J& mat_io,
                         J& mat_oi,
                         J& mat_oo) const
  {
    static_assert(AlwaysFalse<J>{}, "Not implemented, please use numerical jacobian");
    assert(lfsu_di.size() == lfsv_di.size());
    assert(lfsu_do.size() == lfsv_do.size());


    const auto& entity_f = ig.intersection();
    const auto& entity_i = ig.inside();
    const auto& entity_o = ig.outside();

    auto geo_f = entity_f.geometry();
    auto geo_i = entity_i.geometry();
    auto geo_o = entity_o.geometry();
    auto geo_in_i = entity_f.geometryInInside();
    auto geo_in_o = entity_f.geometryInOutside();

    std::size_t components_i = _components[domain_i].size();
    std::size_t components_o = _components[domain_o].size();
    assert(components_i == lfsu_di.degree());
    assert(components_o == lfsu_do.degree());

    auto x_coeff_local_i = [&](const std::size_t& component,
                               const std::size_t& dof) {
      return x_i(lfsu_di.child(component), dof);
    };
    auto x_coeff_local_o = [&](const std::size_t& component,
                               const std::size_t& dof) {
      return x_o(lfsu_do.child(component), dof);
    };

    auto accumulate_ii = [&](const std::size_t& component_i,
                             const std::size_t& dof_i,
                             const std::size_t& component_j,
                             const std::size_t& dof_j,
                             const auto& value) {
      mat_ii.accumulate(lfsu_di.child(component_i),
                        dof_i,
                        lfsu_di.child(component_j),
                        dof_j,
                        value);
    };

    auto accumulate_io = [&](const std::size_t& component_i,
                             const std::size_t& dof_i,
                             const std::size_t& component_j,
                             const std::size_t& dof_j,
                             const auto& value) {
      mat_io.accumulate(lfsu_di.child(component_i),
                        dof_i,
                        lfsu_do.child(component_j),
                        dof_j,
                        value);
    };

    auto accumulate_oi = [&](const std::size_t& component_i,
                             const std::size_t& dof_i,
                             const std::size_t& component_j,
                             const std::size_t& dof_j,
                             const auto& value) {
      mat_oi.accumulate(lfsu_do.child(component_i),
                        dof_i,
                        lfsu_di.child(component_j),
                        dof_j,
                        value);
    };

    auto accumulate_oo = [&](const std::size_t& component_i,
                             const std::size_t& dof_i,
                             const std::size_t& component_j,
                             const std::size_t& dof_j,
                             const auto& value) {
      mat_oo.accumulate(lfsu_do.child(component_i),
                        dof_i,
                        lfsu_do.child(component_j),
                        dof_j,
                        value);
    };

    const auto& local_basis_i = lfsu_di.child(0).finiteElement().localBasis();
    const auto& local_basis_o = lfsu_do.child(0).finiteElement().localBasis();

    using LocalBasis = std::decay_t<decltype(local_basis_i)>;
    using RF = typename LocalBasis::Traits::RangeFieldType;
    using Range = typename LocalBasis::Traits::RangeType;
    using Jacobian = typename LocalBasis::Traits::JacobianType;
    constexpr int dim = Grid::dimension;

    std::vector<Range> phiu_i;
    std::vector<Range> phiu_o;

    std::vector<Jacobian> jacphi_i;
    std::vector<Jacobian> jacphi_o;

    std::vector<FieldVector<RF, dim>> gradphi_i(local_basis_i.size());
    std::vector<FieldVector<RF, dim>> gradphi_o(local_basis_o.size());

    auto normal_f = ig.centerUnitOuterNormal();

    auto dn_i = (geo_f.center() - geo_i.center()).two_norm();
    auto dn_o = (geo_f.center() - geo_o.center()).two_norm();

    for (const auto& it : quadratureRule(geo_f, 3)) {
      const auto& position_f = it.position();

      std::fill(gradphi_i.begin(), gradphi_i.end(), 0.);
      std::fill(gradphi_o.begin(), gradphi_o.end(), 0.);

      // position of quadrature point in local coordinates of elements
      const auto position_i = geo_in_i.global(position_f);
      const auto position_o = geo_in_o.global(position_f);

      // evaluate basis functions
      local_basis_i.evaluateFunction(position_i, phiu_i);
      local_basis_o.evaluateFunction(position_o, phiu_o);

      // evaluate concentrations at quadrature point
      std::fill(_components[domain_i].begin(),_components[domain_i].end(),0.);
      std::fill(_components[domain_o].begin(),_components[domain_o].end(),0.);

      for (std::size_t comp = 0; comp < components_i; comp++)
        for (std::size_t dof = 0; dof < local_basis_i.size(); dof++)
          _components[domain_i][comp] += x_coeff_local_i(comp, dof) * phiu_i[dof];

      for (std::size_t comp = 0; comp < components_o; comp++)
        for (std::size_t dof = 0; dof < local_basis_o.size(); dof++)
          _components[domain_o][comp] += x_coeff_local_o(comp, dof) * phiu_o[dof];

      std::fill(_components_dn[domain_i].begin(),_components_dn[domain_i].end(),0.);
      std::fill(_components_dn[domain_o].begin(),_components_dn[domain_o].end(),0.);

      if (local_basis_i.order() == 0)
      {
        for (std::size_t comp_i = 0; comp_i < components_i; comp_i++)
        {
          auto comp_o_it = _component_offset.find({domain_i,domain_o,comp_i});
          if (comp_o_it != _component_offset.end())
          {
            std::size_t comp_o = comp_o_it->second;
            _components_dn[domain_i][comp_i] = (_components[domain_o][comp_o] - _components[domain_i][comp_i])/dn_i;
          }
        }
      } else {
        local_basis_i.evaluateJacobian(position_i,jacphi_i);
        auto jac_i = geo_i.jacobianInverseTransposed(position_i);
        for (std::size_t i=0; i<gradphi_i.size(); i++)
          jac_i.mv(jacphi_i[i][0],gradphi_i[i]);
        for (std::size_t comp_i = 0; comp_i < components_i; comp_i++)
          for (std::size_t j = 0; j < gradphi_i.size(); j++)
            _components_dn[domain_i][comp_i] += (normal_f * gradphi_i[j]) * x_coeff_local_i(comp_i, j);
      }

      if (local_basis_o.order() == 0)
      {
        for (std::size_t comp_o = 0; comp_o < components_i; comp_o++)
        {
          auto comp_i_it = _component_offset.find({domain_o,domain_i,comp_o});
          if (comp_i_it != _component_offset.end())
          {
            std::size_t comp_i = comp_i_it->second;
            _components_dn[domain_o][comp_o] = - (_components[domain_i][comp_i] - _components[domain_o][comp_o])/dn_o;
          }
        }
      } else {
        local_basis_o.evaluateJacobian(position_o,jacphi_o);
        auto jac_o = geo_o.jacobianInverseTransposed(position_o);
        for (std::size_t i=0; i<gradphi_o.size(); i++)
          jac_o.mv(jacphi_o[i][0],gradphi_o[i]);

        for (std::size_t comp_o = 0; comp_o < components_o; comp_o++)
          for (std::size_t j = 0; j < gradphi_i.size(); j++)
            _components_dn[domain_o][comp_o] = - (normal_f * gradphi_o[j]) * x_coeff_local_o(comp_o, j);
      }



      // integration factor
      auto factor = it.weight() * geo_f.integrationElement(position_f);

      const auto& parser_jac_ii = _outflow_jac_parser.find({domain_i,domain_o})->second;
      const auto& parser_jac_io = _outflow_cross_jac_parser.find({domain_i,domain_o})->second;
      for (std::size_t comp_i = 0; comp_i < components_i; comp_i++)
      {
        for (std::size_t comp_ii = 0; comp_ii < components_i; comp_ii++)
          for (std::size_t i = 0; i < lfsu_di.child(comp_i).size(); i++)
            for (std::size_t j = 0; j < lfsu_di.child(comp_ii).size(); j++)
              accumulate_ii(comp_i, i, comp_ii, j, factor * parser_jac_ii[comp_i*components_i+comp_ii].Eval() * phiu_i[j]);

        for (std::size_t comp_io = 0; comp_io < components_o; comp_io++)
          for (std::size_t i = 0; i < lfsu_di.child(comp_i).size(); i++)
            for (std::size_t j = 0; j < lfsu_do.child(comp_io).size(); j++)
              accumulate_io(comp_i, i, comp_io, j, factor * parser_jac_io[comp_i*components_o+comp_io].Eval() * phiu_o[j]);
      }

      const auto& parser_jac_oo = _outflow_jac_parser.find({domain_o,domain_i})->second;
      const auto& parser_jac_oi = _outflow_cross_jac_parser.find({domain_o,domain_i})->second;
      for (std::size_t comp_o = 0; comp_o < components_o; comp_o++)
      {
        for (std::size_t comp_oo = 0; comp_oo < components_o; comp_oo++)
          for (std::size_t i = 0; i < lfsu_do.child(comp_o).size(); i++)
            for (std::size_t j = 0; j < lfsu_do.child(comp_oo).size(); j++)
              accumulate_oo(comp_o, i, comp_oo, j, factor * parser_jac_oo[comp_o*components_o+comp_oo].Eval() * phiu_o[j]);

        for (std::size_t comp_oi = 0; comp_oi < components_i; comp_oi++)
          for (std::size_t i = 0; i < lfsu_do.child(comp_o).size(); i++)
            for (std::size_t j = 0; j < lfsu_di.child(comp_oi).size(); j++)
              accumulate_oi(comp_o, i, comp_oi, j, factor * parser_jac_oi[comp_o*components_i+comp_oi].Eval() * phiu_i[j]);
      }
    }
  }
};

/**
 * @brief      This class describes a PDELab temporal local operator for multi
 *             domain diffusion reaction systems.
 * @details    This class describre the operatrions for local integrals required
 *             for diffusion reaction system. The operator is only valid for
 *             entities contained in the grid. The local finite element is used
 *             for caching shape function evaluations.And the jacobian method
 *             switches between numerical and analytical jacobians. This local
 *             operator creates internally an individual local operator for
 *             every subdomain in the grid
 * @todo       Add numerical jacobian methods
 *
 * @tparam     Grid                The grid
 * @tparam     LocalFiniteElement  The local finite element
 * @tparam     JM                  The jacobian method
 */
template<class Grid,
         class SubLocalOperator,
         JacobianMethod JM = JacobianMethod::Analytical>
class TemporalLocalOperatorMultiDomainDiffusionReaction
  : public Dune::PDELab::LocalOperatorDefaultFlags
  , public Dune::PDELab::InstationaryLocalOperatorDefaultMethods<double>
{
  static_assert(Concept::isMultiDomainGrid<Grid>());

  using GridView = typename Grid::SubDomainGrid::LeafGridView;
  using SubLOP = SubLocalOperator;

  std::vector<std::shared_ptr<SubLOP>> _lops;

public:
  //! pattern assembly flags
  static constexpr bool doPatternVolume = true;

  //! residual assembly flags
  static constexpr bool doAlphaVolume = true;

  /**
   * @brief      Constructs a new instance.
   *
   * @param[in]  grid            The grid
   * @param[in]  config          The configuration
   * @param[in]  finite_element  The local finite element
   */
  TemporalLocalOperatorMultiDomainDiffusionReaction(
    std::shared_ptr<const Grid> grid,
    const ParameterTree& config)
  {
    const auto& compartments = config.sub("compartments").getValueKeys();
    _lops.resize(compartments.size());

    for (std::size_t i = 0; i < _lops.size(); ++i) {
      const std::string compartement = compartments[i];

      int sub_domain_id =
        config.sub("compartments").template get<int>(compartement);
      GridView sub_grid_view = grid->subDomain(sub_domain_id).leafGridView();

      const auto& sub_config = config.sub(compartments[i]);
      _lops[i] = std::make_shared<SubLOP>(sub_grid_view, sub_config);
    }
  }

  /**
   * @copydoc TemporalLocalOperatorDiffusionReactionCG::pattern_volume
   */
  template<typename LFSU, typename LFSV, typename LocalPattern>
  void pattern_volume(const LFSU& lfsu,
                      const LFSV& lfsv,
                      LocalPattern& pattern) const
  {
    for (std::size_t i = 0; i < _lops.size(); ++i)
      _lops[i]->pattern_volume(lfsu.child(i), lfsv.child(i), pattern);
  }

  /**
   * @copydoc LocalOperatorMultiDomainDiffusionReaction::alpha_volume
   */
  template<typename EG, typename LFSU, typename X, typename LFSV, typename R>
  void alpha_volume(const EG& eg,
                    const LFSU& lfsu,
                    const X& x,
                    const LFSV& lfsv,
                    R& r) const
  {
    for (std::size_t i = 0; i < _lops.size(); ++i) {
      if (lfsu.child(i).size() > 0) {
        const auto& sub_lfsu = lfsu.child(i);
        const auto& sub_lfsv = lfsv.child(i);
        _lops[i]->alpha_volume(eg, sub_lfsu, x, sub_lfsv, r);
      }
    }
  }

  /**
   * @copydoc LocalOperatorMultiDomainDiffusionReaction::jacobian_volume
   */
  template<typename EG, typename LFSU, typename X, typename LFSV, typename Mat>
  void jacobian_volume(const EG& eg,
                       const LFSU& lfsu,
                       const X& x,
                       const LFSV& lfsv,
                       Mat& mat) const
  {
    for (std::size_t i = 0; i < _lops.size(); ++i) {
      if (lfsu.child(i).size() > 0) {
        const auto& sub_lfsu = lfsu.child(i);
        const auto& sub_lfsv = lfsv.child(i);
        _lops[i]->jacobian_volume(eg, sub_lfsu, x, sub_lfsv, mat);
      }
    }
  }

  /**
   * @copydoc LocalOperatorMultiDomainDiffusionReaction::jacobian_apply_volume
   */
  template<typename EG, typename LFSU, typename X, typename LFSV, typename R>
  void jacobian_apply_volume(const EG& eg,
                             const LFSU& lfsu,
                             const X& x,
                             const X& z,
                             const LFSV& lfsv,
                             R& r) const
  {
    for (std::size_t i = 0; i < _lops.size(); ++i) {
      if (lfsu.child(i).size() > 0) {
        const auto& sub_lfsu = lfsu.child(i);
        const auto& sub_lfsv = lfsv.child(i);
        _lops[i]->jacobian_apply_volume(
          eg, sub_lfsu, x, z, sub_lfsv, r);
      }
    }
  }

  /**
   * @copydoc LocalOperatorMultiDomainDiffusionReaction::jacobian_apply_volume
   */
  template<typename EG, typename LFSU, typename X, typename LFSV, typename R>
  void jacobian_apply_volume(const EG& eg,
                             const LFSU& lfsu,
                             const X& x,
                             const LFSV& lfsv,
                             R& r) const
  {
    for (std::size_t i = 0; i < _lops.size(); ++i) {
      if (lfsu.child(i).size() > 0) {
        const auto& sub_lfsu = lfsu.child(i);
        const auto& sub_lfsv = lfsv.child(i);
        _lops[i]->jacobian_apply_volume(eg, sub_lfsu, x, sub_lfsv, r);
      }
    }
  }
};

} // namespace Dune::Copasi

#endif // DUNE_COPASI_LOCAL_OPERATOR_MULTIDOMAIN_DIFFUSION_REACTION_HH
