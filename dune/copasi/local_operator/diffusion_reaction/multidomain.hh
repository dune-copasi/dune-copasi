#ifndef DUNE_COPASI_LOCAL_OPERATOR_MULTIDOMAIN_DIFFUSION_REACTION_HH
#define DUNE_COPASI_LOCAL_OPERATOR_MULTIDOMAIN_DIFFUSION_REACTION_HH

#include <dune/copasi/common/enum.hh>
#include <dune/copasi/parser/factory.hh>
#include <dune/copasi/concepts/grid.hh>

#include <dune/assembler/common/trace.hh>
#include <dune/assembler/concepts/discrete_function_space.hh>

#include <dune/pdelab/localoperator/flags.hh>
#include <dune/pdelab/localoperator/idefault.hh>

#include <dune/logging.hh>

#include <dune/common/parametertree.hh>

#include <algorithm>
#include <map>
#include <vector>
#include <utility>
#include <set>

namespace Dune::Copasi {

/**
 * @brief      This class describes a local operator for multi domain
 *             diffusion reaction systems.
 * @details    This class describre the operatrions for local integrals required
 *             for diffusion reaction system. This local operator maintains an
 *             individual local operator for every subdomain in the grid.
 *
 * @tparam     Space               A multidomain discrete function space
 * @tparam     DomainLocalOperator    Local operator of sub-domain
 */
template<Dune::Assembler::Concept::DiscreteFunctionSpace Space,
         class DomainLocalOperator>
class LocalOperatorMultiDomainDiffusionReaction
  : public Dune::PDELab::LocalOperatorDefaultFlags
  , public Dune::PDELab::InstationaryLocalOperatorDefaultMethods<double>
{

  static constexpr std::size_t unused_domain = ~std::size_t(0);

#if DUNE_COPASI_HAVE_MEMBRANE_SPACE
  static constexpr auto compartments_path = Assembler::multiIndex(Indices::_0);
#else
  static constexpr auto compartments_path = Assembler::multiIndex();
#endif

  using DomainLOP = DomainLocalOperator;

  Logging::Logger _logger;

  Space _space;

  // interior domain local operators
  std::vector<DomainLOP> _lops;

  // same-named components (domain_i,domain_o,component_i) -> component_o
  std::map<std::array<std::size_t, 3>, std::size_t> _component_offset;

  // multidomain cross jacobian indices (domain_i,domain_o,component_i) -> (self jacobian, cross jacobian indices)
  std::map<std::array<std::size_t, 3>, std::array<std::set<std::size_t>,2>> _outflow_jac_map;

  /// map: (domain_i,domain_o) -> vector of outflow_i parsers
  std::map<std::array<std::size_t,2>,std::vector<std::unique_ptr<Parser>>> _outflow_parser;
  std::map<std::array<std::size_t,2>,std::vector<std::unique_ptr<Parser>>> _outflow_self_jac_parser;
  std::map<std::array<std::size_t,2>,std::vector<std::unique_ptr<Parser>>> _outflow_cross_jac_parser;

#if DUNE_COPASI_HAVE_MEMBRANE_SPACE
  // jacobian pattern indices (membrane_id, component_f) -> self_jacobian_indices
  std::map<std::array<std::size_t, 2>, std::set<std::size_t>> _membrane_self_jac_map;

  // cross jacobian pattern indices (membrane_id, component_f, domain) -> cross_jacobian_indices
  std::map<std::array<std::size_t, 3>, std::set<std::size_t>> _membrane_cross_jac_map;

  // jacobian pattern indices (domain_i, domain_o, component) -> outflow/membrane indices
  std::map<std::array<std::size_t, 3>, std::set<std::size_t>> _out_membrane_jac_map;


  std::map<std::array<std::size_t,2>,std::vector<std::unique_ptr<Parser>>> _outflow_mem_jac_parser;

  // /// map: membrane_id -> vector of species parsers at the membrane
  std::vector<std::vector<std::unique_ptr<Parser>>> _diffusion_parser;
  std::vector<std::vector<std::unique_ptr<Parser>>> _reaction_parser;
  std::vector<std::vector<std::unique_ptr<Parser>>> _reaction_self_jac_parser;
  std::map<std::array<std::size_t,2>,std::vector<std::unique_ptr<Parser>>> _reaction_cross_jac_parser;

  using SubDomainIndex = typename Space::EntitySet::Grid::SubDomainIndex;
  // membrane_id -> (domain_i,domain_o)  where domain_i <= domain_o
  // std::vector<std::array<SubDomainIndex,2>> _id2membrane;
  // (domain_i,domain_o) -> membrane_id
  std::map<std::array<std::size_t,2>, std::size_t> _membrane2id;

  // membrane_id, component -> species
  mutable std::vector<std::vector<double>> _u_membrane;
#endif

  // value of trial spaces (for parser evaluation)
  std::vector<std::vector<double>> _components;
  std::vector<std::vector<double>> _components_dn;
  Dune::ParameterTree _config;

  static constexpr auto dim = Space::EntitySet::dimension;
  using BulkLocalBasis = std::decay_t<decltype(TypeTree::child(_space.localView().tree(),compartments_path).child(0).child(0).finiteElement().localBasis())>;
  using BulkRangeField = typename BulkLocalBasis::Traits::RangeFieldType;
  using BulkRange = typename BulkLocalBasis::Traits::RangeType;
  using BulkJacobian = typename BulkLocalBasis::Traits::JacobianType;

  std::vector<BulkRange> _phiu_i, _phiu_o;
  std::vector<BulkJacobian> _jacphi_i, _jacphi_o;
  std::vector<FieldVector<BulkRangeField, dim>> _gradphi_i, _gradphi_o;

#if DUNE_COPASI_HAVE_MEMBRANE_SPACE
  using MembraneLocalBasis = std::decay_t<decltype(_space.localView().tree().child(Indices::_1).child(0).child(0).child(0).finiteElement().localBasis())>;
  using MembraneRangeField = typename MembraneLocalBasis::Traits::RangeFieldType;
  using MembraneRange = typename MembraneLocalBasis::Traits::RangeType;
  using MembraneJacobian = typename MembraneLocalBasis::Traits::JacobianType;

  std::vector<MembraneRange> _phiu_f;
  std::vector<MembraneJacobian> _jacphi_f;
  mutable std::vector<FieldVector<MembraneRangeField, dim>> _gradphi_f;
#endif

public:

  static_assert(not DomainLOP::doSkeletonTwoSided);
  static_assert(not DomainLOP::doSkipEntity);

  //! selective assembly flags
  static constexpr bool doSkipEntity = false;
  static constexpr bool doSkipIntersection = true;

  //! pattern assembly flags
  static constexpr bool doPatternVolume = DomainLOP::doPatternVolume;
  static constexpr bool doPatternSkeleton = true;

  //! residual assembly flags
  static constexpr bool doAlphaVolume = DomainLOP::doAlphaVolume;
  static constexpr bool doAlphaSkeleton = true;


  /**
   * @brief      Constructs a new instance.
   *
   * @param[in]  grid            The grid
   * @param[in]  config          The configuration
   * @param[in]  finite_element  The local finite element
   */
  LocalOperatorMultiDomainDiffusionReaction(
    const Space& space,
    const ParameterTree& config)
    : _logger(Logging::Logging::componentLogger({}, "model"))
    , _space(space)
    , _config(config)
  {
    setup();
  }

  LocalOperatorMultiDomainDiffusionReaction(const LocalOperatorMultiDomainDiffusionReaction& other)
    : LocalOperatorMultiDomainDiffusionReaction(other._space, other._config)
  {
    setup();
  }

  void setup()
  {
    TRACE_EVENT("dune", "LocalOperator::MultiDomain::ParserSetUp");

    _lops.clear();
    auto compartments_space = _space.subSpace(compartments_path);
    for (std::size_t domain = 0; domain != compartments_space.degree(); ++domain) {
      auto domain_space = compartments_space.subSpace(Assembler::multiIndex(domain));
      const auto& domain_config = _config.sub(domain_space.name(), true);
      _lops.emplace_back(domain_space, domain_config);
    }

    _components.resize(compartments_space.degree());
    _components_dn.resize(compartments_space.degree());
    for (std::size_t domain = 0; domain != compartments_space.degree(); ++domain) {
      auto domain_space = compartments_space.subSpace(Assembler::multiIndex(domain));
      _components[domain].assign(domain_space.degree(), std::numeric_limits<double>::quiet_NaN());
      _components_dn[domain].assign(domain_space.degree(), std::numeric_limits<double>::quiet_NaN());
    }

#if DUNE_COPASI_HAVE_MEMBRANE_SPACE
    auto space_ms = _space.subSpace(Assembler::multiIndex(Indices::_1));
    _u_membrane.resize(space_ms.degree());
    _reaction_parser.resize(space_ms.degree());
    _diffusion_parser.resize(space_ms.degree());
    _reaction_self_jac_parser.resize(space_ms.degree());
    for (std::size_t membrane_id = 0; membrane_id != space_ms.degree(); ++membrane_id) {
      auto space_m = space_ms.subSpace(Assembler::multiIndex(membrane_id));
      _u_membrane[membrane_id].assign(space_m.degree(), std::numeric_limits<double>::quiet_NaN());
    }
#endif

   // Helper lambda to register variables into the parser
    auto add_var =
      [&](auto& parser, const auto& expr, const auto& name, auto& var) {
        if (expr.find(name) != std::string::npos) {
          _logger.trace(4, "Adding variable: {}"_fmt, name);
          parser->define_variable(name, &var);
        }
      };

    // Create boundary parsers
    for (std::size_t domain_i = 0; domain_i != compartments_space.degree(); ++domain_i) {
      for (std::size_t domain_o = 0; domain_o != compartments_space.degree(); ++domain_o) {

        if (domain_i == domain_o)
          continue;

        auto space_i = compartments_space.subSpace(Assembler::multiIndex(domain_i));
        auto space_o = compartments_space.subSpace(Assembler::multiIndex(domain_o));

        auto add_domain_vars = [&](auto domain, auto& parser, const auto& expr) {
          auto domain_space = compartments_space.subSpace(Assembler::multiIndex(domain));
          for (std::size_t comp = 0; comp != domain_space.degree(); ++comp) {
            auto space_comp = domain_space.subSpace(Assembler::multiIndex(comp));
            add_var(parser, expr, space_comp.name(), _components[domain][comp]);
            add_var(parser, expr, fmt::format("d{}__dn", space_comp.name()), _components_dn[domain][comp]);
          }
        };

#if DUNE_COPASI_HAVE_MEMBRANE_SPACE
        bool membrane_found = false;

        { // build _membrane2id from space names
          auto membrane_io = fmt::format("{}-{}", space_i.name(), space_o.name());
          auto membrane_oi = fmt::format("{}-{}", space_o.name(), space_i.name());
          for (std::size_t membrane_id = 0; membrane_id != space_ms.degree(); ++membrane_id) {
            auto space_m = space_ms.subSpace(Assembler::multiIndex(membrane_id));
            if (membrane_io == space_m.name() or membrane_oi == space_m.name()) {
              _membrane2id[{domain_i,domain_o}] = _membrane2id[{domain_o,domain_i}] = membrane_id;
              membrane_found = true;
              break;
            }
          }
        }

        if (not membrane_found)
          continue;
        auto membrane_id = _membrane2id.at({domain_i,domain_o});
        auto space_m = _space.subSpace(Assembler::multiIndex(Indices::_1, membrane_id));
        
        auto add_membrane_vars = [&](auto& parser, const auto& expr) {
          for (std::size_t comp = 0; comp != space_m.degree(); ++comp) {
            auto space_comp = space_m.subSpace(Assembler::multiIndex(comp));
            add_var(parser, expr, space_comp.name(), _u_membrane[membrane_id][comp]);
          }
        };
#else
        auto add_membrane_vars = [&](auto&, const auto&) {};
#endif
        _logger.trace("Transmission condition: {} - {}"_fmt, space_i.name(), space_o.name());

        auto add_interface_vars = [&](auto& parser, const auto& expr){
          add_domain_vars(domain_i, parser, expr);
          add_domain_vars(domain_o, parser, expr);
          add_membrane_vars(parser, expr);
        };


#if DUNE_COPASI_HAVE_MEMBRANE_SPACE
        const auto& reaction_config = _config.sub(fmt::format("{}.boundary.{}.reaction", space_i.name(), space_o.name()));
        const auto& diffusion_config = _config.sub(fmt::format("{}.boundary.{}.diffusion", space_i.name(), space_o.name()));

        auto& parser_m_react = _reaction_parser[membrane_id];
        auto& parser_m_diff = _diffusion_parser[membrane_id];
        auto& parser_m_self_jac = _reaction_self_jac_parser[membrane_id];
        auto& parser_m_cross_jac = _reaction_cross_jac_parser[{membrane_id, domain_i}];

        parser_m_react.resize(space_m.degree());
        parser_m_diff.resize(space_m.degree());
        parser_m_self_jac.resize(space_m.degree() * space_m.degree());
        parser_m_cross_jac.resize(space_m.degree() * space_i.degree());

        // Do jacobian
        for (std::size_t reaction_m = 0; reaction_m != space_m.degree(); ++reaction_m) {
          auto space_reaction_m = space_m.subSpace(Assembler::multiIndex(reaction_m));

          if (not reaction_config.hasKey(space_reaction_m.name())) // TODO assert that is in the out?
            continue;

          std::string mem_react_expr = reaction_config.template get<std::string>(space_reaction_m.name());
          std::string mem_diff_expr = diffusion_config.template get<std::string>(space_reaction_m.name());
          _logger.trace(2, "Setup reaction expression ({}): {}"_fmt, space_reaction_m.name(), mem_react_expr);
          _logger.trace(2, "Setup diffusion expression ({}): {}"_fmt, space_reaction_m.name(), mem_diff_expr);

          std::set<std::size_t> mem_self_map, mem_cross_map;

          auto& parser_react = parser_m_react.at(reaction_m);
          parser_react = make_parser();
          add_interface_vars(parser_react, mem_react_expr);
          parser_react->set_expression(mem_react_expr);
          parser_react->compile();

          auto& parser_diff = parser_m_diff.at(reaction_m);
          parser_diff = make_parser();
          parser_diff->set_expression(mem_diff_expr);
          parser_diff->compile();
          mem_self_map.insert(reaction_m);

          auto& reaction_jac_config = reaction_config.sub("jacobian");

          // Do mrmbrane/membrane jacobian
          for (std::size_t reaction_mm = 0; reaction_mm != space_m.degree(); ++reaction_mm) {
            auto space_reaction_mm = space_m.subSpace(Assembler::multiIndex(reaction_mm));
            std::string jac_name = fmt::format("d{}__d{}", space_reaction_m.name(), space_reaction_mm.name());
            std::string expr = reaction_jac_config.get(jac_name, "0");
            _logger.trace(2, "Setup ({}/{}) membrane reaction jacobian expression ({}): {}"_fmt, space_m.name(), space_m.name(), jac_name, expr);
            std::size_t jac_index = reaction_m * space_m.degree() + reaction_mm;
            auto& jac_parser = parser_m_self_jac.at(jac_index);
            jac_parser = make_parser();
            add_interface_vars(jac_parser, expr);
            jac_parser->set_expression(expr);
            jac_parser->compile();

            if (mem_react_expr.find(space_reaction_mm.name()) != std::string::npos)
              mem_self_map.insert(reaction_mm);
          }
          _membrane_self_jac_map[{ membrane_id, reaction_m }] = std::move(mem_self_map);

          // Do membrane/domain jacobian
          for (std::size_t reaction_mi = 0; reaction_mi != space_i.degree(); ++reaction_mi) {
            auto space_reaction_mi = space_i.subSpace(Assembler::multiIndex(reaction_mi));
            std::string jac_name = fmt::format("d{}__d{}", space_reaction_m.name(), space_reaction_mi.name());
            std::string expr = reaction_jac_config.get(jac_name, "0");
            _logger.trace(2, "Setup ({}/{}) membrane reaction jacobian expression ({}): {}"_fmt, space_m.name(), space_i.name(), jac_name, expr);
            std::size_t jac_index = reaction_m * space_i.degree() + reaction_mi;
            auto& jac_parser = parser_m_cross_jac.at(jac_index);
            jac_parser = make_parser();
            add_interface_vars(jac_parser, expr);
            jac_parser->set_expression(expr);
            jac_parser->compile();

            if (mem_react_expr.find(space_reaction_mi.name()) != std::string::npos)
              mem_cross_map.insert(reaction_mi);
          }
          _membrane_cross_jac_map[{ membrane_id, reaction_m, domain_i}] = std::move(mem_cross_map);
        }
#endif

        auto boundary_outflow = fmt::format("{}.boundary.{}.outflow", space_i.name(), space_o.name());
        if (not _config.hasSub(boundary_outflow))
          continue;
        auto& outflow_config = _config.sub(boundary_outflow);

        auto& parser_o = _outflow_parser[{ domain_i, domain_o }];
        parser_o.resize(space_i.degree());
        
        auto& parser_o_self_jac = _outflow_self_jac_parser[{ domain_i, domain_o }];
        parser_o_self_jac.resize(space_i.degree() * space_i.degree());
        
        auto& parser_o_cross_jac = _outflow_cross_jac_parser[{ domain_i, domain_o }];
        parser_o_cross_jac.resize(space_i.degree() * space_o.degree());

#if DUNE_COPASI_HAVE_MEMBRANE_SPACE
        auto& parser_o_mem_jac = _outflow_mem_jac_parser[{ domain_i, domain_o }];
        parser_o_mem_jac.resize(space_i.degree() * space_m.degree());
#endif

        // Do jacobian
        for (std::size_t outflow_i = 0; outflow_i != space_i.degree(); ++outflow_i) {
          auto space_outflow_i = space_i.subSpace(Assembler::multiIndex(outflow_i));

          std::string out_expr = outflow_config.template get<std::string>(space_outflow_i.name());
          _logger.trace(2, "Setup expression ({}): {}"_fmt, space_outflow_i.name(), out_expr);

          std::set<std::size_t> out_self_map, out_cross_map;

          auto& parser = parser_o.at(outflow_i);
          parser = make_parser();
          add_interface_vars(parser, out_expr);
          parser->set_expression(out_expr);
          parser->compile();

          auto& outflow_jac_config = outflow_config.sub("jacobian");

          // Do self jacobian
          for (std::size_t outflow_ii = 0; outflow_ii != space_i.degree(); ++outflow_ii) {
            auto space_outflow_ii = space_i.subSpace(Assembler::multiIndex(outflow_ii));
            std::string jac_name = fmt::format("d{}__d{}", space_outflow_i.name(), space_outflow_ii.name());
            std::string expr = outflow_jac_config.get(jac_name, "0");
            _logger.trace(2, "Setup ({}/{}) outflow jacobian expression ({}): {}"_fmt, space_i.name(), space_i.name(), jac_name, expr);
            std::size_t jac_index = outflow_i * space_i.degree() + outflow_ii;
            auto& jac_parser = parser_o_self_jac.at(jac_index);
            jac_parser = make_parser();
            add_interface_vars(jac_parser, expr);
            jac_parser->set_expression(expr);
            jac_parser->compile();
            if (out_expr.find(space_outflow_ii.name()) != std::string::npos)
              out_self_map.insert(outflow_ii);
          }

          // Do cross jacobian
          for (std::size_t outflow_io = 0; outflow_io != space_o.degree(); ++outflow_io) {
            auto space_outflow_io = space_o.subSpace(Assembler::multiIndex(outflow_io));
            std::string jac_name = fmt::format("d{}__d{}", space_outflow_i.name(), space_outflow_io.name());
            std::string expr = outflow_jac_config.get(jac_name, "0");
            _logger.trace(2, "Setup ({}/{}) outflow jacobian expression ({}): {}"_fmt, space_i.name(), space_o.name(), jac_name, expr);
            std::size_t jac_index = outflow_i * space_o.degree() + outflow_io;
            auto& jac_parser = parser_o_cross_jac.at(jac_index);
            jac_parser = make_parser();
            add_interface_vars(jac_parser, expr);
            jac_parser->set_expression(expr);
            jac_parser->compile();
            if (out_expr.find(space_outflow_io.name()) != std::string::npos)
              out_cross_map.insert(outflow_io);
          }
          _outflow_jac_map[{ domain_i, domain_o, outflow_i }] = {std::move(out_self_map), std::move(out_cross_map)};

#if DUNE_COPASI_HAVE_MEMBRANE_SPACE
          // Do membrane jacobian
          std::set<std::size_t> out_memb_map;
          for (std::size_t outflow_m = 0; outflow_m != space_m.degree(); ++outflow_m) {
            auto space_outflow_m = space_m.subSpace(Assembler::multiIndex(outflow_m));
            std::string jac_name = fmt::format("d{}__d{}", space_outflow_i.name(), space_outflow_m.name());
            std::string expr = outflow_jac_config.get(jac_name, "0");
            _logger.trace(2, "Setup ({}/{}) outflow jacobian expression ({}): {}"_fmt, space_i.name(), space_m.name(), jac_name, expr);
            std::size_t jac_index = outflow_i * space_m.degree() + outflow_m;
            auto& jac_parser = parser_o_mem_jac.at(jac_index);
            jac_parser = make_parser();
            add_interface_vars(jac_parser, expr);
            jac_parser->set_expression(expr);
            jac_parser->compile();
            if (out_expr.find(space_outflow_m.name()) != std::string::npos)
              out_memb_map.insert(outflow_m);
          }
          _out_membrane_jac_map[{ domain_i, domain_o, outflow_i }] = std::move(out_memb_map);
#endif
        }
      }
    }

    {
      // (debug info) log interface pattern
      _logger.debug("Interface jacobian pattern:"_fmt);
      auto cache_interface = std::array{ unused_domain, unused_domain };
      for (auto [key, out_sets] : _outflow_jac_map) {
        auto [domain_i, domain_o, comp_i] = key;
        const auto& [out_self_set, out_cross_set] = out_sets;
        auto space_i = compartments_space.subSpace(Assembler::multiIndex(domain_i));
        auto space_o = compartments_space.subSpace(Assembler::multiIndex(domain_o));
        if (cache_interface != std::array{domain_i, domain_o})
          _logger.debug(2,
                        "{} -> {}"_fmt,
                        space_i.name(),
                        space_o.name());
        cache_interface = {domain_i, domain_o};
        auto space_comp_i = space_i.subSpace(Assembler::multiIndex(comp_i));
        for (auto self : out_self_set) {
          auto space_comp_self = space_i.subSpace(Assembler::multiIndex(self));
          _logger.debug(4,
                        "{} -> {}"_fmt,
                        space_comp_i.name(),
                        space_comp_self.name());
        }
        for (auto cross : out_cross_set) {
          auto space_comp_cross = space_o.subSpace(Assembler::multiIndex(cross));
          _logger.debug(4,
                        "{} -> {}"_fmt,
                        space_comp_i.name(),
                        space_comp_cross.name());
        }
      }
    }
  }

  template<class EG>
  auto subDomain(const EG& entity) const
  {
    auto domain_set = _space.entitySet().indexSet().subDomains(entity);
    assert(domain_set.size() == 1);
    return *(domain_set.begin());
  }

  template<class IG>
  bool skip_intersection(const IG& ig) const
  {
    auto domain_i = subDomain(ig.inside());

    bool skip = true;
    if constexpr (DomainLOP::doSkipIntersection)
      skip &= _lops[domain_i].skip_intersection(ig);

    if (ig.neighbor()) {
      auto domain_o = subDomain(ig.outside());
      skip &= (domain_i == domain_o);
    }
    return skip;
  }

  /**
   * @copydoc LocalOperatorDiffusionReactionCG::pattern_volume
   */
  template<class Entity, typename LFSU, typename LFSV, typename LocalPattern>
  void pattern_volume(const Entity& entity, const LFSU& lfsu,
                      const LFSV& lfsv,
                      LocalPattern& pattern) const
  {
    auto domain = subDomain(entity);
    auto domain_path = join(compartments_path, Assembler::multiIndex(domain));
    _lops[domain].pattern_volume(entity, child(lfsu, domain_path), child(lfsv, domain_path), pattern);
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
   */
  void interface_pattern_skeleton(const auto& ig,
                        std::size_t domain_i, std::size_t domain_o,
                        const auto& lfsu_i,
                        const auto& lfsv_i,
                        const auto& lfsu_o,
                        const auto& lfsv_o,
                        auto& pattern_ii,
                        auto& pattern_io,
                        auto& pattern_oi,
                        auto& pattern_oo) const
  {
    auto domain_i_path = join(compartments_path, Assembler::multiIndex(domain_i));
    auto domain_o_path = join(compartments_path, Assembler::multiIndex(domain_o));
    const auto& lfsu_di = child(lfsu_i, domain_i_path);
    const auto& lfsv_di = child(lfsv_i, domain_i_path);
    const auto& lfsu_do = child(lfsu_o, domain_o_path);
    const auto& lfsv_do = child(lfsv_o, domain_o_path);

#if DUNE_COPASI_HAVE_MEMBRANE_SPACE
    auto membrane_id_it = _membrane2id.find({domain_i, domain_o});
    assert((membrane_id_it != _membrane2id.end()) && "Membrane ID is missing");
    auto membrane_id = membrane_id_it->second;

    auto membrane_path = Assembler::multiIndex(Indices::_1, membrane_id);
    const auto& lfsu_f = child(lfsu_i, membrane_path);
    const auto& lfsv_f = child(lfsv_i, membrane_path);
    auto face_id = ig.indexInInside();
#endif

    for (std::size_t comp_i = 0; comp_i < lfsu_di.degree(); ++comp_i) {
      auto& lfsv_ci = lfsv_di.child(comp_i);
      auto oit = _outflow_jac_map.find({ domain_i, domain_o, comp_i });
      if (oit != _outflow_jac_map.end()) {
        const auto& [self_map, cross_map] = oit->second;
        for (auto&& comp_ii : self_map) {
          const auto& lfsu_cii = lfsu_di.child(comp_ii);
          for (std::size_t dof_i = 0; dof_i != lfsv_ci.size(); ++dof_i)
            for (std::size_t dof_ii = 0; dof_ii != lfsu_cii.size(); ++dof_ii)
              pattern_ii.addLink(lfsv_ci, dof_i, lfsu_cii, dof_ii);
        }
        for (auto&& comp_io : cross_map) {
          const auto& lfsu_co = lfsu_do.child(comp_io);
          for (std::size_t dof_i = 0; dof_i != lfsv_ci.size(); ++dof_i)
            for (std::size_t dof_o = 0; dof_o != lfsu_co.size(); ++dof_o)
              pattern_io.addLink(lfsv_ci, dof_i, lfsu_co, dof_o);
        }
      }
#if DUNE_COPASI_HAVE_MEMBRANE_SPACE
      auto mit = _out_membrane_jac_map.find({ domain_i, domain_o, comp_i });
      if (mit != _out_membrane_jac_map.end()) {
        const auto& memb_map = mit->second;
        for (auto&& comp_if : memb_map) {
          const auto& lfsu_cf = lfsu_f.child(comp_if).child(face_id);
          for (std::size_t dof_i = 0; dof_i != lfsv_ci.size(); ++dof_i)
            for (std::size_t dof_f = 0; dof_f != lfsu_cf.size(); ++dof_f)
              pattern_ii.addLink(lfsv_ci, dof_i, lfsu_cf, dof_f);
        }
      }
#endif
    }

    for (std::size_t comp_o = 0; comp_o != lfsu_do.degree(); ++comp_o) {
      auto& lfsv_co = lfsv_do.child(comp_o);
      auto oit = _outflow_jac_map.find({ domain_o, domain_i, comp_o });
      if (oit != _outflow_jac_map.end()) {
        const auto& [self_map, cross_map] = oit->second;
        for (auto&& comp_oo : self_map) {
          const auto& lfsu_coo = lfsu_do.child(comp_oo);
          for (std::size_t dof_o = 0; dof_o != lfsv_co.size(); ++dof_o)
            for (std::size_t dof_oo = 0; dof_oo != lfsu_coo.size(); ++dof_oo)
              pattern_oo.addLink(lfsv_co, dof_o, lfsu_coo, dof_oo);
        }
        for (auto&& comp_oi : cross_map) {
          const auto& lfsu_ci = lfsu_di.child(comp_oi);
          for (std::size_t dof_o = 0; dof_o != lfsv_co.size(); ++dof_o)
            for (std::size_t dof_i = 0; dof_i != lfsu_ci.size(); ++dof_i)
              pattern_oi.addLink(lfsv_co, dof_o, lfsu_ci, dof_i);
        }
      }
#if DUNE_COPASI_HAVE_MEMBRANE_SPACE
      auto mit = _out_membrane_jac_map.find({ domain_o, domain_i, comp_o });
      if (mit != _out_membrane_jac_map.end()) {
        const auto& memb_map = mit->second;
        for (auto&& comp_of : memb_map) {
          const auto& lfsu_cf = lfsu_f.child(comp_of).child(face_id);
          for (std::size_t dof_o = 0; dof_o != lfsv_co.size(); ++dof_o)
            for (std::size_t dof_f = 0; dof_f != lfsu_cf.size(); ++dof_f)
              pattern_oi.addLink(lfsv_co, dof_o, lfsu_cf, dof_f);
        }
      }
#endif
    }



#if DUNE_COPASI_HAVE_MEMBRANE_SPACE
    for (std::size_t comp_f = 0; comp_f < lfsu_f.degree(); ++comp_f) {
      auto& lfsv_cf = lfsv_f.child(comp_f).child(face_id);
      auto mit = _membrane_self_jac_map.find({ membrane_id, comp_f });
      if (mit != _membrane_self_jac_map.end()) {
        const auto& self_map = mit->second;
        for (auto&& comp_ff : self_map) {
          const auto& lfsu_cff = lfsu_f.child(comp_ff).child(face_id);
          for (std::size_t dof_f = 0; dof_f != lfsv_cf.size(); ++dof_f)
            for (std::size_t dof_ff = 0; dof_ff != lfsu_cff.size(); ++dof_ff)
              pattern_ii.addLink(lfsv_cf, dof_f, lfsu_cff, dof_ff);
        }
      }
      auto cit_i = _membrane_cross_jac_map.find({ membrane_id, comp_f, domain_i });
      if (cit_i != _membrane_cross_jac_map.end()) {
        const auto& cross_map = cit_i->second;
        for (auto&& comp_fi : cross_map) {
          const auto& lfsu_cfi = lfsu_di.child(comp_fi);
          for (std::size_t dof_f = 0; dof_f != lfsv_cf.size(); ++dof_f)
            for (std::size_t dof_fi = 0; dof_fi != lfsu_cfi.size(); ++dof_fi)
              pattern_ii.addLink(lfsv_cf, dof_f, lfsu_cfi, dof_fi);
        }
      }
      auto cit_o = _membrane_cross_jac_map.find({ membrane_id, comp_f, domain_o });
      if (cit_o != _membrane_cross_jac_map.end()) {
        const auto& cross_map = cit_o->second;
        for (auto&& comp_fo : cross_map) {
          const auto& lfsu_cfo = lfsu_do.child(comp_fo);
          for (std::size_t dof_f = 0; dof_f != lfsv_cf.size(); ++dof_f)
            for (std::size_t dof_fo = 0; dof_fo != lfsu_cfo.size(); ++dof_fo)
              pattern_io.addLink(lfsv_cf, dof_f, lfsu_cfo, dof_fo);
        }
      }
    }
#endif
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
  template<class IG, typename LFSU, typename LFSV, typename LocalPattern>
  void pattern_skeleton(const IG& ig,
                        const LFSU& lfsu_i,
                        const LFSV& lfsv_i,
                        const LFSU& lfsu_o,
                        const LFSV& lfsv_o,
                        LocalPattern& pattern_ii,
                        LocalPattern& pattern_io,
                        LocalPattern& pattern_oi,
                        LocalPattern& pattern_oo) const
  {
    auto domain_i = subDomain(ig.inside());
    auto domain_o = subDomain(ig.outside());

    assert(lfsu_i.degree() == lfsv_i.degree());
    assert(lfsu_o.degree() == lfsv_o.degree());

    if (domain_i != domain_o)
      interface_pattern_skeleton(ig,domain_i,domain_o,lfsu_i,lfsv_i,lfsu_o,lfsv_o,pattern_ii,pattern_io,pattern_oi,pattern_oo);
    else if constexpr (DomainLOP::doPatternSkeleton) {
      auto domain_i_path = join(compartments_path, Assembler::multiIndex(domain_i));
      auto domain_o_path = join(compartments_path, Assembler::multiIndex(domain_o));
      _lops[domain_i].pattern_skeleton(ig,
                          child(lfsu_i, domain_i_path),
                          child(lfsv_i, domain_i_path),
                          child(lfsu_o, domain_o_path),
                          child(lfsv_o, domain_o_path),
                          pattern_ii,pattern_io,pattern_oi,pattern_oo);
    }
  }

  /**
   * @brief      Sets the time.
   *
   * @param[in]  t     The new time
   */
  void setTime(double t)
  {
    Dune::PDELab::InstationaryLocalOperatorDefaultMethods<double>::setTime(t);
    for (auto& lp : _lops)
      lp.setTime(t);
  }

  /**
   * @copydoc LocalOperatorDiffusionReactionCG::jacobian_apply_volume
   * @details    This particular operator does a jacobian apply volume for the
   *             LocalOperatorDiffusionReactionCG corresponding to incoming
   * entity
   */
  template<typename EG, typename LFSU, typename X, typename LFSV, typename R>
  void jacobian_apply_volume(const EG& entity,
                             const LFSU& lfsu,
                             const X& x,
                             const X& z,
                             const LFSV& lfsv,
                             R& r)
  {
    auto domain = subDomain(entity);
    auto domain_path = join(compartments_path, Assembler::multiIndex(domain));
    const auto& sub_lfsu = child(lfsu, domain_path);
    const auto& sub_lfsv = child(lfsv, domain_path);
    _lops[domain].jacobian_apply_volume(entity, sub_lfsu, x, z, sub_lfsv, r);
  }

  /**
   * @copydoc LocalOperatorDiffusionReactionCG::jacobian_apply_volume
   * @details    This particular operator does a jacobian apply volume for the
   *             LocalOperatorDiffusionReactionCG corresponding to incoming
   * entity
   */
  template<typename EG, typename LFSU, typename X, typename LFSV, typename R>
  void jacobian_apply_volume(const EG& entity,
                             const LFSU& lfsu,
                             const X& x,
                             const LFSV& lfsv,
                             R& r)
  {
    auto domain = subDomain(entity);
    auto domain_path = join(compartments_path, Assembler::multiIndex(domain));
    const auto& sub_lfsu = child(lfsu, domain_path);
    const auto& sub_lfsv = child(lfsv, domain_path);
    _lops[domain].jacobian_apply_volume(entity, sub_lfsu, x, sub_lfsv, r);
  }

  /**
   * @copydoc LocalOperatorDiffusionReactionCG::jacobian_volume
   * @details    This particular operator does a jacobian volume for the
   *             LocalOperatorDiffusionReactionCG corresponding to incoming
   * entity
   */
  template<typename EG, typename LFSU, typename X, typename LFSV, typename M>
  void jacobian_volume(const EG& entity,
                       const LFSU& lfsu,
                       const X& x,
                       const LFSV& lfsv,
                       M& mat)
  {
    auto domain = subDomain(entity);
    auto domain_path = join(compartments_path, Assembler::multiIndex(domain));
    const auto& sub_lfsu = child(lfsu, domain_path);
    const auto& sub_lfsv = child(lfsv, domain_path);
    _lops[domain].jacobian_volume(entity, sub_lfsu, x, sub_lfsv, mat);
  }

  /**
   * @copydoc LocalOperatorDiffusionReactionCG::alpha_volume
   * @details    This particular operator does a alpha volume for the
   *             LocalOperatorDiffusionReactionCG corresponding to incoming
   * entity
   */
  template<typename EG, typename LFSU, typename X, typename LFSV, typename R>
  void alpha_volume(const EG& entity,
                    const LFSU& lfsu,
                    const X& x,
                    const LFSV& lfsv,
                    R& r)
  {
    auto domain = subDomain(entity);
    auto domain_path = join(compartments_path, Assembler::multiIndex(domain));
    const auto& sub_lfsu = child(lfsu, domain_path);
    const auto& sub_lfsv = child(lfsv, domain_path);
    _lops[domain].alpha_volume(entity, sub_lfsu, x, sub_lfsv, r);
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
                      R& r_o)
  {
    const auto& entity_i = ig.inside();
    const auto& entity_o = ig.outside();

    assert(lfsu_i.degree() == lfsv_i.degree());
    assert(lfsu_o.degree() == lfsv_o.degree());

    auto domain_i = subDomain(entity_i);
    auto domain_o = subDomain(entity_o);

    if (domain_i != domain_o)
      interface_alpha_skeleton(domain_i,domain_o,ig,lfsu_i,x_i,lfsv_i,lfsu_o,x_o,lfsv_o,r_i,r_o);
    else if constexpr (DomainLocalOperator::doAlphaSkeleton) {
      auto domain_i_path = join(compartments_path, Assembler::multiIndex(domain_i));
      auto domain_o_path = join(compartments_path, Assembler::multiIndex(domain_o));
      _lops[domain_i].alpha_skeleton(ig,
                child(lfsu_i, domain_i_path),
                x_i,
                child(lfsv_i, domain_i_path),
                child(lfsu_o, domain_o_path),
                x_o,
                child(lfsv_o, domain_o_path),
                r_i, r_o);
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
  void jacobian_apply_skeleton(const IG& ig,
                      const LFSU& lfsu_i,
                      const X& x_i,
                      const X& z_i,
                      const LFSV& lfsv_i,
                      const LFSU& lfsu_o,
                      const X& x_o,
                      const X& z_o,
                      const LFSV& lfsv_o,
                      R& r_i,
                      R& r_o)
  {
    const auto& entity_i = ig.inside();
    const auto& entity_o = ig.outside();

    assert(_lops.size() >= lfsu_i.degree());
    assert(_lops.size() >= lfsu_o.degree());
    assert(lfsu_i.degree() == lfsv_i.degree());
    assert(lfsu_o.degree() == lfsv_o.degree());

    auto domain_i = subDomain(entity_i);
    auto domain_o = subDomain(entity_o);

    if (domain_i != domain_o)
      interface_jacobian_apply_skeleton(domain_i,domain_o,ig,lfsu_i,x_i,z_i,lfsv_i,lfsu_o,x_o,z_o,lfsv_o,r_i,r_o);
    else if constexpr (DomainLocalOperator::doAlphaSkeleton) {
      auto domain_i_path = join(compartments_path, Assembler::multiIndex(domain_i));
      auto domain_o_path = join(compartments_path, Assembler::multiIndex(domain_o));
      _lops[domain_i].jacobian_apply_skeleton(ig,
                child(lfsu_i, domain_i_path),
                x_i, z_i,
                child(lfsv_i, domain_i_path),
                child(lfsu_o, domain_o_path),
                x_o, z_o,
                child(lfsv_o, domain_o_path),
                r_i, r_o);
    }
  }

  template<class R, class X>
  struct PseudoJacobian {
    void accumulate(const auto& ltest, auto test_dof, const auto& ltrial, auto trial_dof, auto value) {
      _r.accumulate(ltest, test_dof, _z(ltrial, trial_dof) * value);
    }

    R& _r;
    const X& _z;
  };


  template<typename IG, typename LFSU, typename X, typename LFSV, typename R>
  void interface_jacobian_apply_skeleton(std::size_t domain_i, std::size_t domain_o,
                      const IG& ig,
                      const LFSU& lfsu_i,
                      const X& x_i,
                      const X& z_i,
                      const LFSV& lfsv_i,
                      const LFSU& lfsu_o,
                      const X& x_o,
                      const X& z_o,
                      const LFSV& lfsv_o,
                      R& r_i,
                      R& r_o)
  {
    PseudoJacobian<R, X> mat_ii{r_i, z_i}, mat_io{r_i, z_o}, mat_oi{r_o, z_i}, mat_oo{r_o, z_o};
    interface_jacobian_skeleton(domain_i, domain_o, ig, lfsu_i, x_i, lfsv_i, lfsu_o, x_o, lfsv_o, mat_ii, mat_io, mat_oi, mat_oo);
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
                      const LFSU& lfsu_i,
                      const X& x_i,
                      const LFSV& lfsv_i,
                      const LFSU& lfsu_o,
                      const X& x_o,
                      const LFSV& lfsv_o,
                      R& r_i,
                      R& r_o)
  {
    auto domain_i_path = join(compartments_path, Assembler::multiIndex(domain_i));
    auto domain_o_path = join(compartments_path, Assembler::multiIndex(domain_o));
    const auto& lfsu_di = child(lfsu_i, domain_i_path);
    const auto& lfsv_di = child(lfsv_i, domain_i_path);
    const auto& lfsu_do = child(lfsu_o, domain_o_path);
    const auto& lfsv_do = child(lfsv_o, domain_o_path);

#if DUNE_COPASI_HAVE_MEMBRANE_SPACE
    assert(ig.conforming());
    auto membrane_id_it = _membrane2id.find({domain_i, domain_o});
    assert((membrane_id_it != _membrane2id.end()) && "Membrane ID is missing");
    auto membrane_id = membrane_id_it->second;

    auto membrane_path = Assembler::multiIndex(Indices::_1, membrane_id);
    const auto& lfsu_f = child(lfsu_i, membrane_path);
    const auto& lfsv_f = child(lfsv_i, membrane_path);
    auto face_id = ig.indexInInside();

    MembraneLocalBasis const * local_basis_f = nullptr;
    if (lfsv_f.degree() != 0)
      local_basis_f = &lfsv_f.child(0).child(face_id).finiteElement().localBasis();

    if (local_basis_f)
      _gradphi_f.resize(local_basis_f->size());
#endif

    const auto& entity_f = ig;
    const auto& entity_i = ig.inside();
    const auto& entity_o = ig.outside();

    auto geo_f = entity_f.geometry();
    auto geo_i = entity_i.geometry();
    auto geo_o = entity_o.geometry();
    auto geo_in_i = entity_f.geometryInInside();
    auto geo_in_o = entity_f.geometryInOutside();

    assert(_components[domain_i].size() == lfsu_di.degree());
    assert(_components[domain_o].size() == lfsu_do.degree());

    BulkLocalBasis const * local_basis_i = nullptr;
    BulkLocalBasis const * local_basis_o = nullptr;

    if (lfsu_di.degree() != 0)
      local_basis_i = &lfsu_di.child(0).finiteElement().localBasis();
    if (lfsu_do.degree() != 0)
      local_basis_o = &lfsu_do.child(0).finiteElement().localBasis();


    if (local_basis_i)
      _gradphi_i.resize(local_basis_i->size());
    if (local_basis_o)
      _gradphi_o.resize(local_basis_o->size());

    auto normal_f = ig.centerUnitOuterNormal();

    auto dn_i = (geo_f.center() - geo_i.center()).two_norm();
    auto dn_o = (geo_f.center() - geo_o.center()).two_norm();

    for (const auto& it : quadratureRule(geo_f, 3)) {
      const auto& position_f = it.position();

      // position of quadrature point in local coordinates of reference elements
      const auto position_i = geo_in_i.global(position_f);
      const auto position_o = geo_in_o.global(position_f);

      // evaluate basis functions
      if (local_basis_i)
        local_basis_i->evaluateFunction(position_i, _phiu_i);
      if (local_basis_o)
        local_basis_o->evaluateFunction(position_o, _phiu_o);

      // evaluate concentrations at quadrature point
      std::fill(_components[domain_i].begin(),_components[domain_i].end(),0.);
      std::fill(_components[domain_o].begin(),_components[domain_o].end(),0.);
      for (std::size_t comp = 0; comp != lfsu_di.degree(); ++comp)
        for (std::size_t dof = 0; dof != local_basis_i->size(); ++dof)
          _components[domain_i][comp] += x_i(lfsu_di.child(comp), dof) * _phiu_i[dof];

      for (std::size_t comp = 0; comp != lfsu_do.degree(); ++comp)
        for (std::size_t dof = 0; dof != local_basis_o->size(); ++dof)
          _components[domain_o][comp] += x_o(lfsu_do.child(comp), dof) * _phiu_o[dof];

#if DUNE_COPASI_HAVE_MEMBRANE_SPACE
      if (local_basis_f) {
        const auto& geo_if = entity_i.template subEntity<1>(face_id).geometry();
        auto position_if = geo_if.local(geo_f.global(position_f));
        local_basis_f->evaluateFunction(position_if, _phiu_f);

        std::fill(_u_membrane[membrane_id].begin(),_u_membrane[membrane_id].end(),0.);
        for (std::size_t comp = 0; comp != lfsu_f.degree(); ++comp)
          for (std::size_t dof = 0; dof != local_basis_f->size(); ++dof)
            _u_membrane[membrane_id][comp] += x_i(lfsu_f.child(comp).child(face_id), dof) * _phiu_f[dof];

        local_basis_f->evaluateJacobian(position_if, _jacphi_f);
        auto jac_inv_i = geo_if.jacobianInverse(position_if);
        for (std::size_t dof = 0; dof != _jacphi_f.size(); ++dof)
          _gradphi_f[dof] = (_jacphi_f[dof] * jac_inv_i)[0];
      }
#endif

      std::fill(_components_dn[domain_i].begin(),_components_dn[domain_i].end(),0.);
      std::fill(_components_dn[domain_o].begin(),_components_dn[domain_o].end(),0.);

      std::fill(_gradphi_i.begin(), _gradphi_i.end(), 0.);
      std::fill(_gradphi_o.begin(), _gradphi_o.end(), 0.);

      if (local_basis_i) {
        if (local_basis_i->order() == 0) {
          for (std::size_t comp_i = 0; comp_i != lfsu_di.degree(); ++comp_i) {
            _components_dn[domain_i][comp_i] = std::numeric_limits<double>::quiet_NaN();
            auto comp_o_it = _component_offset.find({domain_i,domain_o,comp_i});
            if (comp_o_it != _component_offset.end()) {
              std::size_t comp_o = comp_o_it->second;
              _components_dn[domain_i][comp_i] = (_components[domain_o][comp_o] - _components[domain_i][comp_i])/dn_i;
            }
          }
        } else {
          local_basis_i->evaluateJacobian(position_i,_jacphi_i);
          auto jac_i = geo_i.jacobianInverse(position_i);
          for (std::size_t dof = 0; dof != _gradphi_i.size(); ++dof)
            _gradphi_i[dof] = (_jacphi_i[dof]*jac_i)[0];
          for (std::size_t comp_i = 0; comp_i != lfsu_di.degree(); ++comp_i)
            for (std::size_t dof = 0; dof != _gradphi_i.size(); ++dof)
              _components_dn[domain_i][comp_i] += (normal_f * _gradphi_i[dof]) * x_i(lfsu_di.child(comp_i), dof);
        }
      }

      if (local_basis_o) {
        if (local_basis_o->order() == 0) {
          for (std::size_t comp_o = 0; comp_o != lfsu_do.degree(); ++comp_o) {
            _components_dn[domain_o][comp_o] = std::numeric_limits<double>::quiet_NaN();
            auto comp_i_it = _component_offset.find({domain_o,domain_i,comp_o});
            if (comp_i_it != _component_offset.end()) {
              std::size_t comp_i = comp_i_it->second;
              _components_dn[domain_o][comp_o] = (_components[domain_i][comp_i] - _components[domain_o][comp_o])/dn_o;
            }
          }
        } else {
          local_basis_o->evaluateJacobian(position_o,_jacphi_o);
          auto jac_o = geo_o.jacobianInverse(position_o);
          for (std::size_t dof = 0; dof != _gradphi_o.size(); ++dof)
            _gradphi_o[dof] = (_jacphi_o[dof]*jac_o)[0];
          for (std::size_t comp_o = 0; comp_o != lfsu_do.degree(); ++comp_o)
            for (std::size_t dof = 0; dof != _gradphi_i.size(); ++dof)
              _components_dn[domain_o][comp_o] -= (normal_f * _gradphi_o[dof]) * x_o(lfsu_do.child(comp_o), dof);
        }
      }

      // integration factor
      auto factor = it.weight() * geo_f.integrationElement(position_f);

      if (lfsu_di.degree() != 0) {
        auto outflow_parser_io_it = _outflow_parser.find({domain_i,domain_o});
        const auto& parser_io = outflow_parser_io_it->second;
        for (std::size_t comp_i = 0; comp_i != lfsv_di.degree(); ++comp_i) {
          double value = parser_io[comp_i]->eval();
          for (std::size_t dof = 0; dof != lfsv_di.child(comp_i).size(); ++dof)
            r_i.accumulate(lfsv_di.child(comp_i), dof, factor * value * _phiu_i[dof]);
        }
      }

      if (lfsu_do.degree() != 0) {
        auto outflow_parser_oi_it = _outflow_parser.find({domain_o,domain_i});
        if (outflow_parser_oi_it != _outflow_parser.end()) {
          for (std::size_t comp_o = 0; comp_o != lfsv_do.degree(); ++comp_o) {
            const auto& parser_oi = outflow_parser_oi_it->second;
            double value = parser_oi[comp_o]->eval();
            for (std::size_t dof = 0; dof != lfsv_do.child(comp_o).size(); ++dof)
              r_o.accumulate(lfsv_do.child(comp_o), dof, factor * value * _phiu_o[dof]);
          }
        }
      }

#if DUNE_COPASI_HAVE_MEMBRANE_SPACE
      const auto& react_parser_f = _reaction_parser[membrane_id];
      const auto& diff_parser_f = _diffusion_parser[membrane_id];
      for (std::size_t comp_f = 0; comp_f != lfsv_f.degree(); ++comp_f) {
        const auto& lfsv_fc = lfsv_f.child(comp_f).child(face_id);
        for (std::size_t dof = 0; dof != lfsv_fc.size(); ++dof) {
          double reaction = react_parser_f[comp_f]->eval();
          auto diffusion = diff_parser_f[comp_f]->eval();

          FieldVector<MembraneRangeField, dim> graduh(.0);
          for (std::size_t dof = 0; dof != _gradphi_f.size(); ++dof)
            graduh += _gradphi_f[dof] * x_i(lfsv_fc, dof);

          r_i.accumulate(lfsv_fc, dof, 00.0000*(diffusion * dot(graduh, _gradphi_f[dof]) - reaction * _phiu_f[dof]) * factor);
        }
      }
#endif
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
                         J& mat_oo)
  {
    assert(_lops.size() >= lfsu_i.degree());
    assert(_lops.size() >= lfsu_o.degree());
    assert(lfsu_i.degree() == lfsv_i.degree());
    assert(lfsu_o.degree() == lfsv_o.degree());

    const auto& entity_i = ig.inside();
    const auto& entity_o = ig.outside();

    auto domain_i = subDomain(entity_i);
    auto domain_o = subDomain(entity_o);

    if (domain_i != domain_o)
      interface_jacobian_skeleton(domain_i,domain_o,ig,lfsu_i,x_i,lfsv_i,lfsu_o,x_o,lfsv_o,mat_ii,mat_io,mat_oi,mat_oo);
    else if constexpr (DomainLocalOperator::doAlphaSkeleton) {
      auto domain_i_path = join(compartments_path, Assembler::multiIndex(domain_i));
      auto domain_o_path = join(compartments_path, Assembler::multiIndex(domain_o));
      _lops[domain_i].jacobian_skeleton(ig,
                child(lfsu_i, domain_i_path),
                x_i,
                child(lfsv_i, domain_i_path),
                child(lfsu_o, domain_o_path),
                x_o,
                child(lfsv_o, domain_o_path),
                mat_ii, mat_io, mat_oi, mat_oo);
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
   */
  void interface_jacobian_skeleton(std::size_t domain_i, std::size_t domain_o,
                         const auto& ig,
                         const auto& lfsu_i,
                         const auto& x_i,
                         const auto& lfsv_i,
                         const auto& lfsu_o,
                         const auto& x_o,
                         const auto& lfsv_o,
                         auto& mat_ii,
                         auto& mat_io,
                         auto& mat_oi,
                         auto& mat_oo)
  {
    auto domain_i_path = join(compartments_path, Assembler::multiIndex(domain_i));
    auto domain_o_path = join(compartments_path, Assembler::multiIndex(domain_o));
    const auto& lfsu_di = child(lfsu_i, domain_i_path);
    const auto& lfsv_di = child(lfsv_i, domain_i_path);
    const auto& lfsu_do = child(lfsu_o, domain_o_path);
    const auto& lfsv_do = child(lfsv_o, domain_o_path);

    const auto& entity_f = ig;
    const auto& entity_i = ig.inside();
    const auto& entity_o = ig.outside();

    auto geo_f = entity_f.geometry();
    auto geo_i = entity_i.geometry();
    auto geo_o = entity_o.geometry();
    auto geo_in_i = entity_f.geometryInInside();
    auto geo_in_o = entity_f.geometryInOutside();

    BulkLocalBasis const * local_basis_i = nullptr;
    BulkLocalBasis const * local_basis_o = nullptr;

    if (lfsu_di.degree() != 0)
      local_basis_i = &lfsu_di.child(0).finiteElement().localBasis();
    if (lfsu_do.degree() != 0)
      local_basis_o = &lfsu_do.child(0).finiteElement().localBasis();

    if (local_basis_i)
      _gradphi_i.resize(local_basis_i->size());
    if (local_basis_o)
      _gradphi_o.resize(local_basis_o->size());

#if DUNE_COPASI_HAVE_MEMBRANE_SPACE
    assert(ig.conforming());
    auto membrane_id_it = _membrane2id.find({domain_i, domain_o});
    assert((membrane_id_it != _membrane2id.end()) && "Membrane ID is missing");
    auto membrane_id = membrane_id_it->second;

    auto membrane_path = Assembler::multiIndex(Indices::_1, membrane_id);
    const auto& lfsu_f = child(lfsu_i, membrane_path);
    const auto& lfsv_f = child(lfsv_i, membrane_path);
    auto face_id = ig.indexInInside();

    MembraneLocalBasis const * local_basis_f = nullptr;
    if (lfsv_f.degree() != 0)
      local_basis_f = &lfsv_f.child(0).child(face_id).finiteElement().localBasis();

    if (local_basis_f)
      _gradphi_f.resize(local_basis_f->size());
#endif

    auto normal_f = ig.centerUnitOuterNormal();

    auto dn_i = (geo_f.center() - geo_i.center()).two_norm();
    auto dn_o = (geo_f.center() - geo_o.center()).two_norm();

    for (const auto& it : quadratureRule(geo_f, 3)) {
      const auto& position_f = it.position();

      std::fill(_gradphi_i.begin(), _gradphi_i.end(), 0.);
      std::fill(_gradphi_o.begin(), _gradphi_o.end(), 0.);

      // position of quadrature point in local coordinates of elements
      const auto position_i = geo_in_i.global(position_f);
      const auto position_o = geo_in_o.global(position_f);

      // evaluate basis functions
      if (local_basis_i)
        local_basis_i->evaluateFunction(position_i, _phiu_i);
      if (local_basis_o)
        local_basis_o->evaluateFunction(position_o, _phiu_o);

      // evaluate concentrations at quadrature point
      std::fill(_components[domain_i].begin(),_components[domain_i].end(),0.);
      std::fill(_components[domain_o].begin(),_components[domain_o].end(),0.);

      for (std::size_t comp = 0; comp != lfsu_di.degree(); ++comp)
        for (std::size_t dof = 0; dof != local_basis_i->size(); ++dof)
          _components[domain_i][comp] += x_i(lfsu_di.child(comp), dof) * _phiu_i[dof];

      for (std::size_t comp = 0; comp != lfsu_do.degree(); ++comp)
        for (std::size_t dof = 0; dof != local_basis_o->size(); ++dof)
          _components[domain_o][comp] += x_o(lfsu_do.child(comp), dof) * _phiu_o[dof];

#if DUNE_COPASI_HAVE_MEMBRANE_SPACE
      if (local_basis_f) {
        const auto& geo_if = entity_i.template subEntity<1>(face_id).geometry();
        auto position_if = geo_if.local(geo_f.global(position_f));
        local_basis_f->evaluateFunction(position_if, _phiu_f);

        std::fill(_u_membrane[membrane_id].begin(),_u_membrane[membrane_id].end(),0.);
        for (std::size_t comp = 0; comp != lfsu_f.degree(); ++comp)
          for (std::size_t dof = 0; dof != local_basis_f->size(); ++dof)
            _u_membrane[membrane_id][comp] += x_i(lfsu_f.child(comp).child(face_id), dof) * _phiu_f[dof];

        local_basis_f->evaluateJacobian(position_if, _jacphi_f);
        auto jac_inv_i = geo_if.jacobianInverse(position_if);
        for (std::size_t dof = 0; dof != _jacphi_f.size(); ++dof)
          _gradphi_f[dof] = (_jacphi_f[dof] * jac_inv_i)[0];
      }
#endif

      std::fill(_components_dn[domain_i].begin(),_components_dn[domain_i].end(),0.);
      std::fill(_components_dn[domain_o].begin(),_components_dn[domain_o].end(),0.);

      if (local_basis_i) {
        if (local_basis_i->order() == 0) {
          for (std::size_t comp_i = 0; comp_i != lfsu_di.degree(); ++comp_i) {
            _components_dn[domain_i][comp_i] = std::numeric_limits<double>::quiet_NaN();
            auto comp_o_it = _component_offset.find({domain_i,domain_o,comp_i});
            if (comp_o_it != _component_offset.end()) {
              std::size_t comp_o = comp_o_it->second;
              _components_dn[domain_i][comp_i] = (_components[domain_o][comp_o] - _components[domain_i][comp_i])/dn_i;
            }
          }
        } else {
          local_basis_i->evaluateJacobian(position_i,_jacphi_i);
          auto jac_i = geo_i.jacobianInverse(position_i);
          for (std::size_t dof = 0; dof != _gradphi_i.size(); ++dof)
            _gradphi_i[dof] = (_jacphi_i[dof]*jac_i)[0];
          for (std::size_t comp_i = 0; comp_i != lfsu_di.degree(); ++comp_i)
            for (std::size_t dof = 0; dof != _gradphi_i.size(); ++dof)
              _components_dn[domain_i][comp_i] += x_i(lfsu_di.child(comp_i), dof) * (normal_f * _gradphi_i[dof]);
        }
      }

      if (local_basis_o) {
        if (local_basis_o->order() == 0) {
          for (std::size_t comp_o = 0; comp_o != lfsu_di.degree(); ++comp_o) {
            _components_dn[domain_o][comp_o] = std::numeric_limits<double>::quiet_NaN();
            auto comp_i_it = _component_offset.find({domain_o,domain_i,comp_o});
            if (comp_i_it != _component_offset.end()) {
              std::size_t comp_i = comp_i_it->second;
              _components_dn[domain_o][comp_o] = (_components[domain_i][comp_i] - _components[domain_o][comp_o])/dn_o;
            }
          }
        } else {
          local_basis_o->evaluateJacobian(position_o,_jacphi_o);
          auto jac_o = geo_o.jacobianInverse(position_o);
          for (std::size_t dof = 0; dof != _gradphi_o.size(); ++dof)
            _gradphi_o[dof] = (_jacphi_o[dof]*jac_o)[0];
          for (std::size_t comp_o = 0; comp_o != lfsu_do.degree(); ++comp_o)
            for (std::size_t dof = 0; dof != _gradphi_i.size(); ++dof)
              _components_dn[domain_o][comp_o] -= x_o(lfsu_do.child(comp_o), dof) * (normal_f * _gradphi_o[dof]);
        }
      }

      // integration factor
      auto factor = it.weight() * geo_f.integrationElement(position_f);

      const auto& parser_self_jac_ii = _outflow_self_jac_parser.find({domain_i,domain_o})->second;
      const auto& parser_self_jac_io = _outflow_cross_jac_parser.find({domain_i,domain_o})->second;
      const auto& parser_self_jac_oo = _outflow_self_jac_parser.find({domain_o,domain_i})->second;
      const auto& parser_self_jac_oi = _outflow_cross_jac_parser.find({domain_o,domain_i})->second;

#if DUNE_COPASI_HAVE_MEMBRANE_SPACE
      const auto& parser_self_jac_of = _outflow_mem_jac_parser.find({ domain_o, domain_i })->second;
      const auto& parser_self_jac_if = _outflow_mem_jac_parser.find({ domain_i, domain_o })->second;
#endif

      for (std::size_t comp_i = 0; comp_i != lfsu_di.degree(); ++comp_i) {
        const auto& lfsv_ci = lfsv_di.child(comp_i);
        auto oit = _outflow_jac_map.find({ domain_i, domain_o, comp_i });
        if (oit != _outflow_jac_map.end()) {
          const auto& [out_self_map, out_cross_map] = oit->second;
          for (auto comp_ii : out_self_map) {
            const auto& lfsu_cii = lfsu_di.child(comp_ii);
            auto jac_value = parser_self_jac_ii[comp_i*lfsv_di.degree() + comp_ii]->eval();
            for (std::size_t dof_i = 0; dof_i != lfsv_ci.size(); ++dof_i)
              for (std::size_t dof_ii = 0; dof_ii != lfsu_cii.size(); ++dof_ii)
                mat_ii.accumulate(lfsv_ci,  dof_i,
                                  lfsu_cii, dof_ii,
                                  factor * jac_value * _phiu_i[dof_ii]);
          }
          for (auto comp_io : out_cross_map) {
            auto jac_value = parser_self_jac_io[comp_i*lfsv_do.degree() + comp_io]->eval();
            const auto& lfsu_cio = lfsu_do.child(comp_io);
            for (std::size_t dof_i = 0; dof_i != lfsv_ci.size(); ++dof_i)
              for (std::size_t dof_io = 0; dof_io != lfsu_cio.size(); ++dof_io)
                mat_io.accumulate(lfsv_ci,  dof_i,
                                  lfsu_cio, dof_io,
                                  factor * jac_value * _phiu_o[dof_io]);
          }
        }
#if DUNE_COPASI_HAVE_MEMBRANE_SPACE
        auto mit = _out_membrane_jac_map.find({ domain_i, domain_o, comp_i });
        if (mit != _out_membrane_jac_map.end()) {
          const auto& out_memb_map = mit->second;
          for (auto comp_if : out_memb_map) {
            auto jac_value = parser_self_jac_if[comp_i*lfsv_di.degree() + comp_if]->eval();
            const auto& lfsu_cif = lfsu_f.child(comp_if).child(face_id);
            for (std::size_t dof_i = 0; dof_i != lfsv_di.child(comp_i).size(); ++dof_i)
              for (std::size_t dof_if = 0; dof_if != lfsu_cif.size(); ++dof_if)
                mat_ii.accumulate(lfsv_ci,  dof_i,
                                  lfsu_cif, dof_if,
                                  factor * jac_value * _phiu_f[dof_if]);
          }
        }
#endif
      }

      for (std::size_t comp_o = 0; comp_o != lfsu_do.degree(); ++comp_o) {
        const auto& lfsv_co = lfsv_do.child(comp_o);
        auto oit = _outflow_jac_map.find({ domain_o, domain_i, comp_o });
        if (oit != _outflow_jac_map.end()) {
          const auto& [out_self_map, out_cross_map] = oit->second;
          for (auto comp_oo : out_self_map) {
            const auto& lfsu_coo = lfsu_do.child(comp_oo);
            auto jac_value = parser_self_jac_oo[comp_o*lfsv_do.degree() + comp_oo]->eval();
            for (std::size_t dof_o = 0; dof_o != lfsv_co.size(); ++dof_o)
              for (std::size_t dof_oo = 0; dof_oo != lfsu_coo.size(); ++dof_oo)
                mat_oo.accumulate(lfsv_co,  dof_o,
                                  lfsu_coo, dof_oo,
                                  factor * jac_value * _phiu_o[dof_oo]);
          }
          for (auto comp_oi : out_cross_map) {
            const auto& lfsu_coi = lfsu_di.child(comp_oi);
            auto jac_value = parser_self_jac_oi[comp_o*lfsv_di.degree() + comp_oi]->eval();
            for (std::size_t dof_o = 0; dof_o != lfsv_co.size(); ++dof_o)
              for (std::size_t dof_oi = 0; dof_oi != lfsu_coi.size(); ++dof_oi)
                mat_oi.accumulate(lfsv_co,  dof_o,
                                  lfsu_coi, dof_oi,
                                  factor * jac_value * _phiu_i[dof_oi]);
          }
        }
#if DUNE_COPASI_HAVE_MEMBRANE_SPACE
        auto mit = _out_membrane_jac_map.find({ domain_o, domain_i, comp_o });
        if (mit != _out_membrane_jac_map.end()) {
          const auto& out_memb_map = mit->second;
          for (auto comp_of : out_memb_map) {
            auto jac_value = parser_self_jac_of[comp_o*lfsv_do.degree() + comp_of]->eval();
            const auto& lfsu_cof = lfsu_f.child(comp_of).child(face_id);
            for (std::size_t dof_o = 0; dof_o != lfsv_do.child(comp_o).size(); ++dof_o)
              for (std::size_t dof_of = 0; dof_of != lfsu_cof.size(); ++dof_of)
                mat_oi.accumulate(lfsv_co,  dof_o,
                                  lfsu_cof, dof_of,
                                  factor * jac_value * _phiu_f[dof_of]);
          }
        }
#endif
      }

#if DUNE_COPASI_HAVE_MEMBRANE_SPACE
      const auto& parser_self_jac_f = _reaction_self_jac_parser[membrane_id];
      const auto& parser_cross_jac_fi = _reaction_cross_jac_parser.find({membrane_id,domain_i})->second;
      const auto& parser_cross_jac_fo = _reaction_cross_jac_parser.find({membrane_id,domain_o})->second;
      const auto& diff_parser_f = _diffusion_parser[membrane_id];
      for (std::size_t comp_f = 0; comp_f != lfsu_f.degree(); ++comp_f) {
        const auto& lfsv_cf = lfsv_f.child(comp_f).child(face_id);
        const auto& lfsu_cf = lfsu_f.child(comp_f).child(face_id);

        // diffusion part
        auto diffusion = diff_parser_f[comp_f]->eval();
        for (std::size_t dof_i = 0; dof_i != lfsv_cf.size(); ++dof_i)
          for (std::size_t dof_j = 0; dof_j != lfsu_cf.size(); ++dof_j)
            mat_ii.accumulate(lfsv_cf, dof_i, lfsu_cf, dof_j, 0.0000*diffusion * dot(_gradphi_f[dof_i], _gradphi_f[dof_j]) * factor);

        // reaction part
        auto mit = _membrane_self_jac_map.find({ membrane_id, comp_f });
        if (mit != _membrane_self_jac_map.end()) {
          const auto& mem_self_map = mit->second;
          for (auto comp_ff : mem_self_map) {
            double jac_value = -1.0 * parser_self_jac_f[comp_f*lfsv_f.degree() + comp_ff]->eval();
            const auto& lfsu_cff = lfsu_f.child(comp_ff).child(face_id);
            for (std::size_t dof_f = 0; dof_f != lfsv_cf.size(); ++dof_f)
              for (std::size_t dof_ff = 0; dof_ff != lfsu_cff.size(); ++dof_ff)
                mat_ii.accumulate(lfsv_cf,  dof_f,
                                  lfsu_cff, dof_ff,
                                  0.0000*factor * jac_value * _phiu_f[dof_ff]);
          }
        }
        auto cit_i = _membrane_cross_jac_map.find({ membrane_id, comp_f, domain_i });
        if (cit_i != _membrane_cross_jac_map.end()) {
          const auto& mem_cross_map = cit_i->second;
          for (auto comp_fi : mem_cross_map) {
            double jac_value = -1.0 * parser_cross_jac_fi[comp_f*lfsv_di.degree() + comp_fi]->eval();
            const auto& lfsu_cfi = lfsu_di.child(comp_fi);
            for (std::size_t dof_f = 0; dof_f != lfsv_cf.size(); ++dof_f)
              for (std::size_t dof_fi = 0; dof_fi != lfsu_cfi.size(); ++dof_fi)
                mat_ii.accumulate(lfsv_cf,  dof_f,
                                  lfsu_cfi, dof_fi,
                                  0.0000*factor * jac_value * _phiu_i[dof_fi]);
          }
        }
        auto cit_o = _membrane_cross_jac_map.find({ membrane_id, comp_f, domain_o });
        if (cit_o != _membrane_cross_jac_map.end()) {
          const auto& mem_cross_map = cit_o->second;
          for (auto comp_fo : mem_cross_map) {
            double jac_value = -1.0 * parser_cross_jac_fo[comp_f*lfsv_do.degree() + comp_fo]->eval();
            const auto& lfsu_cfo = lfsu_do.child(comp_fo);
            for (std::size_t dof_f = 0; dof_f != lfsv_cf.size(); ++dof_f)
              for (std::size_t dof_fo = 0; dof_fo != lfsu_cfo.size(); ++dof_fo)
                mat_io.accumulate(lfsv_cf,  dof_f,
                                  lfsu_cfo, dof_fo,
                                  0.0000*factor * jac_value * _phiu_o[dof_fo]);
          }
        }
      }
#endif
    }
  }
};

/**
 * @brief      This class describes a temporal local operator for multi domain
 *             diffusion reaction systems.
 * @details    This class describre the operatrions for local integrals required
 *             for diffusion reaction system. This local operator maintains an
 *             individual local operator for every subdomain in the grid.
 *
 * @tparam     Space               A multidomain discrete function space
 * @tparam     DomainLocalOperator    Local operator of sub-domain
 */
template<Assembler::Concept::DiscreteFunctionSpace Space,
         class DomainLocalOperator>
class TemporalLocalOperatorMultiDomainDiffusionReaction
  : public Dune::PDELab::LocalOperatorDefaultFlags
  , public Dune::PDELab::InstationaryLocalOperatorDefaultMethods<double>
{
  using DomainLOP = DomainLocalOperator;

#if DUNE_COPASI_HAVE_MEMBRANE_SPACE
  static constexpr auto compartments_path = Assembler::multiIndex(Indices::_0);
#else
  static constexpr auto compartments_path = Assembler::multiIndex();
#endif

  std::vector<DomainLOP> _lops;

  Space _space;
  ParameterTree _config;

#if DUNE_COPASI_HAVE_MEMBRANE_SPACE
  static constexpr auto dim = Space::EntitySet::dimension;
  using MembraneLocalBasis = std::decay_t<decltype(_space.localView().tree().child(Indices::_1).child(0).child(0).child(0).finiteElement().localBasis())>;
  using MembraneRange = typename MembraneLocalBasis::Traits::RangeType;

  std::vector<MembraneRange> _phiu_f;
#endif


public:
  //! pattern assembly flags
  static constexpr bool doPatternVolume = true;

  //! residual assembly flags
  static constexpr bool doAlphaVolume = true;

#if DUNE_COPASI_HAVE_MEMBRANE_SPACE
  static constexpr bool doPatternSkeleton = true;
  static constexpr bool doAlphaSkeleton = true;

  // (domain_i,domain_o) -> membrane_id
  std::map<std::array<std::size_t,2>, std::size_t> _membrane2id;  
#endif

  //! selective assembly flags
  static constexpr bool doSkipEntity = false;
  static constexpr bool doSkipIntersection = true;

  /**
   * @brief      Constructs a new instance.
   *
   * @param[in]  space           The discrete function space
   * @param[in]  config          The configuration
   */
  TemporalLocalOperatorMultiDomainDiffusionReaction(
    const Space space,
    const ParameterTree& config)
  : _space(space)
  , _config(config)
  {
    setup_lops();
  }

  void setup_lops()
  {
    TRACE_EVENT("dune", "LocalOperator::MultiDomain::ParserSetUp");

    _lops.clear();
    auto compartments_space = _space.subSpace(compartments_path);
    for (std::size_t domain_i = 0; domain_i != compartments_space.degree(); ++domain_i) {
      auto space_i = compartments_space.subSpace(Assembler::multiIndex(domain_i));
      const auto& domain_config = _config.sub(space_i.name(), true);
      _lops.emplace_back(space_i, domain_config);
      
#if DUNE_COPASI_HAVE_MEMBRANE_SPACE
      auto space_ms = _space.subSpace(Assembler::multiIndex(Indices::_1));
      for (std::size_t domain_o = 0; domain_o != compartments_space.degree(); ++domain_o) {
        if (domain_i == domain_o)
          continue;
        auto space_o = compartments_space.subSpace(Assembler::multiIndex(domain_o));
        { // build _membrane2id from space names
          auto membrane_io = fmt::format("{}-{}", space_i.name(), space_o.name());
          auto membrane_oi = fmt::format("{}-{}", space_o.name(), space_i.name());
          for (std::size_t membrane_id = 0; membrane_id != space_ms.degree(); ++membrane_id) {
            auto space_m = space_ms.subSpace(Assembler::multiIndex(membrane_id));
            if (membrane_io == space_m.name() or membrane_oi == space_m.name()) {
              _membrane2id[{domain_i,domain_o}] = _membrane2id[{domain_o,domain_i}] = membrane_id;
              break;
            }
          }
        }
      }
#endif
    }
  }

  template<class EG>
  auto subDomain(const EG& entity) const
  {
    auto domain_set = _space.entitySet().indexSet().subDomains(entity);
    assert(domain_set.size() == 1);
    return *(domain_set.begin());
  }

  template<class IG>
  bool skip_intersection(const IG& ig) const
  {
    auto domain_i = subDomain(ig.inside());

    bool skip = true;
    if constexpr (DomainLOP::doSkipIntersection)
      skip &= _lops[domain_i].skip_intersection(ig);

    if (ig.neighbor()) {
      auto domain_o = subDomain(ig.outside());
      skip &= (domain_i == domain_o);
    }
    return skip;
  }

  /**
   * @copydoc TemporalLocalOperatorDiffusionReactionCG::pattern_volume
   */
  template<class Entity, typename LFSU, typename LFSV, typename LocalPattern>
  void pattern_volume(const Entity entity, const LFSU& lfsu,
                      const LFSV& lfsv,
                      LocalPattern& pattern) const
  {
    auto domain = subDomain(entity);
    auto domain_path = join(compartments_path, Assembler::multiIndex(domain));
    const auto& sub_lfsu = child(lfsu, domain_path);
    const auto& sub_lfsv = child(lfsv, domain_path);
    _lops[domain].pattern_volume(entity, sub_lfsu, sub_lfsv, pattern);
  }

#if DUNE_COPASI_HAVE_MEMBRANE_SPACE
  template<class IG, typename LFSU, typename LFSV, typename LocalPattern>
  void pattern_skeleton(const IG& ig,
                        const LFSU& lfsu_i,
                        const LFSV& lfsv_i,
                        const LFSU& lfsu_o,
                        const LFSV& lfsv_o,
                        LocalPattern& pattern_ii,
                        LocalPattern& pattern_io,
                        LocalPattern& pattern_oi,
                        LocalPattern& pattern_oo) const
  {
    std::size_t domain_i = subDomain(ig.inside());
    std::size_t domain_o = subDomain(ig.outside());

    assert(lfsu_i.degree() == lfsv_i.degree());
    assert(lfsu_o.degree() == lfsv_o.degree());

    auto membrane_id_it = _membrane2id.find({domain_i, domain_o});
    assert((membrane_id_it != _membrane2id.end()) && "Membrane ID is missing");
    auto membrane_id = membrane_id_it->second;

    auto membrane_path = Assembler::multiIndex(Indices::_1, membrane_id);
    const auto& lfsu_f = child(lfsu_i, membrane_path);
    const auto& lfsv_f = child(lfsv_i, membrane_path);
    const auto face_id = ig.indexInInside();

    assert(lfsv_f.degree() == lfsu_f.degree());
    for (std::size_t k = 0; k != lfsv_f.degree(); ++k) {
      const auto& lfsu_cf = lfsu_f.child(k).child(face_id);
      const auto& lfsv_cf = lfsv_f.child(k).child(face_id);
      for (std::size_t i = 0; i != lfsu_cf.size(); ++i)
        for (std::size_t j = 0; j != lfsv_cf.size(); ++j)
          pattern_ii.addLink(lfsv_cf, i, lfsu_cf, j);
    }
  }
#endif

  /**
   * @copydoc LocalOperatorMultiDomainDiffusionReaction::alpha_volume
   */
  template<typename EG, typename LFSU, typename X, typename LFSV, typename R>
  void alpha_volume(const EG& entity,
                    const LFSU& lfsu,
                    const X& x,
                    const LFSV& lfsv,
                    R& r)
  {
    auto domain = subDomain(entity);
    auto domain_path = join(compartments_path, Assembler::multiIndex(domain));
    const auto& sub_lfsu = child(lfsu, domain_path);
    const auto& sub_lfsv = child(lfsv, domain_path);
    _lops[domain].alpha_volume(entity, sub_lfsu, x, sub_lfsv, r);
  }

  /**
   * @copydoc LocalOperatorMultiDomainDiffusionReaction::jacobian_volume
   */
  template<typename EG, typename LFSU, typename X, typename LFSV, typename Mat>
  void jacobian_volume(const EG& entity,
                       const LFSU& lfsu,
                       const X& x,
                       const LFSV& lfsv,
                       Mat& mat)
  {
    auto domain = subDomain(entity);
    auto domain_path = join(compartments_path, Assembler::multiIndex(domain));
    const auto& sub_lfsu = child(lfsu, domain_path);
    const auto& sub_lfsv = child(lfsv, domain_path);
    _lops[domain].jacobian_volume(entity, sub_lfsu, x, sub_lfsv, mat);
  }

  /**
   * @copydoc LocalOperatorMultiDomainDiffusionReaction::jacobian_apply_volume
   */
  template<typename EG, typename LFSU, typename X, typename LFSV, typename R>
  void jacobian_apply_volume(const EG& entity,
                             const LFSU& lfsu,
                             const X& x,
                             const X& z,
                             const LFSV& lfsv,
                             R& r)
  {
    auto domain = subDomain(entity);
    auto domain_path = join(compartments_path, Assembler::multiIndex(domain));
    const auto& sub_lfsu = child(lfsu, domain_path);
    const auto& sub_lfsv = child(lfsv, domain_path);
    _lops[domain].jacobian_apply_volume(entity, sub_lfsu, x, z, sub_lfsv, r);
  }

#if DUNE_COPASI_HAVE_MEMBRANE_SPACE
  template<typename IG, typename LFSU, typename X, typename LFSV, typename R>
  void alpha_skeleton(const IG& ig,
                      const LFSU& lfsu_i,
                      const X& x_i,
                      const LFSV& lfsv_i,
                      const LFSU& lfsu_o,
                      const X& x_o,
                      const LFSV& lfsv_o,
                      R& r_i,
                      R& r_o)
  {
    const auto& entity_f = ig;
    const auto& entity_i = ig.inside();
    const auto& entity_o = ig.outside();

    auto geo_f = entity_f.geometry();
    auto geo_i = entity_i.geometry();
    auto geo_in_i = entity_f.geometryInInside();

    assert(lfsu_i.degree() == lfsv_i.degree());
    assert(lfsu_o.degree() == lfsv_o.degree());

    std::size_t domain_i = subDomain(entity_i);
    std::size_t domain_o = subDomain(entity_o);

    assert(ig.conforming());
    auto membrane_id_it = _membrane2id.find({domain_i, domain_o});
    assert((membrane_id_it != _membrane2id.end()) && "Membrane ID is missing");
    auto membrane_id = membrane_id_it->second;

    auto membrane_path = Assembler::multiIndex(Indices::_1, membrane_id);
    const auto& lfsu_f = child(lfsu_i, membrane_path);
    const auto& lfsv_f = child(lfsv_i, membrane_path);
    auto face_id = ig.indexInInside();

    MembraneLocalBasis const * local_basis_f = nullptr;
    if (lfsv_f.degree() == 0)
      return;
    else
      local_basis_f = &lfsv_f.child(0).child(face_id).finiteElement().localBasis();

    for (const auto& it : quadratureRule(geo_f, 3)) {
      const auto& position_f = it.position();

      const auto& geo_if = entity_i.template subEntity<1>(face_id).geometry();
      local_basis_f->evaluateFunction(geo_if.local(geo_f.global(position_f)), _phiu_f);

      auto factor = it.weight() * geo_f.integrationElement(position_f);

      // contribution for each component
      for (std::size_t comp = 0; comp != lfsu_f.degree(); ++comp){
        const auto& lfsu_cf = lfsu_f.child(comp).child(face_id);
        const auto& lfsv_cf = lfsv_f.child(comp).child(face_id);

        double u = 0.0;
        for (std::size_t dof = 0; dof != lfsu_cf.size(); ++dof)
          u += x_i(lfsu_cf, dof) * _phiu_f[dof];

        for (std::size_t dof = 0; dof != lfsv_cf.size(); ++dof)
          r_i.accumulate(lfsv_cf, dof, 0.0000*u * _phiu_f[dof] * factor);
      }
    }
  }

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
                         J& mat_oo)
  {
    const auto& entity_f = ig;
    const auto& entity_i = ig.inside();
    const auto& entity_o = ig.outside();

    auto geo_f = entity_f.geometry();
    auto geo_in_i = entity_f.geometryInInside();

    assert(lfsu_i.degree() == lfsv_i.degree());
    assert(lfsu_o.degree() == lfsv_o.degree());

    std::size_t domain_i = subDomain(entity_i);
    std::size_t domain_o = subDomain(entity_o);

    assert(ig.conforming());
    auto membrane_id_it = _membrane2id.find({domain_i, domain_o});
    assert((membrane_id_it != _membrane2id.end()) && "Membrane ID is missing");
    auto membrane_id = membrane_id_it->second;

    auto membrane_path = Assembler::multiIndex(Indices::_1, membrane_id);
    const auto& lfsu_f = child(lfsu_i, membrane_path);
    const auto& lfsv_f = child(lfsv_i, membrane_path);
    auto face_id = ig.indexInInside();

    MembraneLocalBasis const * local_basis_f = nullptr;
    if (lfsv_f.degree() == 0)
      return;
    else
      local_basis_f = &lfsv_f.child(0).child(face_id).finiteElement().localBasis();

    for (const auto& it : quadratureRule(geo_f, 3)) {
      const auto& position_f = it.position();

      const auto& geo_if = entity_i.template subEntity<1>(face_id).geometry();
      local_basis_f->evaluateFunction(geo_if.local(geo_f.global(position_f)), _phiu_f);

      auto factor = it.weight() * geo_f.integrationElement(position_f);

      // contribution for each component
      for (std::size_t comp = 0; comp != lfsu_f.degree(); ++comp){
        const auto& lfsu_cf = lfsu_f.child(comp).child(face_id);
        const auto& lfsv_cf = lfsv_f.child(comp).child(face_id);

        for (std::size_t dof_i = 0; dof_i != lfsv_cf.size(); ++dof_i)
          for (std::size_t dof_j = 0; dof_j != lfsu_cf.size(); ++dof_j)
            mat_ii.accumulate(lfsv_cf, dof_i, lfsu_cf, dof_j, 0.0000*_phiu_f[dof_i] * _phiu_f[dof_j] * factor);
      }
    }
  }

  template<class R, class X>
  struct PseudoJacobian {
    void accumulate(const auto& ltest, auto test_dof, const auto& ltrial, auto trial_dof, auto value) {
      _r.accumulate(ltest, test_dof, _z(ltrial, trial_dof) * value);
    }

    R& _r;
    const X& _z;
  };

  template<typename IG, typename LFSU, typename X, typename LFSV, typename R>
  void jacobian_apply_skeleton(const IG& ig,
                      const LFSU& lfsu_i,
                      const X& x_i,
                      const X& z_i,
                      const LFSV& lfsv_i,
                      const LFSU& lfsu_o,
                      const X& x_o,
                      const X& z_o,
                      const LFSV& lfsv_o,
                      R& r_i,
                      R& r_o)
  {
    PseudoJacobian<R, X> mat_ii{r_i, z_i}, mat_io{r_i, z_o}, mat_oi{r_o, z_i}, mat_oo{r_o, z_o};
    jacobian_skeleton(ig, lfsu_i, x_i, lfsv_i, lfsu_o, x_o, lfsv_o, mat_ii, mat_io, mat_oi, mat_oo);
  }
#endif
};

} // namespace Dune::Copasi

#endif // DUNE_COPASI_LOCAL_OPERATOR_MULTIDOMAIN_DIFFUSION_REACTION_HH
