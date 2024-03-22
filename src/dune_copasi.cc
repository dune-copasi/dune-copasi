
#ifndef DUNE_COPASI_GRID_DIMENSIONS
// comma separated list of dimensions to compile
#define DUNE_COPASI_GRID_DIMENSIONS 2
#endif

#include <dune-copasi-config.hh>

#include <dune/copasi/common/stepper.hh>
#include <dune/copasi/common/exceptions.hh>
#include <dune/copasi/common/fmt_style.hh>

#if HAVE_PERFETTO
#include <dune/copasi/common/ostream_to_spdlog.hh>
#endif

#include <dune/copasi/grid/make_multi_domain_grid.hh>
#include <dune/copasi/model/factory.hh>
#include <dune/copasi/model/local_equations/functor_factory_parser.hh>
#include <dune/copasi/model/model.hh>
#include <dune/copasi/parser/context.hh>
#include <dune/copasi/parser/factory.hh>
#include <dune/copasi/parser/parser.hh>

#include <dune/pdelab/common/trace.hh>

#include <dune/grid/multidomaingrid/mdgridtraits.hh>
#include <dune/grid/multidomaingrid/multidomaingrid.hh>

#include <dune/grid/uggrid.hh>
#include <dune/grid/yaspgrid.hh>
#include <dune/grid/yaspgrid/coordinates.hh>

#include <dune/common/exceptions.hh>
#include <dune/common/parallel/mpihelper.hh>
#include <dune/common/parametertree.hh>
#include <dune/common/parametertreeparser.hh>
#include <dune/common/hybridutilities.hh>

#include <spdlog/spdlog.h>

#include <fmt/core.h>
#include <fmt/ranges.h>

#include <nlohmann/json.hpp>

// includes a string with json contents generated a configuration time
#include "config_opts.hh"

#include <algorithm>
#include <array>
#include <cstddef>
#include <exception>
#include <functional>
#include <filesystem>
#include <iostream>
#include <memory>
#include <optional>
#include <type_traits>
#include <sstream>
#include <string>
#include <string_view>
#include <utility>
#include <vector>

void collapse_subtypes(const auto& json, std::set<std::string_view>& set) {
  if (json.contains(".type"))
    set.insert(json.at(".type").template get<std::string_view>());
  if (json.contains(".options")) {
    for (auto& [key, values] : json.at(".options").items())
      collapse_subtypes(values, set);
  }
}

void print_config_node(const auto& config_node, unsigned verbosity, unsigned indent_level) {
  auto& [key, json] = config_node;
  std::string indent(indent_level * 2, ' ');
  fmt::print("{}{}", indent, DUNE_COPASI_FMT_STYLED_BOLD(key));
  std::set<std::string_view> type_set;
  collapse_subtypes(json, type_set);
  if (not type_set.empty()) {
    std::string types_str = fmt::format("{}", fmt::join(type_set, "|"));
    fmt::print("={}\n", fmt::format("<{}>",DUNE_COPASI_FMT_STYLED_ITALIC(types_str)));
  }
  else
    fmt::print("\n");

  indent += "    ";
  if (json.contains(".brief"))
    fmt::print("{}{}\n", indent, DUNE_COPASI_FMT_STYLED_DARK_GRAY(json.at(".brief").template get<std::string_view>()));
  if (json.contains(".default"))
    fmt::print("{}{}{}\n", indent, DUNE_COPASI_FMT_STYLED_ITALIC("Default: "), DUNE_COPASI_FMT_STYLED_ITALIC(json.at(".default").template get<std::string_view>()));


  if (verbosity == 0)
    return;

  indent += "  ";
  if (json.contains(".details")) {
    fmt::print("\n");
    for (auto& details_j : json.at(".details")) {
      auto details = details_j.template get<std::string_view>();
      if (indent.size() >= 120)
        throw Dune::Copasi::format_exception(Dune::IOError{}, "Maximum indentation reached!");
      std::size_t max_width = 120u - indent.size();
      std::size_t l = 0;
      // find lines to print
      while (l < details.size()) {
        // split lines by max text width
        auto line = details.substr(l, max_width);
        // split lines by new lines
        line = line.substr(0, line.find('\n'));
        // break last word into the next line
        if (line.size() == max_width and not line.ends_with(' ') and
            not details.substr(l + line.size()).starts_with(' ') and
            line.rfind(' ') != line.find(' '))
          line = line.substr(0, line.rfind(' '));
        l += line.size();
        // trim spaces
        line.remove_prefix(std::min(line.find_first_not_of(' '), line.size()));
        fmt::print("{}{}\n", indent, DUNE_COPASI_FMT_STYLED_DARK_GRAY(line));
      }
    }
    fmt::print("\n");
  }

  if (json.contains(".options")) {
    for (auto& node : json.at(".options").items())
      print_config_node(node, std::max(verbosity,1u) - 1, indent.size() / 2);
  }
}

void
program_help(std::string_view prog_name, bool long_help)
{
  fmt::print(
    "USAGE: {} [options]\n\n"
    "Numerical simulator for diffusion-reaction systems on single or multiple compartments\n\n"
    "Website: <https://dune-copasi.netlify.app>\n\n"
    "OPTIONS:\n\n"
    "  -h / --help          - Display this help.\n"
    "  --help-full          - Display this help with long descriptions.\n"
    "  --version            - Display the version of this program.\n"
    "  --parser-default     - Display the default parser.\n"
    "  --parser-list        - Display the parsers available for this program.\n"
    "  --dimension-list     - Display the grid dimensions available for this program.\n"
    "  --dump-config        - Dumps configuration in the INI format to stdout.\n"
    "  --config=<string>    - Specifies a config file in INI format. See Configuration Options.\n"
    "  --{{key}}={{value}}      - Overrides key=value sections of the config file. See Configuration Options.\n\n",
    prog_name);

  auto data = nlohmann::json::parse(config_opts_json);
  if (not data.empty()) {
    fmt::print("Configuration Options:\n\n");
    for(auto& node : data.items())
      print_config_node(node, long_help ? 1000 : 0, 1);
  }

  std::cout << std::endl;
  if (long_help) {
    fmt::print(
      "EXAMPLE:\n"
      "  Step diffusion:\n"
      "    The following command will solve the laplace equation for a\n"
      "    1D scalar field variable named 'u' in a compartment named 'domain'.\n"
      "    The initial condition is an expression which has an step when the\n"
      "    coordinate 'x' is 0.5. The storage term is set to 1.0 because the\n"
      "    problem is transient and the diffusion coefficient is set to 0.001\n"
      "    to blend the step over time. Finally, the solution will be written\n"
      "    in the vtk format using the keyword 'step_diffusion' used to format\n"
      "    the output files.\n\n"
      "      {} \\\n"
      "        --grid.dimension=1 \\\n"
      "        --grid.refinement_level=5 \\\n"
      "        --compartments.domain.type=expression \\\n"
      "        --compartments.domain.expression=1 \\\n"
      "        --model.scalar_field.u.compartment=domain \\\n"
      "        --model.scalar_field.u.initial.expression=\"(position_x > 0.5)\" \\\n"
      "        --model.scalar_field.u.storage.expression=1 \\\n"
      "        --model.scalar_field.u.cross_diffusion.u.expression=0.001 \\\n"
      "        --model.writer.vtk.path=step_diffusion\n\n",
      prog_name);
  }
}

int
main(int argc, char** argv)
{
  const std::filesystem::path prog_path = argv[0];

  int end_code = 0;
  Dune::ParameterTree config;

  auto log_error = [](std::string message){
    std::istringstream stream(message);
    std::string line;
    while (std::getline(stream, line)) {
      spdlog::error("  {}", line);
    }
  };

  std::vector<char*> cmd_args(argv, argv + argc);
  auto is_help = [](std::string_view opt) { return opt.starts_with("--help") or opt == "-h"; };
  if (auto hlp = std::ranges::find_if(cmd_args, is_help); hlp != cmd_args.end()) {
    program_help(prog_path.filename().string(), std::string_view{ *hlp } == "--help-full");
    return 0;
  }

  auto is_version = [](std::string_view opt) { return opt == "--version"; };
  if (std::ranges::find_if(cmd_args, is_version) != cmd_args.end()) {
    fmt::print("{}\n", DUNE_COPASI_VERSION);
    return 0;
  }

  auto is_parser_default = [](std::string_view opt) { return opt == "--parser-default"; };
  if (std::ranges::find_if(cmd_args, is_parser_default) != cmd_args.end()) {
    fmt::print("{}\n", Dune::Copasi::default_parser_str);
    return 0;
  }

  auto is_parser_list = [](std::string_view opt) { return opt == "--parser-list"; };
  if (std::ranges::find_if(cmd_args, is_parser_list) != cmd_args.end()) {
    std::unique_ptr<Dune::Copasi::Parser> parser;
    std::vector<std::string> parsers;
    for (auto [parser_str, _] : Dune::Copasi::string2parser)
      try {
        parser = Dune::Copasi::make_parser(parser_str);
        parsers.emplace_back(parser_str);
      } catch (...) {}
    fmt::print("{}\n", parsers);
    return 0;
  }

  auto is_dimension_list = [](std::string_view opt) { return opt == "--dimension-list"; };
  if (std::ranges::find_if(cmd_args, is_dimension_list) != cmd_args.end()) {
    fmt::print("{}\n", std::array{DUNE_COPASI_GRID_DIMENSIONS});
    return 0;
  }

  bool dump_config = false;
  auto is_dump_config = [](std::string_view opt) { return opt == "--dump-config"; };
  if (auto dmp = std::ranges::find_if(cmd_args, is_dump_config); dmp != cmd_args.end()) {
    cmd_args.erase(dmp);
    dump_config = true;
  }

  try {
    auto is_config = [](std::string_view opt) { return opt.starts_with("--config="); };
    if (auto cfg_it = std::ranges::find_if(cmd_args, is_config); cfg_it != cmd_args.end()) {
      auto cfg_file = std::string{ *cfg_it }.substr(9);
      if (not exists(std::filesystem::path{cfg_file})) {
        throw Dune::Copasi::format_exception(Dune::IOError{}, "Configuration file '{}' does not exist", cfg_file);
      }
      if (not dump_config) {
        spdlog::info("Reading configuration file '{}'", cfg_file);
      }
      Dune::ParameterTreeParser::readINITree(cfg_file, config);
    }
    Dune::ParameterTreeParser::readNamedOptions(cmd_args.size(), cmd_args.data(), config, {});
  } catch (Dune::Exception& e) {
    spdlog::error("Invalid arguments!");
    log_error(e.what());
    spdlog::error("dune-copasi finished with some errors :(");
    return 1;
  }
  if (dump_config) {
    config.report(std::cout);
    return 0;
  }

  spdlog::info("Starting dune-copasi (version: {})", DUNE_COPASI_VERSION);
  auto trace_path = config.get("trace.path", "");
#if HAVE_PERFETTO
  std::unique_ptr<Dune::PDELab::TracingSession> tracing_session;
  if (not trace_path.empty()) {
    tracing_session = [trace_path]() {
      [[maybe_unused]] auto ostream_guard = Dune::Copasi::ostream2spdlog();
      return std::make_unique<Dune::PDELab::TracingSession>(trace_path);
    }();
  }
#else
  if (not trace_path.empty()) {
    spdlog::warn("This executable cannot generate traces. The 'trace.path' argument will be ignored.");
  }
#endif
  {
    TRACE_EVENT("dune", "MPI::Init");
    Dune::MPIHelper::instance(argc, argv);
  }

  try {
    using namespace Dune::Copasi;

    // find grid dimension
    std::size_t config_dim;
    if (config.hasKey("grid.dimension")) {
      config_dim = config.template get<std::size_t>("grid.dimension");
    } else if (config.hasKey("grid.extensions")) {
      config_dim = config.template get<std::vector<double>>("grid.extensions").size();
    } else {
      config_dim = 2;
    }

    if (config.hasKey("grid.axis_names"))
      axis_names = config.template get<std::vector<std::string>>("grid.axis_names");
    spdlog::info("Axis names are set to: {}", axis_names);
    if (axis_names.size() < config_dim)
      spdlog::warn("There are less axis names ({}) than number dimensions ({})", axis_names.size(), config_dim);
    if (not std::ranges::unique(axis_names).empty())
      throw format_exception(Dune::IOError{}, "Axis names are not unique!");

    auto parser_context = std::make_shared<ParserContext>(config.sub("parser_context"));

    // unroll dimensions and select the one case defined at run-time
    constexpr auto dims = std::index_sequence<DUNE_COPASI_GRID_DIMENSIONS>{};
    Dune::Hybrid::switchCases(
      dims,
      config_dim,
      [&](auto dim) {
        // get a pointer to the grid
        const auto max_subdomains = 64;
        auto [md_grid_ptr, cell_data] = [&] {
          using MDGTraits = Dune::mdgrid::FewSubDomainsTraits<dim, max_subdomains>;
          if constexpr (dim < 2) {
            using MDGrid = Dune::mdgrid::MultiDomainGrid<Dune::YaspGrid<dim, Dune::EquidistantOffsetCoordinates<double,dim>>, MDGTraits>;
            return make_multi_domain_grid_with_cell_data<MDGrid>(config, parser_context);
          } else {
            using MDGrid = Dune::mdgrid::MultiDomainGrid<Dune::UGGrid<dim>, MDGTraits>;
            return make_multi_domain_grid_with_cell_data<MDGrid>(config, parser_context);
          }
        }();

        const auto& model_config = config.sub("model");
        if (model_config.hasSub("compartments"))
          spdlog::warn(
            "The section '[model.compartments]' will be ignored, use '[compartments]' instead");
        config.sub("model.compartments") = config.sub("compartments");

        if (auto comp_size = model_config.sub("compartments").getSubKeys().size();
            comp_size > max_subdomains)
          spdlog::error("The model can have at most '{}' compartments but '{}' were defined. To "
                        "have more compartments recompile with appropriate grid triats!",
                        max_subdomains,
                        comp_size);

        using MDGrid = std::decay_t<decltype(*md_grid_ptr)>;
        using SDGridView = typename MDGrid::SubDomainGrid::Traits::LeafGridView;
        using SpeciesQuantity = double;
        using TimeQuantity = double;
        using DurationQuantity = double;
        using Model = Model<MDGrid, SDGridView, SpeciesQuantity, TimeQuantity>;

        auto parser_type = string2parser.at(config.get("model.parser_type", default_parser_str));
        auto functor_factory =
          std::make_shared<FunctorFactoryParser<dim>>(parser_type, std::move(parser_context));
        std::shared_ptr model = make_model<Model>(model_config, functor_factory, std::move(cell_data));

        // create time stepper
        const auto& time_config = model_config.sub("time_step_operator");
        auto decrease_factor = time_config.get("time_step_decrease_factor", 0.5);
        auto increase_factor = time_config.get("time_step_increase_factor", 1.1);
        auto begin_time = time_config.get("time_begin", TimeQuantity{ 0. });
        auto end_time = time_config.get("time_end", TimeQuantity{ 1. });
        std::optional<DurationQuantity> min_step;
        std::optional<DurationQuantity> max_step;
        if (time_config.hasKey("time_step_min")) {
          min_step = time_config.template get<DurationQuantity>("time_step_min");
        }
        if (time_config.hasKey("time_step_max")) {
          max_step = time_config.template get<DurationQuantity>("time_step_max");
        }
        auto initial_step = time_config.get("time_step_initial", max_step.value_or(.1));

        using State = typename Model::State;
        auto stepper = SimpleAdaptiveStepper<State, TimeQuantity, DurationQuantity>{
          decrease_factor, increase_factor, min_step, max_step
        };

        // setup functor for each step
        std::function<void(const State&)> on_each_step;
        // setup the vtk writer
        auto vtk_path = model_config.get("writer.vtk.path", "");
        // define each time step to write out (standard value = 0. implies every step is written out)
        auto time_step = model_config.get("writer.time_step", TimeQuantity{ 0. });
        auto time_write = begin_time;

        on_each_step = [model_config, vtk_path, &time_write, time_step, model](const auto& state) {
          // write to vtk if requested
          if (not vtk_path.empty() ) {
            if( state.time >= time_write - 1e-9 * time_step){ // small epsilon correction
              model->write_vtk(state, vtk_path, true);
              time_write += time_step;
            }
          }
          // evaluate transform/reduce operations on the model state
          if (model_config.hasSub("reduce")) {
            model->reduce(state, model_config.sub("reduce"));
          }
        };

        auto in = model->make_state(std::move(md_grid_ptr), model_config);
        in->time = begin_time;

        // interpolate initial conditions into the model state
        model->interpolate(*in, model->make_initial(*(in->grid), model_config));
        if (on_each_step) {
          on_each_step(*in); // write initial condition
        }

        auto step_operator = model->make_step_operator(*in, model_config);

        stepper.evolve(*step_operator, *in, *in, initial_step, end_time, on_each_step).or_throw();
      },
      [config_dim]() {
        throw Dune::Copasi::format_exception(
          Dune::IOError{},
          "\n  This executable cannot simulate in {}D grids!\n  Available dimensions are: {}.",
          config_dim,
          std::array{DUNE_COPASI_GRID_DIMENSIONS});
      });
  } catch (Dune::NotImplemented& e) {
    spdlog::error("Feature is not implemented:");
    log_error(e.what());
    end_code = 1;
  } catch (Dune::Exception& e) {
    spdlog::error("Dune reported error:");
    log_error(e.what());
    end_code = 1;
  } catch (std::exception& e) {
    spdlog::error("C++ reported error:");
    log_error(e.what());
    end_code = 1;
  } catch (...) {
    spdlog::error("Unknown exception thrown!");
    end_code = 1;
  }

  if (end_code != 0) {
    spdlog::error("dune-copasi finished with some errors :(");
  } else {
    spdlog::info("dune-copasi successfully finished :)");
  }

#if HAVE_PERFETTO
  {
    auto ostream_guard = Dune::Copasi::ostream2spdlog();
    tracing_session.reset();
  }
#endif

  return end_code;
}
