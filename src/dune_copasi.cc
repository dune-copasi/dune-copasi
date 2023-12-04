
#ifndef DUNE_COPASI_GRID_DIMENSIONS
// comma separated list of dimensions to compile
#define DUNE_COPASI_GRID_DIMENSIONS 2
#endif

#include <dune/copasi/common/stepper.hh>
#include <dune/copasi/grid/make_multi_domain_grid.hh>
#include <dune/copasi/model/factory.hh>
#include <dune/copasi/model/local_equations/functor_factory_parser.hh>
#include <dune/copasi/model/model.hh>
#include <dune/copasi/parser/context.hh>
#include <dune/copasi/parser/factory.hh>

#include <dune/pdelab/common/trace.hh>

#include <dune/grid/multidomaingrid/mdgridtraits.hh>
#include <dune/grid/multidomaingrid/multidomaingrid.hh>

#include <dune/grid/uggrid.hh>
#include <dune/grid/yaspgrid.hh>

#include <dune/common/exceptions.hh>
#include <dune/common/float_cmp.hh>
#include <dune/common/fvector.hh>
#include <dune/common/indices.hh>
#include <dune/common/parallel/mpihelper.hh>
#include <dune/common/parametertree.hh>
#include <dune/common/parametertreeparser.hh>

#include <spdlog/spdlog.h>

#include <fmt/color.h>
#include <fmt/core.h>

// Generated options for the done-copasi configuration ini-file
#if __has_include("dune-copasi-config-file-options.hh")
#include "dune-copasi-config-file-options.hh"
#else
static const std::vector<std::array<std::string, 4>> config_file_opts;
#endif

#include <algorithm>
#include <array>
#include <cassert>
#include <chrono>
#include <cstddef>
#include <exception>
#include <functional>
#include <iostream>
#include <memory>
#include <optional>
#include <set>
#include <sstream>
#include <string>
#include <utility>
#include <vector>


void
program_help(std::string_view prog_name, bool long_help)
{
  fmt::print(
    "USAGE: {} [options]\n\n"
    "Numerical simulator for diffusion-reaction systems on single or multiple compartments\n\n"
    "Website: <https://dune-copasi.netlify.app>\n\n"
    "OPTIONS:\n\n"
    "  -h / --help          - Display this help\n"
    "  --help-full          - Display this help with long descriptions\n"
    "  --version            - Display the version of this program\n"
    "  --config=<string>    - Specifies a config file in INI format. See Configuration Options\n"
    "  --dump-config        - Dumps configuration in the INI format to stdout\n"
    "  --{{key}}={{value}}      - Overrides key=value sections of the config file\n\n",
    prog_name);

  if (not config_file_opts.empty()) {
    fmt::print("Configuration Options:\n\n");
    for (auto [key, type, short_doc, long_doc] : config_file_opts) {
#if FMT_VERSION >= 90000
      auto sty_dash = fmt::styled("--" + key, fmt::emphasis::bold);
      auto sty_time = fmt::styled(type, fmt::emphasis::italic);
      auto sty_short_doc = fmt::styled(short_doc, fmt::fg(fmt::color::dark_gray));
#else
      auto sty_dash = "--";
      auto sty_time = type;
      auto sty_short_doc = short_doc;
#endif
      fmt::print("  {}={}\n     {}\n", sty_dash, sty_time, sty_short_doc);
      if (long_help and not long_doc.empty()) {
        std::istringstream iss(long_doc);
        for (std::string line; std::getline(iss, line);) {
#if FMT_VERSION >= 90000
          auto sty_line = fmt::styled(line, fmt::fg(fmt::color::dark_gray));
#else
          auto sty_line = line;
#endif
          fmt::print("       {}\n", sty_line);
        }
      }
    }
  }
  std::cout << std::endl;
}

int
main(int argc, char** argv)
{
  const Dune::Copasi::fs::path prog_path = argv[0];

  int end_code = 0;
  Dune::ParameterTree config;

  std::vector<char*> cmd_args(argv, argv + argc);
  auto is_help = [](std::string_view opt) { return opt.starts_with("--help") or opt == "-h"; };
  if (auto hlp = std::ranges::find_if(cmd_args, is_help); hlp != cmd_args.end()) {
    program_help(prog_path.filename().string(), std::string_view{ *hlp } == "--help-full");
    return 0;
  }

  auto is_version = [](std::string_view opt) { return opt == "--version"; };
  if (std::ranges::find_if(cmd_args, is_version) != cmd_args.end()) {
    std::cout << DUNE_COPASI_VERSION << std::endl;
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
      if (not exists(Dune::Copasi::fs::path{cfg_file})) {
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
    spdlog::error("{}", e.what());
    spdlog::error("dune-copasi finished with some errors :(");
    return 1;
  }
  if (dump_config) {
    config.report(std::cout);
    return 0;
  }

  spdlog::info("Starting dune-copasi");
  auto trace_path = config.get("trace.path", "");
#if HAVE_PERFETTO
  std::unique_ptr<Dune::PDELab::TracingSession> tracing_session;
  if (not trace_path.empty()) {
    tracing_session = [trace_path]() {
      auto ostream_guard = Dune::Copasi::ostream2spdlog();
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

    auto parser_context = std::make_shared<ParserContext>(config.sub("parser_context"));

    // find grid dimension
    std::size_t config_dim;
    if (config.hasKey("grid.dimension")) {
      config_dim = config.template get<std::size_t>("grid.dimension");
    } else if (config.hasKey("grid.extensions")) {
      config_dim = config.template get<std::vector<double>>("grid.extensions").size();
    } else {
      config_dim = 2;
    }

    // unroll dimensions and select the one case defined at run-time
    constexpr auto dims = std::index_sequence<DUNE_COPASI_GRID_DIMENSIONS>{};
    Dune::Hybrid::switchCases(
      dims,
      config_dim,
      [&](auto dim) {
        // get a pointer to the grid
        auto md_grid_ptr = [&] {
          using MDGTraits = Dune::mdgrid::DynamicSubDomainCountTraits<dim, 10>;
          if constexpr (dim < 2) {
            using MDGrid = Dune::mdgrid::MultiDomainGrid<Dune::YaspGrid<dim, Dune::EquidistantOffsetCoordinates<double,dim>>, MDGTraits>;
            return make_multi_domain_grid<MDGrid>(config, parser_context);
          } else {
            using MDGrid = Dune::mdgrid::MultiDomainGrid<Dune::UGGrid<dim>, MDGTraits>;
            return make_multi_domain_grid<MDGrid>(config, parser_context);
          }
        }();

        const auto& model_config = config.sub("model");
        if (model_config.hasSub("compartments"))
          spdlog::warn(
            "The section '[model.compartments]' will be ignored, use '[compartments]' instead");
        config.sub("model.compartments") = config.sub("compartments");

        using MDGrid = std::decay_t<decltype(*md_grid_ptr)>;
        using SDGridView = typename MDGrid::SubDomainGrid::Traits::LeafGridView;
        using SpeciesQuantity = double;
        using TimeQuantity = double;
        using DurationQuantity = double;
        using Model = Model<MDGrid, SDGridView, SpeciesQuantity, TimeQuantity>;

        auto parser_type = string2parser.at(config.get("model.parser_type", default_parser_str));
        auto functor_factory =
          std::make_shared<FunctorFactoryParser<MDGrid>>(parser_type, std::move(parser_context));
        std::shared_ptr model = make_model<Model>(model_config, functor_factory);

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
        Dune::IOError excep;
        excep.message(
          fmt::format("This executable cannot simulate in {}D grids!", std::size_t{ config_dim }));
        throw excep;
      });
  } catch (Dune::NotImplemented& e) {
    spdlog::error("Feature is not implemented:");
    spdlog::error("{}", e.what());
    end_code = 1;
  } catch (Dune::Exception& e) {
    spdlog::error("Dune reported error:");
    spdlog::error("{}", e.what());
    end_code = 1;
  } catch (std::exception& e) {
    spdlog::error("C++ reported error:");
    spdlog::error("{}", e.what());
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
