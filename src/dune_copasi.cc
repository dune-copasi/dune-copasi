
#ifndef DUNE_COPASI_GRID_DIMENSION
#define DUNE_COPASI_GRID_DIMENSION 2
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

#if DUNE_COPASI_GRID_DIMENSION < 2
#include <dune/grid/yaspgrid.hh>
#else
#include <dune/grid/uggrid.hh>
#endif

#include <dune/common/exceptions.hh>
#include <dune/common/float_cmp.hh>
#include <dune/common/fvector.hh>
#include <dune/common/indices.hh>
#include <dune/common/parallel/mpihelper.hh>
#include <dune/common/parametertree.hh>
#include <dune/common/parametertreeparser.hh>

#include <spdlog/spdlog.h>

#include <fmt/format.h>

#include <algorithm>
#include <array>
#include <cassert>
#include <chrono>
#include <cstddef>
#include <exception>
#include <filesystem>
#include <fmt/core.h>
#include <functional>
#include <iostream>
#include <memory>
#include <optional>
#include <set>
#include <sstream>
#include <string>
#include <utility>
#include <vector>

constexpr auto exec_dimension = Dune::index_constant<DUNE_COPASI_GRID_DIMENSION>{};

#if DUNE_COPASI_GRID_DIMENSION < 2
using HostGrid = Dune::YaspGrid<exec_dimension>;
#else
using HostGrid = Dune::UGGrid<exec_dimension>;
#endif

using MDGTraits = Dune::mdgrid::DynamicSubDomainCountTraits<exec_dimension, 10>;
using MDGrid = Dune::mdgrid::MultiDomainGrid<HostGrid, MDGTraits>;

int
main(int argc, char** argv)
{
  spdlog::info("Starting dune-copasi");
  const fs::path prog_name = argv[0];
  const std::string string_help =
    fmt::format("Usage: {} [options] CONFIG_FILE\n\n"
                "Execute numerical simulation of reaction-diffusion systems on single\n"
                "compartments. The CONFIG_FILE is a DUNE INI file with the\n"
                "parameters to perform the simulation.\n",
                prog_name.filename().string());

  Dune::ParameterTree config;

  int end_code = 0;
  std::vector<std::string> cmd_line_args(argv, argv + argc);

  if (argc == 1) {
    std::cerr << string_help;
    return 1;
  }

  if (cmd_line_args[1] == "--help" or cmd_line_args[1] == "-h") {
    std::cout << string_help;
    return 0;
  }

  // Read and parse ini file
  spdlog::info("Reading configuration file: '{}'", cmd_line_args.back());
  Dune::ParameterTreeParser::readINITree(cmd_line_args.back(), config);

#if HAVE_PERFETTO
  auto tracing_session = Dune::PDELab::TracingSession{ prog_name.filename().string() + ".pftrace" };
#endif

  {
    TRACE_EVENT("dune", "MPI::Init");
    Dune::MPIHelper::instance(argc, argv);
  }

  try {
    // detailed report of the input parameter tree
    std::stringstream ss;
    config.report(ss);
    spdlog::info("----");
    for (std::string line; std::getline(ss, line);) {
      spdlog::info("{}", line);
    }
    spdlog::info("----");

    using namespace Dune::Copasi;

    auto parser_context = std::make_shared<ParserContext>(config.sub("parser_context"));

    std::shared_ptr md_grid_ptr = make_multi_domain_grid<MDGrid>(config, parser_context);
    config.sub("model.compartments") = config.sub("compartments");

    using SDGridView = typename MDGrid::SubDomainGrid::Traits::LeafGridView;
    using SpeciesQuantity = double;
    using TimeQuantity = double;
    using DurationQuantity = double;
    using Model = Model<MDGrid, SDGridView, SpeciesQuantity, TimeQuantity>;

    const auto& model_config = config.sub("model", true);
    auto parser_type = string2parser.at(config.get("parser_type", default_parser_str));
    auto functor_factory = std::make_shared<FunctorFactoryParser<exec_dimension>>(
      parser_type, std::move(parser_context));
    std::unique_ptr<Model> model = make_model<Model>(model_config, functor_factory);

    // create time stepper
    const auto& time_config = model_config.sub("time_step_operator", true);
    auto decrease_factor = time_config.get("time_step_decrease_factor", 0.5);
    auto increase_factor = time_config.get("time_step_increase_factor", 1.5);
    auto begin_time = time_config.template get<TimeQuantity>("time_begin");
    auto end_time = time_config.template get<TimeQuantity>("time_end");
    auto initial_step = time_config.template get<DurationQuantity>("time_step_initial");
    std::optional<DurationQuantity> min_step;
    std::optional<DurationQuantity> max_step;
    if (time_config.hasKey("min_step")) {
      min_step = time_config.template get<DurationQuantity>("time_step_min");
    }
    if (time_config.hasKey("max_step")) {
      max_step = time_config.template get<DurationQuantity>("time_step_max");
    }

    using State = typename Model::State;
    auto stepper = SimpleAdaptiveStepper<State, TimeQuantity, DurationQuantity>{
      decrease_factor, increase_factor, min_step, max_step
    };

    // setup writer
    std::function<void(const State&)> output_writter;
    const auto& writer_config = model_config.sub("writer");
    const auto& writer_type = writer_config.get("type", "vtk");
    if (writer_type == "vtk") {
      auto file = model_config.get("writer.file_path", "");
      output_writter = [&](const auto& state) {
        if (not file.empty()) {
          model->write(state, file, true);
        }
      };
    }

    auto in = model->make_state(md_grid_ptr, model_config);
    in->time = begin_time;

    // interpolate initial conditions into the model state
    model->interpolate(*in, model->make_initial(*md_grid_ptr, model_config));
    output_writter(*in); // write initial condition

    auto step_operator = model->make_step_operator(*in, model_config);

    stepper.evolve(*step_operator, *in, *in, initial_step, end_time, output_writter).or_throw();

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
  return end_code;
}
