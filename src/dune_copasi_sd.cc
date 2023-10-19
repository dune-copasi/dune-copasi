#ifdef HAVE_DUNE_COPASI_CONFIG_H
#include <dune/copasi/config.h>
#endif

#include <dune/copasi/common/stepper.hh>
#include <dune/copasi/grid/mark_stripes.hh>
#include <dune/copasi/grid/multidomain_gmsh_reader.hh>
#include <dune/copasi/model/diffusion_reaction.hh>
#ifndef DUNE_COPASI_SD_LIBRARY
#include <dune/copasi/model/diffusion_reaction.cc>
#endif

#include <dune/assembler/common/entityset_threadpool.hh>

#if HAVE_UNITS
#include <dune/assembler/quantity_grid/grid.hh>
#endif

#include <dune/grid/multidomaingrid.hh>

#include <dune/grid/io/file/gmshreader.hh>
#include <dune/grid/uggrid.hh>

#include <dune/logging/logging.hh>

#include <dune/common/exceptions.hh>
#include <dune/common/hybridutilities.hh>
#include <dune/common/parallel/mpihelper.hh>
#include <dune/common/parametertree.hh>
#include <dune/common/parametertreeparser.hh>

#if HAVE_UNITS
#include <units/isq/si/length.h>
#endif

#include <vector>
#include <string>
#include <iostream>
#include <chrono>


std::unique_ptr<perfetto::TracingSession> StartTracing() {
  // The trace config defines which types of data sources are enabled for
  // recording. In this example we just need the "track_event" data source,
  // which corresponds to the TRACE_EVENT trace points.
  perfetto::TraceConfig cfg;
  cfg.add_buffers()->set_size_kb(1024*100);
  auto* ds_cfg = cfg.add_data_sources()->mutable_config();
  ds_cfg->set_name("track_event");

  auto tracing_session = perfetto::Tracing::NewTrace();
  tracing_session->Setup(cfg);
  tracing_session->StartBlocking();
  return tracing_session;
}

void StopTracing(std::unique_ptr<perfetto::TracingSession> tracing_session) {
  // Make sure the last event is closed for this example.
  perfetto::TrackEvent::Flush();

  // Stop tracing and read the trace data.
  tracing_session->StopBlocking();
  std::vector<char> trace_data(tracing_session->ReadTraceBlocking());

  // Write the result into a file.
  // Note: To save memory with longer traces, you can tell Perfetto to write
  // directly into a file by passing a file descriptor into Setup() above.
  std::ofstream output;
  output.open("dune-copasi-sd.pftrace", std::ios::out | std::ios::binary);
  output.write(&trace_data[0], std::streamsize(trace_data.size()));
  output.close();
  PERFETTO_LOG("Trace written in dune-copasi-sd.pftrace file.");
}


static std::string string_help =
  "Usage: dune-copasi-sd CONFIG_FILE\n"
  "\n"
  "Execute numerical simulation of reaction-diffusion systems on single\n"
  "compartments. The CONFIG_FILE is a DUNE INI file with the\n"
  "parameters to perform the simulation.\n";

int
main(int argc, char** argv)
{
  perfetto::TracingInitArgs args;
  args.backends |= perfetto::kInProcessBackend;
  perfetto::Tracing::Initialize(args);
  perfetto::TrackEvent::Register();

  auto tracing_session = StartTracing();

  // Give a custom name for the traced process.
  perfetto::ProcessTrack process_track = perfetto::ProcessTrack::Current();
  perfetto::protos::gen::TrackDescriptor desc = process_track.Serialize();
  desc.mutable_process()->set_process_name("Main");
  perfetto::TrackEvent::SetTrackDescriptor(process_track, desc);

  TRACE_EVENT_BEGIN("dune", "Main");

  std::vector<std::string> cmd_line_args(argv, argv+argc);

  if (cmd_line_args.size() != 2) {
    std::cerr << string_help;
    return 1;
  }

  if (cmd_line_args[1] == "--help" or cmd_line_args[1] == "-h") {
    std::cout << string_help;
    return 0;
  }

  auto stime_c = std::chrono::system_clock::now();
  int end_code = 0;

  // initialize mpi
  auto& mpi_helper = Dune::MPIHelper::instance(argc, argv);

  // Read and parse ini file
  const std::string config_filename = cmd_line_args[1];
  Dune::ParameterTree config;
  Dune::ParameterTreeParser::readINITree(config_filename, config);

  // initialize loggers
  auto comm = mpi_helper.getCommunication();
  Dune::Logging::Logging::init(comm, config.sub("logging"));
  auto log = Dune::Logging::Logging::logger(config);
  using namespace Dune::Literals;

  try {
    std::time_t stime_t = std::chrono::system_clock::to_time_t(stime_c);
    std::tm stime_tm;
    std::memcpy(&stime_tm, std::localtime(&stime_t), sizeof(std::tm));
    auto stime_s = fmt::format("{:%a %F %T %Z}", stime_tm);
    log.notice("Starting dune-copasi(sd) at {}"_fmt, stime_s);
    log.info("Reading configuration file: '{}'"_fmt, config_filename);

    // detailed report of the input paramter tree
    std::stringstream ss;
    config.report(ss);
    log.detail(2, "----"_fmt);
    for (std::string line; std::getline(ss, line);)
      log.detail(2, "{}"_fmt, line);
    log.detail(2, "----"_fmt);

    std::size_t order = config.template get<std::size_t>("model.order");
    std::size_t dim = config.template get<std::size_t>("grid.dimension");

    if (dim != 2 and dim != 3)
      DUNE_THROW(Dune::IOError, "Only 2D and 3D grids are alloed!");

    // lambda that instantiates and evolves a multidomain model
    auto evolve_model = [](const auto& config, auto dim, auto order) {
      using namespace Dune::Copasi;

      using HostGrid = Dune::UGGrid<dim>;
#if HAVE_UNITS
      using DistanceRepresentation = double;
      using Distance = units::isq::si::length<units::isq::si::metre, DistanceRepresentation>;
      using QuantityGrid = Dune::Assembler::QuantityGrid<HostGrid,Distance>;
#else
      using QuantityGrid = HostGrid;
#endif
      using MDGTraits = Dune::mdgrid::DynamicSubDomainCountTraits<dim, 1>;
      using MDGrid = Dune::mdgrid::MultiDomainGrid<QuantityGrid, MDGTraits>;

      auto& grid_config = config.sub("grid", true);
      auto level = grid_config.template get<int>("initial_level");
      auto grid_file = grid_config.template get<std::string>("file");

      auto [md_grid_ptr, host_grid_ptr] =
        MultiDomainGmshReader<MDGrid>::read(grid_file);

      auto grid_log = Dune::Logging::Logging::componentLogger({}, "grid");
      grid_log.detail("Applying refinement of level: {}"_fmt, level);

      // for (int i = 1; i <= level; i++) {
      //   mark_stripes(*host_grid_ptr);
      //   md_grid_ptr->preAdapt();
      //   md_grid_ptr->adapt();
      //   md_grid_ptr->postAdapt();
      // }

      auto& model_config = config.sub("model", true);
      auto& compartments_map = model_config.sub("compartments", true);
      if (compartments_map.getValueKeys().size() != 1)
        DUNE_THROW(
          Dune::NotImplemented,
          "Multiple compartments per model are not allowed in this executable");

      // get subdomain grid as a shared pointer
      auto compartment = compartments_map.getValueKeys().front();
      int domain = compartments_map.template get<int>(compartment);
      auto grid_ptr =
        Dune::stackobject_to_shared_ptr(md_grid_ptr->subDomain(domain));

      using Grid = typename MDGrid::SubDomainGrid;
      using GridView = typename Grid::Traits::LeafGridView;
      using EntitySet = GridView;
      using Traits = ModelDiffusionPkReactionTraits<Grid, EntitySet, /*order*/ 1, true>;
      ModelDiffusionReaction<Traits> model(grid_ptr, model_config, EntitySet{grid_ptr->leafGridView()});

      // setup writer
      auto file = model_config.get("writer.file_path", "");
      auto write_output = [=](const auto& state) {
        if (not file.empty())
          state.write(file, true);
      };
      write_output(model.state()); // write initial condition

      // create time stepper
      using TimeQuantity = typename ModelDiffusionReaction<Traits>::TimeQuantity;
      auto timestep_config = model_config.sub("time_stepping", true);
      auto end_time = TimeQuantity{timestep_config.template get<double>("end")};
      auto initial_step = TimeQuantity{timestep_config.template get<double>("initial_step")};
      auto stepper = Dune::Copasi::make_default_stepper<TimeQuantity>(timestep_config);

      auto in = model.state();
      auto step_operator = model.make_step_operator(in, {});

      auto out = in;
      stepper.evolve(*step_operator, in, out, initial_step, end_time, write_output);
    };

    // maximum polynomial order instantiated
    constexpr auto max_order = 2;
    if (order > max_order)
      DUNE_THROW(Dune::IOError,
                 "Finite element order " << order
                                         << " is not supported by dune-copasi");

    // static loop that instantiates 2D and 3D polynomial orders but runs the
    // dynamic ones
    auto order_range = Dune::range(Dune::index_constant<max_order + 1>{});
    Dune::Hybrid::forEach(order_range, [&](auto static_order) {
      // instantiate models with 2D grids
      if (dim == 2 and order == static_order)
        evolve_model(config, Dune::Indices::_2, static_order);

      // instantiate models with 3D grids
      if (dim == 3 and order == static_order)
#ifdef DUNE_COPASI_COMPILE_3D
        evolve_model(config, Dune::Indices::_3, static_order);
#else
        DUNE_THROW(Dune::IOError,"dune_copasi was not built with 3D support");
#endif
    });

  } catch (Dune::Exception& e) {
    log.error("Dune reported error:"_fmt);
    log.error(2, "{}"_fmt, e.what());
    end_code = 1;
  } catch (std::exception& e) {
    log.error("C++ reported error:"_fmt);
    log.error(2, "{}"_fmt, e.what());
    end_code = 1;
  } catch (...) {
    log.error("Unknown exception thrown!"_fmt);
    end_code = 1;
  }

  if (end_code)
    log.notice("dune-copasi(sd) finished with some errors :("_fmt);
  else
    log.notice("dune-copasi(sd) successfully finished :)"_fmt);
  TRACE_EVENT_END("dune");
  StopTracing(std::move(tracing_session));
  return end_code;
}
