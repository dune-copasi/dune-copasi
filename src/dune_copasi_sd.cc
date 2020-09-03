#ifdef HAVE_CONFIG_H
#include "config.h"
#endif

#include <dune/copasi/common/enum.hh>
#include <dune/copasi/common/stepper.hh>
#include <dune/copasi/grid/mark_stripes.hh>
#include <dune/copasi/grid/multidomain_gmsh_reader.hh>
#include <dune/copasi/model/diffusion_reaction.hh>

#include <dune/grid/multidomaingrid.hh>

#include <dune/grid/io/file/gmshreader.hh>
#include <dune/grid/uggrid.hh>

#include <dune/logging/logging.hh>

#include <dune/common/exceptions.hh>
#include <dune/common/parallel/mpihelper.hh>
#include <dune/common/parametertree.hh>
#include <dune/common/parametertreeparser.hh>

#include <ctime>

int
main(int argc, char** argv)
{
  auto stime_c = std::chrono::system_clock::now();
  int end_code = 0;

  // initialize mpi
  auto& mpi_helper = Dune::MPIHelper::instance(argc, argv);

  // Read and parse ini file
  if (argc != 2)
    DUNE_THROW(Dune::IOError, "Wrong number of arguments");
  const std::string config_filename = argv[1];

  Dune::ParameterTree config;
  Dune::ParameterTreeParser::readINITree(config_filename, config);

  // initialize loggers
  auto comm = mpi_helper.getCollectiveCommunication();
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

    // Grid setup
    constexpr int dim = 2;
    using HostGrid = Dune::UGGrid<dim>;
    using MDGTraits = Dune::mdgrid::DynamicSubDomainCountTraits<dim, 1>;
    using MDGrid = Dune::mdgrid::MultiDomainGrid<HostGrid, MDGTraits>;

    auto& grid_config = config.sub("grid",true);
    auto level = grid_config.get<int>("initial_level", 0);

    auto grid_file = grid_config.get<std::string>("file");

    auto [md_grid_ptr, host_grid_ptr] =
      Dune::Copasi::MultiDomainGmshReader<MDGrid>::read(grid_file);

    using Grid = typename MDGrid::SubDomainGrid;
    using GridView = typename Grid::Traits::LeafGridView;

    auto grid_log = Dune::Logging::Logging::componentLogger({}, "grid");
    grid_log.detail("Applying refinement of level: {}"_fmt, level);

    for (int i = 1; i <= level; i++) {
      Dune::Copasi::mark_stripes(*host_grid_ptr);
      md_grid_ptr->preAdapt();
      md_grid_ptr->adapt();
      md_grid_ptr->postAdapt();
    }

    auto& model_config = config.sub("model",true);
    auto& compartments_map = model_config.sub("compartments",true);
    int order = model_config.get<int>("order");

    // TODO: Use OS for different domains when is ready
    if (compartments_map.getValueKeys().size() != 1)
      DUNE_THROW(
        Dune::NotImplemented,
        "Multiple compartments per model are not allowed in this executable");

    // get subdomain grid as a shared pointer
    auto compartment = compartments_map.getValueKeys().front();
    int domain = compartments_map.template get<int>(compartment);
    std::shared_ptr<Grid> grid_ptr =
      Dune::stackobject_to_shared_ptr(md_grid_ptr->subDomain(domain));

    // create time stepper
    auto timestep_config = model_config.sub("time_stepping",true);
    auto end_time = timestep_config.template get<double>("end");
    auto initial_step = timestep_config.template get<double>("initial_step");
    auto stepper = Dune::Copasi::make_default_stepper(timestep_config);

    auto file = model_config.get("writer.file_path", "");
    auto write_output = [=](const auto& state) {
      if (not file.empty())
        state.write(file,true);
    };

    if (order == 0) {
      constexpr int Order = 0;
      using ModelTraits =
        Dune::Copasi::ModelPkDiffusionReactionTraits<Grid, GridView, Order>;
      Dune::Copasi::ModelDiffusionReaction<ModelTraits> model(grid_ptr,
                                                              model_config);
      write_output(model.state()); // write initial condition
      stepper.evolve(model, initial_step, end_time, write_output);
    } else if (order == 1) {
      constexpr int Order = 1;
      using ModelTraits =
        Dune::Copasi::ModelP0PkDiffusionReactionTraits<Grid, GridView, Order>;
      Dune::Copasi::ModelDiffusionReaction<ModelTraits> model(grid_ptr,
                                                              model_config);
      write_output(model.state()); // write initial condition
      stepper.evolve(model, initial_step, end_time, write_output);
    } else if (order == 2) {
      constexpr int Order = 2;
      using ModelTraits =
        Dune::Copasi::ModelP0PkDiffusionReactionTraits<Grid, GridView, Order>;
      Dune::Copasi::ModelDiffusionReaction<ModelTraits> model(grid_ptr,
                                                              model_config);
      write_output(model.state()); // write initial condition
      stepper.evolve(model, initial_step, end_time, write_output);
    } else {
      DUNE_THROW(Dune::IOError,
                 "Finite element order "
                   << order << " is not supported by dune-copasi(sd)");
    }
  } catch (Dune::Exception& e) {
    log.error("Dune reported error:"_fmt);
    log.error(2,"{}"_fmt, e.what());
    end_code = 1;
  } catch (std::exception& e) {
    log.error("C++ reported error:"_fmt);
    log.error(2,"{}"_fmt, e.what());
    end_code = 1;
  } catch (...) {
    log.error("Unknown exception thrown!"_fmt);
    end_code = 1;
  }

  if (end_code)
    log.notice("dune-copasi(sd) finished with some errors :("_fmt);
  else
    log.notice("dune-copasi(sd) successfully finished :)"_fmt);
  return end_code;
}
