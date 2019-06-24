#ifdef HAVE_CONFIG_H
# include "config.h"
#endif

#include <dune/copasi/model_diffusion_reaction.hh>
#include <dune/copasi/model_diffusion_reaction.cc>

#include <dune/grid/utility/structuredgridfactory.hh>

#include <dune/logging/logging.hh>

#include <dune/common/exceptions.hh>
#include <dune/common/parametertree.hh>
#include <dune/common/parametertreeparser.hh>
#include <dune/common/parallel/mpihelper.hh>

#include <iostream>

// diffusion-reaction problem description
template<typename Number>
class Parameterization
{
  Number time;
  Number lambda;
public:
  //! Export the number of components
  enum {components=3};

  //! export the number type given as template parameter
  typedef Number value_type;

  //! return type for vector with components entries
  typedef std::vector<Number> Vector;

  //! Type for entries of sparse Jacobian in point block format
  struct JacobianEntry {
    int row, col;
    Number value;
    JacobianEntry () {}
    JacobianEntry (int r, int c, Number v)
      : row(r), col(c), value(v)
    {}
  };

  //! Export the number of nonzero components in the Jacobian
  enum {nonzeroes=6};

  //! return type for Jacobian matrix
  typedef std::vector<JacobianEntry> Jacobian;

  //! Constructor without arg sets nonlinear term to zero
  Parameterization () : lambda(1.0) {}

  //! store time for subsequent evaluations
  void setTime (Number t_)
  {
    time = t_;
  }

  //! reaction mechanism
  template<typename E, typename X, typename U>
  Vector f (const E& e, const X& x, U& u) const
  {
    Vector rv(components);
    rv[0] = +u[0]*u[1]*lambda;
    rv[1] = +u[0]*u[1]*lambda;
    rv[2] = -u[0]*u[1]*lambda;
    return rv;
  }

  //! Jacobian of reaction mechanism
  template<typename E, typename X, typename U>
  Jacobian nablaf (const E& e, const X& x, U& u) const
  {
    Jacobian rv(nonzeroes);
    rv[0] = JacobianEntry(0,0,+u[1]*lambda);
    rv[1] = JacobianEntry(0,1,+u[0]*lambda);
    rv[2] = JacobianEntry(1,0,+u[1]*lambda);
    rv[3] = JacobianEntry(1,1,+u[0]*lambda);
    rv[4] = JacobianEntry(2,0,-u[1]*lambda);
    rv[5] = JacobianEntry(2,1,-u[0]*lambda);
    return rv;
  }

  //! Diffusion coefficients
  template<typename E, typename X>
  Vector D (const E& e, const X& x) const
  {
    Vector rv(components);
    rv[0] = 1.0;
    rv[1] = 2.0;
    rv[2] = 0.02;
    return rv;
  }

  //! boundary condition type function (true = Dirichlet)
  template<typename I, typename X>
  bool b (const I& i, const X& x) const
  {
    return false;
  }

  //! Dirichlet extension and initial condition
  template<typename E, typename X>
  Vector g (const E& e, const X& x) const
  {
    auto global = e.geometry().global(x);
    auto r = global.two_norm();
    Vector rv(components);
    if (r>1.5)
      rv[0] = 1.0;
    else if (r<1.1)
      rv[0] = 0.0;
    else
      rv[0] = (r-1.1)/0.4*(r-1.1)/0.4;

    if (r>1.0)
      rv[1] = 0.0;
    else if (r<0.6)
      rv[1] = 1.0;
    else
      rv[1] = (1.0-r)/0.4*(1.0-r)/0.4;

    rv[2] = 0.0;
    return rv;
  }
};

int main(int argc, char** argv)
{

  try{
    // initialize mpi
    auto& mpi_helper = Dune::MPIHelper::instance(argc, argv);
    auto comm = mpi_helper.getCollectiveCommunication();

    // Read and parse ini file
    if (argc!=2)
      DUNE_THROW(Dune::IOError, "Wrong number of arguments");
    const std::string config_filename = argv[1];

    Dune::ParameterTree config;
    Dune::ParameterTreeParser ptreeparser;
    ptreeparser.readINITree(config_filename, config);

    // initialize loggers
    Dune::Logging::Logging::init(comm,config.sub("logging"));

    using namespace Dune::Literals;
    auto log = Dune::Logging::Logging::logger(config);
    log.notice("Starting dune-copasi"_fmt);


    // test the current code... 

    // create a grid
    constexpr int dim = 2;
    using Grid = Dune::UGGrid<dim>;
    using Domain = Dune::FieldVector<double,2>;

    auto& grid_config = config.sub("grid");
    auto level = grid_config.get<int>("initial_level",0);
    auto upper_right = grid_config.get<Domain>("extensions",{1.,1.});
    auto elements = grid_config.get<std::array<uint, 2>>("cells",{10,10});

    log.info("Creating a rectangular grid in {}D"_fmt, dim);
    log.debug("Grid extensions: {}"_fmt, upper_right);
    log.debug("Grid cells: {}"_fmt, elements);

    Domain origin(.0);

    std::shared_ptr<Grid> grid;
    grid = Dune::StructuredGridFactory<Grid>::createCubeGrid(origin,
                                                             upper_right,
                                                             elements);

    log.debug("Applying global refinement of level: {}"_fmt, level);
    grid->globalRefine(level);

    // instantiate a model
    auto& model_config = config.sub("model");
    using Param = Parameterization<double>;
    Dune::Copasi::ModelDiffusionReaction<3,Param> model(grid,model_config);

    model.run();

    return 0;
  }
  catch (Dune::Exception &e){
    std::cerr << "Dune reported error: " << e << std::endl;
  }
  catch (...){
    std::cerr << "Unknown exception thrown!" << std::endl;
  }
}
