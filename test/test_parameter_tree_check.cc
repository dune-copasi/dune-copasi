#ifdef HAVE_CONFIG_H
#include "config.h"
#endif


#include <dune/copasi/parameter_tree_check.hh>

#include <dune/common/exceptions.hh>
#include <dune/common/parametertreeparser.hh>

int
main(int argc, char** argv)
{
  bool failed = false;

  try {

    Dune::ParameterTree param_a, param_b, param_c;

    param_a["a"] = "1";
    param_a["A.a"] = "11";
    param_b["b"] = "2";
    param_c["C.c.a"] = "333";

    if (not Dune::Copasi::eq(param_a,param_a))
      DUNE_THROW(Dune::RangeError,"eq operator is wrong");

    if (Dune::Copasi::eq(param_a,param_b))
      DUNE_THROW(Dune::RangeError,"eq operator is wrong");

    Dune::Copasi::add(param_a,param_b);

    param_b["a"] = "1";
    if (Dune::Copasi::eq(param_a,param_b))
      DUNE_THROW(Dune::RangeError,"eq operator is wrong");

    param_b["A.a"] = "11";
    if (not Dune::Copasi::eq(param_a,param_b))
      DUNE_THROW(Dune::RangeError,"eq operator is wrong");

    Dune::Copasi::diff(param_a,param_b);
    if (not Dune::Copasi::eq(param_a,{}))
      DUNE_THROW(Dune::RangeError,"diff operator is wrong");

    Dune::ParameterTree param_c_copy(param_c);
    Dune::Copasi::add(param_c,param_b);
    Dune::Copasi::add(param_a,param_b);
    Dune::Copasi::diff(param_c,param_a);

    if (not Dune::Copasi::eq(param_c,param_c_copy))
      DUNE_THROW(Dune::RangeError,"diff operator is wrong");

    const std::string config_filename = argv[1];
    Dune::ParameterTree config;
    Dune::ParameterTreeParser ptreeparser;
    ptreeparser.readINITree(config_filename, config);

    auto config_model = Dune::Copasi::get_model(config);
    std::cout << "*******************Valid configuration model*******************" << std::endl;
    config_model.report();

    std::cout << "***************************Diff***************************" << std::endl;
    Dune::Copasi::diff(config,config_model);
    config.report();

    return failed;
  } catch (Dune::Exception& e) {
    std::cerr << "Dune reported error: " << e << std::endl;
  } catch (...) {
    std::cerr << "Unknown exception thrown!" << std::endl;
  }
}
