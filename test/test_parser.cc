#ifdef HAVE_DUNE_COPASI_CONFIG_H
#include <dune/copasi/config.h>
#endif

#include <dune/copasi/parser/factory.hh>

#include <dune/common/exceptions.hh>

#include <cassert>

auto test_parser(Dune::Copasi::Parser& parser){
  double x = 20, y = 30, z = 10;
  parser.define_variable("x", &x);
  parser.define_variable("y", &y);
  parser.define_variable("z", &z);
  parser.define_constant("c",10);
  int i = 2;
  parser.define_function("f0",[=](){return i;});
  parser.define_function("f1",[](auto x){return x*2;});
  parser.define_function("f2",[&](auto x, auto y){return std::pow(x,i)+std::pow(y,i);});
  i = 3;

  parser.compile();

  return parser.eval();
}

auto test_parsers(std::string expression){

  if(HAVE_MUPARSER){
    auto parser = make_parser(expression, Dune::Copasi::ParserType::MuParser);
    test_parser(*parser);
  }

  if(HAVE_EXPRTK){
    auto parser = make_parser(expression, Dune::Copasi::ParserType::ExprTk);
    test_parser(*parser);
  }

  if(HAVE_SYMENGINE){
    auto parser = make_parser(expression, Dune::Copasi::ParserType::SymEngine);
    test_parser(*parser);

    auto sbml_parser = make_parser(expression, Dune::Copasi::ParserType::SymEngineSBML);
    test_parser(*sbml_parser);
  }
}

int
main(int argc, char** argv)
{
  try {
    bool passed = true;

    using namespace Dune::Copasi;

    test_parsers("x*y + z*c + f0()*f1(x) + f2(2*x, 2*y)");

    return passed;
  } catch (Dune::Exception& e) {
    std::cerr << "Dune reported error: " << e << std::endl;
  } catch (...) {
    std::cerr << "Unknown exception thrown!" << std::endl;
  }
}
