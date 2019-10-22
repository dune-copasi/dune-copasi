#ifdef HAVE_CONFIG_H
#include "config.h"
#endif

#include <dune/copasi/model/base.hh>

#include <dune/common/exceptions.hh>
#include <dune/common/float_cmp.hh>

namespace Dune::Copasi {

ModelBase::ModelBase(const Dune::ParameterTree& config)
  : _logger(Logging::Logging::componentLogger(config, "model"))
  , _adapt_policy(AdaptivityPolicy::None)
  , _begin_time(config.template get<double>("begin_time"))
  , _end_time(config.template get<double>("end_time"))
  , _current_time(config.template get<double>("current_time", 0.))
{
  _logger.debug("ModelBase constructed"_fmt);
}

ModelBase::~ModelBase()
{
  _logger.debug("ModelBase deconstructed"_fmt);
}

void
ModelBase::set_policy(AdaptivityPolicy adapt_policy)
{
  _adapt_policy = adapt_policy;
}

AdaptivityPolicy
ModelBase::adaptivity_policy() const
{
  return _adapt_policy;
}

void
ModelBase::mark_grid()
{
  DUNE_THROW(Dune::NotImplemented, "'mark_grid' not implemented!");
}

void
ModelBase::pre_adapt_grid()
{}

void
ModelBase::adapt_grid()
{
  if (_adapt_policy != AdaptivityPolicy::None) {
    DUNE_THROW(Dune::NotImplemented, "'adapt_grid' not implemented");
  } else {
    DUNE_THROW(Dune::InvalidStateException, "Invalid adaptation policy");
  }
}

void
ModelBase::post_adapt_grid()
{}

const double&
ModelBase::begin_time() const
{
  return _end_time;
}

double&
ModelBase::begin_time()
{
  return _end_time;
}

const double&
ModelBase::end_time() const
{
  return _end_time;
}

double&
ModelBase::end_time()
{
  return _end_time;
}

const double&
ModelBase::current_time() const
{
  return _current_time;
}

double&
ModelBase::current_time()
{
  return _current_time;
}

void
ModelBase::run()
{
  auto do_step = [&]() {
    return Dune::FloatCmp::lt(current_time(), end_time());
  };

  while (do_step()) {
    step();

    if (adaptivity_policy() != AdaptivityPolicy::None)
      if (do_step()) {
        mark_grid();
        pre_adapt_grid();
        adapt_grid();
        post_adapt_grid();
      }
  }
}

} // namespace Dune::Copasi