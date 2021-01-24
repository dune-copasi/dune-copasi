#ifdef HAVE_CONFIG_H
#include "config.h"
#endif

#include <dune/copasi/model/base.hh>

#include <dune/common/exceptions.hh>
#include <dune/common/float_cmp.hh>

namespace Dune::Copasi {

ModelBase::ModelBase(const Dune::ParameterTree& config)
  : _logger(Logging::Logging::componentLogger({}, "model"))
  , _adapt_policy(AdaptivityPolicy::None)
{
  _logger.trace("ModelBase constructed"_fmt);
}

ModelBase::~ModelBase()
{
  _logger.trace("ModelBase deconstructed"_fmt);
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

} // namespace Dune::Copasi
