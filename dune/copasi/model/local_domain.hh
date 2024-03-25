#ifndef DUNE_COPASI_MODEL_LOCAL_DOMAIN_HH
#define DUNE_COPASI_MODEL_LOCAL_DOMAIN_HH

#include <dune/common/fvector.hh>

#include <functional>
#include <string_view>
#include <string>
#include <vector>

namespace Dune::Copasi {

template<std::size_t dim>
struct LocalDomain
{
  LocalDomain() = default;
  LocalDomain(const LocalDomain&) = delete;
  LocalDomain(LocalDomain&&) = delete;

  LocalDomain& operator=(const LocalDomain&) = delete;
  LocalDomain& operator=(LocalDomain&&) = delete;

  virtual ~LocalDomain() {}

  virtual void forEachValue(Codim<0>,
                            const std::function<void(std::string_view,
                                                     const Dune::FieldVector<double, 1>&,
                                                     const Dune::FieldVector<double, dim>&)>&) const
  {
  }
  virtual void forEachValue(Codim<1>,
                            const std::function<void(std::string_view,
                                                     const FieldVector<double, 1>&,
                                                     const FieldVector<double, dim - 1>&)>&) const
  {
  }

  // local data related to general position metrics
  FieldVector<double, dim> position;
  FieldVector<double, dim> normal;
  double time = 0.;
  double entity_volume = 0.;
  double in_volume = 0;
  double in_boundary = 0;
  double in_skeleton = 0;

  // local data related to the cell element and intended to be used by the CellData
  std::vector<double> cell_values;  // data addresses used by the parser
  std::vector<bool> cell_mask; // masks whether data is valid
  std::vector<std::string> cell_keys;  // data member names
};

} // namespace Dune::Copasi

#endif // DUNE_COPASI_MODEL_LOCAL_DOMAIN_HH
