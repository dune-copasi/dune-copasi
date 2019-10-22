#ifndef DUNE_COPASI_PARAMETER_PARSER_HH
#define DUNE_COPASI_PARAMETER_PARSER_HH

#include <functional>

#include <dune/common/parametertree.hh>

namespace Dune::Copasi {

bool
eq(const ParameterTree& config_l, const ParameterTree& config_r)
{
  // absolutely not obtimal for big trees!
  for (auto&& key : config_r.getValueKeys())
    if (not config_l.hasKey(key))
      return false;

  for (auto&& key : config_l.getValueKeys())
    if (not config_r.hasKey(key))
      return false;

  for (auto&& section : config_r.getSubKeys())
    if (not eq(config_l.sub(section), config_r.sub(section)))
      return false;

  for (auto&& section : config_l.getSubKeys())
    if (not eq(config_l.sub(section), config_r.sub(section)))
      return false;

  return true;
}

// merge parameter trees
// config_l += config_r
void
add(ParameterTree& config_l,
    const ParameterTree& config_r,
    bool trow_on_override = true)
{
  for (auto&& key : config_r.getValueKeys()) {
    if (config_l.hasKey(key))
      DUNE_THROW(RangeError,
                 "config file addition failed due to duplucated key: " << key);
    config_l[key] = config_r[key];
  }
  for (auto&& section : config_r.getSubKeys())
    add(config_l.sub(section), config_r.sub(section));
}

// diff parameter trees
// config_l -= config_r
void
diff(ParameterTree& config_l, const ParameterTree& config_r)
{
  ParameterTree config_l_copy(config_l);
  config_l = {};

  for (auto&& key : config_l_copy.getValueKeys())
    if (not config_r.hasKey(key))
      config_l[key] = config_l_copy[key];

  for (auto&& section : config_l_copy.getSubKeys()) {
    diff(config_l_copy.sub(section), config_r.sub(section));
    if (config_l_copy.sub(section).getValueKeys().size() > 0 or
        config_l_copy.sub(section).getSubKeys().size() >
          0) // remove empty sections
      config_l.sub(section) = config_l_copy.sub(section);
  }
}

// get keys in vector
ParameterTree
get_keys(const ParameterTree& config, const std::vector<std::string>& keys)
{
  ParameterTree new_config;
  for (auto&& key : keys)
    new_config[key] = config[key];
  return new_config;
}

// get all keys (no section)
ParameterTree
get_keys(const ParameterTree& config)
{
  return get_keys(config, config.getValueKeys());
}

// get sections in vector
ParameterTree
get_sections(
  const ParameterTree& config,
  const std::vector<std::function<ParameterTree(const ParameterTree&)>>&
    sections)
{
  ParameterTree new_config;
  for (auto&& section : sections)
    add(new_config, section(config));
  return new_config;
}

ParameterTree
get_grid(const ParameterTree& config)
{
  std::vector<std::string> grid_keys{ "file", "initial_level" };
  auto grid_config = get_keys(config.sub("grid"), grid_keys);
  // todo : check file is valid
  // todo : check initial_level is > 0

  ParameterTree new_config;
  new_config.sub("grid") = grid_config;
  return new_config;
}

ParameterTree
get_initial(const ParameterTree& config)
{
  auto initial_config = get_keys(config.sub("initial"));
  // todo : check that they are math expressions
  // todo : keys are ordered

  ParameterTree new_config;
  new_config.sub("initial") = initial_config;
  return new_config;
}

// check_var_consistency with respect to initial sections
ParameterTree
get_diffusion(const ParameterTree& config, bool check_var_consistency = true)
{
  ParameterTree diffusion_config;
  if (check_var_consistency) {
    auto initial_config = get_keys(config.sub("initial"));
    diffusion_config =
      get_keys(config.sub("diffusion"), initial_config.getValueKeys());
  } else {
    diffusion_config = get_keys(config.sub("diffusion"));
  }
  // todo : check that they are math expressions
  // todo : keys are ordered

  ParameterTree new_config;
  new_config.sub("diffusion") = diffusion_config;
  return new_config;
}

// check_var_consistency with respect to initial sections
ParameterTree
get_operator(const ParameterTree& config, bool check_var_consistency = true)
{
  ParameterTree operator_config;
  if (check_var_consistency) {
    auto initial_config = get_keys(config.sub("initial"));
    operator_config =
      get_keys(config.sub("operator"), initial_config.getValueKeys());
  } else {
    operator_config = get_keys(config.sub("operator"));
  }
  // todo : check that they are signed integers
  // todo : keys are ordered

  ParameterTree new_config;
  new_config.sub("operator") = operator_config;
  return new_config;
}

ParameterTree
get_jacobian(const ParameterTree& config,
             std::vector<std::string> base_variables)
{
  auto jacobian_config = get_keys(config.sub("jacobian"));
  std::size_t jac_size = jacobian_config.getValueKeys().size();
  std::size_t base_size = base_variables.size();

  // it's important that jacobian has the right size and is ordered
  if (jac_size != base_size * base_size)
    DUNE_THROW(RangeError, "jacobian section has wrong size");

  std::size_t count(0);
  for (auto&& var_i : base_variables) {
    for (auto&& var_j : base_variables) {
      std::string jac_key = jacobian_config.getValueKeys()[count];
      auto found_i = jac_key.find(var_i);
      if (found_i == std::string::npos)
        DUNE_THROW(RangeError,
                   "Jacobian key '" << jac_key
                                    << "' does not contain its base key i '"
                                    << var_i << "'");
      auto found_j = jac_key.find(var_j);
      if (found_j == std::string::npos)
        DUNE_THROW(RangeError,
                   "Jacobian key '" << jac_key
                                    << "' does not contain its base key j '"
                                    << var_j << "'");
      count++;
    }
  }

  // todo : check that they are math expressions
  // todo : keys are ordered

  ParameterTree new_config;
  new_config.sub("jacobian") = jacobian_config;
  return new_config;
}

ParameterTree
get_reaction(const ParameterTree& config, bool check_var_consistency = true)
{
  ParameterTree reaction_config;
  if (check_var_consistency) {
    auto initial_config = get_keys(config.sub("initial"));
    reaction_config =
      get_keys(config.sub("reaction"), initial_config.getValueKeys());
  } else {
    reaction_config = get_keys(config.sub("reaction"));
  }

  auto jacobian_config =
    get_jacobian(config.sub("reaction"), reaction_config.getValueKeys());
  add(reaction_config, jacobian_config);
  // todo : check that they are math expressions
  // todo : keys are ordered

  ParameterTree new_config;
  new_config.sub("reaction") = reaction_config;
  return new_config;
}

ParameterTree
get_compartment_i(const ParameterTree& config, const std::string& name)
{
  std::vector<std::function<ParameterTree(const ParameterTree&)>> sections;
  sections.push_back(get_initial);
  sections.push_back([](auto i) { return get_diffusion(i); });
  sections.push_back([](auto i) { return get_reaction(i); });
  sections.push_back([](auto i) { return get_operator(i); });
  auto compart_config = get_sections(config.sub(name), sections);

  ParameterTree new_config;
  new_config.sub(name) = compart_config;
  return new_config;
}

ParameterTree
get_compartments(const ParameterTree& config)
{
  auto compart_config = get_keys(config.sub("compartments"));
  // todo : check that they are valid keys

  ParameterTree new_config;
  new_config.sub("compartments") = compart_config;
  return new_config;
}

ParameterTree
get_model(const ParameterTree& config)
{
  auto compartments = get_compartments(config.sub("model"));

  std::vector<std::function<ParameterTree(const ParameterTree&)>> sections;
  for (auto&& compartment : compartments.sub("compartments").getValueKeys())
    sections.push_back(
      [&](auto ini) { return get_compartment_i(ini, compartment); });

  auto compartment_i = get_sections(config.sub("model"), sections);

  add(compartments, compartment_i);

  ParameterTree new_config;
  new_config.sub("model") = compartments;

  // bug in parameter tree. If keys go first, sub section is not well set!
  std::vector<std::string> model_keys{ "begin_time", "end_time", "time_step" };
  add(new_config.sub("model"), get_keys(config.sub("model"), model_keys));

  // todo : check that they are valid keys
  return new_config;
}

} // namespace Dune::Copasi

#endif // DUNE_COPASI_PARAMETER_PARSER_HH