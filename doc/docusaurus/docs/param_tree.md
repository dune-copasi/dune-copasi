---
id: param_tree
title: Configuration Options
sidebar_label: Config Options
---

import ConfigOpts from './ConfigOpts.js'


## Key-Value Pairs

The simulator may be configured by passing a collection of `key=value` pair that change the behavior of the solver.
These options may be passed via a configuration file with the [`DUNE` ini file convention](ini_file.md) or by the command line `--{key}={value}`. Note that command line arguments take precedence over the config file parameters.

::::info Placeholders
Whenever a key in this documentation contains curly braces `{...}` it means that it is a placeholder for a user defined key.

:::tip Example
The option `compartments.domain.expression` means that the `{cmp}` placeholder in `compartments.{cmp}.expression` is being replaced with the `domain` keyword and is a valid option.
:::
::::


::::info Multiple keys
* Whenever a key in this documentation contains square brackets `[...]` it means that any of the options expanded by the contents separated by `|` is a valid option.

:::tip Example
The option `model.time_[begin|end]` means that both `model.time_begin` and `model.time_end` are valid options.
:::
::::

Below you can find an expandable list of the parameters that may be used in the solver:

<ConfigOpts/>

<!-- ## Example

:::tip Example
An example of an ini file with all the parameters required by
`dune-copasi` for the [Gray-Scott model](http://mrob.com/pub/comp/xmorphia/F420/F420-k610.html)
is

```ini title="config.ini"
# TODO
```
::: -->
