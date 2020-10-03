---
id: ini_file
title: Ini File
sidebar_label: Ini File
---

import Tabs from '@theme/Tabs';
import TabItem from '@theme/TabItem';

The [INI file](https://en.wikipedia.org/wiki/INI_file) creates a
[parameter tree](param_tree.md) to initializate and run the models.
The one we use follows the `DUNE` convention. In short, the data is composed
of a *key–value pairs* on the form `key = value`.

### Comments

Lines starting with a hash `#` will be ignored

```ini
# This is a comment
```

### Value

The `value` on the `key=value` pair may be of different kinds. The following is
list of possible accepted types

| Type | Description |
| ---- | ----------- |
| `string`      | Sequence of characters terminated an [end of line](https://en.wikipedia.org/wiki/Newline)
| `float`       | Decimal number in floating point representation
| `integer`     | Integer number
| `path`        | String that represent filesystem paths absolute or relative to a program
| `math_expr`   | [Mathematical Expression](math_expr.md)

Examples

```ini
string = I am a good looking string
float = 1e-2
integer = 1
path = /path/to/directory/
math_expr = pi*sin(x)*sin(y)
```

### Group

A `key` is a `string` which may be contained on a common
group of parameters by preceding the line with the group name between
brackets, or/and by preceding the `key` by the group name and a dot. Notice that
the two forms may be combined.

As an example, the three following ini files, are equivalent
<Tabs
  defaultValue="form1"
  values={[
      {label: 'No grouping', value: 'form1', },
      {label: 'Grouping', value: 'form2', },
      {label: 'Nested grouping', value: 'form3', },
    ]
  }>

  <TabItem value="form1">

```ini
# No preceding section
section.subsection.first = 1.0
section.subsection.second = 2.0
```

  </TabItem>
  <TabItem value="form2">

```ini
[section]
subsection.first = 1.0
subsection.second = 2.0
```

  </TabItem>


  <TabItem value="form3">

```ini
[section.subsection]
first = 1.0
second = 2.0
```

  </TabItem>
</Tabs>
