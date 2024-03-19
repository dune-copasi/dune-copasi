---
id: input_data
title: Input Data
sidebar_label: Input Data
---

## TIFF Images

Each `key=value` pair on the `[model.data]` section defines a
[Context-Dependet-Function](math_expr.md) with 2 arguments named as the `key`
identifier. The corresponding `value` should contain the directory path for a
16-bit grayscale TIFF file. These function are permitted to be used on the
initialization of spatial variables for the compartments
(i.e. `[model.<compartment>.initial]`). For example, on the initialization of
variables for the `nucleus` compartment, the following math expressions for `u`
and `v` are valid.

```ini
[model.data]
fnc_1 = file_1.tiff
fnc_2 = file_2.tiff

[model.nucleus.initial]
u = fnc_1(x, y) + fnc_2( -x, -y)
v = 30 * fnc_1(x, t)
```

:::caution Out of domain behavior
If the resulting function is evaluated outside the domain of the TIFF image, the
arguments will be [clamped](https://en.wikipedia.org/wiki/Clamping_(graphics)) to
the nearest valid point in the domain.
:::
