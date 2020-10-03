---
id: input_data
title: Input Data
sidebar_label: Input Data
---

## TIFF Images

Each `key=value` pair on the `[model.data]` section defines a
[Context-Dependet-Function](math_expr.md) with 2 arguments named as the `key`
indentifier. The corresponding `value` should contain the directory path for a
16-bit grayscale TIFF file. These function are permitted to be used on the
initialization of spacial variables for the compartments
(i.e. `[mode.<compartment>.initial]`). For example, on the initialization of
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

If the resulting function is evaluated outside the domain of the TIFF image, the
function will evaluate to `0`.
