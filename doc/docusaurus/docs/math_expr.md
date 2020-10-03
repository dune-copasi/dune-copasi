---
id: math_expr
title: Mathematical Expressions
sidebar_label: Math Expressions
---

Mathematical expressions are permitted on different sections of the
[Parameter Tree](param_tree.md). Such expressions may contain typical math
operations such as `+`, `-`, `*`, `^`, `sin`, `tan`, etc. The complete list of
possibilities is found [here](https://beltoforion.de/en/muparser/features.php).
That is, the `muParser` built in operators.

Furthermore, we provide extra variables/functions that may be included on the
math expressions.

## Variables
### Context-Independent

The following list of variables are unconditionally available when defining
expressions:

| Variable | Description |
| -------- | ----------- |
| `t`      | Time
| `x`      | First dimension in the cartesian coordinate system
| `y`      | Second dimension in the cartesian coordinate system
| `z`      | Third dimension in the cartesian coordinate system (only in 3D)
| `dim`    | Dimension of the euclidian space
| `pi`     | Approximation to the number π

### Context-Dependent

Additionally, there are other variables that are conditionally available to
define the expression. These variables depend on the context of its definition.
An example of it, is the reaction expressions within a compartment. If the
section `[model.nucleous.reaction]` defines expresions `u`, `v`, and `w`,
then, they all will be available to each other. E.g.

```ini
[model.nucleous.reaction]
u = u*v*w
v = u*v*w
w = u*v*w
```

## Functions

### Context-Dependent

There are also functions that may be conditionally available depending on the
context. In such case, the function will have a custom name and fixed number of
arguments and will be available for use on the mathematical expression.

Currently we only provide these for defining TIFF input images as 2D functions
for initial conditions (See [Input Data](input_data.md) section).
