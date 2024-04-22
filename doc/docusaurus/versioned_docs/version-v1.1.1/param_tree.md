---
id: param_tree
title: Parameter Tree
sidebar_label: Parameter Tree
---

## Reference Keys

A parameter tree is a collection of *key–value pairs* that may be grouped by
sections. We call it a *tree* because sections may contain also other sections.
Here, we follow the [`DUNE` ini file convention](ini_file.md) to refer to
sections and keywords. In addition, keywords between angle brackets (i.e.
`<key>`) are references to **multiple** user-defined keywords of the same kind.

### `[grid]`

The grid is entred by a [GMSH file][GMSH] and should be formed by triangles and
rectangles in 2D or tetraedra and hexahedra in 3D. It must define
*physical group*s, each referring to a different compartment.
Each *physical group* can only be formed the same collection of geometric types.

| Key | Type | Description |
| ----------------| -------------- | ----- |
| `dimension`     | `integer` | The dimension of the grid |
| `file`          | `path` | [GMSH v2 file][GMSH] directory absolute or relative to the executable |
| `initial_level` | `integer` | The refinement level to start the simulation with |

:::caution Bug in GMSH reader
The GMSH file should not contain the `$PhysicalGroup` section at the
beginning of the file (this is a known bug in the `dune-grid` GMSH reader).
:::

### `[model]`

This section is the one that contains most of the aspects for simulation.

| Key | Type | Description |
| -----------|-----| -------------- |
| `order`           | `integer` | Finite Element polynomial order. if `order==0` FV else CG |
| `[data]` | subsection | List of input files |
| `[writer]` | subsection | Options for vtk writer |
| `[time_stepping]` | subsection | Options for time-stepping |
| `[compartment]`   | subsection | List of compartmentes to simulate |
| `[<comp_k>]` | subsection | Options for equations per each compartment |

Whenever a compartment is composed by cubes (2D) or hexahedra (3D), then such
compartment is modeled using Finite Volumes. These compartments are
suitable to model membranes processes.

The subsection `[<comp_k>]` is a placeholder for the name of all
compartment given in `[compartment]`. See `[model.compartment]` subsection for
more information.

#### `[model.data]`

This section will contain all the input files for the initial conditions of the
model.

| Key | Type | Description |
| -----------|-----| -------------- |
| `<data_name>` | `path` | 16-bit grayscale TIFF file |

See [Input Data](input_data.md) section for more information about usage.

#### `[model.output]`

The model output is controlled by the following parameters:

| Key | Type | Description |
| -----------|-----| -------------- |
| `file_path` | `string` | Name used to output VTK files |

#### `[model.time_stepping]`

This section sets the settings for the evolution in time of the forward model.

| Key | Type | Description |
| -----------|-----| -------------- |
| `rk_method` | `rk_id` | Runge-Kutta method |
| `begin` | `float` | Initial time of the forward simulation |
| `end`   | `float` | Final time of the forward simulation |
| `initial_step` | `float` | Time step at the beginning of the simulation |
| `min_step` | `float` | Maximum time step allowed |
| `max_step` | `float` | Minimum time step allowed |
| `decrease_factor` | `float` | Time step decrease factor when time step fails |
| `increase_factor` | `float` | Time step increase factor when time step succeeds |
| `[newton]` | subsection | Parameters for the Newton method at each timestep |

The possible `rk_id`s can be found at the end of this document.

#### `[model.newton]`

| Key | Type | Description |
| -----------|-----| -------------- |
| `reduction` | `float` | Threshold for termination, relative to first linear defect
| `use_max_norm` | `bool` | Use the maximum norm as a stopping criterion. This helps loosen the tolerance when solving for stationary solutions of nonlinear time-dependent problems.
| `absolute_limit` | `float` | Threshold for termination, absolute to linear defect
| `min_linear_reduction` | `float` | The linear reduction will be determined as minimum of this and the one needed to achieve second order newton convergence
| `fixed_linear_reduction`   | `bool` | Whenever `true`, the linear reduction rate will always be fixed to `min_linear_reduction`
| `max_iterations` | `integer` | Maximum number of linear searches allowed
| `reassemble_threshold` | `float` | If the reduction drops below this the linear operator is reassembled to ensure convergence
| `keep_matrix` | `bool` | Whether the jacobian matrix should be kept across time steps
| `force_iteration` | `bool` | Enforce first iteration even if required reduction is achieved
| `[linear_search]` | subsection | Parameters for linear search |

#### `[model.newton.linear_search]`

| Key | Type | Description |
| -----------|-----| -------------- |
| `strategy` | `string` | `noLineSearch` or `hackbuschReusken` or `hackbuschReuskenAcceptBest`
| `max_iterations` | `integer` | Maximum number of line search iterations
| `damping_factor`   | `float` | Multiplier to line search parameter after each iteration
| `accept_best`   | `bool` | Accept the best line search parameter if there was any improvement, even if the convergence criterion was not

#### `[model.compartment]`

The compartment section is filled with `key=value` pairs that assign a physical
group id (`value`) to an compartment name (`key`).

| Key | Type | Description |
| -----------|-----| -------------- |
| `<comp_k>` | `integer` | Physical group id assigned to `<comp_k>` subdomain |

Each compartment id maps to a *physical group* in the gmsh identifiers.
Although the gmsh format allows you to name such physical groups, we have no
support for it and identifiers must be used. Notice that `dune-copasi` uses
0-based indices while `gmsh` uses 1-based indices. In other words,
`<gmsh_physical_group> = <dune_copasi_compartment> + 1`.

:::tip Example
Let's say that there is two *physical groups* in our gmsh file
and we are going to name them as `nucleus` and `cytoplasm` compartments:

```ini
[model.compartments]
 # nucleus corresponds to the physical group 1 in the gmsh file
nucleus  = 0
 # cytoplasm corresponds to the physical group 2 in the gmsh file
cytoplasm = 1

# They then become available as new compartment subsections
[model.nucleus]
# Parameters for the nucleus compartment

[model.cytoplasm]
# Parameters for the nucleus compartment
```
:::

#### `[model.<comp_k>]`

Each compartment will define its own initial conditions,
its diffusion-reaction system, etc. A set of *variables* is be assigned to each
compartment. The amount *variables* may be different on each compartment, but
the same namings must be used on all the `[model.<comp_k>]` subsections.

| Key | Type | Description |
| -----------|-----| -------------- |
| `[initial]` | subsection | List of *variables* and their initial conditions |
| `[diffusion]` | subsection | List of *variables* and their diffusion coefficient |
| `[reaction]` | subsection | List of *variables* and their reaction networks |
| `[boundary]` | subsection | Definition of boundary fluxes on this `<comp_k>` |

##### `[model.<comp_k>.initial]`

The `initial` subsection allows the initialization for each of the *variables*.
Expressions in this section may contain [Input Data](input_data.md) functions.

| Key | Type | Description |
| -----------|-----| -------------- |
| `<var>` | `math_expr` | Math expression to initialize variable `<var>` |

##### `[model.<comp_k>.diffusion]`

The `diffusion` subsection defines the self-diffusion coefficient associated
with each *variable*.

| Key | Type | Description |
| -----------|-----| -------------- |
| `<var>` | `math_expr` | Diffusion math expression assigned to `<var>` |

##### `[model.<comp_k>.reaction]`

This subsection defines the reaction network associated with each *variable*.

| Key | Type | Description |
| -----------|-----| -------------- |
| `<var>` | `math_expr` | Reaction network expression assigned to `<var>` |
| `[jacobian]` | subsection | Jacobian expressions for compartment reaction networks |

*Variables* (`<var>`s) in this section differ from other math expressions in that
all *variables* names on this `<comp_k>` are available to form the
expression.

:::tip Example
The variables `u` and `v` are available to each other for a subsection
`[model.nucleus.reaction]`
```ini
[model.nucleus.reaction]
u = u*v
v = u*v
```
:::

##### `[model.<comp_k>.reaction.jacobian]`
This subsection lists the
[jacobian matrix](https://en.wikipedia.org/wiki/Jacobian_matrix_and_determinant)
entries for the `reaction` section. It must follow the syntax of
`d<var_i>__d<var_j>`, which reads as the *partial derivative of `<var_i>` with
respect to `<var_j>`*.

| Key | Type | Description |
| -----------|-----| -------------- |
| `d<var_i>__d<var_j>` | `math_expr` | Reaction network jacobian entries |

Similar to `[model.<comp_k>.reaction]`, math expressions here allow all
compartmenet *variables* to be used.

:::tip Example
The section `[model.nucleus]` for a [Gray-Scott model with F=0.0420 and
k=0.0610](http://mrob.com/pub/comp/xmorphia/F420/F420-k610.html) may look like
this:

```ini
[model.nucleus.initial]
u = 0.5
v = (x>0) && (x<0.5) && (y>0.) && (y<0.5) ? 1 : 0

[model.nucleus.diffusion]
u = 2e-5
v = 2e-5/2

[model.nucleus.reaction]
u = -u*v^2+0.0420*(1-u)
v = u*v^2-(0.0420+0.0610)*v

[model.nucleus.reaction.jacobian]
du__du = -v^2-0.0420
du__dv = -2*u*v
dv__du = v^2
dv__dv = 2*u*v-(0.0420+0.0610)
```
:::

##### `[model.<comp_k>.boundary]`

The boundary section encapsulates the information for the transmmission conditios
$\mathcal{T}^{kl}_\star$. It contains subsections for every other compartment
`<comp_l>` making reference to the membrane defined by $\Gamma^{kl}$, where
`<comp_k>` corresponds to the interior $k$-th domain $\Omega^k$ and `<comp_l>`
to the exterior $l$-th domain $\Omega^l$.

| Key | Type | Description |
| -----------|-----| -------------- |
| `[<comp_l>.outflow]` | subsection | Outflow for membranes between `<comp_k>` and `<comp_l>` |

##### `[model.<comp_k>.boundary.<comp_l>.outflow]`

This section contains the transmmission condition $\mathcal{T}^{kl}_i$ for each
*variable*, `<var_i>` in `<comp_k>`.

| Key | Type | Description |
| -----------|-----| -------------- |
| `<var_i>` | `math_expr` | Transmission condition for `<var_i>` |
| `[jacobian]` | subsection | Outflow jacobian entries |

The math expressions in this subsection are allowed to use all *variables* from
both `<comp_k>` and `<comp_l>` compartments. Because they may have the same
name, the parser distinguish them by assigin a `_i` (inside) suffix for variables
in `<comp_k>` and `_o` (outside) for variables in `<comp_l>`. Additionally,
their normal gradients, $\nabla u_i^k\cdot \mathbf{n}^k$ and
$\nabla u_j^l\cdot \mathbf{n}^l$, are also available as `d<var_i>__dn`. Here,
$\mathbf{n}^k$ is the unit outer normal vector of $\Omega^k$ on the
intersection.

:::tip Example
Let's say that our configuration file defines two compartments, `nucleus` ($n$)
and `cytoplasm` ($c$), and both of them contain a variable `u` ($u^n$ and
$u^c$). Then, a boundary condition $\mathcal{T}^{nc}(u^n,u^c) =
u^n-u^c$ for the nucleus and $\mathcal{T}^{cn}(u^c,u^n) =
-\mathsf{D}^n\nabla u^n\cdot \mathbf{n}_T$ for the
cytoplasm is
```ini
[model.nucleus.boundary.cytoplasm.outflow]
# Flux proportional to the difference of species
u = u_i-u_o

[model.cytoplasm.boundary.nucleus.outflow]
# Flux in cytoplasm equal to the outflow of u in nucleus (assume D = 0.1)
u = 0.1 * du_o__dn
```
:::

:::caution
For compartments solved with Finite Volume method, the normal gradient
`d<var_i>__dn` is only available for variables that exist on both compartments.
Otherwise, it will be evaluated to
[NaN](https://en.wikipedia.org/wiki/NaN#Quiet_NaN).
:::

##### `[model.<comp_k>.boundary.<comp_l>.outflow.jacobian]`

The jacobian for the outflow works similarly as the one for the reaction. That
is, all *variables* from inside and outside as well as their normal gradients
are available. Additionally, you have to distinguish between jacobians with
respect to inside and outside *variables*.

| Key | Type | Description |
| -----------|-----| -------------- |
| `d<var_i>__d<var_j>_i` | `math_expr` | Self-domain jacobian entries |
| `d<var_i>__d<var_j>_o` | `math_expr` | Cross-domain jacobian entries |

Now, since in most cases these jacobians are zero, there is only need to specify
the non-zero entries. The rest will be assumed to be zero.

:::tip Example
Following the example for the outflow, its jacobian would be

```ini
[model.nucleus.boundary.cytoplasm.outflow.jacobian]
du__du_i = 1
du__du_o = -1

# [model.cytoplasm.boundary.nucleus.outflow.jacobian]
# du__du_i = 0
# du__du_o = 0
```
:::

### `[logging]`

:::caution Work in progress
:::

The logging settings are directly forwarded to the `dune-logging` module. Please
check its doxygen documentation for detailed information.

:::tip Example
A simple configuration is the following:
```ini
[logging]
# possible levels: off, critical, error, warn, notice, info, debug, trace, all
default.level = info

[logging.sinks.stdout]
pattern = [{reldays:0>2}-{reltime:8%T}] [{backend}] {msg}
```
:::

## Runge-Kutta methods (`rk_id`)

The following are the accepted Runge-Kutta methods in our time stepping schemes

|  Method |
| ----------------------- |
| `explicit_euler`        |
| `implicit_euler`        |
| `heun`                  |
| `shu_3`                 |
| `runge_kutta_4`         |
| `alexander_2`           |
| `fractional_step_theta` |
| `alexander_3`           |



## Example

:::tip Example
An example of an ini file with all the parameters required by
`dune-copasi` for the [Gray-Scott model](http://mrob.com/pub/comp/xmorphia/F420/F420-k610.html)
is

```ini title="config.ini"
[grid]
dimension = 2
file = path/to/grid.msh
initial_level = 1

[model]
order = 1

[model.data]
input_data = path/to/image.tiff

[model.time_stepping]
rk_method = alexander_2
begin = 0
end = 1
initial_step = 0.3
min_step = 1e-3
max_step = 0.35
decrease_factor = 0.5
increase_factor = 1.5

[model.time_stepping.newton]
reduction = 1e-8
min_linear_reduction = 1e-3
fixed_linear_reduction = false
max_iterations = 40
absolute_limit = 1e-12
reassemble_threshold = 0.0
keep_matrix = true
force_iteration = false

[model.time_stepping.newton.linear_search]
strategy = hackbuschReusken
max_iterations = 10
damping_factor = 0.5

[model.compartments]
domain = 0

[model.domain.initial]
u = 0.5 + input_data(x,y)
v = (x>0) && (x<0.5) && (y>0.) && (y<0.5) ? 1 : 0

[model.domain.diffusion]
u = 2e-5
v = 2e-5/2

[model.domain.reaction]
u = -u*v^2+0.0420*(1-u)
v = u*v^2-(0.0420+0.0610)*v

[model.domain.reaction.jacobian]
du__du = -v^2-0.0420
du__dv = -2*u*v
dv__du = v^2
dv__dv = 2*u*v-(0.0420+0.0610)

[model.writer]
file_path = gray_scott
```
:::

[GMSH]: http://gmsh.info/
