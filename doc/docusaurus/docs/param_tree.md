---
id: param_tree
title: Parameter Tree
sidebar_label: Parameter Tree
---

## Reference Keys

A parameter tree is a collection of *key–value pairs* that may be grouped by
sections. We call it a *tree* because sections may contain also other sections.
Here, we follow the [`DUNE` ini file convention](ini_file.md) to refer to
sections and keyords. In addition, keywords between angle brackets (i.e.
`<key>`) refer to **multiple** user defined keywords.

### `[grid]`

The grid is entred by a [GMSH file][GMSH] and should be formed by triangles and
recangles in 2D, tetraedra and hexahedra in 3D. It must define
*physical group*s each refering to a different compartment.
Each *physical group* can only be formed the same collection of geometric types.

:::caution Bug in GMSH reader
The GMSH file should not contain the `$PhysicalGroup` section at the
begining of the file (this is a knwon bug in the `dune-grid` GMSH reader).
:::

| Key | Type |Description |
| ----------------| -------------- | ----- |
| `dimension`     | `integer` | The dimension of the grid |
| `file`          | `path` | [GMSH v2 file][GMSH] directory absolute or relative to the executable |
| `initial_level` | `integer` | The refinement level to start the simulation with |


### `[model]`

This section is the one that contains most of the aspects for simulation.

| Key | Type | Description |
| -----------|-----| -------------- |
| `order`           | `integer` | Finite Element polynomial order. if `order==0` FV else CG |
| `[data]` | subsection | List of input files |
| `[writer]` | subsection | Options for vtk writer |
| `[time_stepping]` | subsection | Options for time-stepping |
| `[compartment]`   | subsection | List of compartmentes to simulate |
| `[<compartment>]` | subsection | Options for equations per each compartment |

Whenever a compartment is coposed by cubes (2D) or hexahedra (3D), then such
compartment is modeled using Finite Volumes. These compartments are
suitable to model membranes processes.

The subsection `[<compartment>]` is a placeholder for the name of all
compartment given in `[compartment]`. See `[model.compartment]` subsection for
more information.

#### `[model.data]`

This section will contain all the input files for the initial conditions of the
model.

| Key | Type | Description |
| -----------|-----| -------------- |
| `[<data_name>]` | `path` | 16-bit grayscale TIFF file |

See [Input Data](input_data.md) section for more information about usage.

#### `[model.output]`

The model ouput is controlled by the following parameters:

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
| `min_step` | `float` | Maxumum time step allowed |
| `max_step` | `float` | Minimum time step allowed |
| `decrease_factor` | `float` | Time step decrease factor when time step fails |
| `increase_factor` | `float` | Time step increase factor when time step succeds |

The possible `rk_id`s can be found at the end of this document.

#### `[model.compartment]`

The compartment section is filled with `key=value` pairs that assign a phisical
group id (`value`) to an compartment name (`key`).

| Key | Type | Description |
| -----------|-----| -------------- |
| `<compartment>` | `integer` | Phisical group id assigned to `<compartment>` subdomain |

Each compartment id maps to a *physical group* in the gmsh identifiers.
Although the gmsh format allows you to name such physical groups.
Notice that `dune-copasi` uses 0-based
 indices while `gmsh` uses 1-based indices. In other words,
`<gmsh_physical_group> = <dune_copasi_compartment> + 1`.

For example, let's say that there is two *physical groups* in our gmsh file
and we are going to name them as `nucleous` and `cytoplasm` compartments:


```ini
[model.compartments]
 # nucleous corresponds to the physical group 1 in the gmsh file
nucleous  = 0
 # cytoplasm corresponds to the physical group 2 in the gmsh file
cytoplasm = 1

# They then become available as new compartment subsections
[model.nucleous]
# Parameters for the nucleous compartment

[model.cytoplasm]
# Parameters for the nucleous compartment
```

#### `[model.<compartment>]`

Each compartment will define its own initial conditions,
its diffusion-reaction system, etc. A set of *variables* is be assigned to each
compartment. The amount *variables* may be different on each compartment, but
the same namings must be used on all the `[model.<compartment>]` subsections.

| Key | Type | Description |
| -----------|-----| -------------- |
| `[initial]` | subsection | List of *variables* and their initial conditions |
| `[diffusion]` | subsection | List of *variables* and their diffusion coefficient |
| `[reaction]` | subsection | List of *variables* and their reaction networks |
| `[outflux]` | subsection | Definition of output fluxes on this `<compartment>` |

##### `[model.<compartment>.initial]`

The `initial` subsection allows the initialization for each of the *variables*.
Expressions in this section may contain [Input Data](input_data.md) functions.

| Key | Type | Description |
| -----------|-----| -------------- |
| `<var>` | `math_expr` | Math expression to initializate variable `<var>` |

##### `[model.<compartment>.diffusion]`

The `diffusion` subsection defines the diffusion associated with each
*variable*.

| Key | Type | Description |
| -----------|-----| -------------- |
| `<var>` | `math_expr` | Diffusion math expression assigned to `<var>` |

##### `[model.<compartment>.reaction]`

This subsection defines the reaction network associated with each *variable*.

| Key | Type | Description |
| -----------|-----| -------------- |
| `<var>` | `math_expr` | Reaction network expression assigned to `<var>` |
| `[jacobian]` | subsection | Jacobian expressions for compartment reaction networks |

*Variables* (`<var>`s) in this section differ from other math expressions in that
all *variables* names on this `<compartment>` are available to form the
expression. E.g.

```ini
[model.nucleous.reaction]
u = u*v
v = u*v
```

##### `[model.<compartment>.reaction.jacobian]`
This subsection lists the
[jacobian matrix](https://en.wikipedia.org/wiki/Jacobian_matrix_and_determinant)
entries for the `reaction` section. It must follow the syntax of
`d<var_i>__d<var_j>`, which reads as the *partial derivative of `<var_i>` with
respect to `<var_j>`*.

| Key | Type | Description |
| -----------|-----| -------------- |
| `d<var_i>__d<var_j>` | `math_expr` | Reaction network jacobian entry |

Similar to `[model.<compartment>.reaction]`, math expressions here allow all
compartmenet *variables* to be used. For example, the section `[model.nucleous]`
for a [Gray-Scott model with `F=0.0420` and `k=0.0610`](http://mrob.com/pub/comp/xmorphia/F420/F420-k610.html)
may look like this:

```ini
[model.nucleous.initial]
u = 0.5
v = (x>0) && (x<0.5) && (y>0.) && (y<0.5) ? 1 : 0

[model.nucleous.diffusion]
u = 2e-5
v = 2e-5/2

[model.nucleous.reaction]
u = -u*v^2+0.0420*(1-u)
v = u*v^2-(0.0420+0.0610)*v

[model.nucleous.reaction.jacobian]
du__du = -v^2-0.0420
du__dv = -2*u*v
dv__du = v^2
dv__dv = 2*u*v-(0.0420+0.0610)
```

##### `[model.<compartment>.outflux]`

The fluxes are set automatically to [dirichlet-dirichlet](https://en.wikipedia.org/wiki/Dirichlet_boundary_condition)
boundary conditions if the *variable* is shared between the two
intersecting domains. No key is required for this.

:::tip Comming feature
Custom fluxes will be available in upcoming releases.
:::

### `[logging]`

:::caution Work in progress
:::

The logging settings are directly forwarded to the `dune-logging` module. Please
check its doxygen documentation for detailed information. A simple configuration
is the following:

```ini
[logging]
# possible levels: off, critical, error, waring, notice, info, debug, trace, all
default.level = trace

[logging.sinks.stdout]
pattern = [{reldays:0>2}-{reltime:8%T}] [{backend}] {msg}
```

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

A example of an ini file with all the parameters required by
`DuneCopasi` for the [Gray-Scott model](http://mrob.com/pub/comp/xmorphia/F420/F420-k610.html)
is

```ini
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

[GMSH]: http://gmsh.info/
