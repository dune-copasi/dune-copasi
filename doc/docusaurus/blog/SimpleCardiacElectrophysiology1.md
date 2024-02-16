---
title: Cardiac Electrophysiology simulations with Dune Copasi - I
description: A tutorial for cardiac electrophysiology simulations using Dune Copasi.
slug: cardiacEP-simulation-v1
authors:
  - name: Dylan Vermoortele
    title: Postdoctoral Researcher, Cardiovascular Imaging and Dynamics, KU Leuven, Belgium
tags: [Electrophysiology, Dune Copasi]
hide_table_of_contents: false
---

import Tabs from '@theme/Tabs';
import TabItem from '@theme/TabItem';

In this blog I will show you can configure the Dune Copasi simulation environment to perform a very simple cardiac electrophysiology simulation. I will show you how you can simulate the Mitchell-Schaeffer cardiomyocyte electrophysiology model in a 2D environment. First, I will show you how you can obtain a working version of Dune Copasi within a docker environment. Thereafter, I will introduce the electrophysiology model and how to implement it.

![alt text](images/EP_simulation1_activation.png)

<!-- truncate -->

We will simulate the Mitchell-Schaefer model where we stimulate the bottom left corner using a Gaussian point source. This will result in the activation and repolarization of our 2D square grid. Lets get started.  

## Running Dune Copasi in Docker

The first step is to obtain the docker image of Dune Copasi. This can simply be done by pulling it from the repository.
``` sh
docker pull registry.dune-project.org/copasi/dune-copasi/dune-copasi:debian--master
```
We also create a working directory where we will save the simulation files.
``` sh
mkdir -m o+rw workdir && cd workdir
```
Having fetched the docker container we are now set to test the Dune Copasi solver. To do this we will create a simple initialization file. We write the code of the minimal example to the file `test.ini`:  
``` sh
echo "[grid]
dimension = 2
extensions = 2 2
origin = -1 -1
refinement_level = 5

[compartments.domain]
type = expression
expression = 1

[parser_context]
diffusion.type = constant
diffusion.value = 0.005

gauss.type = function
gauss.expression = x, y, z: exp(-(x^2+y^2+z^2)/(4*diffusion)) / (4*3.14159265359*diffusion)

[model.scalar_field.u]
compartment = domain
cross_diffusion.u.expression = diffusion
initial.expression = gauss(position_x, position_y, position_z)
storage.expression = 1

[model.time_step_operator]
time_begin = 1
time_end = 4
" >> test.ini
```
By parsing this initialization file to the Dune Copasi solver it will exactly know how to construct the problem and the simulation will be started. This can be done by:
```
docker run -v $PWD:/dunecopasi dune-copasi --config=test.ini
```
The Dune Copasi solver will now construct the problem and start running. Once it finished you should get the message `dune-copasi successfully finished :)`. You now just solved the heat equation with an initial Gaussian distribution (we did not save it though ;)). Although this example was rather simple it does nicely illustrate the use of the Dune Copasi solver. You define the problem that you want to solve in your initialization file and pass this to Dune Copasi which does all the rest! For now do not bother too much about the syntax of the initialization file, we will explain everything once we construct our electrophysiology model. However, before implementing the initialization file of our cardiac electrophysiology simulation we will describe, in the next section, the equations that we want to implement.

## The mathematical model of cardiac electrophysiology

### The Mitchell-Schaeffer model of cell electrophysiology

Cardiomyocytes maintain an ion concentration gradient between the intracellular and extracellular space. This ion concentration gradient results in a potential across the cell membrane (i.e. the transmembrane potential). A typical resting membrane potential is around $-85 mV$. A disturbance of this resting state can result in a cardiac action potential. This is a brief change in membrane potential of heart cells resulting in a typical transmembrane potential waveform. This typical waveform is determined by the dynamics of the ion channels in the cell membrane. In the very simple Mitchell-Schaeffer model the action potential is described by 2 ion channel currents. The governing equations are given by

$$
\begin{aligned}
  C_m \partial_t V_m = &J_{in}(V_m,h) + J_{out}(V_m),\\
  &J_{in}(V_m,h) = \frac{z {V_m}^2 (1-V_m)}{\tau_{in}},\\
  &J_{out}(V_m)  = \frac{V_m}{\tau_{out}},
\end{aligned}
$$

where $\tau_{in}$ and $\tau_{out}$ are time-scale variables, $J_{in}$ and $J_{out}$ describe the respective inward and outward current, $V_m$ is the normalized transmembrane potential and $h$ is a gating variable. The dynamics of the gating variable are given by

$$
\frac{\partial z}{\partial t} =  
\begin{cases}
  \frac{z-1}{\tau_{open}},& \text{if } u\leq u_0,\\
  \frac{z}{\tau_{close}},& \text{if }  u > u_0.
\end{cases}
$$
where $\tau_{in}$ and $\tau_{out}$ are time-scale variables for the gating variable and $u_0$ defines the critical threshold for depolarization.

### The bi-domain equation of cardiac electrophysiology

In the section above we described a model capturing the basic electrophysiology of cardiomyocytes. Now we want to extend this by coupling multiple cardiomyocytes yielding a spatial description of the cardiac electrophysiology. To obtain this we will depart from Ohm's law stating

$$
  J = \sigma E.
$$  

In our specific case we have two compartments. We will use an intracellular and extracellular current density and potential. We obtain

$$
\begin{aligned}
  J_{i} &= \sigma_i E_i = \sigma_i \nabla \varphi_i,\\
  J_{e} &= \sigma_e E_e = \sigma_e \nabla \varphi_e. \\
\end{aligned}
$$

Furthermore, we assume that any current leaving the intracellular space flows into the extracellular (and vice-versa) and that there is no charge accumulation. Therefore, the change in current density in the extracellular space has to be equal to the current flowing into the intracellular space. Using Gauss' law (i.e. the divergence theorem) this can be written as

$$  
  \nabla \cdot J_{i} + \nabla \cdot J_{e} = 0 = \nabla \cdot \sigma_i \nabla \varphi_i + \nabla \cdot \sigma_e \nabla \varphi_e.
$$

The current flowing from the intracellular into the extracellular space is nothing else than the transmembrane current $I_m$ (per surface area). The conservation of current can be rewritten as

$$
\begin{aligned}
  I_m &= \nabla \cdot \sigma_i \nabla \varphi_i = -\nabla \cdot \sigma_e \nabla \varphi_e.\\
\end{aligned}
$$

The transmembrane current $I_m$ can be modeled using the cable equation where we combine the ion channel current $I_{ion}$ with a capacitative component

$$
\begin{aligned}
  I_m &= A_m \Big( C_m \frac{\partial V_m }{\partial t} + I_{ion} \Big)\\
\end{aligned}
$$

where $V_m = \varphi_i - \varphi_e$ is the transmembrane potential and $A_m$ is the surface density. Rewriting these equations (by eliminating $\varphi_i$) yields the standard form of the bidomain equation of cardiac electrophysiology:

$$
\begin{aligned}
  A_m \Big( C_m \frac{\partial V_m }{\partial t} + I_{ion} \Big) &= \nabla \cdot \sigma_i \nabla V_m + \nabla \cdot \sigma_i \nabla \varphi_e ,\\
  0 &= \nabla \cdot \sigma_i \nabla V_m + \nabla \cdot (\sigma_i + \sigma_e) \nabla \varphi_e . \\
\end{aligned}
$$

## EP simulations using Dune Copasi

### The Dune Copasi Ini file

Now that we have discussed the mathematical problem we can tend our attention to implementing this so that we can simulate it using Dune Copasi. To do this we will write a configuration file (typically we use the .ini extension) that entails all the information that is needed for the Dune Copasi solver to construct the problem at hand. The configuration file consists of a parameter tree containing the parameters. The nodes at the top level can be divided in 4 categories: grid, parser_context, compartments, model.

| Category | Description  |
|     ---- | -----------  |
| `grid`                  | The nodes under the grid header describe the basic properties of the grid. For example the dimension, extension, origin, etc.
| `parser_context`        | Under the parser_context node we are able to define constant and functions that can be used to define the mathematical model.
| `compartments`          | In this category we define the properties related to division of the problem to multiple compartments.
| `model`                 | This is the heart of the mathematical problem where we define the reaction terms, diffusion, etc.

The ini files we use follow the DUNE convention. In short, the data is composed of a key–value pairs on the form ```key = value```. For example, to define a property of the category grid we write

```ini
grid.property = value
```

In essence the ini file is nothing else than a parameter tree of key-value pairs. For convenience one can use sections to avoid rewriting the leading part and provide more structure. This is done by using ``` [section.subsection] ``` syntax. The keys below will now be within the designated subsection. As an example, the three following ini files, are equivalent

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


### Setting up the Mitchell-Schaefer model

Let us get started with creating the ini file for our electrophysiology simulation. Open your favorite text editor and you can start building the configuration file. First we will define our grid.

#### The grid properties

The grid is defined as follows:

 ```ini
 [grid]
 dimension = 2
 extensions = 3 3
 refinement_level = 8
 ```

 This will create a square grid existing of 2 triangles with 3 by 3 extension. Thereafter the grid will be refined 8 times. The solver can also load Gmsh meshes but we will come to that at a later point. For the moment this will be all that we need to define our grid. The next step is the  ```parser_context```.

#### The parser context

 In the parser context we will define the constants and functions that we will use in our model. This is easiest seen by looking at an example:

 ```ini
 [parser_context]
 tau_in.type = constant
 tau_in.value = 0.1 # ms

 tau_out.type = constant
 tau_out.value = 1 # ms

 tau_open.type = constant
 tau_open.value = 150 # ms

 tau_close.type = constant
 tau_close.value = 120 # ms

 gamma_i.type = constant
 gamma_i.value = 3e-3 # 1/(\Omega*cm)

 A_m.type = constant
 A_m.value = 20e-1

 C_m.type = constant
 C_m.value = 1.0

 u_0.type = constant
 u_0.value = 0.13

 gauss.type = function
 gauss.expression = x, y, z: exp(-(x^2+y^2+z^2)/(4*0.0001)) / (4*3.14159265359*0.0001)/400

 pulse.type = function
 pulse.expression = t, t0, dt: sqrt((t-t0)^2)< dt ? 1 : 0
 ```

 To create a constant variable lets say ```tau_in``` we simply write ```tau_in.type = constant```)). The value of the constant then assigned by writing ```tau_in.value = 0.1```. In a similar fashion a function is defined. For our Gaussian function which we name ```gauss``` we simply write ```gauss.type = function```. The function itself is then defined by typing ```gauss.expression = x,y,z : exp(...``` where the ```x,y,z``` before the colon denote the function arguments and the function itself is defined on the right of the colon. Note that standard mathematical functions such as ```exp```, ```sin```, etc. are available from the underlying parser allowing great freedom to define mathematical expressions. At this point we already have a grid and defined some constants and functions, now we will start defining the problem.

#### Defining the compartments

Although this example is fairly simple and we only need 1 compartment (which is the complete square), the concept to define multiple compartments is similar and therefore we will already discuss this briefly. Defining a compartment is similar to defining a defining a function. For our example we have one compartment which we will call the domain. To define this compartment type under the compartments section ```domain.type = expression```. This indicates that we will provide an expression which will be evaluated for each entity and if this expression evaluates to be non-zero then it will be added to the compartment. Here we only have one compartment which is all of the space so we set this expression equal to the constant 1.

   ```ini
   [compartments]
   domain.type = expression
   domain.expression = 1
    ```

#### The model properties

The first step in defining the model is to chose the parser that will be used to define the problem. Here we will go for [ExprTk](https://github.com/ArashPartow/exprtk).

  ```ini
  [model]
  parser_type = ExprTk
  ```
Once that is selected we will start defining the first variable of our problem. In this case we will start with the transmembrane potential which we refer to as $$V_m$$. To do this we will create a subsection ```model.scalar_field.V``` with the properties of our problem.

```ini
[model.scalar_field.V]
compartment = domain
storage.expression = 1
cross_diffusion.V.expression = gamma_i / (A_m * C_m)
cross_diffusion.phi.expression = gamma_i / (A_m * C_m)
reaction.expression = (z*V^2*(1-V)/tau_in - V/tau_out)/C_m +  0.85*pulse(time, 10, 20)*gauss(position_x - 0.15, position_y - 0.15, position_z)
reaction.jacobian.V.expression = (z*V*(2-3*V)/tau_in-1/tau_out)/C_m
reaction.jacobian.z.expression = ((1-V)*V^2/tau_in)/C_m
initial.expression = 0.01
```

The ```compartment``` property will define the compartment on which the scalar_field $V$ is defined. Here we only had one compartment ```domain``` that corresponds to the complete grid. Since the equation determining $V$ is time-dependent we have to store the previous time-step for use in the time-stepping scheme. We do this by setting ```storage.expression``` to 1. The following 2 lines define the cross diffusion of respectively $V_m$ and $\varphi_e$. Then we define the reaction term of the problem and the associated Jacobian entries of the reaction term. Note that we defined a function in the reaction term explicitly using the position $$(x,y,z)$$ and time $$t$$ to assemble the problem. At last we provide a function defining the initial state of the problem. In a similar fashion we define the extracellular potential field $\varphi$.

```ini
[model.scalar_field.phi]
compartment = domain
storage.expression = 0
cross_diffusion.V.expression = gamma_i / (A_m * C_m)             # (gamma_i)/(A_m * C_m)
cross_diffusion.phi.expression = 1.25 * gamma_i / (A_m * C_m)    # (gamma_i + gamma_e)/(A_m * C_m)
reaction.expression = 0
constrain.boundary.expression = ((position_x < 0.025) & (position_y < 0.025)) ? 0.0 : no_value
initial.expression = 0.00
```
The $\varphi$ equation is not time-dependent so we set ```storage.expression```to zero.  Further, we add the cross diffusion components and set the reaction term to zero. Note that explicitly setting the reaction term to zero is not necessary as the solver sets all non-defined properties to zero. Thereafter, we add a constraint to a small boundary region in the bottom left corner. Note that adding this constraint is technically not necessary, but since the $\varphi_e$ variable only appears as gradients the solution is only determined up to a constant. So by setting the bottom left corner to zero we consider the potentials to be with respect to the potential in the bottom left corner. At last we provide the initial value. Similarly we define the gating variables $z$.

```ini
[model.scalar_field.z]
compartment = domain
storage.expression = 1
reaction.expression = (V <= u_0) ? ( (1 - z)/(tau_open) ): (-z/tau_close)
reaction.jacobian.z.expression = (V <= u_0) ? ( -1/(tau_open))  : (-1/tau_close)
initial.expression = 1
```

Since this equations has no diffusion components, we simply don't define them rendering the coefficients of the diffusion and cross-diffusion zero. By defining the problem as time-dependent with a reaction term and no diffusion this is in fact just a ordinary differential equation. So at this point we have defined the mathematical properties of the problem. The next step is to define the underlying time-step approach to solve these equations. This is done by setting the ```model.time_step_operator``` properties.

```ini
[model.time_step_operator]
type = ImplicitEuler
time_step_initial = 0.500
time_step_max = 0.500
time_end = 800
linear_solver.type = SuperLU
```
For this problem we will use an implicit-Euler time-stepping scheme and then we define the initial and maximal time-step. The simulation will step from time-point ```time_start``` to ```time_end```. In our case this will be from $0$ (default value) to $800$. To solve the linear problem we will use the SuperLU direct solver. The last step is to write the solution to a file. This is done by initializing the writer.

```ini
[model.writer]
time_step = 1.000
vtk.path = 2D_Mitchell_Schaefer_model
```

At this point we have defined all the properties of the problem and we are ready to simulate the model! Let us save the configuration file to ``EP_simulation1.ini``.

### Running the Mitchell-Schaefer model

Running the code is the same as before, simply call:

```sh
docker run -v $PWD:/dunecopasi dune-copasi --config=EP_simulation1.ini
```

The simulation will start and write the ```*.vtk```files to the specified directory. You can easily inspect these files by opening them in [Paraview](https://www.paraview.org/).

![alt text](images/EP_simulation1_activation.png)

At this point we have successfully simulated the activation of a 2D square grid using the Mitchell-Schaefer cell model. Although this simulation model is very simple it does already illustrate the capabilities of the Dune Copasi solver. In the next post we will dig deeper into the available options and explore true multi-compartmental problems.
