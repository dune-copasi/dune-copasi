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

<!-- truncate -->

## Running Dune Copasi in Docker

The first step is to obtain the source of Dune Copasi. This can simply be done by pulling it from the repository.
```
git clone https://gitlab.dune-project.org/copasi/dune-copasi.git
```
Go into the Dune Copasi source code and build the docker image. We also create a working directory which we will mount to the docker image.
```
cd dune-copasi
docker build -t dune-copasi .
mkdir -m o+rw workdir && cd workdir
```
Building the docker container and compiling the solver can take some time. Once the docker container is build we are now set to test the Dune Copasi solver. To do this we will create a simple initialization file. We write the code of the minimal example to the file `test.ini`:  
```
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
gauss.expression = x, y, z: exp(-(x^2+y^2+z^2)/(4*t*diffusion)) / (4*3.14159265359*diffusion)

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
  &J_{in}(V_m,h) = \frac{h {V_m}^2 (1-V_m)}{\tau_{in}},\\
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

Let us get started with creating the ini file for our electrophysiology simulation. First we will define our grid:

 ```ini
 [grid]
 dimension = 2
 extensions = 3 3
 refinement_level = 8
 ```

 This will create a square grid existing of 2 triangles with 3 by 3 extension. Thereafter the grid will be refined 8 times. The solver can also load Gmsh meshes but we will come to that at a later point. For the moment this will be all that we need to define our grid. The next step is the  ```parser_context```. In the parser context we will define the constants and functions that we will use in our model.

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

 u_0.type = constant
 u_0.value = 0.13

 gauss.type = function
 gauss.expression = x, y, z: exp(-(x^2+y^2+z^2)/(4*0.0001)) / (4*3.14159265359*0.0001)/400

 periodic.type = function
 periodic.expression = t : cos(t) > 0.95 ? 1 : 0
 ```
