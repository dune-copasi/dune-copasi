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

In this blog I will show you can configure the Dune Copasi simulation environment to perform a very simple cardiac electrophysiology simulation. I will show you how you can simulate the Mitchell-Schaeffer cardiomyocyte electrophysiology model in a 2D environment. First I will show you how you can obtain a working version of Dune Copasi within a docker environment. Thereafter, I will introduce the electrophysiology model and how to implement it.

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
  \partial_t V_m = &J_{in}(V_m,h) + J_{out}(V_m),\\
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

WORK IN PROGRESS ... 
