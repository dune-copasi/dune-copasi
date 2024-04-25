---
title: "DuneCopasi: A multi-compartment reaction-diffusion simulator for systems biology"
tags:
  - Diffusion
  - Reaction Networks
  - Multiple Compartments
  - Continuous Galerkin Method
  - Partial Differential Equations
  - Finite Element Method
  - Systems Biology
  - C++
  - WebAssembly
  - Docker
  - DUNE
authors:
  - firstname: Santiago
    family: Ospina De Los Ríos
    orcid: 0000-0003-0814-9670
    affiliation: "1, 2"
  - firstname: Peter
    family: Bastian
    affiliation: "1"
  - firstname: Liam
    family: Keegan
    orcid: 0000-0002-0654-4979
    affiliation: "3"
  - firstname: Sven
    family: Sahle
    orcid: 0000-0002-5458-7404
    affiliation: "4"
  - firstname: Dylan
    family: Vermoortele
    orcid: 0000-0001-8769-783X
    affiliation: "6"
  - firstname: Lilija
    family: Wehling
    orcid: 0000-0002-8697-5348
    affiliation: "4, 5"

affiliations:
  - name: Interdisciplinary Center for Scientific Computing (IWR), Heidelberg University, Germany
    index: 1
  - name: Heidelberg Graduate School of Mathematical and Computational Methods for the Sciences (HGS MathComp), Heidelberg University, Germany
    index: 2
  - name: Scientific Software Center (SSC), Heidelberg University, Germany
    index: 3
  - name: BioQuant, Centre for Organismal Studies (COS), Heidelberg University, Germany
    index: 4
  - name: Institute of Pathology, University Hospital Heidelberg, Germany
    index: 5
  - name: Cardiovascular Imaging and Dynamics, KU Leuven, Belgium
    index: 6
# date: 25 January 2024
bibliography: paper.bib
---

# Summary
DuneCopasi is a C++ library providing a numerical simulator for systems of reaction-diffusion equations on single and multiple compartments using the **D**istributed and **U**nified **N**umerics **E**nvironment DUNE [@bastian_concepts_2021]. It supports multi-threaded computations with continuous finite elements on unstructured grids in 1, 2, and 3 dimensions. An application program is made available as a Docker Image, as a Web Assembly binary, or as a GUI through the Spatial Model Editor [@liam_keegan_2023] for the manipulation of spatial models of bio-chemical reactions.

# Model Problem

In short, our program solves the system of partial differential equations (PDE) formulated in the form of

$$
\begin{aligned}
  \partial_t (\phi_{ik} u_{ik}) &= \nabla \cdot \vec{\mathcal{j}}_{ik}+\mathcal{R}_{ik}(\mathbf{u})\qquad &\text{in }\Omega_k\subset \mathbb{R}^d,\\
  u_{ik} &= u^{(0)}_{ik} \qquad &\text{on }\Gamma_k^D\subseteq\partial\Omega_k,\\
  \vec{\mathcal{j}}_{ik}\cdot \mathbf{n}_k &= \sum_{l\in T_k}\mathcal{T}_{ikl}(\mathbf{u}) \qquad &\text{on }\partial\Omega_k \setminus \Gamma^D_{k},\\
  \vec{\mathcal{j}}_{ik} &:= \sum_{j\in N_k}\mathsf{D}_{ijk}\nabla u_{jk}\\
\end{aligned}
$$

where $i\in N_k$ is the subscript for each mass balance equation in the $k$-th compartment and with $k\in K$, $T_k:=\{l\in K : \Gamma_{kl}\neq\varnothing\}$, $\mathbf{u}_k:=(u_{1k},\ldots,u_{\text{dim}(N_k)k})$, $\mathbf{u} := (\mathbf{u}_1, \ldots, \mathbf{u}_{\text{dim}(K)})$, $\mathbf{n}_k$ the unit outer normal on $\Omega_k$, $\overline{\Omega}:=\bigcup_{k\in K}\overline{\Omega}_k$, and $\Gamma_{kl}$ is $\partial\Omega_k\cap\overline{\Omega}_l$ if $k \neq l$ or $\partial\Omega_k\cap\overline{\Omega}$ otherwise.


All of the parameterization for the equations are easily configurable at run-time using the command line or a configuration file. That is, the underlying grid, the partition for each sub-domain $\Omega_k$ on the grid, the reaction operator $\mathcal{R}_{ik}(\mathbf{u})$, the non-linear transmission conditions $\mathcal{T}_{ikl}(\mathbf{u})$ on $\Gamma_{kl}$, the dirichlet boundary conditions $u^{(0)}_{ik}$ on $\Gamma^D_{k}$, the storage terms $\phi_{ik}$, and the cross-diffusion terms $\mathsf{D}_{ijk}$.

# Research and Use Cases

The principal purpose of our project has been to bridge the gap between systems biology and scientific computing by providing researchers with an accessible and reproducible spatial simulator. An illustration of this is the study presented by @wehling_2022, depicted in \autoref{fig:wehling_2022}, where simulations with the Spatial Model Editor aided in assessing and better understanding the mechanisms of the YAP and TAZ proteins that regulate cell proliferation in liver cells. @elias_2014 present further examples that fit our model problem in this context. Additionally, DuneCopasi stands as a versatile tool capable of accommodating diverse computational needs beyond its initial focus on systems biology. For example, the general purpose nature allows for seamless integration into electrophysiology simulations, without requiring any modifications. The emerging concept of identifying patients based on personalized cardiac electrophysiology simulations [@arevalo_2016] underscores the demand for easy-to-use and efficient simulators. \autoref{fig:cardiac_ep} demonstrates how DuneCopasi emerges as a flexible solver poised to meet these evolving needs.

![Mathematical modelling predicts that nuclear phosphorylation controls spatial localization of Hippo signalling pathway components YAP and TAZ. Here, we compare two model topologies - "Model 1" (canonical) and "Model 2" (alternative model) - concerning the intracellular distribution of YAP and TAZ proteins (Pr) and their phosphorylated counterparts (pPr). If the phosphorylation of YAP/TAZ takes place exclusively outside the nucleus, as shown in "Model 1", PDE simulation indicates low spatial accordance with the experimentally measured subcellular localization of YAP/TAZ.  Whereas, "Model 2" describes YAP/TAZ protein phosphorylation and dephosphorylation in the nucleus. The simulation of "Model 2" agrees with experimentally measured subcellular distribution of YAP/TAZ proteins. \label{fig:wehling_2022}](wehling_2022.pdf)

![A DuneCopasi based pipeline for cardiac electrophysiology simulation of activation and repolarization of the ventricles. A personalized biventricular geometry is constructed based on cardiac magnetic resonance imaging of long and short-axis images. The bi-domain equation of cardiac electrophysiology using the Mitchell-Schaefer cell model is solved on the biventricular geometry allowing to simulate the cardiac activation and repolarization. \label{fig:cardiac_ep}](cardiac_ep.png)

# Capabilities

Many features of DuneCopasi have been designed and developed with the specific case of systems biology in mind. Thus, a substantial effort has been made to have a library that is interoperable with systems biology data assets and requirements needed by their practitioners. Among others, include:

* a single executable configurable at run-time,
* a run-time mathematical expression parser that understands the Systems Biology Markup Language SBML [@sbml_2003] specification,
* the input of custom grid data and image files in the TIFF format for initial spatial concentrations and other parameters,
* a powerful, yet simple, boundary/transmission condition specification for each compartment to account for generic trans-membrane fluxes,
* a (non-linear) cross-diffusion specification that allows any species to cross-diffuse into other mass balance equations,
* a specification to compare results with user-defined objective functions,
* an in-place function interpolator that reduces the computational cost of evaluating common expensive reaction operations, and
* an embedded random field generator by @kempf_2023 to represent statistical spatial variations on the domain.

Furthermore, DuneCopasi is a stand-alone multi-compartment reaction-diffusion solver that may easily be used for many other fields of research and engineering fitting our model problem.

# Statement of Need

<!--- taken from proposal --->
In the context of cell biology, computational modelling has become an established and important technique to understand and discover biological processes. Comparing model simulations and experimental data is used to set up models, validate models, and test hypotheses about biological processes. Up to now, the majority of biological systems studied in this way have been spatially homogeneous, i.e., experimental data are typically time-resolved concentrations or amounts of biochemical species, and the mathematical description is based on ordinary differential equations. Technological development in live-cell imaging during the last years led to the availability of an increasing amount of spatio-temporal data which calls for a corresponding shift towards spatially resolved models.

One of the natural extensions of spatially homogeneous models to account for spatio-temporal dynamics is its formulation in terms multi-compartment reaction-diffusion system of equations. Finite element frameworks like DUNE [@bastian_grid_2008], Deal.II [@dealII_1995], Necktar++ [@nektar_2015], or FreeFem++ [@freefem_2012] are too generic by design and do not address the specific requirements of the computational modelling practices in systems biology out of the box. In particular, most finite element frameworks do not consider the multi-compartment systems in their design resulting, if at all possible, in inefficient simulations or obscure tricks to force this feature. Here, we extended DUNE-PDELab [@bastian_pdelab_2010] and the work by @muthing_2015 with efficient data structures especially tailored for this task.

In the space of systems biology, two well-known software packages also provide similar features to DuneCopasi, namely Morpheus [@morpheus_2014] and VCell [@vcell_1997]. Morpheus is a computational tool for multi-cellular systems which follows a cellular automata design to set rules of interactions between individual cells. It focuses on modelling approaches like cellular Potts models combined with gradients modelled by PDEs, rather than multi-compartment PDE models. Furthermore, due to its design, it falls short for moderate and big PDE computations since explicit solvers are severely limited by the required time steps. On the other hand, VCell provides a fully implicit finite volume solver on structured grids that can very well manage the multi-compartment case in addition to membrane unknowns. In comparison, our solution aims to resolve the geometry and the transmission conditions directly in the weak formulation of the problem while also being designed to provide implicit and monolithic solvers as well as tailored preconditioners for the underlying linear solvers.

# Acknowledgements

We want to thank all the contributions that aided in the development and deployment of the package. This work has been funded and supported by the German Federal Ministry of Education and Research (BMBF) FKZ 031L0158.

# References
