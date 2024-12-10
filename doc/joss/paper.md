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

DuneCopasi is a C++ library designed to simulate how chemical reactions and diffusion processes occur in space, which is crucial for understanding many biological systems. It allows users to model these processes in both simple and complex environments, using grids that represent one, two, or three dimensions. Optimized for modern computers, DuneCopasi can be used as a standalone command-line application, accessed through tools like Docker or a web-based terminal, as a C++ library, or integrated into a graphical interface. One such interface is the **S**patial **M**odel **E**ditor (SME) [@keegan_2023], which helps users create and manipulate biological reaction models. The project was initiated to provide the numerical foundation for a spatial modelling tool for biochemical reaction networks, complementing the **CO**mplex **PA**thway **SI**mulator (COPASI) software [@copasi_2006]. This development resulted from a collaboration between the COPASI team and the **D**istributed and **U**nified **N**umerics **E**nvironment (DUNE) team [@bastian_concepts_2021].

# Background

In the context of cell biology, computational modelling has become an essential technique for understanding and discovering biological processes. By comparing model simulations with experimental data, researchers can set up models, validate them, and test hypotheses about these processes. Traditionally, most biological systems studied this way have been spatially homogeneous. That is, the data are typically time-resolved concentrations or quantities of biochemical species, modelled using **O**rdinary **D**ifferential **E**quations (ODEs).

However, recent advances in live-cell imaging technology have made more detailed spatio-temporal data available, leading to a growing need for models that capture not only the time dynamics but also the spatial distribution of these biochemical species. This has prompted a shift toward spatially resolved models.

# Model Problem

A natural extension of spatially homogeneous ODE models is to incorporate spatio-temporal dynamics by formulating the problem as a system of **P**artial **D**ifferential **E**quations (PDEs) across multiple compartments. In such models, each compartment’s PDE represents the reaction-diffusion processes of biochemical species in a specific physical domain (e.g., the cytosol), while boundary conditions between compartments (e.g., membrane fluxes) govern the interactions between them.

Specifically, our program solves the mass balance equation for the species $u_{ik}$ in the $k$-th compartment $\Omega_k$ for every $k\in K$ and $i \in N_k$. Each mass balance equation is given by

$$
\begin{aligned}[c]
  \partial_t (\phi_{ik} u_{ik}) &= \nabla \cdot \mathbf{j}_{ik}+\mathcal{R}_{ik}(\mathbf{u})\quad &\text{in }\Omega_k\subset \mathbb{R}^d,\\
  u_{ik} &= u^{(0)}_{ik} \quad &\text{on }\Gamma_k^D\subseteq\partial\Omega_k,\\
  \mathbf{j}_{ik}\cdot \mathbf{n}_k &= \sum_{l\in T_k}\mathcal{T}_{ikl}(\mathbf{u}) \quad &\text{on }\partial\Omega_k \setminus \Gamma^D_{k},\\
\end{aligned}
\qquad\text{with}\quad
\begin{aligned}[c]
\mathbf{j}_{ik} &:= \sum_{j\in N_k}\mathsf{D}_{ijk}\nabla u_{jk}.
\end{aligned}
$$

Here, $\mathbf{u}_k := (u_{1k}, \ldots, u_{\dim(N_k)k})$ represents the vector of species concentrations in compartment $k$, while the full vector of species concentrations across all compartments is denoted as $\mathbf{u} := (\mathbf{u}_1, \ldots, \mathbf{u}_{\dim(K)})$. The unit outer normal vector on the boundary $\partial \Omega_k$ is represented by $\mathbf{n}_k$, and the set of neighbouring compartments to $k$ is given by $T_k := \{ l \in K : \partial \Omega_k \cap \overline{\Omega}_l \neq \varnothing \}$, indicating the compartments that share a boundary with $\Omega_k$.

The reaction operator $\mathcal{R}_{ik}(\mathbf{u})$ governs the local reaction dynamics within $\Omega_k$, while the storage terms $\phi_{ik}$ account for species accumulation in the compartment. Cross-diffusion terms $\mathsf{D}_{ijk}$ describe how species diffuse between different species within the same compartment. The non-linear transmission conditions $\mathcal{T}_{ikl}(\mathbf{u})$ represent the outflow of species $u_{ik}$ from compartment $\Omega_k$, where the outflow can either move to a neighbouring compartment $l$ (if $l \neq k$) or exit the system. Likewise, Dirichlet boundary conditions $u^{(0)}_{ik}$ are imposed on the subset $\Gamma^D_k$ of the boundary $\partial \Omega_k$, specifying the fixed concentrations of species on that portion of the boundary.

The parameters governing these equations are fully configurable at run-time, either through the command line or via a configuration file. By allowing full control over these parameters, users can adapt the software to simulate a wide variety of biochemical processes, from simple reactions in homogeneous environments to complex, multi-compartmental systems with intricate boundary conditions and interactions.

# Capabilities

Many features of DuneCopasi have been designed and developed with the specific case of systems biology in mind. Thus, a substantial effort has been made to have a library that is interoperable with systems biology data assets and requirements needed by their practitioners. Among others, include:

* a single executable configurable at run-time,
* a run-time mathematical expression parser that understands the Systems Biology Markup Language SBML [@sbml_2003] specification,
* the input of custom grid data and image files in the TIFF format for initial spatial concentrations and other parameters,
* a powerful, yet simple, boundary/transmission condition specification for each compartment to account for generic trans-membrane fluxes,
* a (non-linear) cross-diffusion specification that allows any species to cross-diffuse into other mass balance equations,
* a specification to compare results with user-defined objective functions,
* an in-place function interpolator that reduces the computational cost of evaluating common expensive reaction operations, and
* an embedded random field generator [@kempf_2023] to represent statistical spatial variations on the domain.

Furthermore, DuneCopasi is a stand-alone multi-compartment reaction-diffusion solver that may easily be used for many other fields of research and engineering fitting our model problem (e.g. \autoref{fig:cardiac_ep}).

![DuneCopasi is a highly flexible simulator that is run-time configurable using dedicated configuration files. DuneCopasi can autonomysly interpret mathematical equations, assemble the computational problems and solve these equations based on the configuration file. Here we illustrate the integration of DuneCopasi within a computational electrophysiology pipeline allowing to solve the electrophysiological mono-domain equations with Mitchell-Schaefer cell model using a medical imaging derived from biventricular geometry. \label{fig:cardiac_ep}](cardiac_ep.png)

# Statement of Need

Finite element frameworks like DUNE [@bastian_grid_2008], Deal.II [@dealII_1995], Necktar++ [@nektar_2015], or FreeFem++ [@freefem_2012] are too generic by design and don't address the specific requirements of the computational modelling practices in systems biology out of the box. In particular, most finite element frameworks don't consider the multi-compartment systems in their design, resulting, if at all possible, in inefficient simulations or obscure tricks to force this feature. Here, we extended DUNE-PDELab [@bastian_pdelab_2010;@muthing_2015] with efficient data structures especially tailored for this task.

In the space of systems biology, two well-known software packages also provide similar features to DuneCopasi, namely Morpheus [@morpheus_2014] and VCell [@vcell_1997]. Morpheus is a computational tool for multi-cellular systems which follows a cellular automata design to set rules of interactions between individual cells. It focuses on modelling approaches like cellular Potts models combined with gradients modelled by PDEs, rather than multi-compartment PDE models. Furthermore, due to its design, it falls short for moderate and big PDE computations since explicit solvers are severely limited by the required time steps. On the other hand, VCell provides a fully implicit finite volume solver on structured grids that can very well manage the multi-compartment case in addition to membrane unknowns. In comparison, our solution aims to resolve the geometry and the transmission conditions directly in the weak formulation of the problem while also being designed to provide implicit and monolithic solvers as well as tailored preconditioners for the underlying linear solvers.

# Research and Use Cases

The principal purpose of our project has been to bridge the gap between systems biology and scientific computing by providing researchers with an accessible and reproducible spatial simulator. An illustration of this is the study depicted in \autoref{fig:wehling_2022} [@wehling_2022], where simulations with the SME [@keegan_2023] aided in assessing and better understanding the mechanisms of the YAP and TAZ proteins that regulate cell proliferation in liver cells. Further examples that fit our model problem in this context have been presented by others [@elias_2014]. Additionally, DuneCopasi stands as a versatile tool capable of accommodating diverse computational needs beyond its initial focus on systems biology. For example, the general purpose nature allows for seamless integration into electrophysiology simulations, without requiring any modifications. The emerging concept of identifying patients based on personalized cardiac electrophysiology simulations [@arevalo_2016] underscores the demand for easy-to-use and efficient simulators. \autoref{fig:cardiac_ep} demonstrates how DuneCopasi emerges as a flexible solver poised to meet these evolving needs.

![Mathematical modelling predicts that nuclear phosphorylation controls spatial localization of Hippo signalling pathway components YAP and TAZ. Here, we compare two model topologies - "Model 1" (canonical) and "Model 2" (alternative model) - concerning the intracellular distribution of YAP and TAZ proteins (Pr) and their phosphorylated counterparts (pPr). If the phosphorylation of YAP/TAZ takes place exclusively outside the nucleus, as shown in "Model 1", PDE simulation indicates low spatial accordance with the experimentally measured subcellular localization of YAP/TAZ.  Whereas, "Model 2" describes YAP/TAZ protein phosphorylation and dephosphorylation in the nucleus. The simulation of "Model 2" agrees with experimentally measured subcellular distribution of YAP/TAZ proteins, as reported in [@wehling_2022]. \label{fig:wehling_2022}](wehling_2022.pdf)

# Acknowledgements

We want to thank all the contributions that aided in the development and deployment of the package. Special thanks to Ursula Kummer from the COPASI team for valuable guidance. This work has been funded and supported by the German Federal Ministry of Education and Research (BMBF) FKZ 031L0158.

# References
