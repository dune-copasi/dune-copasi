---
id: math_model
title: Mathematical Model
sidebar_label: Math Model
---

In the following, we decribe (in rather mathemathical terms) the actual
model problem that we want to solve.

## Continuous Problem

In principle, we want to describe processes that are *continuously* observable in
time and space (e.g.
[continuum mechanics](https://en.wikipedia.org/wiki/Continuum_mechanics#Concept_of_a_continuum)).
Additionally, we are interested in systems of reaction-diffusion equations that
compartmentalized. That is, where different compartments have different
reaction-diffusion processes and the interaction between compartments happens
*through* their compartment intesection (e.g. a
[living cell](https://en.wikipedia.org/wiki/Cell_(biology))).

### Compartments

We say that each compartments, denoted by $\Omega^k$, is an open, bounded,
connected subset of Euclidean space $\mathbb{R}^d$, with Lipschitz boundary
$\partial\Omega^k$ and space dimension $d=2$ or $d=3$. Moreover, we require
that the compartments are the non-overlapping decomposition of an open, bounded,
and connected domain $\Omega\subset\mathbb{R}^d$, such that,
$$
  \overline{\Omega}=\bigcup_{k=1}^{K}\overline{\Omega}^k\quad\text{and}\quad\Omega^k\cap\Omega^l=\emptyset\quad \text{with}\quad k\ne l,
$$
where $K$ is the total number of compartments.

### Species
We will distinguish species within the same compartment with a subscript,
usually $i$, and species on different compartments with a superscript, usually
$k$. For example, $u_i^k$ is the $i$-th species on the $k$-th compartment.
Furthermore, we will denote bold letters to refer to all $N^k$ species in the
$k$-th compartment, i.e, $\bm{u}^k:=[u^k_1,\dots,u_{N^k}^k]^T$.

### Transport

We then characterize transport with a **diffusive flux operator**,
$\mathcal{D}_i^k$, for the $i$-th species by
$$
\mathcal{D}^k_i\left(\bm{u}^k\right):=-\mathsf{D}^k_{i\star} \nabla \bm{u}^k,
$$
where $\mathsf{D}^k_{i\star}$ represents the $i$-th row of the self and
cross-dispersion tensor $\mathsf{D}^k$, in other words,
$\mathsf{D}^k_{i\star}\nabla\bm{u}^k:=\sum_{i=1}^{N^k}\mathsf{D}^k_{ij}\nabla u^k_j$.
Ultimately, applying the law of conservation of mass we get

### Reaction Network

A (bio-)chemical reaction is the transformation of one species to another, i.e,
$u_i^k\to u_j^k$. There are many natural processes that present a rich reaction
network. Particularly, we will represent the *deterministic* rate of change of
$u^k_i$ caused by chemical reactions with the **reaction operator**
$\mathcal{R}_i^k\left(\bm{u}^k\right)$. Such an operator is required to be a
Lipschitz function and be mass conservative.

### Membrane

We call the boundaries of the compartments as *membranes*. In particular, we
designate its geometry to be a $(d-1)$-manifold defined with respect to the
boundary of the compartments, i.e.,
$$
\Gamma^{kl}:=
  \left\{
  \begin{matrix}
    \partial\Omega_k\cap\partial\Omega_l, \quad& \text{if }k\ne l\quad& \text{interior boundaries}\\
    \partial\Omega_k\cap\partial\Omega,   \quad& \text{if }k =  l\quad& \text{exterior boundaries}.
  \end{matrix}
  \right.
$$
where $\partial\Omega$ represents the boundaries of a domain $\Omega$.


### Transmission Conditions

The flux rate at which compartment species are transformed and moved across the
membrane depends on its concentration and the dispersion coefficients at
which species can move on the surroundings of the membrane. In our model, we say
that such a transport is equal to the outer flux $\mathcal{D}_i^k$ of the
species $i$ leaving the compartment $k$ and is defined by a general
**transmission condition**, $\mathcal{T}^{kl}_i$, that may take the form of
typical Dirichlet, Neumann, and Robin boundary condition as well as complex
chemical reaction networks, i.e.,
$$
\mathcal{D}_i^{k}\left(\bm{u}^k\right)\cdot\mathbf{n}^k = \mathcal{T}_i^{kl}\left(\bm{u}^k,\bm{u}^l\right)\qquad \text{on }\Gamma^{kl},
$$
where $\mathbf{n}^k$ is the outer normal vector on $\Omega^k$. As is natural,
the transmission conditions must be mass conservative.

### Strong Formulation

Joining all definitions from above, we obtain a *boundary value problem* (BVP)
which it reads as follow:

Given an initial condition ${u_{(0)}}_i^k$ and a final time $T$, find $u_i^k$
such that

$$
\begin{aligned}
\partial_t u_i^k = -\nabla\cdot\mathcal{D}_i^k\left(\bm{u}^k\right) + \mathcal{R}_i^k\left(\bm{u}^k\right) &\qquad \text{in }\Omega^k\times(0,T),\\
\mathcal{D}_i^k\left(\bm{u}^k\right)\cdot \mathbf{n}^k = \mathcal{T}^{kl}_i(\bm{u}^k,\bm{u}^l) &\qquad \text{on }\Gamma^{kl}\times(0,T),\\
u_i^k = {u_{(0)}}_i^k &\qquad \text{in }\Omega^k,
\end{aligned}
$$
for every $k,l=1,\ldots,K$, and $i=1,\ldots,N^k$.
