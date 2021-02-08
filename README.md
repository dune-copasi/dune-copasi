[![Build Status](https://gitlab.dune-project.org/copasi/dune-copasi/badges/master/pipeline.svg)](https://gitlab.dune-project.org/copasi/dune-copasi/pipelines)
[![Build Status](https://github.com/dune-copasi/dune-copasi/workflows/CI%20Builds/badge.svg?branch=master)](https://github.com/dune-copasi/dune-copasi/actions?query=branch%3Amaster+)
[![Netlify Status](https://api.netlify.com/api/v1/badges/6fc6d371-87df-49b5-8e72-e1873fa5d54b/deploy-status)](https://app.netlify.com/sites/dune-copasi/deploys)

# DuneCopasi

Solver for reaction-diffusion systems in multiple compartments

 * Solve a reaction-diffusion system for each comartment
 * Each compartment may have different system with different number of variables
 * Neumann flux at the interface of compartments for variables with
   the same name on the two compartments
 * Easy to modify configuration file
 * Initial conditions can be a TIFF file or/and a math expression
 * Using the finite element
 * Output in the VTK format

Get started [here](https://dune-copasi.netlify.app/docs/).
