[![Build Status](https://gitlab.dune-project.org/copasi/dune-copasi/badges/master/pipeline.svg)](https://gitlab.dune-project.org/copasi/dune-copasi/pipelines)
[![Build Status](https://travis-ci.org/dune-copasi/dune-copasi.svg?branch=master)](https://travis-ci.org/dune-copasi/dune-copasi)
[![Build status](https://ci.appveyor.com/api/projects/status/e7w7u5dt50kue5sb/branch/master?svg=true)](https://ci.appveyor.com/project/SoilRos/dune-copasi-3gv18/branch/master)
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
