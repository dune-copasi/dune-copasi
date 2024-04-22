---
id: output
title: Output
sidebar_label: Output
---

The output of `dune_copasi_md` is a set of `vtk` files for the timesteps and the
different compartments. The compendium of all timesteps is written in a `pvd`
file on the same directory where the executables was called, while the actual
data files for each timestep is written on folders following the
`[model.writer.file_path]` value.

The `pvd` files can be opened by [ParaView][paraview] or [VisiIt][visit] to
inspect the results.

[paraview]: https://www.paraview.org/
[visit]: https://wci.llnl.gov/simulation/computer-codes/visit/
