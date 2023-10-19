Membranes:
  * Local operator coupling between compartment-membrane
    - The mass residual for membranes seems to be wrong. Strange values!
  * Set up boundary membranes

Pattern generation seems to be broken if estimate is too low
Parameters for solvers not ported yet
Use a preconditioner
Change logger to spdlog
Homogenize use of 'compartment' and 'domain' variable names
Use UMFPack for smal problems
Add block preconditioner (maybe store LU decpmposition or copy sparse to dense and back solve)
Check resize (pattern creation) on newton derivative -> is because step operator is created every time anew
