
set(DUNE_COPASI_GRID_DIMENSIONS 1 2 3 CACHE STRING "List of grid dimension for models")
set(DUNE_COPASI_PRECOMPILED_MODE TRUE CACHE BOOL "If dune-copasi should be pre-compile object libraries for model classes")

foreach(dim ${DUNE_COPASI_GRID_DIMENSIONS})
  set(DUNE_COPASI_${dim}D_FEM_ORDERS 1 CACHE STRING "List of polynomial order Pk of the finite elements for the ${dim}D library")
  if(DUNE_COPASI_PRECOMPILED_MODE)
    foreach(order ${DUNE_COPASI_${dim}D_FEM_ORDERS})
      set(target dune-copasi-lib-${dim}D-P${order})
      dune_add_library(${target} STATIC EXPORT_NAME Copasi::ModelLibraryP${order}${dim}D)
      target_link_libraries(${target}
        PRIVATE
          $<TARGET_NAME_IF_EXISTS:SuiteSparse::SuiteSparse>
          $<TARGET_NAME_IF_EXISTS:SuperLU::SuperLU>
      )
      target_compile_definitions(${target}
        PRIVATE
          $<$<TARGET_EXISTS:SuiteSparse::SuiteSparse>:ENABLE_SUITESPARSE=1>
          $<$<TARGET_EXISTS:SuperLU::SuperLU>:ENABLE_SUPERLU=1>
      )
      target_precompile_headers(${target}
        PUBLIC
          $<BUILD_INTERFACE:dune/copasi/model/diffusion_reaction_sc.hh>
          $<BUILD_INTERFACE:dune/copasi/model/diffusion_reaction_mc.hh>
      )
      target_sources(${target} PRIVATE diffusion_reaction_sc.cc)
      target_sources(${target} PRIVATE diffusion_reaction_mc.cc)
      target_compile_definitions(${target}
        INTERFACE
          DUNE_COPASI_PRECOMPILED_MODE
          PRIVATE
          DUNE_COPASI_GRID_DIMENSION=${dim}
          DUNE_COPASI_FEM_ORDER=${order})
      target_link_libraries(${target}
        PRIVATE Dune::Copasi::Common
        PUBLIC  ${DUNE_LIBS})
      target_link_libraries(dune-copasi INTERFACE ${target})
    endforeach()
  endif()
endforeach()
