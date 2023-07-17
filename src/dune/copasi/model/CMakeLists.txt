
set(DUNE_COPASI_GRID_DIMENSIONS 2 CACHE STRING "List of grid dimension for models")
set(DUNE_COPASI_PRECOMPILED_MODE TRUE CACHE BOOL "If dune-copasi should be pre-compile object libraries for model classes")

foreach(dim ${DUNE_COPASI_GRID_DIMENSIONS})
  set(DUNE_COPASI_${dim}D_MAX_FEM_ORDER 1 CACHE STRING "Maximum polynomial order Pk of the finite elements for the ${dim}D library")
  if(DUNE_COPASI_PRECOMPILED_MODE)
    foreach(order RANGE 1 ${DUNE_COPASI_${dim}D_MAX_FEM_ORDER})
      set(target dune-copasi-lib-${dim}D-P${order})
      dune_add_library(${target} STATIC EXPORT_NAME Copasi::ModelLibraryP${order}${dim}D)
      target_precompile_headers(${target}
        PUBLIC
          $<BUILD_INTERFACE:dune/copasi/model/diffusion_reaction_sc.hh>
          $<BUILD_INTERFACE:dune/copasi/model/diffusion_reaction_mc.hh>
      )
      target_sources(${target} PRIVATE diffusion_reaction_sc.cc)
      # target_sources(${target} PRIVATE diffusion_reaction_mcm.cc)
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
