
add_library(dune-copasi-model INTERFACE)
add_library(Dune::Copasi::Model ALIAS dune-copasi-model)
set_target_properties(dune-copasi-model PROPERTIES EXPORT_NAME Copasi::Model)

target_link_libraries(dune-copasi-model
  INTERFACE
    $<TARGET_NAME_IF_EXISTS:SuiteSparse::SuiteSparse>
    $<TARGET_NAME_IF_EXISTS:SuperLU::SuperLU>
)
target_compile_definitions(dune-copasi-model
  INTERFACE
    $<$<TARGET_EXISTS:SuiteSparse::SuiteSparse>:ENABLE_SUITESPARSE=1>
    $<$<TARGET_EXISTS:SuperLU::SuperLU>:ENABLE_SUPERLU=1>
)

set(DUNE_COPASI_GRID_DIMENSIONS 1 2 3 CACHE STRING "List of grid dimension for models")
set(DUNE_COPASI_PRECOMPILED_MODE OFF CACHE BOOL "If dune-copasi should be pre-compile object libraries for model classes")

foreach(dim ${DUNE_COPASI_GRID_DIMENSIONS})
  set(DUNE_COPASI_${dim}D_FEM_ORDERS 1 CACHE STRING "List of polynomial order Pk of the finite elements for the ${dim}D library")
  if(DUNE_COPASI_PRECOMPILED_MODE)
    foreach(order ${DUNE_COPASI_${dim}D_FEM_ORDERS})
      set(target dune-copasi-lib-${dim}D-P${order})
      add_library(${target} STATIC)
      add_library(Dune::Copasi::ModelLibraryP${order}${dim}D ALIAS ${target})
      set_target_properties(${target} PROPERTIES EXPORT_NAME Copasi::ModelLibraryP${order}${dim}D)
      add_dune_suitesparse_flags(${target})
      add_dune_superlu_flags(${target})
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

      target_link_libraries(dune-copasi-model INTERFACE ${target})
      install(TARGETS ${target} EXPORT dune-copasi-targets ARCHIVE COMPONENT Development)
    endforeach()
  endif()
endforeach()

# target configuration for downstream projects
install(TARGETS dune-copasi-model EXPORT dune-copasi-targets ARCHIVE COMPONENT Development)
