
add_library(dune-copasi-model INTERFACE)
add_library(Dune::Copasi::Model ALIAS dune-copasi-model)
set_target_properties(dune-copasi-model PROPERTIES EXPORT_NAME Copasi::Model)

find_package(METIS)

set(DUNE_COPASI_GRID_DIMENSIONS "2" CACHE STRING "List of grid dimension for models")
option(DUNE_COPASI_CONCURRENT_ASSEMBLY "Whether dune-copasi should be compiled with support for concurrent for assembly" ON)
option(DUNE_COPASI_PRECOMPILED_MODE "Whether dune-copasi should be pre-compile object libraries for model classes" ON)

set(model_libraries
  $<TARGET_NAME_IF_EXISTS:SuiteSparse::UMFPACK>
  $<TARGET_NAME_IF_EXISTS:SuperLU::SuperLU>
  $<TARGET_NAME_IF_EXISTS:METIS::METIS>
  $<TARGET_NAME_IF_EXISTS:PTScotch::Scotch>
)

set(model_definitions
  $<$<TARGET_EXISTS:SuiteSparse::UMFPACK>:HAVE_SUITESPARSE_UMFPACK=1>
  $<$<TARGET_EXISTS:SuperLU::SuperLU>:ENABLE_SUPERLU=1>
  $<$<TARGET_EXISTS:METIS::METIS>:HAVE_METIS=1>
  $<$<TARGET_EXISTS:PTScotch::Scotch>:HAVE_SCOTCH_METIS=1>
  $<$<BOOL:${DUNE_COPASI_CONCURRENT_ASSEMBLY}>:DUNE_COPASI_CONCURRENT_ASSEMBLY=1>
)

target_link_libraries(dune-copasi-model INTERFACE ${model_libraries})
target_compile_definitions(dune-copasi-model INTERFACE ${model_definitions})

foreach(dim ${DUNE_COPASI_GRID_DIMENSIONS})
  set(DUNE_COPASI_${dim}D_DIFFUSION_REACTION_FEM_ORDERS 1 CACHE STRING "List of polynomial order Pk of the finite elements for the ${dim}D library")
  if(DUNE_COPASI_PRECOMPILED_MODE)
    foreach(order ${DUNE_COPASI_${dim}D_DIFFUSION_REACTION_FEM_ORDERS})
      set(target dune-copasi-lib-${dim}D-P${order})
      add_library(${target} STATIC)
      add_library(Dune::Copasi::ModelLibraryP${order}${dim}D ALIAS ${target})
      set_target_properties(${target} PROPERTIES EXPORT_NAME Copasi::ModelLibraryP${order}${dim}D)
      target_sources(${target} PRIVATE model_single_compartment.cc)
      target_sources(${target} PRIVATE model_multi_compartment.cc)
      target_compile_definitions(${target}
        INTERFACE
          DUNE_COPASI_PRECOMPILED_MODE
        PUBLIC
          ${model_definitions}
        PRIVATE
          DUNE_COPASI_GRID_DIMENSION=${dim}
          DUNE_COPASI_FEM_ORDER=${order}
      )
      target_link_libraries(${target}
        PRIVATE
          ${model_libraries}
          Dune::Copasi::Common
        PUBLIC
          ${DUNE_LIBS}
      )

      target_link_libraries(dune-copasi-model INTERFACE ${target})
      install(TARGETS ${target} EXPORT dune-copasi-targets ARCHIVE COMPONENT Development)
    endforeach()

    list(JOIN DUNE_COPASI_${dim}D_DIFFUSION_REACTION_FEM_ORDERS "," orders)
    target_compile_definitions(dune-copasi-model
      INTERFACE
        DUNE_COPASI_${dim}D_DIFFUSION_REACTION_FEM_ORDERS=${orders}
    )
  endif()
endforeach()

# target configuration for downstream projects
install(TARGETS dune-copasi-model EXPORT dune-copasi-targets ARCHIVE COMPONENT Development)
