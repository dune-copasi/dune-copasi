# packaging rules
include(InstallRequiredSystemLibraries)

set(CPACK_GENERATOR "STGZ;TGZ;TZ")
set(CPACK_PACKAGE_DESCRIPTION_SUMMARY "Solver for reaction-diffusion systems in multiple compartments")
set(CPACK_PACKAGE_VENDOR "IWR, Universität Heidelberg")
set(CPACK_PACKAGE_CONTACT "santiago@dune-project.org")
set(CPACK_RESOURCE_FILE_LICENSE "${CMAKE_CURRENT_SOURCE_DIR}/LICENSE.md")
set(CPACK_PACKAGE_VERSION_MAJOR ${DUNE_VERSION_MAJOR})
set(CPACK_PACKAGE_VERSION_MINOR ${DUNE_VERSION_MINOR})
set(CPACK_PACKAGE_VERSION_PATCH ${DUNE_VERSION_REVISION})
set(CPACK_PACKAGE_DIRECTORY "${CMAKE_CURRENT_BINARY_DIR}/packages")

set(CPACK_DEB_COMPONENT_INSTALL ON)
set(CPACK_DEBIAN_PACKAGE_SHLIBDEPS ON)
set(CPACK_DEBIAN_PACKAGE_SECTION science)
set(CPACK_DEBIAN_PACKAGE_HOMEPAGE "https://gitlab.dune-project.org/copasi/dune-copasi/")
if(${MPI_FOUND})
  set(CPACK_DEBIAN_PACKAGE_DEPENDS "openmpi-bin")
endif()

include(CPack)
