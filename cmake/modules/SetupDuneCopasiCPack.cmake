# packaging rules
include(InstallRequiredSystemLibraries)

set(CPACK_GENERATOR "STGZ;TGZ;TZ")
set(CPACK_PACKAGE_DESCRIPTION_SUMMARY "Solver for reaction-diffusion systems in multiple compartments")
set(CPACK_PACKAGE_VENDOR "IWR, Universität Heidelberg")
set(CPACK_PACKAGE_CONTACT "santiago.ospina@iwr.uni-heidelberg.de")
# set(CPACK_PACKAGE_DESCRIPTION_FILE "${CMAKE_CURRENT_SOURCE_DIR}/README.md")
set(CPACK_RESOURCE_FILE_LICENSE "${CMAKE_CURRENT_SOURCE_DIR}/LICENSE.md")
set(CPACK_PACKAGE_VERSION_MAJOR ${DUNE_VERSION_MAJOR})
set(CPACK_PACKAGE_VERSION_MINOR ${DUNE_VERSION_MINOR})
set(CPACK_PACKAGE_VERSION_PATCH ${DUNE_VERSION_REVISION})
set(CPACK_PACKAGE_DIRECTORY "${CMAKE_CURRENT_BINARY_DIR}/packages")

set(CPACK_DEB_COMPONENT_INSTALL ON)
set(CPACK_DEBIAN_PACKAGE_SHLIBDEPS ON)
set(CPACK_DEBIAN_PACKAGE_SECTION science)
set(CPACK_DEBIAN_PACKAGE_HOMEPAGE "https://gitlab.dune-project.org/copasi/dune-copasi/")
set(CPACK_DEBIAN_PACKAGE_DEPENDS "openmpi-bin")
# set(CPACK_DEBIAN_Runtime_PACKAGE_DEPENDS "libscotchparmetis-dev,libldl2,libspqr2,libumfpack5,libarpack++2c2a,libsuperlu5,libgmpxx4ldbl,libopenblas-base,libtiff5,libmuparser2v5")

include(CPack)
