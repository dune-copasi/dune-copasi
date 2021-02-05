# This settings define the installed cmake config file. It adds a minimal setup
# for the dune requirements without cluttering the whole project with old CMake

set(DUNE_CUSTOM_PKG_CONFIG_SECTION
"
include(CMakeFindDependencyMacro)
find_dependency(dune-common REQUIRED)
list(APPEND CMAKE_MODULE_PATH \${dune-common_MODULE_PATH} \${dune-copasi_MODULE_PATH})
include(DuneMacros)

include(CheckCXXFeatures)
include(DuneCxaDemangle)
include(DuneMPI)

set(DUNE_PYTHON_VIRTUALENV_SETUP ${DUNE_PYTHON_VIRTUALENV_SETUP})

set(ProjectName dune-copasi)

find_file(\${ProjectName}_MODULE_DIR dune.module
  HINTS
    \"\${\${ProjectName}_PREFIX}\"
    \"\${\${ProjectName}_PREFIX}/lib/dunecontrol/\${ProjectName}\"
    \"\${\${ProjectName}_PREFIX}/lib64/dunecontrol/\${ProjectName}\"
  REQUIRED
  NO_DEFAULT_PATH
)
string(REPLACE dune.module \"\" \${ProjectName}_MODULE_DIR \${\${ProjectName}_MODULE_DIR})

dune_module_information(\${\${ProjectName}_MODULE_DIR})
dune_create_dependency_tree()
dune_process_dependency_macros()
unset(ProjectName)

# get dune-copasi dependencies
find_dependency(muparser)
find_dependency(TIFF)
find_dependency(Filesystem)
if(${USE_FALLBACK_FILESYSTEM})
  find_package(ghc_filesystem)
endif()

#import the target
get_filename_component(_dir \"\${CMAKE_CURRENT_LIST_FILE}\" PATH)
include(\"\${_dir}/dune-copasi-targets.cmake\")
")
