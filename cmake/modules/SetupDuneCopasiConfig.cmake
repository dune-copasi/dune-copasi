# This settings define the installed cmake config file. It adds a minimal setup
# for the dune requirements without cluttering the whole project with old CMake

set_property(GLOBAL PROPERTY DUNE_MODULE_LIBRARIES dune-copasi)

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
define_property(GLOBAL PROPERTY DUNE_MODULE_LIBRARIES
  BRIEF_DOCS \"List of libraries of the module. DO NOT EDIT!\"
  FULL_DOCS \"List of libraries of the module. Used to generate CMake's package configuration files. DO NOT EDIT!\")

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
")
