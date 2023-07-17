/* begin dune-copasi
   put the definitions for config.h specific to
   your project here. Everything above will be
   overwritten
*/

/* Defines the version of dune-copasi */
#define DUNE_COPASI_VERSION "@DUNE_COPASI_VERSION@"

/* Defines the major version of dune-copasi */
#define DUNE_COPASI_VERSION_MAJOR @DUNE_COPASI_VERSION_MAJOR@

/* Defines the minor version of dune-copasi */
#define DUNE_COPASI_VERSION_MINOR @DUNE_COPASI_VERSION_MINOR@

/* Defines the revision of dune-copasi */
#define DUNE_COPASI_VERSION_REVISION @DUNE_COPASI_VERSION_REVISION@

// std::move_only_function with noexcept attribute seems to be broken in clang
#if defined(__clang__)
#define DUNE_COPASI_FUNCTOR_NOEXCEPT
#else
#define DUNE_COPASI_FUNCTOR_NOEXCEPT noexcept
#endif

/* end dune-copasi
   Everything below here will be overwritten
*/
