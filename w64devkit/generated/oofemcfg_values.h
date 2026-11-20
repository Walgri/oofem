/* Version macros that CMake would pass on the command line (used by oofemcfg.C only). */
#define __OOFEM_MAJOR_VERSION "2"
#define __OOFEM_MINOR_VERSION "5"
#define __OOFEM_VERSION "2.5-companionfea"
#define __OOFEM_COPYRIGHT "Copyright (C) 1994-2025 Borek Patzak"
#define __OOFEM_GIT_HASH "companionfea"
#define __OOFEM_GIT_REPOURL "https://github.com/Walgri/oofem.git"
#define __OOFEM_GIT_BRANCH "companionfea"
#if defined(_WIN32)
#define __HOST_TYPE "x86_64-Windows"
#else
#define __HOST_TYPE "x86_64-Linux"
#endif
#define __HOST_NAME "w64devkit"
#define __MODULE_LIST "fm tm sm mpm dss"
