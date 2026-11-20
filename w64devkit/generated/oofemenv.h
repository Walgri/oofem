/* Configuration for the build without CMake (w64devkit / gcc + make). */
#ifndef oofemenv_h
#define oofemenv_h
#include <cstddef>
#if !defined(_WIN32)
#define HAVE_EXECINFO_H 1      /* backtrace: glibc only, not available in MinGW */
#endif
#define HAVE_CBRT 1
typedef std::size_t indexType;
#include "oofem_export.h"
#endif
