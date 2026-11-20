/* Static build into a single executable: no DLL exports. */
#ifndef OOFEM_EXPORT_H
#define OOFEM_EXPORT_H
#define OOFEM_EXPORT
#define OOFEM_NO_EXPORT
#define OOFEM_DEPRECATED __attribute__ ((__deprecated__))
#define OOFEM_DEPRECATED_EXPORT OOFEM_DEPRECATED
#define OOFEM_DEPRECATED_NO_EXPORT OOFEM_DEPRECATED
#endif
