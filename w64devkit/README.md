# Building OOFEM on Windows with w64devkit (no CMake)

This folder builds a single static `oofem.exe` with the modules of the default CMake configuration
(`sm`, `tm`, `fm`, `mpm`, XML input) plus `dss` and `iml`, using only
[w64devkit](https://github.com/skeeto/w64devkit) (GCC + make). No CMake, Visual Studio, MSYS2 or
external libraries are needed. The same Makefile works on Linux.

1. Download `w64devkit-x64-<version>.7z.exe` (not the "mini" variant, which has no C++) and
   extract it, e.g. to `C:\w64devkit`.
2. Put the OOFEM sources in a folder without spaces in the path, e.g. `C:\fem\oofem`. The XML
   reader needs the `extern/pugixml` submodule: clone with `git clone --recursive`, or run
   `git submodule update --init extern/pugixml` in an existing clone.
3. Start `C:\w64devkit\w64devkit.exe` and run

       cd C:/fem/oofem/w64devkit
       make -j8

   About 700 files, 20 minutes of CPU time. The result, `oofem.exe`, is static: copy it anywhere.
4. Run the tests:

       make test          # tests/regression, as CTest (except the catch2 unit tests)
       make benchmark     # the same, plus the larger models in benchmark/ (about 64000 DOFs)

   Both must end with `failed: 0`. The list of failed tests is in `tests_failed.txt`, the log of
   each test in `tests_scratch/`.

Files:

| File | Content |
|---|---|
| `Makefile` | build rules; all objects are linked directly (static registration of classes) |
| `sources.mk` | list of compiled sources (same set as CMake with default modules plus DSS) |
| `generated/` | headers that CMake would generate (`oofemenv.h`, `oofem_export.h`, version macros) |
| `run_tests.sh` | test runner (busybox sh or Linux sh) |
| `benchmark/` | larger regression and performance models, not part of the regular test suite |

Troubleshooting:

- a compile error with a newer GCC (missing implicit includes): `make -j8 OPT="-O2 -fpermissive"`
  as a workaround, and report the message;
- a very slow build: the antivirus scans every created file; exclude the folder;
- `make clean` removes all objects.
