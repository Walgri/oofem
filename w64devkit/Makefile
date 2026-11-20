# OOFEM build without CMake: a single executable, no external dependencies.
# Modules as in the default CMake configuration (sm, tm, fm, mpm, XML input) plus DSS and IML.
#
#   Windows, w64devkit:   open w64devkit.exe, enter this folder and run
#                             make -j8
#                         result: oofem.exe (static: copy it anywhere, no DLLs needed)
#   Linux:                make -j8      (same Makefile, result: oofem)
#
#   make test       runs the regression tests (run_tests.sh): upstream tests and the tests of this fork
#   make benchmark  same, plus the larger models in benchmark/
#   make clean      removes the object files
#
# Notes:
#  - all object files are linked directly into the executable (no .a library): OOFEM elements,
#    materials and solvers register themselves through static objects, which the linker would
#    drop from a library;
#  - linking uses a response file (objects.txt) because the command line with 600 files
#    exceeds the Windows limit.

S       := ..
OBJ     := obj
CXX     := g++
CC      := gcc
OPT     := -O2

ifeq ($(OS),Windows_NT)
  EXE    := oofem.exe
  LDFLAGS := -static -s
else
  EXE    := oofem
  LDFLAGS :=
endif

DEF      := -DNDEBUG -D__SM_MODULE -D__TM_MODULE -D__FM_MODULE -D__MPM_MODULE -D_USE_XML -DOOFEM_STATIC_DEFINE
INC      := -Igenerated -I$(S)/src -I$(S)/src/core -I$(S)/src/core/T3d/src -I$(S)/src/fm -I$(S)/src/tm \
            -I$(S)/src/sm -I$(S)/src/sm/Elements -I$(S)/src/sm/Materials -I$(S)/src/dss -I$(S)/src/mpm \
            -I$(S)/extern/pugixml/src -I$(S)
CXXFLAGS := -std=c++17 $(OPT) -w -MMD -MP $(DEF) $(INC)
CFLAGS   := -std=c99 $(OPT) -w -MMD -MP $(DEF) $(INC)

include sources.mk
OBJECTS := $(addprefix $(OBJ)/,$(addsuffix .o,$(basename $(SOURCES))))

$(EXE): $(OBJECTS)
	$(file >objects.txt,$(OBJECTS))
	$(CXX) -o $@ @objects.txt $(LDFLAGS)
	@echo "done: $@"

$(OBJ)/%.o: $(S)/%.C
	@mkdir -p $(dir $@)
	$(CXX) $(CXXFLAGS) -c $< -o $@

$(OBJ)/%.o: $(S)/%.cpp
	@mkdir -p $(dir $@)
	$(CXX) $(CXXFLAGS) -c $< -o $@

$(OBJ)/%.o: $(S)/%.c
	@mkdir -p $(dir $@)
	$(CC) $(CFLAGS) -c $< -o $@

# version macros live in a file, not on the command line (no quotes to escape)
$(OBJ)/src/core/oofemcfg.o: CXXFLAGS += -include generated/oofemcfg_values.h

.PHONY: clean test benchmark
clean:
	rm -rf $(OBJ) objects.txt $(EXE)

test: $(EXE)
	sh run_tests.sh ./$(EXE)

benchmark: $(EXE)
	sh run_tests.sh --benchmark ./$(EXE)

-include $(OBJECTS:.o=.d)
