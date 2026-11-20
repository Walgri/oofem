#!/bin/sh
# Runs the OOFEM regression tests for the modules built by this Makefile, following the rules of
# tests/CMakeLists.txt (CTest): for each group, tests/regression/<group>/*.in and, with XML support,
# *.xml; tests/regression/sm/*.sh (workflow tests); the special mpm tests. The catch2 unit tests in
# tests/unit are not built by this Makefile (CTest runs them).
#
# The benchmarks in tests/regression/benchmark/<group> are included, as CTest registers them too.
#
# A .in/.xml test passes when "oofem -Werror -f <test>" exits with status 0: the #%BEGIN_CHECK% block
# in the input is verified by OOFEM itself (errorcheck module). A .sh test receives the executable as
# $1 and -Werror as $2, and passes when it exits with status 0.
#
#   sh run_tests.sh [--benchmark] [oofem executable] [parallel jobs]
#
#   --benchmark   also runs the larger models in w64devkit/benchmark (about 64000 DOFs, ~30 s)
#
# The tests run in a scratch copy (the sources stay clean). Results:
#   tests_passed.txt, tests_failed.txt, and the log of every test in tests_scratch/<test>.log
#
# Works in w64devkit (busybox sh) and on Linux.

HERE=$(cd "$(dirname "$0")" && pwd)
BENCH=0
if [ "$1" = "--benchmark" ]; then BENCH=1; shift; fi
EXE=${1:-$HERE/oofem}
[ -x "$EXE" ] || EXE=$HERE/oofem.exe
case "$EXE" in /*|[A-Za-z]:*) ;; *) EXE=$(pwd)/$EXE ;; esac
JOBS=${2:-2}
TIMEOUT=${TEST_TIMEOUT:-600}
SCRATCH=$HERE/tests_scratch
REG=$HERE/../tests/regression

# groups built by this Makefile (see DEF in the Makefile), as in tests/CMakeLists.txt
GROUPS="sm sm/NAFEMS fm tm tmsm tmfm mpm tmsmdssiml"
# tests that CTest runs without -Werror (BARE_CMD)
BARE=" adapt02.in "
# tests registered separately by CTest (special mpm tests below)
EXCLUDE=" mpms_err_dualto.in mpms_warn_dofquantity.in "

rm -rf "$SCRATCH"
mkdir -p "$SCRATCH/benchmark"
for g in $GROUPS; do
    case "$g" in */*) ;; *) [ -d "$REG/$g" ] && cp -r "$REG/$g" "$SCRATCH/$g" ;; esac
done
for g in $GROUPS; do
    case "$g" in */*) ;; *) [ -d "$REG/benchmark/$g" ] && cp -r "$REG/benchmark/$g" "$SCRATCH/benchmark/$g" ;; esac
done
if [ $BENCH -eq 1 ]; then
    cp -r "$HERE/benchmark" "$SCRATCH/benchmark/w64devkit"
fi

run_one() {
    # $1 = path relative to $SCRATCH, e.g. sm/beam2d_1.in
    dir=$(dirname "$1"); name=$(basename "$1")
    log="$SCRATCH/$1.log"
    case "$name" in
        *.sh)
            ( cd "$SCRATCH/$dir" && timeout "$TIMEOUT" sh "$name" "$EXE" -Werror ) > "$log" 2>&1 ;;
        mpms_err_dualto.in)
            ( cd "$SCRATCH/$dir" && timeout "$TIMEOUT" "$EXE" -f "$name" 2>&1 | grep -q 'supplied by more than one unknown' ) > "$log" 2>&1 ;;
        mpms_warn_dofquantity.in)
            ( cd "$SCRATCH/$dir" && timeout "$TIMEOUT" "$EXE" -f "$name" 2>&1 | tee "$log.out" | grep -q 'conventionally carries' ) > "$log" 2>&1 ;;
        *)
            werror=-Werror
            case "$BARE" in *" $name "*) werror= ;; esac
            ( cd "$SCRATCH/$dir" && timeout "$TIMEOUT" "$EXE" $werror -f "$name" ) > "$log" 2>&1 ;;
    esac
    if [ $? -eq 0 ]; then echo "PASS $1"; else echo "FAIL $1"; fi
}

LIST=$SCRATCH/list.txt
: > "$LIST"
cd "$SCRATCH" || exit 1
for g in $GROUPS; do
    [ -d "$g" ] || continue
    for f in "$g"/*.in "$g"/*.xml; do
        [ -f "$f" ] || continue
        case "$f" in *.fragment.xml|*.frag.xml) continue ;; esac
        case "$EXCLUDE" in *" $(basename "$f") "*) continue ;; esac
        echo "$f" >> "$LIST"
    done
done
ls sm/*.sh >> "$LIST" 2>/dev/null
[ -d mpm ] && echo mpm/mpms_err_dualto.in >> "$LIST" && echo mpm/mpms_warn_dofquantity.in >> "$LIST"
for g in $GROUPS; do
    ls benchmark/"$g"/*.in 2>/dev/null >> "$LIST"
done
if [ $BENCH -eq 1 ]; then
    ls benchmark/w64devkit/*.in >> "$LIST"
fi
cd "$HERE" || exit 1
N=$(wc -l < "$LIST")
echo "running $N tests with $JOBS jobs..."

# simple parallel loop (xargs -P is not available everywhere)
RESULTS=$SCRATCH/results.txt
: > "$RESULTS"
i=0
while read -r t; do
    run_one "$t" >> "$RESULTS" &
    i=$((i + 1))
    if [ $((i % JOBS)) -eq 0 ]; then wait; fi
done < "$LIST"
wait

sort -k2 "$RESULTS" | grep '^PASS' | cut -c6- > "$HERE/tests_passed.txt"
sort -k2 "$RESULTS" | grep '^FAIL' | cut -c6- > "$HERE/tests_failed.txt"
P=$(wc -l < "$HERE/tests_passed.txt"); F=$(wc -l < "$HERE/tests_failed.txt")
echo "passed: $P   failed: $F   (of $N)"
if [ "$F" -gt 0 ]; then
    echo "failed tests (log in $SCRATCH/<test>.log):"
    sed 's/^/  /' "$HERE/tests_failed.txt"
fi
[ "$F" -eq 0 ]
