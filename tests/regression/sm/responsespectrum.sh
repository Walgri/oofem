#
# Response spectrum analysis (ResponseSpectrum) of a single degree of freedom oscillator:
# cantilever column, L = 3, EI = 3e7 * 6.75e-4, lumped mass m = 10 on top, constant spectrum Sa = 2.
# Exact values: base shear m Sa = 20, top displacement Sa / omega^2 = 2 / 225 (plus shear deformation),
# base moment m Sa L = 60. With P-Delta (axial load 1000 on top) the base moment grows to 60 + P d.
# The combined results are printed as RSA lines in the .out file, so they are checked here with awk.
#
OOFEM=$1
shift
failed=0

# value of field n on the RSA line "<kind> [<id>]" in direction <dir>
rsa() {     # file dir kind id n
    awk -v dir="$2" -v kind="$3" -v id="$4" -v n="$5" '
        $1 == "RSA" && $2 == "direction" { d = $3; next }
        d == dir && $1 == "RSA" && $2 == kind && (id == "" || $3 == id) { print $n; exit }' "$1"
}
check() {   # description value expected relative_tolerance
    r=$(awk -v a="$2" -v b="$3" -v t="$4" 'BEGIN { if (a == "") { print "missing"; exit }
        e = (a - b) / (b == 0 ? 1 : b); if (e < 0) e = -e; print (e <= t ? "ok" : "FAILED") }')
    echo "$1: $2 (expected $3) $r"
    [ "$r" = ok ] || failed=1
}

$OOFEM "$@" -f responsespectrum01.in.0 || exit 1
check "base shear" "$(rsa responsespectrum01.out X base_shear '' 3)" 20.0 1e-6
check "top displacement" "$(rsa responsespectrum01.out X node 11 4)" 8.8888936889e-03 1e-6
check "base moment" "$(rsa responsespectrum01.out X forces 1 10)" -60.0 1e-6
check "base moment, direction Y" "$(rsa responsespectrum01.out Y forces 1 9)" 60.0 1e-6

$OOFEM "$@" -f responsespectrum02.in.0 || exit 1
check "base moment with P-Delta" "$(rsa responsespectrum02.out X forces 1 10)" -70.841192237 1e-6

exit $failed
