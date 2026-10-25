#!/bin/bash
# Broken-certificate tests for check_lift: every mutation of a valid stage-2 log must be rejected (FAILURES > 0 or a
# parse error), and the unmodified log must pass.  Usage: check_lift_tests.sh CHECK_LIFT LOG   (LOG with N, F, P, B leaves)
set -u
CL=$1; LOG=$2; T=$(mktemp -d); trap 'rm -rf "$T"' EXIT
pass=0; fail=0
expect() {    # expect NAME accept|reject FILE
    out=$("$CL" "$3" 2>&1); rc=$?
    if [ "$2" = accept ]; then ok=$([ $rc = 0 ] && echo y); else ok=$([ $rc != 0 ] && echo y); fi
    if [ "$ok" = y ]; then pass=$((pass + 1)); echo "ok    $1 ($2): $(echo "$out" | tail -1 | cut -c1-90)"
    else fail=$((fail + 1)); echo "WRONG $1 (expected $2): $(echo "$out" | tail -1 | cut -c1-90)"; fi
}
nth() { grep -n "$1" "$LOG" | sed -n "${2}p" | cut -d: -f1; }      # line number of the k-th match
expect original accept "$LOG"
# tree structure
l=$(nth '^S ' 5); sed "${l}d" "$LOG" > $T/a; expect "split line deleted" reject $T/a
l=$(nth '^S ' 5); sed "${l}s/^S [0-9]*/S 13/" "$LOG" > $T/a; expect "split coordinate changed" reject $T/a
grep -v '^END' "$LOG" > $T/a; expect "END missing" reject $T/a
sed 's/^END complete/END incomplete/' "$LOG" > $T/a; expect "END incomplete" reject $T/a
{ cat "$LOG"; echo "F 1 0 0 0 0 0 1 0 0 0 0 0 1 0 0 0 0 0 1 0 0 0 0 0 1"; } > $T/a; expect "trailing leaf" reject $T/a
# root fibre
rr() { python3 -c "import sys; L = open(sys.argv[1]).read().split('\n'); lo, hi = (float.fromhex(x) for x in L[1].split()); i = int(sys.argv[2]); d = float(sys.argv[3])
lo, hi = (lo + d, hi) if i == 0 else (lo, hi + d); L[1] = '%s %s' % (lo.hex(), hi.hex()); open(sys.argv[4], 'w').write('\n'.join(L))" "$LOG" "$@"; }
rr 0 1e-3 $T/a; expect "r-range starts above 2^(2/5)" reject $T/a
rr 1 -1e-3 $T/a; expect "r-range ends below 8^(1/5)" reject $T/a
rr 0 2e-16 $T/a; expect "r-range starts 1 ulp below 2^(2/5) (still covers)" accept $T/a
rr 0 8.9e-16 $T/a; expect "r-range starts 2 ulps above 2^(2/5)" reject $T/a
awk 'NR == 3 {print "0x0p+0 0x1.f5c28f5c28f5cp-2"; next} {print}' "$LOG" > $T/a; expect "phi_1 range [0, 0.49]" reject $T/a
awk 'NR == 5 {print "-0x1p-1 0x1.f5c28f5c28f5cp-2"; next} {print}' "$LOG" > $T/a; expect "phi_3 range [-0.5, 0.49]" reject $T/a
# N leaves
l=$(nth '^N ' 3); awk -v l=$l 'NR == l {$2 = $2 + 1} {print}' "$LOG" > $T/a; expect "N: k off by one" reject $T/a
l=$(nth '^N ' 3); sed "${l}s/^N \(-*[0-9]*\) .*/N \1 0 0 0 0/" "$LOG" > $T/a; expect "N: c = 0" reject $T/a
# F leaves
l=$(nth '^F ' 2); sed "${l}s/.*/F 2 0 0 0 0 0 1 0 0 0 0 0 1 0 0 0 0 0 1 0 0 0 0 0 1/" "$LOG" > $T/a; expect "F: not unimodular" reject $T/a
sed 's/^F .*/F 1 0 0 0 0 0 1 0 0 0 0 0 1 0 0 0 0 0 1 0 0 0 0 0 1/' "$LOG" > $T/a; expect "F: every flag the identity" reject $T/a
l=$(nth '^F ' 2); sed "${l}s/^F .*/N 0 1 0 0 0/" "$LOG" > $T/a; expect "F leaf replaced by an N leaf" reject $T/a
# P leaves
l=$(nth '^P ' 2); awk -v l=$l 'NR == l + 1 {$2 = "-" $2} {print}' "$LOG" > $T/a; expect "P: negative weight" reject $T/a
l=$(nth '^P ' 2); awk -v l=$l 'NR > l && NR <= l + 9 && /^[0-4] / {$2 = "0x0p+0"} {print}' "$LOG" > $T/a
expect "P: all weights zero" reject $T/a
awk '/^P / {p = $2; print; next} p > 0 {p--; print $1, $2, "1 0 0 0 0 0 1 0 0 0 0 0 1 0 0 0 0 0 1 0 0 0 0 0 1"; next} {print}' "$LOG" > $T/a; expect "P: every flag the identity" reject $T/a
# B leaves
if grep -q '^B ' "$LOG"; then
    l=$(nth '^B ' 1); sed "${l}s/^B \([01]\) .*/B \1 1 0 0 0 0 0 1 0 0 0 0 0 1 0 0 0 0 0 1 0 0 0 0 0 1/" "$LOG" > $T/a; expect "B: W the identity" reject $T/a
    l=$(nth '^B ' 1); awk -v l=$l 'NR == l {$2 = 1 - $2} {print}' "$LOG" > $T/a; expect "B: wrong target (D5 <-> D5*)" reject $T/a
    l=$(nth '^B ' 1); awk -v l=$l 'NR == l {$3 = 2 * $3 + ($3 == 0)} {print}' "$LOG" > $T/a; expect "B: W not unimodular" reject $T/a
fi
echo "passed $pass, wrong $fail"
[ $fail = 0 ]
