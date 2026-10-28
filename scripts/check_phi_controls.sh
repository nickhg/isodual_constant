#!/bin/bash
# Negative controls for check_phi at the two points where the bound is attained: M* (Phi_5 = 2^16/(3^4 5^3) = P) and
# pi_v D5 (Phi_5 = 6).  For each, a one-leaf log is built whose root box has half-width 1e-7 about a Psi-minimal
# basis of the lattice (the first basis in inputs/mstar_U.txt, resp. inputs/pivd5_U.txt), and checked with the
# threshold equal to the maximum, 10^-6 below it, and 10^-3 above it:
#   - without the balls (MSTAR/PVD5 unset), the leaf must be REJECTED at the maximum and below it, although the checker
#     may bisect it 24 more times: no certificate of Phi < threshold exists on a box containing the maximizer;
#   - 10^-3 above the maximum, without the balls, it must be ACCEPTED (by the kappa/eta bound), so the rejection
#     comes from the value of Phi there, not from the box;
#   - with the balls, it must be accepted at the maximum and below it (Proposition [M*], resp. [pi_v D5], proves
#     Phi <= P, resp. <= 6, there), and the checker must refuse a threshold above the maximum.
# Usage: scripts/check_phi_controls.sh CHECK_PHI        (about a second; needs python3)
set -u
CP=$1; R=$(cd "$(dirname "$0")/.." && pwd); T=$(mktemp -d); trap 'rm -rf "$T"' EXIT
python3 - "$R/inputs" "$T" <<'EOF'
import sys, math
from fractions import Fraction as Fr
inp, out = sys.argv[1:]
fams = {'mstar': ([[18, -9, 0, 0], [-9, 18, -9, -9], [0, -9, 16, 1], [0, -9, 1, 16]], 'mstar_U.txt', Fr(65536, 10125)),
        'pvd5': ([[4, -2, 0, 0], [-2, 3, -2, 0], [0, -2, 4, -2], [0, 0, -2, 4]], 'pivd5_U.txt', Fr(6))}
W = 1e-7
for name, (P, fn, phi) in fams.items():
    u = [int(x) for x in open('%s/%s' % (inp, fn)).readline().split()]
    U = [u[4 * i:4 * i + 4] for i in range(4)]
    Q = [[sum(U[i][a] * P[a][b] * U[j][b] for a in range(4) for b in range(4)) for j in range(4)] for i in range(4)]
    mu = [[Fr(0)] * 4 for _ in range(4)]; r = [Fr(0)] * 4          # Gram-Schmidt, exactly
    for i in range(4):
        for j in range(i): mu[i][j] = (Q[i][j] - sum(mu[i][k] * mu[j][k] * r[k] for k in range(j))) / r[j]
        r[i] = Q[i][i] - sum(mu[i][k] ** 2 * r[k] for k in range(i))
    det = r[0] * r[1] * r[2] * r[3]
    y = [math.log(r[0]) - math.log(det) / 4, math.log(r[0] * r[1]) - math.log(det) / 2,
         math.log(r[0] * r[1] * r[2]) - 3 * math.log(det) / 4]       # det-1 normalisation, y_k = log G_k
    c = y + [float(mu[i][j]) for i, j in ((1, 0), (2, 0), (2, 1), (3, 0), (3, 1), (3, 2))]
    for tag, f in (('eq', Fr(1)), ('low', 1 - Fr(1, 10 ** 6)), ('high', 1 + Fr(1, 10 ** 3))):
        with open('%s/%s_%s.log' % (out, name, tag), 'w') as o:
            o.write('ROOT n=4 dim=9 phi=%s c=1.8204 ver=negative-control\n' % float(phi * f).hex())
            for v in c: o.write('%s %s\n' % ((v - W).hex(), (v + W).hex()))
            o.write('L 0\nEND complete\n')
EOF
pass=0; fail=0
expect() {    # expect NAME accept|reject|refuse ENV... : accept = exit 0; reject = FAILURES 1; refuse = header refused
    local name=$1 want=$2; shift 2
    out=$(env EXTRA=24 "$@" 2>&1); rc=$?
    case $want in
        accept) ok=$([ $rc = 0 ] && echo "$out" | grep -q 'FAILURES 0' && echo y) ;;
        reject) ok=$([ $rc = 4 ] && echo "$out" | grep -q 'FAILURES 1' && echo y) ;;
        refuse) ok=$([ $rc = 1 ] && echo "$out" | grep -q 'needs PHI' && echo y) ;;
    esac
    if [ "$ok" = y ]; then pass=$((pass + 1)); echo "ok    $name ($want)"
    else fail=$((fail + 1)); echo "WRONG $name (expected $want): $(echo "$out" | tail -1 | cut -c1-90)"; fi
}
for f in mstar:MSTAR:mstar_U.txt pvd5:PVD5:pivd5_U.txt; do
    IFS=: read -r name var balls <<< "$f"
    expect "$name, threshold = max, no balls"        reject "$CP" $T/${name}_eq.log
    expect "$name, threshold = max(1-1e-6), no balls" reject "$CP" $T/${name}_low.log
    expect "$name, threshold = max(1+1e-3), no balls" accept "$CP" $T/${name}_high.log
    expect "$name, threshold = max, balls"           accept "$var=$R/inputs/$balls" "$CP" $T/${name}_eq.log
    expect "$name, threshold = max(1-1e-6), balls"   accept "$var=$R/inputs/$balls" "$CP" $T/${name}_low.log
    expect "$name, threshold = max(1+1e-3), balls"   refuse "$var=$R/inputs/$balls" "$CP" $T/${name}_high.log
done
echo "passed $pass, wrong $fail"
[ $fail = 0 ]
