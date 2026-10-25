#!/bin/bash
# Test of the vendored CORE-MATH exp/log (third_party/core-math): N random arguments over the checkers' ranges, results
# in round-upward and round-to-nearest compared with the correctly rounded values from mpmath at 60 digits.
# Usage: scripts/test_crexplog.sh [N]   (2026-10-24: 200,000 samples, 0 mismatches in all four cases)
set -eu
N=${1:-200000}; T=$(mktemp -d); trap 'rm -rf "$T"' EXIT; C=src/third_party/core-math
gcc -O2 -mfma -ffp-contract=off -c $C/exp.c -o $T/exp.o
gcc -O2 -mfma -ffp-contract=off -I$C -c $C/log.c -o $T/log.o
gcc -O2 -frounding-math -ffp-contract=off scripts/crexplog/gen.c $T/exp.o $T/log.o -lm -o $T/gen
$T/gen "$N" > $T/out.txt
python3 - "$T/out.txt" <<'PY'
import sys, math, mpmath as mp
mp.mp.dps = 60
def up(x):
    d = float(x); return d if mp.mpf(d) >= x else math.nextafter(d, math.inf)
bad = dict(exp_up=0, exp_rn=0, log_up=0, log_rn=0); n = 0
for l in open(sys.argv[1]):
    xe, eu, en, xl, lu, ln = (float.fromhex(t) for t in l.split()); n += 1
    te, tl = mp.exp(mp.mpf(xe)), mp.log(mp.mpf(xl))
    bad["exp_up"] += eu != up(te); bad["exp_rn"] += en != float(te)
    bad["log_up"] += lu != up(tl); bad["log_rn"] += ln != float(tl)
print(n, "samples; mismatches against correct rounding:", bad)
sys.exit(1 if any(bad.values()) else 0)
PY
