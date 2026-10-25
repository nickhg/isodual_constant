#!/usr/bin/env python3
"""Check that the union of the certify roots (as printed by cover_check.py) contains the box B of the appendix.

B, for 4-dimensional M with det 1, Phi_5(M) >= 7, and a Psi-minimal, size-reduced, sign-normalised basis:
  y_k <= ybar_k, where ybar_2 = log (4/3)^2 (LLL bound) and ybar_1 = ybar_3 = log mu_4 = (1/2) log 2 (window bound,
          sharper than the LLL bound (3/2) log(4/3));
  y_k >= log c - sum_{j != k} ybar_j, with c = 63/32 (kappa^2 <= 32/9, so Phi >= 7 forces eta >= 63/32);
  mu_{i,i-1} in [0, 1/2], the other mu_ij in [-1/2, 1/2].
The endpoints are computed to 50 digits; the union must contain B, compared exactly with the hex floats.
Usage: cover_check.py logs... | root_box_check.py
"""
import sys, re
from fractions import Fraction
from mpmath import mp, mpf, log

mp.dps = 50
ybar = [log(2) / 2, 2 * log(mpf(4) / 3), log(2) / 2]
logc = log(mpf(63) / 32)
B = [(logc - (sum(ybar) - ybar[k]), ybar[k]) for k in range(3)]
B += [(mpf(0), mpf(1) / 2) if sub else (mpf(-1) / 2, mpf(1) / 2) for sub in (1, 0, 1, 0, 0, 1)]   # mu10 mu20 mu21 mu30 mu31 mu32

text = sys.stdin.read()
m = re.search(r"tile y1 in \[(\S+), (\S+)\]", text)
U = [(float.fromhex(m.group(1)), float.fromhex(m.group(2)))]
U += [(float.fromhex(a), float.fromhex(b)) for a, b in re.findall(r"coordinate \d+: \[(\S+), (\S+)\]", text)]
assert len(U) == 9, "expected 9 coordinates"
ok = True
for t, ((ulo, uhi), (blo, bhi)) in enumerate(zip(U, B)):
    good = mpf(Fraction(ulo).numerator) / Fraction(ulo).denominator <= blo and \
           mpf(Fraction(uhi).numerator) / Fraction(uhi).denominator >= bhi
    ok &= good
    print("coordinate %d: union [%.12f, %.12f]  B [%.12f, %.12f]  %s" % (t + 1, ulo, uhi, blo, bhi, "ok" if good else "NOT CONTAINED"))
print("union contains B" if ok else "FAILED: union does not contain B")
sys.exit(0 if ok else 1)
