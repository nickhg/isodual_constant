#!/usr/bin/env python3
"""Check, in exact rational arithmetic, that the union of the certify roots (as printed by cover_check.py) contains
the box B of Appendix C for the threshold's eta cut c.

B, for 4-dimensional M with det 1, Phi_5(M) >= PHI, and a Psi-minimal, size-reduced, sign-normalised basis:
  y_k <= ybar_k, where ybar_2 = log (4/3)^2 (LLL bound) and ybar_1 = ybar_3 = log mu_4 = (1/2) log 2 (window bound,
          sharper than the LLL bound (3/2) log(4/3));
  y_k >= log c - sum_{j != k} ybar_j, with c = 9 PHI / 32 (kappa^2 <= 32/9, so Phi >= PHI forces eta >= 9 PHI / 32):
          c = 63/32 for Phi_5 < 7 (App. C), 2048/1125 for Phi_5 <= P (App. D), 27/16 for stage 1 of eta_5 = 6 (App. E);
  mu_{i,i-1} in [0, 1/2], the other mu_ij in [-1/2, 1/2].
Every comparison is exact. Each endpoint u of the union is a double, hence a rational; the logarithmic endpoints of B
are compared through exp: u <= (1/2) log Q  iff  exp(2u) <= Q, with Q rational, and exp(2u) is enclosed between
rational bounds (Taylor series with an explicit remainder, after halving the argument). A comparison that the
enclosure cannot decide counts as a failure.
Usage: cover_check.py logs... | root_box_check.py [c]        (c a rational such as 63/32; default 63/32)
"""
import sys, re
from fractions import Fraction as F


def exp_bounds(x, terms=40):
    """Rational lo <= exp(x) <= hi for rational x."""
    k = 0
    while abs(x) > F(1, 4):
        x /= 2; k += 1
    s, t = F(0), F(1)
    for n in range(terms):
        s += t; t = t * x / (n + 1)
    rem = 2 * abs(t)              # |x| <= 1/4: the tail after `terms` terms is at most twice its first term
    lo, hi = s - rem, s + rem
    for _ in range(k):            # exp(x 2^k) = exp(x)^(2^k); squaring is monotone on positive bounds
        lo, hi = lo * lo, hi * hi
    assert 0 < lo <= hi
    return lo, hi


def le_half_log(u, Q):          # u <= (1/2) log Q ?   True only if proved
    return exp_bounds(2 * F(u))[1] <= Q


def ge_half_log(u, Q):          # u >= (1/2) log Q ?
    return exp_bounds(2 * F(u))[0] >= Q


c = F(sys.argv[1]) if len(sys.argv) > 1 else F(63, 32)
# (1/2) log of: the lower ends exp(2 blo_k) and the upper ends exp(2 ybar_k)
Qlo = [81 * c * c / 512, c * c / 4, 81 * c * c / 512]      # k = 1, 3: c (9/16) / sqrt 2;  k = 2: c / 2
Qhi = [F(2), F(256, 81), F(2)]                            # k = 1, 3: sqrt 2;  k = 2: 16/9
MU = [(F(0), F(1, 2)) if sub else (F(-1, 2), F(1, 2)) for sub in (1, 0, 1, 0, 0, 1)]   # mu10 mu20 mu21 mu30 mu31 mu32

text = sys.stdin.read()
m = re.search(r"tile y1 in \[(\S+), (\S+)\]", text)
U = [(float.fromhex(m.group(1)), float.fromhex(m.group(2)))]
U += [(float.fromhex(a), float.fromhex(b)) for a, b in re.findall(r"coordinate \d+: \[(\S+), (\S+)\]", text)]
assert len(U) == 9, "expected 9 coordinates"
ok = True
for t, (ulo, uhi) in enumerate(U):
    if t < 3:
        good = le_half_log(ulo, Qlo[t]) and ge_half_log(uhi, Qhi[t])
        desc = "B = [(1/2) log %s, (1/2) log %s]" % (Qlo[t], Qhi[t])
    else:
        blo, bhi = MU[t - 3]
        good = F(ulo) <= blo and F(uhi) >= bhi
        desc = "B = [%s, %s]" % (blo, bhi)
    ok &= good
    print("coordinate %d: union [%.12f, %.12f]  %s  %s" % (t + 1, ulo, uhi, desc, "ok" if good else "NOT CONTAINED"))
print("c = %s: union contains B" % c if ok else "FAILED: union does not contain B (c = %s)" % c)
sys.exit(0 if ok else 1)
