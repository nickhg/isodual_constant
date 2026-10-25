"""Exact check that D4 and A2 (x) A2 lie in the isodual family F through M_* = {G : U (c G^{-1}) U^T = G} (Prop.
mstarprops (4)): from the numerical minimisers of family.py, recover an integral basis change B with B G B^T ~ R, and
check in exact rational arithmetic that H = B^{-1} R B^{-T} satisfies U (c' H^{-1}) U^T = H for the scale c' = det(H)^{1/2}
(the equation is homogeneous: G = lambda H is in F iff c' = 135 / lambda^2)."""
import numpy as np, itertools, os, sys
from fractions import Fraction as Fr
os.environ["CONFIRM_ONLY"] = "1"
sys.argv = ["x"]
exec(open(os.path.join(os.path.dirname(os.path.abspath(__file__)), "family.py")).read().split("# M_* on the geodesic")[0])
def best_basis(G, R, tol=0.35):
    G = n1(G); Rn = n1(R); rmax = max(np.diag(Rn)) * (1 + tol); Gi = np.linalg.inv(G)
    rad = [int(np.sqrt(rmax * Gi[i, i])) for i in range(4)]
    ks = np.array([k for k in itertools.product(*[range(-r, r + 1) for r in rad]) if any(k)])
    nv = np.einsum('ij,jk,ik->i', ks, G, ks); m = nv <= rmax; ks, nv = ks[m], nv[m]; ip = ks @ G @ ks.T
    cand = [np.where(np.abs(nv - Rn[i, i]) <= tol * Rn[i, i])[0] for i in range(4)]; best = [9.0, None]
    def rec(i, ch):
        if i == 4:
            B = ks[ch]
            if abs(round(np.linalg.det(B))) == 1:
                d = logdist(Rn, B @ G @ B.T)
                if d < best[0]: best[0], best[1] = d, B.copy()
            return
        for c in cand[i]:
            if all(abs(ip[c, ch[j]] - Rn[i, j]) <= tol * np.sqrt(Rn[i, i] * Rn[j, j]) for j in range(i)): rec(i + 1, ch + [c])
    rec(0, []); return best
def frmat(A): return [[Fr(int(round(x))) for x in r] for r in A]
def mul(A, B): return [[sum(A[i][k] * B[k][j] for k in range(len(B))) for j in range(len(B[0]))] for i in range(len(A))]
def tr(A): return [list(r) for r in zip(*A)]
def inv(A):
    n = len(A); M = [list(A[i]) + [Fr(int(i == j)) for j in range(n)] for i in range(n)]
    for c in range(n):
        p = next(r for r in range(c, n) if M[r][c] != 0); M[c], M[p] = M[p], M[c]
        d = M[c][c]; M[c] = [x / d for x in M[c]]
        for r in range(n):
            if r != c and M[r][c] != 0: f = M[r][c]; M[r] = [a - f * b for a, b in zip(M[r], M[c])]
    return [r[n:] for r in M]
def det(A):
    import sympy; return Fr(str(sympy.Matrix([[sympy.Rational(x.numerator, x.denominator) for x in r] for r in A]).det()))
Uf = frmat(U)
for nm, R, x in (("D4", D4, [-0.06686362, -0.06537841, 0.12987742]), ("A2xA2", K, [0.18733742, 0.18317619, -0.36388849])):
    d, B = best_basis(fam(x), R)
    Bf = frmat(B); Bi = inv(Bf); H = mul(mul(Bi, frmat(R)), tr(Bi))            # the copy of R in M_*'s coordinates
    dh = det(H); import math
    c2 = dh                                                                  # c'^2 = det H (n = 4)
    from sympy import sqrt, Rational, nsimplify
    cp = Fr(str(nsimplify(sqrt(Rational(dh.numerator, dh.denominator)))))
    lhs = mul(mul(Uf, [[cp * v for v in r] for r in inv(H)]), tr(Uf))
    print("%s: numerical distance %.2e; integral B found; H rational, det %s, c' = %s; U (c' H^-1) U^T == H exactly: %s"
          % (nm, d, dh, cp, lhs == H))

# Is M_* exactly on the geodesic from the D4 copy H_D to the A2 x A2 copy H_A?  The geodesic is
# H_D^{1/2} (H_D^{-1/2} H_A H_D^{-1/2})^t H_D^{1/2} (after scaling); so M_* (P0) is on it iff X = H_D^{-1} P0 and
# Y = H_D^{-1} H_A commute and, simultaneously diagonalised, the normalised eigenvalues satisfy mu_i = lambda_i^t.
import sympy as sp
Hs = {}
for nm, R, x in (("D4", D4, [-0.06686362, -0.06537841, 0.12987742]), ("A2xA2", K, [0.18733742, 0.18317619, -0.36388849])):
    d, B = best_basis(fam(x), R); Bi = inv(frmat(B)); Hs[nm] = sp.Matrix(mul(mul(Bi, frmat(R)), tr(Bi))).applyfunc(lambda v: sp.Rational(v.numerator, v.denominator))
HD, HA, PP = Hs["D4"], Hs["A2xA2"], sp.Matrix(P0.astype(int).tolist())
X = HD.inv() * PP; Y = HD.inv() * HA
print("X Y == Y X exactly:", (X * Y - Y * X).is_zero_matrix)
print("H_D =", HD.tolist()); print("H_A =", HA.tolist())
