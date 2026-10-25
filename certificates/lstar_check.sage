# L_* (Remark rem:phi5): build it from P_0, phi_0 = (0,0,1/3,1/3), r = 18; check WC's integral form, shells, the
# flag with prefix determinants 6, 27, 108, 405, 1350, and compute eta(L_*) exactly with psi.sage.
import os
load("psi.sage")
from sage.all import QQ, ZZ, matrix, vector, pari

P0 = matrix(QQ, [[18, -9, 0, 0], [-9, 18, -9, -9], [0, -9, 16, 1], [0, -9, 1, 16]])
phi = vector(QQ, [0, 0, QQ(1) / 3, QQ(1) / 3]); r = QQ(18)
G = matrix(QQ, 5, 5)
G[0, 0] = r
for i in range(4):
    G[0, i + 1] = G[i + 1, 0] = r * phi[i]
    for j in range(4): G[i + 1, j + 1] = P0[i, j] + r * phi[i] * phi[j]
print("min of P0:", pari(P0).qfminim(None, 0)[1], " r(1 - 1/9) =", r * (1 - QQ(1) / 9))
Gs = G / 3
print("G/3 integral:", all(x in ZZ for x in Gs.list()), " det", Gs.det(), "=", factor(Gs.det()))
W = matrix(ZZ, [[6, 0, 0, 2, 2], [0, 6, -3, 0, 0], [0, -3, 6, -3, -3], [2, 0, -3, 6, 1], [2, 0, -3, 1, 6]])
print("WC form det", W.det())
# equivalence of G/3 and WC's form: same genus is not enough; test isometry with pari qfauto/qfisom
iso = pari(matrix(ZZ, Gs)).qfisom(pari(W))
print("isometric to WC's form:", iso != 0)
v = pari(W).qfminim(14)
from collections import Counter
cnt = Counter()
for x in matrix(ZZ, v[2]).columns(): cnt[x * W * x] += 2
print("shells (norm: vectors):", sorted(cnt.items()))
F = matrix(ZZ, [[0, 0, 0, 0, 1], [0, 0, 1, 0, 0], [0, 1, 0, 0, 0], [0, 0, 0, 1, 0], [1, -1, -2, -1, -1]])
print("flag unimodular:", F.det(), " prefix dets:", [(F[:k] * W * F[:k].transpose()).det() for k in range(1, 6)])
pd = [(F[:k] * W * F[:k].transpose()).det() for k in range(1, 5)]
print("flag bound:", prod(pd) / W.det() ** 2)
e, fl = eta(W) if "NOETA" not in os.environ else (None, None)
print("eta(L_*) exact =", e, "=", float(e))
