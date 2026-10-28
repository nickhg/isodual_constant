# Independent check of "E6 is a strict local maximum of eta" (B. Moulton's e6_theorem9.py, 26 Oct), by a different
# route.  Flags: chains L1 < ... < L5 with det L_k = 2,3,4,4,4 (the Hermite bounds for min 2, gamma_1..gamma_5), so
# L_k is the critical lattice A1, A2, A3, D4, D5.  Completeness with ROOTS only: L_{k+1} = L_k + Z w with w reduced
# mod L_k, |w|^2 <= det L_{k+1}/det L_k + mu(L_k)^2 = 3/2+1/2, 4/3+2/3, 1+1, 1+1 = 2 (mu = covering radius), so w
# is a root.  Then: design identity for the first-order functionals (Sym^2_0 irreducible), which gives a radius.
import json, itertools, time
t0 = time.time()
n = 6
G0 = matrix(ZZ, [[2,0,-1,0,0,0],[0,2,0,-1,0,0],[-1,0,2,-1,0,0],[0,-1,-1,2,-1,0],[0,0,0,-1,2,-1],[0,0,0,0,-1,2]])
assert G0.det() == 3
# roots: all x in Z^6 with x G0 x^T = 2 (Fincke-Pohst via pari: qfminim returns one of each +-pair)
R = pari(G0).qfminim(2)
roots = matrix(ZZ, R[2]).columns()
roots = roots + [-r for r in roots]
assert len(roots) == 72 and all(r * G0 * r == 2 for r in roots)
TARGET = [2, 3, 4, 4, 4]
def key(B): return tuple(matrix(ZZ, B).echelon_form().list())
# DFS over chains; children of L_k are the distinct L_k + Z w (w root) with the target determinant
flags = []
memo = {}
def children(B):
    k0 = key(B)
    if k0 in memo: return memo[k0]
    out = {}
    k = B.nrows()
    for w in roots:
        C = matrix(ZZ, B.rows() + [w])
        if C.rank() == k + 1 and (C * G0 * C.T).det() == TARGET[k]:
            kk = key(C)
            if kk not in out: out[kk] = C
    memo[k0] = list(out.values()); return memo[k0]
lines = {}
for r in roots:
    kk = key([r])
    if kk not in lines: lines[kk] = matrix(ZZ, [r])
def dfs(chain):
    if len(chain) == 5:
        flags.append(chain); return
    for C in children(chain[-1]): dfs(chain + [C])
for B in lines.values(): dfs([B])
print("lines", len(lines), "flags", len(flags), "%.0fs" % (time.time() - t0))
psi2 = 2 * 3 * 4 * 4 * 4
eta = psi2 / sqrt(QQ(G0.det()))**(n - 1)   # eta = prod_{k<n} det L_k / det^{(n-1)/2} (det L_n factor absorbed)
print("eta(E6) =", eta, "=", eta.n(), " 128 sqrt3/9 =", (128*sqrt(3)/9).n())
# first-order functionals on traceless Y (tr(G0^-1 Y) = 0), Lam = G0^{-1/2} Y G0^{-1/2}
G0q = G0.change_ring(QQ); G0i = G0q.inverse()
Ys = []
for i in range(n):
    for j in range(i, n):
        E = matrix(QQ, n, n); E[i, j] = 1; E[j, i] = 1
        Ys.append(E - (G0i * E).trace() / n * G0q)
Ys = [Ys[p] for p in matrix(QQ, [Y.list() for Y in Ys]).pivot_rows()]
assert len(Ys) == 20
H = matrix(QQ, 20, 20, lambda a, b: (G0i * Ys[a] * G0i * Ys[b]).trace())   # <Lam_a, Lam_b>_F
def Qf(chain):
    W = matrix(QQ, n, n)
    for U in chain:
        U = U.change_ring(QQ); W += U.T * (U * G0q * U.T).inverse() * U
    return W
Qs = [Qf(c) for c in flags]
assert all((G0q * Q).trace() == 15 for Q in Qs)
assert all(((G0q * Q)**2).trace() == 55 for Q in Qs)          # |C_f|_F^2 = 55 - 15^2/6 = 35/2
A = matrix(QQ, [[(Q * Y).trace() for Y in Ys] for Q in Qs])   # N x 20, row f = L_f in Y-coords
N = A.nrows()
print("rank", A.rank(), "distinct", len(set(tuple(r) for r in A.rows())))
s1 = sum(A.rows())
print("sum of L_f = 0:", s1 == 0)
# design: sum_f L_f(Lam)^2 = c |Lam|^2.  In Y-coords L_f(Lam) = a_f . y and |Lam|^2 = y^T H y, so the dual
# statement is A^T A = c H, with c = N (35/2) / 20
S = A.T * A
c = QQ(N) * QQ(35) / 2 / 20
print("design A^T A == c H:", S == c * H, " c =", c)
Hi = H.inverse()
print("dual norms all 35/2:", all(a * Hi * a == QQ(35)/2 for a in A.rows()))
json.dump({"flags": [[[int(x) for x in row] for row in C.rows()] for C in [ch[-1] for ch in flags]],
           "chains": [[list(map(int, ch[k].row(k))) for k in range(5)] for ch in flags]},
          open("e6_flags.json", "w"))
save((A, H, Ys), "e6_A.sobj")
print("done %.0fs" % (time.time() - t0))
