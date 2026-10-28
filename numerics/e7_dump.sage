# Dump the C_f vectors of all 544,320 Psi-minimal flags of E7 (chains A1<A2<A3<D4<D5<E6 from roots, as in
# ../certificates/e7_design.sage) in orthonormal coordinates for the trace form on traceless Lambda, deduplicated.
# Also the fundamental-weight directions I - (n/k) P_V for all 127 subsets V, in the same coordinates.
import time, itertools, numpy as np; t0 = time.time()
n = 7
G0 = matrix(ZZ, CartanMatrix(['E', 7])); G0q = G0.change_ring(QQ); G0i = G0q.inverse()
R = pari(G0).qfminim(2); roots = matrix(ZZ, R[2]).columns(); roots += [-r for r in roots]
TARGET = [2, 3, 4, 4, 4, 3]
Ys = []
for i in range(n):
    for j in range(i, n):
        E = matrix(QQ, n, n); E[i, j] = 1; E[j, i] = 1
        Ys.append(E - (G0i * E).trace() / n * G0q)
Ys = [Ys[p] for p in matrix(QQ, [Y.list() for Y in Ys]).pivot_rows()]
D = len(Ys); assert D == 27
H = matrix(QQ, D, D, lambda a, b: (G0i * Ys[a] * G0i * Ys[b]).trace())
def key(B): return tuple(B.echelon_form().list())
def wvec(B):
    U = B.change_ring(QQ); W = U.T * (U * G0q * U.T).inverse() * U
    return vector(QQ, [(W * Y).trace() for Y in Ys])
kids = {}; W = {}
def children(B):
    k0 = key(B)
    if k0 in kids: return kids[k0]
    out = {}; k = B.nrows()
    for w in roots:
        C = matrix(ZZ, B.rows() + [w])
        if C.rank() == k + 1 and (C * G0 * C.T).det() == TARGET[k]:
            out.setdefault(key(C), C)
    kids[k0] = list(out.items()); return kids[k0]
# Cholesky of H: coordinates z with <a, y> = (L^{-1} a) . (L^T y)
Hn = np.array(H.change_ring(RDF)); Lc = np.linalg.cholesky(Hn); Li = np.linalg.inv(Lc)
def wf(kk, B):
    if kk not in W: W[kk] = Li @ np.array(wvec(B).change_ring(RDF))
    return W[kk]
lines = {}
for r in roots:
    B = matrix(ZZ, [r]); lines.setdefault(key(B), B)
out = []
def dfs(kk, B, acc):
    a = acc + wf(kk, B)
    if B.nrows() == 6: out.append(a); return
    for k2, C in children(B): dfs(k2, C, a)
for kk, B in lines.items(): dfs(kk, B, np.zeros(D))
F = np.array(out); print("flags", len(F), "%.0fs" % (time.time() - t0))
Fu = np.unique(np.round(F, 10), axis=0); print("distinct C_f", len(Fu))
print("norms^2", np.unique(np.round((Fu * Fu).sum(1), 8)))
# fundamental weights: columns of G0^{-1} in root coords; Lambda = I - (n/k) P_V, as a form: G0 - (n/k) G0 P_V (traceless wrt G0)
dirs = []
for k in range(1, n + 1):
    for S in itertools.combinations(range(n), k):
        if k == n: continue
        Wt = G0i.matrix_from_columns(list(S))                 # weights in root coordinates
        PV = Wt * (Wt.T * G0q * Wt).inverse() * Wt.T * G0q    # G0-orthogonal projection, root coords
        Lam = G0q * (identity_matrix(QQ, n) - QQ(n) / k * PV) # symmetric form
        y = matrix(QQ, [Y.list() for Y in Ys]).T.solve_right(vector(QQ, Lam.list()))
        z = Lc.T @ np.array(y.change_ring(RDF)); dirs.append((S, z / np.linalg.norm(z)))
np.savez("e7_flags.npz", F=Fu, dirs=np.array([d for _, d in dirs]), subsets=np.array([str(s) for s, _ in dirs]))
print("saved %.0fs" % (time.time() - t0))
