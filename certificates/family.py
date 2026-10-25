"""Is A2 (x) A2 on the isodual family through M_*?  The family F = {G : U (c G^{-1}) U^T = G} (U of Prop. mstarprops (4),
c = 135 at det 135^2) is the fixed set of an isometry of the invariant metric, hence totally geodesic:
F = {P0^{1/2} exp(L) P0^{1/2} : L in T}, T = the +1 eigenspace (dimension 3) of L -> -P0^{1/2} U P0^{-1} ... (computed below
numerically from D tau).  Minimise the invariant distance from F to the copies of A2 (x) A2 (det-normalised) by
multi-start Nelder-Mead; also report the distance from F to D4 for comparison.  Floating point."""
import numpy as np, itertools, os
from scipy.optimize import minimize
P0 = np.array([[18,-9,0,0],[-9,18,-9,-9],[0,-9,16,1],[0,-9,1,16]], float)
U = np.array([[0,0,-1,0],[0,0,0,1],[-1,1,0,-1],[0,-1,1,0]], float)
assert np.allclose(U @ (135 * np.linalg.inv(P0)) @ U.T, P0)
def msqrt(A): w, V = np.linalg.eigh(A); return V @ np.diag(np.sqrt(w)) @ V.T
def mexp(A): w, V = np.linalg.eigh(A); return V @ np.diag(np.exp(w)) @ V.T
S = msqrt(P0); Si = np.linalg.inv(S)
# exactly: with Q = sqrt(135) S^{-1} U^T S^{-1} (orthogonal, since Q^T Q = S^{-1} tau(P0) S^{-1} = I),
# tau(S e^L S) = S e^{-Q^T L Q} S, so F = S exp(T) S with T = {L symmetric : Q^T L Q = -L} (real null space, by SVD)
def tau(G): return U @ (135 * np.linalg.inv(G)) @ U.T
Q = np.sqrt(135) * Si @ U.T @ Si; assert np.allclose(Q.T @ Q, np.eye(4))
sym = []
for i in range(4):
    for j in range(i, 4):
        E = np.zeros((4, 4)); E[i, j] = E[j, i] = 1; sym.append(E)
Mlin = np.array([(Q.T @ E @ Q + E).ravel() for E in sym]).T
_, sv, Vt = np.linalg.svd(Mlin); null = Vt[np.sum(sv > 1e-9):]
T = [sum(c * E for c, E in zip(v, sym)) for v in null]; print("tangent dimension of F:", len(T))
def fam(x): return S @ mexp(sum(a * B for a, B in zip(x, T))) @ S
def n1(G): return G / np.linalg.det(G) ** 0.25
def logdist(R, G):
    w = np.linalg.eigvals(np.linalg.solve(R, G)).real; return float(np.sqrt(np.sum(np.log(w) ** 2)))
def dist(G, R, tol=0.35):
    G = n1(G); R = n1(R); rmax = max(np.diag(R)) * (1 + tol); Gi = np.linalg.inv(G)
    rad = [int(np.sqrt(rmax * Gi[i, i])) for i in range(4)]
    ks = np.array([k for k in itertools.product(*[range(-r, r + 1) for r in rad]) if any(k)])
    nv = np.einsum('ij,jk,ik->i', ks, G, ks); m = nv <= rmax; ks, nv = ks[m], nv[m]; ip = ks @ G @ ks.T
    cand = [np.where(np.abs(nv - R[i, i]) <= tol * R[i, i])[0] for i in range(4)]; best = [9.0]
    def rec(i, ch):
        if i == 4:
            B = ks[ch]
            if abs(round(np.linalg.det(B))) == 1: best[0] = min(best[0], logdist(R, B @ G @ B.T))
            return
        for c in cand[i]:
            if all(abs(ip[c, ch[j]] - R[i, j]) <= tol * np.sqrt(R[i, i] * R[j, j]) for j in range(i)): rec(i + 1, ch + [c])
    rec(0, []); return best[0]
A2 = np.array([[2,-1],[-1,2]], float); K = np.kron(A2, A2)
D4 = np.array([[2,-1,0,0],[-1,2,-1,-1],[0,-1,2,0],[0,-1,0,2]], float)
assert np.allclose(n1(tau(fam([0.3, -0.2, 0.1]))), n1(fam([0.3, -0.2, 0.1]))), "family not invariant"
print("check: tau fixes the family (3 random points):", all(np.allclose(n1(tau(fam(x))), n1(fam(x))) for x in np.random.default_rng(0).normal(size=(3, 3))))
for nm, R in (() if os.environ.get("CONFIRM_ONLY") else (("A2xA2", K), ("D4", D4))):
    best = (9, None); rng = np.random.default_rng(1)
    for s in range(40):
        r = minimize(lambda x: dist(fam(x), R), rng.normal(scale=0.6, size=3), method="Nelder-Mead", options=dict(xatol=1e-9, fatol=1e-12, maxiter=3000))
        if r.fun < best[0]: best = (r.fun, r.x)
    print("min distance from the isodual family through M* to %s: %.3e at x = %s (|x| = %.4f = distance from M* along F)" % (nm, best[0], np.round(best[1], 5), np.linalg.norm(best[1])))

# M_* on the geodesic: fit the two minimisers on one line through 0 (M_*), and evaluate along it
from scipy.optimize import minimize as _min
xd = minimize(lambda x: dist(fam(x), D4), np.array([-0.06686, -0.06538, 0.12988]), method="Nelder-Mead", options=dict(xatol=1e-12, fatol=1e-14, maxiter=4000)).x
xa = minimize(lambda x: dist(fam(x), K), np.array([0.18734, 0.18318, -0.36389]), method="Nelder-Mead", options=dict(xatol=1e-12, fatol=1e-14, maxiter=4000)).x
print("x(D4) =", xd, " x(A2xA2) =", xa, " ratio", xa / xd)
u = xd / np.linalg.norm(xd)
print("collinearity: |x(A2xA2) - (x(A2xA2).u) u| = %.2e" % np.linalg.norm(xa - (xa @ u) * u))
print("along the geodesic t -> fam(t u): t(D4) = %.6f, t(M*) = 0, t(A2xA2) = %.6f" % (xd @ u, xa @ u))
