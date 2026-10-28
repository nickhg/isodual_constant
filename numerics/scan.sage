# Phi_n along the geodesic G(t) = G1^{1/2} (G1^{-1/2} C G1^{-1/2})^t G1^{1/2} (det 1) between G1 and a copy C of the other
# peak (from nearest_local.py): kappa by branch and bound (kappa_bb, honest upper end max(u, lo (1 + tol))), eta by exact
# Psi-minimisation of the Gram matrix rationalised to 1e-13.  Phi_n = kappa^{(n-1)/2} eta, n - 1 = rank.  Numerical.
#   (run from numerics/) sage scan.sage COPIESFILE INDEX NPTS [T0 T1]   (points t = T0 .. T1, default 0 .. 1)
import sys, numpy as np
sys.path.insert(0, ".")
import kappa_bb as K
load("../certificates/psi.sage")
f, idx, npts = sys.argv[1], int(sys.argv[2]), int(sys.argv[3])
L = open(f).read().splitlines()
G1 = np.array(eval(L[0].split(":", 1)[1])); C = np.array(eval([l for l in L if l.startswith("copy")][idx].split(":", 1)[1]))
r = len(G1)
w, V = np.linalg.eigh(G1); S = V @ np.diag(np.sqrt(w)) @ V.T; Si = np.linalg.inv(S)
w2, V2 = np.linalg.eigh(Si @ C @ Si)
print("# relative eigenvalues (det-normalised):", np.round(np.sort(w2), 8), " d =", np.sqrt(np.sum(np.log(w2) ** 2)))
def G(t): M = S @ V2 @ np.diag(w2 ** t) @ V2.T @ S; M = (M + M.T) / 2; return M / np.linalg.det(M) ** (1 / r)
t0, t1 = (float(sys.argv[4]), float(sys.argv[5])) if len(sys.argv) > 5 else (0.0, 1.0)
for t in np.linspace(t0, t1, npts):
    g = G(t); lo, u, cnt = K.kappa(g, tol=1e-8); hi = max(u, lo * (1 + 1e-8))
    GQ = matrix(QQ, r, r, lambda i, j: RR(g[i, j]).nearby_rational(max_error=1e-13)); GQ = (GQ + GQ.transpose()) / 2
    p2, _ = psi2_min(GQ); e = float(p2 / GQ.det() ** ((r + 1) / 2))
    print("%.5f %.8f %.8f %.8f %.8f" % (t, lo ** ((r) / 2) * e, hi ** (r / 2) * e, lo, e), flush=True)
