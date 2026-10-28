# Phi_5 along the isodual geodesic G(t) = H_D (H_D^{-1} H_A)^t (det-normalised) from D4 (t = 0) through M_* (t = log_2(6/5))
# to A2 (x) A2 (t = 1), for the figure: kappa by branch and bound (kappa_bb, honest interval, tol 1e-8, box cap 2e6),
# eta by exact Psi-minimisation (psi.sage) of the Gram matrix rationalised to 1e-13.  Floating point (numerical).
#   (run from numerics/) sage geodesic_phi5.sage > geodesic_phi5.txt
import sys, numpy as np
sys.path.insert(0, ".")
import kappa_bb as K
load("../certificates/psi.sage")
HD = np.array([[2,-1,0,0],[-1,2,-1,-1],[0,-1,2,0],[0,-1,0,2]], float)
HA = np.array([[6,-3,0,0],[-3,6,-3,-3],[0,-3,4,1],[0,-3,1,4]], float)
w, V = np.linalg.eig(np.linalg.solve(HD, HA)); w = w.real; V = V.real; Vi = np.linalg.inv(V)
def G(t):
    M = HD @ (V @ np.diag(w ** t) @ Vi); M = (M + M.T) / 2; return M / np.linalg.det(M) ** 0.25
for t in sorted(set([round(x, 4) for x in np.linspace(-0.3, 1.3, 65)] + [0.0, float(log(6/5) / log(2)), 1.0])):
    g = G(t)
    lo, u, cnt = K.kappa(g, tol=1e-8); hi = max(u, lo * (1 + 1e-8))
    GQ = matrix(QQ, 4, 4, lambda i, j: RR(g[i, j]).nearby_rational(max_error=1e-13)); GQ = (GQ + GQ.transpose()) / 2
    p2, _ = psi2_min(GQ); e = float(p2 / GQ.det() ** (5/2))
    print("%.6f %.8f %.8f %.8f" % (t, lo ** 2 * e, hi ** 2 * e, e), flush=True)
