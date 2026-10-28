# Floating-point landscape of the E6 first-order radius r = min_{|z|=1} max_f g_f.z (orthonormal coordinates),
# over the C_f of the 19,440 Psi-minimal flags (e6_A.sobj).  NUMERICAL ONLY.  Run with sage -python e6_radius.py [starts].
import numpy as np, sys
from sage.all import load, RDF
from scipy.optimize import linprog
A, H, Ys = load("../certificates/e6_A.sobj")     # written by ../certificates/e6_flags.sage
A = np.array(A.change_ring(RDF)); H = np.array(H.change_ring(RDF))
L = np.linalg.cholesky(H)                    # H = L L^T, z = L^T y
G = A @ np.linalg.inv(L).T                   # g_f in orthonormal coords
print("norms", np.unique(np.round((G*G).sum(1), 9)))
rng = np.random.default_rng(1)
Z = rng.standard_normal((20000, 20)); Z /= np.linalg.norm(Z, axis=1, keepdims=True)
v = (Z @ G.T).max(1); print("random dirs: min %.4f  1%% %.4f  median %.4f" % (v.min(), np.percentile(v, 1), np.median(v)))
best = []
for s in range(int(sys.argv[1]) if len(sys.argv) > 1 else 300):
    c = rng.standard_normal(20); x = None
    for it in range(30):
        res = linprog(-c, A_ub=G, b_ub=np.ones(len(G)), bounds=[(None, None)] * 20, method="highs")
        if res.status != 0: break
        if x is not None and np.allclose(res.x, x): break
        x = res.x; c = x
    if x is not None: best.append(np.linalg.norm(x))
best = np.sort(best)[::-1]
print("polar vertex search, top |x|^2:", np.round(best[:8] ** 2, 6), " -> r =", 1 / best[0])
vals, cnt = np.unique(np.round(best ** 2, 5), return_counts=True)
print("local maxima |x|^2 (count):", list(zip(vals[::-1][:10], cnt[::-1][:10])))
