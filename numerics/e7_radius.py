# Floating-point landscape of the E7 first-order radius r = min_{|z|=1} max_f g_f.z (orthonormal coordinates),
# over the distinct C_f of all 544,320 Psi-minimal flags (e7_flags.npz, 118 MB, written by e7_dump.sage in about
# 30 s).  NUMERICAL ONLY.
#   python3 e7_radius.py <starts> <hours> [seed]   (the paper: 300 starts, about 4 hours)
import numpy as np, sys, time
from scipy.optimize import linprog
d = np.load("e7_flags.npz"); G = d["F"]; dirs = d["dirs"]; subs = d["subsets"]
N, D = G.shape; print("distinct flags", N, "dim", D, flush=True)
starts = int(sys.argv[1]) if len(sys.argv) > 1 else 300; hours = float(sys.argv[2]) if len(sys.argv) > 2 else 10
t0 = time.time(); rng = np.random.default_rng(int(sys.argv[3]) if len(sys.argv) > 3 else 7)
target = 5 / 6 * np.sqrt(7 / 6)
# 1. fundamental-weight directions: slope = max_f g_f.z for unit z (r <= slope); Lambda = I - cP gives -min, use -z too
s = np.array([max((G @ z).max(), 0) for z in dirs]); s2 = np.array([(G @ -z).max() for z in dirs])
print("weight directions: best max_f g.z over +z: %.8f (%s), over -z: %.8f (%s); target %.8f"
      % (s.min(), subs[s.argmin()], s2.min(), subs[s2.argmin()], target), flush=True)
order = np.argsort(s2)[:6]; print("  best -z:", [(subs[i], round(s2[i], 8)) for i in order], flush=True)
# 2. random directions
Z = rng.standard_normal((10000, D)); Z /= np.linalg.norm(Z, axis=1, keepdims=True)
v = np.concatenate([(G @ Z[i:i + 500].T).max(0) for i in range(0, 10000, 500)])
print("random dirs: min %.4f  1%% %.4f  median %.4f  (%.0fs)" % (v.min(), np.percentile(v, 1), np.median(v), time.time() - t0), flush=True)
# 3. polar-vertex ascents on P = {x : g_f.x <= 1}: r = 1/max|x|
best = []
for st in range(starts):
    if time.time() - t0 > hours * 3600: print("time box reached after", st, "starts", flush=True); break
    c = rng.standard_normal(D); x = None
    for it in range(30):
        res = linprog(-c, A_ub=G, b_ub=np.ones(N), bounds=[(None, None)] * D, method="highs")
        if res.status != 0: break
        if x is not None and np.allclose(res.x, x, atol=1e-9): break
        x = res.x; c = x
    if x is not None:
        nx = np.linalg.norm(x); act = int((np.abs(G @ x - 1) < 1e-7).sum()); best.append((nx, act))
        print("start %d: |x|^2 %.8f  r %.8f  active %d  (%.0fs)" % (st, nx * nx, 1 / nx, act, time.time() - t0), flush=True)
b = sorted(best, reverse=True)
print("top |x|^2:", [round(nx * nx, 6) for nx, _ in b[:8]], " r =", 1 / b[0][0], " target r =", target)
vals, cnt = np.unique(np.round([nx * nx for nx, _ in b], 5), return_counts=True)
print("local maxima |x|^2 (count):", list(zip(vals[::-1][:10], cnt[::-1][:10])))
print("runs reaching target (|r - target| < 1e-6):", sum(abs(1 / nx - target) < 1e-6 for nx, _ in b), "of", len(b))
