"""Nearest copies of R to G1 (det 1, invariant metric) by local search over bases: from random unimodular starts, apply
elementary moves b_i += s b_j (s = +-1) and row swaps while d(R, B G1 B^T) decreases; collect the distinct local minima
C = B^{-1} R B^{-T} (copies of R in G1's coordinates), nearest first.  Usage: nearest_local.py G1 R [starts] [k]"""
import sys, numpy as np
sys.path.insert(0, __import__("os").path.dirname(__import__("os").path.abspath(__file__)))
from nearest_copies import LAT, n1, logdist
def search(G, R, starts=300, k=6, seed=0):
    rng = np.random.default_rng(seed); n = len(G); G = n1(G); R = n1(R); found = {}
    for s in range(starts):
        B = np.eye(n, dtype=int)
        for _ in range(rng.integers(0, 3 * n)):
            i, j = rng.choice(n, 2, replace=False); B[i] += rng.choice((-1, 1)) * B[j]
        d = logdist(R, B @ G @ B.T)
        while True:
            best = (d, None)
            for i in range(n):
                for j in range(n):
                    if i == j: continue
                    for sg in (1, -1):
                        C = B.copy(); C[i] += sg * C[j]; dd = logdist(R, C @ G @ C.T)
                        if dd < best[0] - 1e-12: best = (dd, C)
                    if i < j:
                        C = B.copy(); C[[i, j]] = C[[j, i]]; dd = logdist(R, C @ G @ C.T)
                        if dd < best[0] - 1e-12: best = (dd, C)
            if best[1] is None: break
            d, B = best
        Bi = np.linalg.inv(B); Cp = Bi @ R @ Bi.T; key = tuple(np.round(Cp, 6).ravel())
        if key not in found or d < found[key][0]: found[key] = (d, Cp)
    return G, sorted(found.values(), key=lambda x: x[0])[:k]
if __name__ == "__main__":
    a, b = sys.argv[1], sys.argv[2]; starts = int(sys.argv[3]) if len(sys.argv) > 3 else 300; k = int(sys.argv[4]) if len(sys.argv) > 4 else 6
    G, cs = search(LAT[a](), LAT[b](), starts, k)
    print("G1 (%s, det 1):" % a, repr(G.tolist()))
    for d, C in cs: print("copy of %s at d = %.6f:" % (b, d), repr(C.tolist()))
