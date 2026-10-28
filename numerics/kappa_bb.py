#!/usr/bin/env python3
"""kappa(M) = sup_phi min_{s != 0} |s|^2 / (1 - ||phi(s)||^2) by best-first branch and bound over the phi-torus.

Only s with |s|^2 <= (4/3) lambda_1^2 can attain the minimum (the shortest vector gives at most that), and they are
enumerated completely (box |s_i| <= sqrt(B (G^-1)_ii)).  A phi-box gets the upper bound min_s |s|^2/(1 - dmax_s^2),
with dmax_s the largest distance to Z of phi(s) over the box; its centre gives a lower bound.  The domain is
phi_1 in [0, 1/2] (phi -> -phi), phi_i in [0, 1] otherwise.  Floating point: a numerical value, not a certificate.
Usage: kappa_bb.py  (tests pi_v E6, D4, A4)
"""
import heapq, itertools, numpy as np

def relevant_vectors(G):
    n = len(G)
    Gi = np.linalg.inv(G)
    l1 = min(G[i, i] for i in range(n))
    # lambda_1 from a small search, then all s with |s|^2 <= (4/3) lambda_1^2
    for _ in range(2):
        B = l1 * 4 / 3 * (1 + 1e-9)
        R = [int(np.floor(np.sqrt(B * Gi[i, i]) + 1e-9)) for i in range(n)]
        ks = np.array([k for k in itertools.product(*[range(-r, r + 1) for r in R])
                       if any(k) and next(x for x in k if x) > 0])
        nv = np.einsum('ij,jk,ik->i', ks, G, ks)
        l1 = nv.min()
    m = nv <= l1 * 4 / 3 * (1 + 1e-9)
    return ks[m], nv[m]

def kappa(G, tol=1e-6, maxboxes=2_000_000):
    ks, nv = relevant_vectors(G)
    n = len(G)
    kpos, kneg = np.clip(ks, 0, None), np.clip(ks, None, 0)
    def ub(lo, hi):
        tlo = kpos @ lo + kneg @ hi
        thi = kpos @ hi + kneg @ lo
        half = np.ceil(tlo - 0.5) + 0.5
        contains_half = (half <= thi) | (thi - tlo >= 1)
        d = lambda t: np.abs(t - np.round(t))
        dmax = np.where(contains_half, 0.5, np.maximum(d(tlo), d(thi)))
        return np.min(nv / (1 - dmax ** 2))
    def val(phi):
        t = ks @ phi
        dd = np.abs(t - np.round(t))
        return np.min(nv / (1 - dd ** 2))
    lo0 = np.zeros(n); hi0 = np.ones(n); hi0[0] = 0.5
    best = val((lo0 + hi0) / 2)
    heap = [(-ub(lo0, hi0), 0, lo0, hi0)]
    cnt = 1
    while heap:
        u, _, lo, hi = heapq.heappop(heap)
        u = -u
        if u <= best * (1 + tol) or cnt > maxboxes:
            return best, u, cnt
        d = int(np.argmax(hi - lo)); mid = (lo[d] + hi[d]) / 2
        for a, b in ((lo, np.where(np.arange(n) == d, mid, hi)), (np.where(np.arange(n) == d, mid, lo), hi)):
            c = val((a + b) / 2)
            if c > best: best = c
            uu = ub(a, b)
            if uu > best * (1 + tol):
                cnt += 1
                heapq.heappush(heap, (-uu, cnt, a, b))
    return best, best, cnt

if __name__ == '__main__':
    import sys, os, time
    sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
    import phi6_local as P
    E6 = P.E6
    piE6 = P.norm1(E6[1:, 1:] - np.outer(E6[0, 1:], E6[0, 1:]) / E6[0, 0])
    D4 = P.norm1(np.array([[2, -1, 0, 0], [-1, 2, -1, -1], [0, -1, 2, 0], [0, -1, 0, 2]], float))
    A4 = P.norm1(2 * np.eye(4) - np.eye(4, k=1) - np.eye(4, k=-1))
    for name, G, exact in (("D4", D4, 36 * np.sqrt(2) / 35), ("A4", A4, 25 / 24 * 2 * 5 ** -0.25),
                           ("piE6", piE6, (64 / 3) ** 0.2)):
        t = time.time(); k, u, c = kappa(G)
        print("%-5s kappa in [%.6f, %.6f]  (expected %.5f)  %d boxes  %.1fs" % (name, k, u, exact, c, time.time() - t))
