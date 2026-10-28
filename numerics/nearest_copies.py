"""Nearest copies of one lattice to another (both det 1, invariant metric d of Theorem D5): all bases B of G1's lattice
whose Gram matrix B G1 B^T matches R within a tolerance, ranked by d(R, B G1 B^T); prints the k nearest distinct copies
C = B^{-1} R B^{-T} (in G1's coordinates).  Used for the geodesic scans between peaks of Phi_n.
Usage: nearest_copies.py n G1 R [k] [tol]   (G1, R: names pivE6, D5s, D5, D4, A2A2, Mstar, pivD5, bcc, fcc, pivA4)"""
import sys, itertools, numpy as np
def cartan(t, n):
    G = 2 * np.eye(n) - np.eye(n, k=1) - np.eye(n, k=-1)
    if t == 'D': G[n-1, n-2] = G[n-2, n-1] = 0; G[n-1, n-3] = G[n-3, n-1] = -1
    if t == 'E': G[n-1, n-2] = G[n-2, n-1] = 0; G[n-1, 2] = G[2, n-1] = -1
    return G
def proj(H, k=0):
    keep = [i for i in range(len(H)) if i != k]
    return H[np.ix_(keep, keep)] - np.outer(H[keep, k], H[k, keep]) / H[k, k]
def n1(G): return G / np.linalg.det(G) ** (1 / len(G))
A2 = cartan('A', 2)
LAT = {"pivE6": lambda: proj(cartan('E', 6)), "D5s": lambda: np.linalg.inv(cartan('D', 5)), "D5": lambda: cartan('D', 5),
       "D4": lambda: cartan('D', 4), "A2A2": lambda: np.kron(A2, A2), "pivD5": lambda: proj(cartan('D', 5), 3),
       "Mstar": lambda: np.array([[18,-9,0,0],[-9,18,-9,-9],[0,-9,16,1],[0,-9,1,16]], float),
       "bcc": lambda: np.linalg.inv(cartan('A', 3)), "fcc": lambda: cartan('A', 3), "pivA4": lambda: np.eye(3) + 0.5}
def logdist(R, S):
    w = np.linalg.eigvals(np.linalg.solve(R, S)).real; return float(np.sqrt(np.sum(np.log(w) ** 2)))
def copies(G, R, tol=0.6, k=5):
    n = len(G); G = n1(G); R = n1(R); rmax = max(np.diag(R)) * (1 + tol); Gi = np.linalg.inv(G)
    rad = [int(np.sqrt(rmax * Gi[i, i])) for i in range(n)]
    ks = np.array([v for v in itertools.product(*[range(-r, r + 1) for r in rad]) if any(v)])
    nv = np.einsum('ij,jk,ik->i', ks, G, ks); m = nv <= rmax; ks, nv = ks[m], nv[m]; ip = ks @ G @ ks.T
    cand = [np.where(np.abs(nv - R[i, i]) <= tol * R[i, i])[0] for i in range(n)]; found = {}
    def rec(i, ch):
        if i == n:
            B = ks[ch]
            if abs(round(np.linalg.det(B))) != 1: return
            d = logdist(R, B @ G @ B.T); Bi = np.linalg.inv(B); C = Bi @ R @ Bi.T
            key = tuple(np.round(C, 6).ravel())
            if key not in found or d < found[key][0]: found[key] = (d, C)
            return
        for c in cand[i]:
            if all(abs(ip[c, ch[j]] - R[i, j]) <= tol * np.sqrt(R[i, i] * R[j, j]) for j in range(i)): rec(i + 1, ch + [c])
    rec(0, [])
    return G, sorted(found.values(), key=lambda x: x[0])[:k]
if __name__ == "__main__":
    n, a, b = int(sys.argv[1]), sys.argv[2], sys.argv[3]; k = int(sys.argv[4]) if len(sys.argv) > 4 else 5
    tol = float(sys.argv[5]) if len(sys.argv) > 5 else 0.6
    G, cs = copies(LAT[a](), LAT[b](), tol, k)
    np.set_printoptions(precision=10, suppress=True, linewidth=200)
    print("G1 (%s, det 1):" % a, repr(G.tolist()))
    for d, C in cs: print("copy of %s at d = %.6f:" % (b, d), repr(C.tolist()))
