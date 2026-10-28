# E7: eta(E7) = 144 = H7, and the exact design identity  sum_f L_f(Lam)^2 = N (28/27) |Lam|^2,  sum_f L_f = 0,
# over all Psi-minimal flags (chains A1<A2<A3<D4<D5<E6 from roots; see e7_count.sage for completeness).
# Sums over the 544,320 chains are done by dynamic programming over the distinct sublattices.
import time; t0 = time.time()
n = 7
G0 = matrix(ZZ, CartanMatrix(['E', 7])); G0q = G0.change_ring(QQ); G0i = G0q.inverse()
R = pari(G0).qfminim(2); roots = matrix(ZZ, R[2]).columns(); roots += [-r for r in roots]
TARGET = [2, 3, 4, 4, 4, 3]
psi2 = prod(TARGET) * G0.det()
print("Psi^2 =", psi2, " eta(E7) = Psi^2/det^4 =", psi2 / G0.det()**((n + 1) / 2))
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
kids = {}
def children(B):
    k0 = key(B)
    if k0 in kids: return kids[k0]
    out = {}; k = B.nrows()
    for w in roots:
        C = matrix(ZZ, B.rows() + [w])
        if C.rank() == k + 1 and (C * G0 * C.T).det() == TARGET[k]:
            out.setdefault(key(C), C)
    kids[k0] = list(out.items()); return kids[k0]
W = {}
def w_of(kk, B):
    if kk not in W: W[kk] = wvec(B)
    return W[kk]
dp = {}
Z = zero_vector(QQ, D); ZM = zero_matrix(QQ, D, D)
def suffix(kk, B):
    # (number of completions, sum of suffix vectors, sum of suffix outer products) below node B
    if B.nrows() == 6: return (1, Z, ZM)
    if kk in dp: return dp[kk]
    c, s, P = 0, Z, ZM
    for k2, C in children(B):
        cc, sc, Pc = suffix(k2, C); w = w_of(k2, C)
        c += cc; s = s + cc * w + sc
        P = P + cc * w.outer_product(w) + w.outer_product(sc) + sc.outer_product(w) + Pc
    dp[kk] = (c, s, P); return dp[kk]
lines = {}
for r in roots:
    B = matrix(ZZ, [r]); lines.setdefault(key(B), B)
N, S1, S2 = 0, Z, ZM
for kk, B in lines.items():
    cc, sc, Pc = suffix(kk, B); w = w_of(kk, B)
    N += cc; S1 = S1 + cc * w + sc; S2 = S2 + cc * w.outer_product(w) + w.outer_product(sc) + sc.outer_product(w) + Pc
print("flags", N, "%.0fs" % (time.time() - t0))
print("sum L_f = 0:", S1 == 0)
c = QQ(N) * 28 / 27
print("design sum L_f L_f^T == N(28/27) H:", S2 == c * H, " c =", c)
# one explicit attaining flag (root flag of the Cartan basis is not a chain of these types in general; take a DFS one)
B = list(lines.values())[0]; chain = [B]
while chain[-1].nrows() < 6: chain.append(children(chain[-1])[0][1])
Hi = H.inverse()
a = sum(w_of(key(C), C) for C in chain)
print("explicit flag, added roots:", [list(chain[k].row(k)) for k in range(6)])
print("its prefix dets:", [(C * G0 * C.T).det() for C in chain], " dual norm^2 =", a * Hi * a)
print("done %.0fs" % (time.time() - t0))
# upper bounds on r along Lam = I - (n/k) P_V, V spanned by k fundamental weights: min over chains by DP
import itertools
Bm = matrix(QQ, [Y.list() for Y in Ys])
def minchain(y):
    memo2 = {}
    def m(kk, B):
        if B.nrows() == 6: return 0
        if kk in memo2: return memo2[kk]
        v = min(w_of(k2, C) * y + m(k2, C) for k2, C in children(B)); memo2[kk] = v; return v
    return min(w_of(kk, B) * y + m(kk, B) for kk, B in lines.items())
res = []
for k in range(1, 4):
    for Sset in itertools.combinations(range(n), k):
        Es = matrix(QQ, n, k, lambda a, b: 1 if a == Sset[b] else 0)
        Y = G0q - QQ(n) / k * Es * G0i.matrix_from_rows_and_columns(Sset, Sset).inverse() * Es.T
        y = Bm.solve_left(vector(QQ, Y.list()))
        mm = minchain(y); res.append((mm**2 / (y * H * y), Sset, mm, y * H * y))
res.sort()
for r2, Sset, mm, nr in res[:4]:
    print("weights", Sset, " min L_f =", mm, " |Lam|^2 =", nr, " ratio =", sqrt(r2), "= %.6f" % sqrt(r2).n())
for name, r in [("design", RIF(28).sqrt() / 27)]:
    for t in [QQ(1)/20, QQ(3)/50, QQ(1)/16, QQ(13)/200]:
        T = RIF(t); v = -r * T + 3 * T**2 * (RIF(6/7).sqrt() * T).exp()
        print(name, t, "neg" if v.upper() < 0 else "NOT", v.upper())
