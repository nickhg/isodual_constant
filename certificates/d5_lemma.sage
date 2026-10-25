r"""
Theorem D5 in the invariant metric: eta(G) <= 6 whenever d(G, G0) = |log(G0^{-1/2} G G0^{-1/2})|_F <= rho, for G0 a
Gram matrix of D5 and G of the same determinant, with equality only at G0.  By isoduality (G -> J G^{-1} J is an
isometry of d) the same holds about every Gram matrix of D5^*.
Exact arithmetic (QQ) for the polytope; RIF for the final scalar inequality.

Notation.  G0 = Cartan matrix of D5 (det 4), G = G0^{1/2} e^Lam G0^{1/2}, Lam traceless, t = |Lam|_F,
Y = G0^{1/2} Lam G0^{1/2} (linear functionals of Lam are rational functionals of Y).  For a flag f = (U_1 < ... < U_4),
h_f(G) = prod_{k<=4} det(U_k G U_k^T) and eta(G) = min_f h_f(G) / det(G)^3, det(G) = det(G0).
  log h_f(G) <= log h_f(G0) + L_f(Y) + 4 tr E      (log det is concave in G, so it lies below its tangent at G0;
                                                    G - G0 = G0^{1/2}(Lam + E)G0^{1/2}, E = e^Lam - I - Lam >= 0, and
                                                    each prefix contributes tr(P_k E) <= tr E, P_k a projection)
  tr E <= t^2/2 exp(sqrt(4/5) t)                   (each eigenvalue of a traceless rank-5 Lam has |.| <= sqrt(4/5) t)
  min_f L_f(Y) <= -r t                             (r = inradius of Gamma = conv{L_f} in the dual norm; exact facets)
Hence log eta(G) - log 6 <= -r t + 2 t^2 exp(sqrt(4/5) t) < 0 for 0 < t <= rho whenever r > 2 rho exp(sqrt(4/5) rho).
The same argument for D4 (n = 4, 144 minimal flags, 9 directions, remainder 3 tr E, eigenvalues <= sqrt(3/4) t)
is selected by a second argument.  In general log eta(G) - log eta(G0) <= -r t + (n-1) t^2/2 exp(sqrt((n-1)/n) t),
so the lemma needs r > (n-1)/2 rho exp(sqrt((n-1)/n) rho).
Usage: sage d5_lemma.sage rho [n]    (rho rational, e.g. 1/10; n = 5 (D5, default) or 4 (D4))
"""
import sys, time
load("psi.sage")
rho = QQ(sys.argv[1]) if len(sys.argv) > 1 else None
n = int(sys.argv[2]) if len(sys.argv) > 2 else 5
G0 = root_gram('D', n)
assert G0.det() == 4
p0 = psi2_min(G0)[0]
eta0 = p0 / G0.det() ** (QQ(n + 1) / 2)
assert eta0 == {5: 6, 4: 3}[n], "eta(D5) = 6, eta(D4) = 3"
dimY = n * (n + 1) // 2 - 1
t0 = time.time()
bases = psi2_all(G0, p0)
print("Psi-minimal bases: %d (%.0fs)" % (len(bases), time.time() - t0))
G0i = G0.inverse()
# traceless directions: Y = G0^{1/2} Lam G0^{1/2} with tr(G0^{-1} Y) = 0
Ys = []
for i in range(n):
    for j in range(i, n):
        E = matrix(QQ, n, n); E[i, j] = 1; E[j, i] = 1
        Ys.append(E - (G0i * E).trace() / n * G0)
Ys = [Ys[p] for p in matrix(QQ, [Y.list() for Y in Ys]).pivot_rows()]; assert len(Ys) == dimY
N = matrix(QQ, dimY, dimY, lambda a, b: (G0i * Ys[a] * G0i * Ys[b]).trace())   # |Lam|_F^2 = y^T N y
fun = lambda W: vector(QQ, [(W * Y).trace() for Y in Ys])
def Lflag(U):
    W = matrix(QQ, n, n)
    for k in range(1, n):
        Uk = U.matrix_from_rows(range(k))
        W += Uk.T * (Uk * G0 * Uk.T).inverse() * Uk
    return fun(W)
LF = sorted(set(tuple(Lflag(U)) for _, U in bases))
print("distinct flag gradients: %d" % len(LF))
# Inradius lower bound without facet enumeration (the facets of a 14-dimensional polytope with 720 vertices are too
# many): choose a basis e_1..e_14 of the gradient space, orthogonal for the dual norm g^T N^{-1} g; for each signed axis
# an exact LP finds the largest t with t (+-e_j) in Gamma = conv{L_f}.  Gamma then contains the cross-polytope with
# half-axes a_j = min(t_j^+, t_j^-) |e_j|_*, hence the dual-norm ball of radius r = (sum_j a_j^{-2})^{-1/2}.
Ni = N.inverse()
E = []
for k in range(dimY):
    v = vector(QQ, [1 if i == k else 0 for i in range(dimY)])
    for w in E: v -= (v * Ni * w) / (w * Ni * w) * w
    E.append(v)
LFv = [vector(QQ, g) for g in LF]
def axis_t(d):
    P = MixedIntegerLinearProgram(maximization=True, solver="PPL")
    lam = P.new_variable(nonnegative=True); t = P.new_variable()
    P.add_constraint(sum(lam[f] for f in range(len(LFv))) == 1)
    for c in range(dimY):
        P.add_constraint(sum(lam[f] * LFv[f][c] for f in range(len(LFv))) == t[0] * d[c])
    P.set_objective(t[0])
    return QQ(P.solve())
inv2 = 0
for k, e in enumerate(E):
    tp, tm = axis_t(e), axis_t(-e)
    assert tp > 0 and tm > 0, "0 is not interior along an axis: no first-order decrease"
    a2 = min(tp, tm) ** 2 * (e * Ni * e)          # a_j^2, exact
    inv2 += 1 / a2
    print("axis %2d: t+ = %s, t- = %s, a_j = %.5f  (%.0fs)" % (k, tp, tm, float(sqrt(a2)), time.time() - t0))
r = 1 / sqrt(RIF(inv2))
print("inradius >= r = %s (cross-polytope bound)" % r)
# the largest admissible rho: 2 rho exp(sqrt(4/5) rho) < r
cf = lambda x: (n - 1) * x / 2 * exp(sqrt(RR(n - 1) / n) * x) - r.lower()
lo, hi = RR(0), RR(1)
for _ in range(60):
    mid = (lo + hi) / 2
    if cf(mid) < 0: lo = mid
    else: hi = mid
print("rho_max ~ %.6f" % lo)
if rho is not None:
    R = RIF; rh = R(rho)
    margin = r - R(n - 1) / 2 * rh * exp(sqrt(R(n - 1) / n) * rh)
    print("rho = %s: margin r - (n-1)/2 rho exp(sqrt((n-1)/n) rho) = %s" % (rho, margin))
    assert margin > 0
    print("LEMMA HOLDS: eta(G) <= %s for d(G, G0) <= %s about D%d (equality only at G0); likewise about its dual by isoduality" % (eta0, rho, n))
