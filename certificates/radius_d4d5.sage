r"""
The first-order radius r (inradius of Gamma = conv{L_f} in the dual norm; notation of d5_lemma.sage) at D4 and D5.
  D4: r = 2/3 exactly, by enumerating every facet of Gamma (9-dimensional, 144 vertices) with Normaliz in exact
      integer arithmetic; with it, eta(G) <= 3 for d(G, G0) <= 33/100 about D4 (and D4^* = D4 up to scaling).
  D5: r <= 4/(3 sqrt5), from one direction Lam_0 = I - 5P, P the orthogonal projection onto a spinor fundamental
      weight: Gamma contains a dual-norm ball of radius r only if min_f L_f(Lam_0) <= -r |Lam_0|_F.
  D4 also: the minimum is attained along Lam_0 = 3(I - 2P), P the projection onto the plane of the two half-spin
      weights (any two of the three minuscule weights, by triality).
In Y-coordinates (Y = G0^{1/2} Lam G0^{1/2}), Lam = I - (n/k)P with P onto the span of the fundamental weights in S
(|S| = k) is Y = G0 - (n/k) E_S ((G0^{-1})_{SS})^{-1} E_S^T, with |Lam|_F^2 = tr(G0^{-1} Y G0^{-1} Y), and
L_f(Lam) = tr(Q_f Y), Q_f = sum_{k<=n-1} U_k^T (U_k G0 U_k^T)^{-1} U_k.
Needs the normaliz binary on the PATH (the D4 facets take about 30 s).
Usage: sage radius_d4d5.sage
"""
import itertools, os, subprocess, tempfile, time
load("psi.sage")
t0 = time.time()

def setup(n):
    G0 = root_gram('D', n); G0i = G0.inverse()
    p0 = psi2_min(G0)[0]
    Qs = []
    for _, U in psi2_all(G0, p0):
        Q = matrix(QQ, n, n)
        for k in range(1, n):
            Uk = U.matrix_from_rows(range(k)); Q += Uk.T * (Uk * G0 * Uk.T).inverse() * Uk
        Qs.append(Q)
    return G0, G0i, Qs

def ray_bound(n, G0, G0i, Qs, S):
    """Upper bound on r from Lam = I - (n/k)P, P onto the fundamental weights in S: (min, #flags at min, |Lam|^2)."""
    ES = matrix(QQ, n, len(S), lambda i, j: 1 if i == S[j] else 0)
    Y = G0 - QQ(n) / len(S) * ES * (ES.T * G0i * ES).inverse() * ES.T
    assert (G0i * Y).trace() == 0
    vals = [(Q * Y).trace() for Q in Qs]
    m = min(vals)
    return m, vals.count(m), (G0i * Y * G0i * Y).trace()

# ---- D5: upper bound along a spinor weight (nodes 3, 4 of root_gram('D', 5) have (G0^{-1})_ss = 5/4)
G0, G0i, Qs = setup(5)
print("D5: %d Psi-minimal bases, (G0^-1)_ii = %s" % (len(Qs), G0i.diagonal()))
for s in [3, 4]:
    m, c, nn = ray_bound(5, G0, G0i, Qs, (s,))
    assert (m, c, nn) == (-8/3, 120, 20)
    print("  spinor node %d: min_f L_f(Lam_0) = %s (%d flags), |Lam_0|^2 = %s  =>  r <= %s = %.6f"
          % (s, m, c, nn, -m / sqrt(nn), float(-m / sqrt(nn))))
print("D5: r <= 4/(3 sqrt5)")

# ---- D4: every facet of Gamma, exactly
n = 4
G0, G0i, Qs = setup(n)
dimY = n * (n + 1) // 2 - 1
Ys = []
for i in range(n):
    for j in range(i, n):
        E = matrix(QQ, n, n); E[i, j] = 1; E[j, i] = 1
        Ys.append(E - (G0i * E).trace() / n * G0)
Ys = [Ys[p] for p in matrix(QQ, [Y.list() for Y in Ys]).pivot_rows()]; assert len(Ys) == dimY
N = matrix(QQ, dimY, dimY, lambda a, b: (G0i * Ys[a] * G0i * Ys[b]).trace())   # |Lam|_F^2 = y^T N y
LF = sorted(set(tuple((Q * Y).trace() for Y in Ys) for Q in Qs))
print("D4: %d Psi-minimal bases, %d distinct gradients" % (len(Qs), len(LF)))
assert len(LF) == 144 and matrix(QQ, LF).rank() == dimY
den = lcm([QQ(x).denominator() for g in LF for x in g])
wd = tempfile.mkdtemp()
with open(os.path.join(wd, "g.in"), "w") as fh:
    fh.write("amb_space %d\npolytope %d\n" % (dimY + 1, len(LF)))
    for g in LF: fh.write(" ".join(str(ZZ(den * x)) for x in g) + "\n")
    fh.write("SupportHyperplanes\n")
subprocess.run(["normaliz", "-c", "-x=1", os.path.join(wd, "g.in")], check=True, stdout=subprocess.DEVNULL)
lines = open(os.path.join(wd, "g.out")).read().split("\n")
k = next(i for i, l in enumerate(lines) if l.strip().endswith("support hyperplanes:") and l.strip()[0].isdigit())
rows = [list(map(ZZ, lines[k + 1 + i].split())) for i in range(int(lines[k].split()[0]))]
# each row (a, c) means a . (den g) + c >= 0 on Gamma; check it against every vertex, independently of Normaliz
V = [vector(ZZ, [den * x for x in g]) for g in LF]
for r_ in rows:
    a, c = vector(ZZ, r_[:-1]), r_[-1]
    sl = [a * v + c for v in V]
    assert min(sl) == 0 and c > 0
    assert matrix(QQ, [list(v) for v, s in zip(V, sl) if s == 0]).rank() == dimY   # vertices on it span a facet (0 not on it)
d2 = [QQ(r_[-1]) ** 2 / (den ** 2 * (vector(QQ, r_[:-1]) * N * vector(QQ, r_[:-1]))) for r_ in rows]
print("D4: %d facets, each re-checked against all 144 vertices; nearest d^2 = %s on %d facets; next d = %s"
      % (len(rows), min(d2), d2.count(min(d2)), ", ".join("%.6f" % sqrt(x) for x in sorted(set(d2))[1:4])))
assert min(d2) == 4/9
for S in itertools.combinations([0, 2, 3], 2):     # the three minuscule nodes of root_gram('D', 4)
    m, c, nn = ray_bound(n, G0, G0i, Qs, S)
    assert -m / sqrt(nn) == 2/3
    print("  Lam_0 = 3(I - 2P), P onto weights %s: min_f L_f = %s (%d flags), |Lam|^2 = %s" % (S, m, c, nn))
r = QQ(2/3); rho = QQ(33/100); R = RIF
margin = R(r) - R(n - 1) / 2 * R(rho) * exp(sqrt(R(n - 1) / n) * R(rho))
print("D4: r = 2/3 exactly; rho = %s: margin r - (3/2) rho exp(sqrt(3/4) rho) = %s" % (rho, margin))
assert margin > 0
print("LEMMA HOLDS: eta(G) <= 3 for d(G, G0) <= 33/100 about D4 (equality only at G0)  (%.0fs)" % (time.time() - t0))
