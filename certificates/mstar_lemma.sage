r"""
Local lemma at M*: Phi_5(G) <= Phi_5(M*) = 2^16/(3^4 5^3) whenever d(G, G0) = |log(G0^{-1/2} G G0^{-1/2})|_F <= rho
(G, G0 of equal det), with equality only at G0.  Exact arithmetic (QQ) for the polytope and the torus cover; RIF for
the final scalar inequalities.  Template: pivd5_lemma.sage.

Notation.  G0 = P = 18 M* (primitive integral, det 135^2, basis of mstar_data.json), phi0 = (0, 0, 1/3, 1/3),
G = G0^{1/2} e^Lam G0^{1/2}, Lam traceless, t = |Lam|_F, Y = G0^{1/2} Lam G0^{1/2}.
Phi~(G) = kappa_raw(G)^2 det(G)^{-1/2} eta~(G), eta~ = min_f h_f / det^{3/2}, h_f = prod_{k<=3} det G_k^f.
Active at (G0, phi0), k_s = q_s/(1 - ||s.phi0||^2) = 18:  9 'third' s (s.phi0 = sigma_s/3 mod Z, q_s = 16) and
3 'integer' s (s.phi0 in Z, q_s = 18); 54 flags.  psi(tau) = -log(1 - tau^2), a = psi'(1/3) = 3/4,
psi'' increasing on [0, 1), psi''(tau) = 2(1 + tau^2)/(1 - tau^2)^2.
  q_s:    log q_s(G) <= log q_s(G0) + l_s(Y) + t^2/4                 (Hoeffding, as at pi_v D5)
  flags:  log h_f(G) <= log h_f(G0) + L_f(Y) + 3/2 t^2 e^{sqrt(3/4) t}  (log det concave)
  U = c t.  Exterior: some third s has ||s.phi|| <= 1/3 - U.  Then psi(||s.phi||) <= psi(1/3) - a U + psi''(1/3) U^2/2
     and log(Phi~/Phi~0) <= B t + rem(t) - 2 a U + psi''(1/3) U^2,  B >= max_{|Lam|=1} [2 l_s + min_f L_f].
  Far: phi in a cell where some short s has q_s(G0)/(1 - dmax_s^2) <= 18 e^{-gamma} (gap absorbs the ball, below).
  Interior: every third s has ||s.phi|| > 1/3 - U, and phi not far.  A rigorous torus cover (below) shows
     phi = +-phi0 + Delta mod Z^4 with |s.Delta| < 1/6 for every active s, where ||s.phi|| = 1/3 + sigma_s s.Delta (third) and |s.Delta| (integer) exactly;
     then sigma_s s.Delta > -U for all thirds, so |s.Delta| <= C_s U for every active s (C_s from the vertices of the
     bounded polytope {Delta : sigma_s s.Delta >= -1}).  Hence (psi convex, psi(u) <= u^2/(1 - u^2))
       third:   psi(||s.phi||) <= psi(1/3) + a sigma_s s.Delta + psi''(1/3 + C_s U) (C_s U)^2 / 2
       integer: psi(|s.Delta|) <= (C_s U)^2 / (1 - (C_s U)^2)
     and max_Delta min_s [l_s(Y) + a_s sigma_s s.Delta] = min_{g in Gamma_kappa} g(Y)  (LP duality), with
       Gamma_kappa = { sum_s lam_s L_s : lam >= 0, sum lam = 2, sum_{third} lam_s sigma_s s = 0 }  (factor 2: kappa^2),
     so log(Phi~/Phi~0) <= -r t + rem(t) + 2 e(t),  r = inradius of Gamma = Gamma_kappa + conv(L_f) (dual norm).
  rem(t) = t^2/2 + 3/2 t^2 e^{sqrt(3/4) t};  the lemma needs  r > rem/t + 2e/t  and  2 a c > B + rem/t + psi''(1/3) c^2 t
  at t = rho (all terms / t are increasing in t).
Data: mstar_data.json, written by mstar_sharp.py --dump (every property used is re-verified here exactly).
Usage: sage mstar_lemma.sage rho c      (rationals, e.g. 1/20 1)
"""
import json, sys, time, os, itertools
from fractions import Fraction as Fr
D = json.load(open(os.path.join(os.path.dirname(os.path.abspath(sys.argv[0])), "mstar_data.json")))
rho = QQ(sys.argv[1]); c = QQ(sys.argv[2])
G0 = matrix(QQ, D["G0"]); G0i = G0.inverse()
assert G0 == matrix(QQ, [[18, -9, 0, 0], [-9, 18, -9, -9], [0, -9, 16, 1], [0, -9, 1, 16]]) and G0.det() == 135^2
assert (135 * G0i).denominator() == 1                    # 135-modular (level divides 135)
phi0 = vector(QQ, [0, 0, 1/3, 1/3])
T = [vector(ZZ, s) for s in D["third"]]; I = [vector(ZZ, s) for s in D["int"]]
F = [[vector(ZZ, v) for v in f] for f in D["flags"]]
q = lambda s: s * G0 * s
dZ = lambda x: abs(x - round(x))
assert len(T) == 9 and len(I) == 3 and all(q(s) == 16 and dZ(s * phi0) == 1/3 for s in T) and all(q(s) == 18 and s * phi0 in ZZ for s in I)
sig = [1 if (s * phi0 - 1/3) in ZZ else -1 for s in T]
assert all((s * phi0 - sg / 3) in ZZ for s, sg in zip(T, sig))
# kappa_raw(G0, phi0) = 18 is attained exactly by these 12 (no other s below 18 at phi0: |s|^2 >= 16 = min, and
# k_s >= q_s >= 18 unless q_s < 18; the q_s = 16 vectors are listed below and all are 'third')
assert all(18 * G0i[i, i] < 25 for i in range(4))     # q_s < 18 forces |s_i| <= sqrt(18 (G0^-1)_ii) < 5
mins = [s for s in itertools.product(range(-4, 5), repeat=4) if any(s) and next(x for x in s if x) > 0 and q(vector(QQ, s)) < 18]
assert all(q(vector(QQ, s)) == 16 for s in mins) and sorted(map(tuple, T)) == sorted(map(tuple, map(vector, mins)))
print("k_s(G0, phi0) = 18 for 9 third and 3 integer vectors; the 9 minimal vectors are the thirds")
for f in F:
    assert len(f) == 3, "flags must have 3 vectors (the minors(3) primitivity test assumes it)"
    assert gcd((matrix(ZZ, f)).minors(3)) == 1
hv = set(prod((matrix(QQ, f[:k]) * G0 * matrix(QQ, f[:k]).T).det() for k in range(1, 4)) for f in F)
assert len(hv) == 1; h0 = hv.pop()
Phi0 = 18^2 / sqrt(G0.det()) * h0 / G0.det()^(3/2)
print("%d flags, prod_k det G_k = %s; Phi~(G0) = %s" % (len(F), h0, Phi0)); assert Phi0 == 2^16 / (3^4 * 5^3)

# traceless Y basis, norm Gram N, functionals (as in pivd5_lemma.sage)
Ys = []
for i in range(4):
    for j in range(i, 4):
        E = matrix(QQ, 4, 4); E[i, j] = 1; E[j, i] = 1
        Ys.append(E - (G0i * E).trace() / 4 * G0)
Ys = [Ys[p] for p in matrix(QQ, [Y.list() for Y in Ys]).pivot_rows()]; assert len(Ys) == 9
N = matrix(QQ, 9, 9, lambda a, b: (G0i * Ys[a] * G0i * Ys[b]).trace()); Ni = N.inverse()
fun = lambda W: vector(QQ, [(W * Y).trace() for Y in Ys])
ell = lambda s: fun(matrix(QQ, 4, 1, list(s)) * matrix(QQ, 1, 4, list(s)) / q(s))
LT = [ell(s) for s in T]; LI = [ell(s) for s in I]
def Lflag(f):
    W = matrix(QQ, 4, 4)
    for k in range(1, 4):
        U = matrix(QQ, f[:k]); W += U.T * (U * G0 * U.T).inverse() * U
    return fun(W)
LF = [Lflag(f) for f in F]
dn2 = lambda g: g * Ni * g

# C_s: |s.Delta| <= C_s for Delta in the polytope {sigma_j t_j.Delta >= -1}
Dpol = Polyhedron(ieqs=[[1] + list(sg * s) for s, sg in zip(T, sig)], base_ring=QQ)
assert Dpol.is_compact(), "the third vectors are not positively dependent: interior not bounded"
Cs = [max(abs(s * vector(QQ, v)) for v in Dpol.vertices_list()) for s in T + I]
CT, CI = max(Cs[:9]), max(Cs[9:])
print("interior polytope: %d vertices; C_third = %s, C_int = %s" % (Dpol.n_vertices(), CT, CI))

# torus cover.  For t <= rho the interior {||t_j.phi|| > 1/3 - U} lies in the U0-interior (U0 = c rho), which is covered
# by cells that are either near +-phi0 (|phi -+ phi0|_inf <= delta0 mod Z^4) or 'far': some short s (q_s <= 40) has
# K = q_s(G0)/(1 - dmax_s^2) <= 18 e^{-gamma} on the cell, with
#   2 gamma = (2 max_s |L_s|_* + |L_f1|_*) rho + rem(rho),
# since then for every G in the ball log(Phi~/Phi~0) <= 2 log(k_s(G, phi)/18) + L_f1(Y) + rem(t)
#   <= -2 gamma + 2 |L_s|_* t + 2 t^2/4 + |L_f1|_* t + 3/2 t^2 e^{..} <= 0     (any one flag f1, Hoeffding for q_s).
SF = [vector(ZZ, s) for s in itertools.product(range(-3, 4), repeat=4) if any(s) and next(x for x in s if x) > 0 and q(vector(QQ, s)) <= 40]
assert all(40 * G0i[i, i] < 16 for i in range(4))     # q_s <= 40 forces |s_i| <= 3: SF is complete (not needed for soundness)
LSn = max(sqrt(RIF(dn2(ell(s)))).upper() for s in SF); LF1 = min(sqrt(RIF(dn2(g))).upper() for g in LF)
Rr = RIF(rho)
gamma2 = (2 * LSn + LF1) * Rr + Rr^2 / 2 + RIF(3) / 2 * Rr^2 * exp(sqrt(RIF(3) / 4) * Rr)
Kthr = (18 * exp(-gamma2 / 2)).lower().exact_rational()
print("%d short vectors for the far test; max |L_s|_* %.4f, |L_f1|_* %.4f; far threshold K <= %.4f" % (len(SF), LSn, LF1, float(Kthr)))
U0 = c * rho
def dmax(lo_, hi_):        # max over [lo_, hi_] of the distance to Z
    if floor(lo_ - 1/2) != floor(hi_ - 1/2) or hi_ - lo_ >= 1: return QQ(1) / 2
    return max(dZ(lo_), dZ(hi_))
def srange(s, lo, hi):
    return (sum(s[i] * (lo[i] if s[i] > 0 else hi[i]) for i in range(4)), sum(s[i] * (hi[i] if s[i] > 0 else lo[i]) for i in range(4)))
def near(lo, hi):          # cell - (+-phi0 + m) has |s.Delta| < 1/6 for every active s (m: nearest integer shift)
    for p0 in (phi0, -phi0):
        m = vector(QQ, [round((lo[i] + hi[i]) / 2 - p0[i]) for i in range(4)])
        lo_ = [lo[i] - p0[i] - m[i] for i in range(4)]; hi_ = [hi[i] - p0[i] - m[i] for i in range(4)]
        if all(-QQ(1) / 6 < srange(s, lo_, hi_)[0] and srange(s, lo_, hi_)[1] < QQ(1) / 6 for s in T + I): return True
    return False
t0 = time.time(); stack = [([QQ(0)] * 4, [QQ(1)] * 4)]; ncell = nacc = next_ = nfar = 0
while stack:
    lo, hi = stack.pop(); ncell += 1
    if any(dmax(*srange(s, lo, hi)) <= 1/3 - U0 for s in T): next_ += 1; continue
    if near(lo, hi): nacc += 1; continue
    if any(q(s) / (1 - dmax(*srange(s, lo, hi))^2) <= Kthr for s in SF): nfar += 1; continue
    d = max(range(4), key=lambda i: hi[i] - lo[i]); m = (lo[d] + hi[d]) / 2
    assert hi[d] - lo[d] > QQ(1) / 2^14, "torus cover does not converge near %s" % [float(x) for x in lo]
    h1 = list(hi); h1[d] = m; l2 = list(lo); l2[d] = m
    stack += [(lo, h1), (l2, hi)]
nsplit = ncell - next_ - nfar - nacc; assert ncell == 2 * nsplit + 1     # a binary tree: leaves = splits + 1
print("torus cover of [0,1]^4: %d leaf cells (%d exterior, %d far, %d near +-phi0; U0 = %s), from %d splits (%d nodes) (%.0fs)"
      % (next_ + nfar + nacc, next_, nfar, nacc, U0, nsplit, ncell, time.time() - t0))
assert CT * U0 < QQ(1) / 6 and CI * U0 < QQ(1) / 6

# exterior constant B (third vectors only), rational mu as in pivd5_lemma.sage (seeded Powell starts)
import numpy as np
from scipy.optimize import minimize
LFn = np.array([[float(x) for x in g] for g in LF]); Nin = np.array(Ni, dtype=float)
Bs = []
for h in range(9):
    g2 = np.array([float(x) for x in 2 * LT[h]])
    obj = lambda m: (lambda v: v @ Nin @ v)(g2 + (np.abs(m) / np.abs(m).sum()) @ LFn)
    m = min((minimize(obj, np.random.default_rng(k).random(len(F)), method="Powell", options={"maxiter": 60000}) for k in range(4)), key=lambda r: r.fun).x
    mu = [QQ(Fr(abs(float(x))).limit_denominator(int(1000000))) for x in m]; sm = sum(mu); mu = [x / sm for x in mu]
    Bs.append(sqrt(RIF(dn2(2 * LT[h] + sum(mu[f] * LF[f] for f in range(len(F)))))).upper())
Bmax = RIF(max(Bs))
print("B <=", Bmax.upper())

# Gamma_kappa: lam in R^12_{>=0}, sum = 2, sum_third lam_s sigma_s s = 0
t0 = time.time()
ieqs = [[0] * (1 + i) + [1] + [0] * (11 - i) for i in range(12)]
eqns = [[-2] + [1] * 12] + [[0] + [sig[j] * T[j][col] for j in range(9)] + [0] * 3 for col in range(4)]
Pl = Polyhedron(ieqs=ieqs, eqns=eqns, base_ring=QQ)
Lall = LT + LI
Gk = Polyhedron(vertices=[sum(v[j] * Lall[j] for j in range(12)) for v in Pl.vertices_list()], base_ring=QQ)
Gm = Gk + Polyhedron(vertices=LF, base_ring=QQ)
assert Gm.dim() == 9
IE = Gm.inequalities_list()
assert all(ie[0] > 0 for ie in IE), "0 is not in the interior of Gamma: no first-order decrease"
r2 = min(ie[0]^2 / (vector(QQ, ie[1:]) * N * vector(QQ, ie[1:])) for ie in IE)
r = sqrt(RIF(r2))
print("Gamma_kappa: %d vertices; Gamma: %d vertices, %d facets (%.0fs); inradius^2 = %s, r = %s"
      % (Gk.n_vertices(), Gm.n_vertices(), len(IE), time.time() - t0, r2, r))

R = RIF; rh = R(rho); cc = R(c); a = R(3) / 4
psi2 = lambda tau: 2 * (1 + tau^2) / (1 - tau^2)^2
remt = rh / 2 + R(3) / 2 * rh * exp(sqrt(R(3) / 4) * rh)               # rem(t)/t at t = rho
U = cc * rh
e_third = psi2(R(1) / 3 + CT * U) * (CT * U)^2 / 2
e_int = (CI * U)^2 / (1 - (CI * U)^2)
interior = r - remt - 2 * max(e_third.upper(), e_int.upper()) / rh
exterior = 2 * a * cc - Bmax - remt - psi2(R(1) / 3) * cc^2 * rh
print("interior margin %s   exterior margin %s" % (interior, exterior))
assert interior > 0 and exterior > 0
print("LEMMA HOLDS: Phi_5(G) <= 2^16/(3^4 5^3) for d(G, G0) <= %s (equality only at G0)" % rho)
