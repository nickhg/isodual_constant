r"""
Local lemma at pi_v D5: Phi_5(G) <= 6 whenever d(G, G0) = |log(G0^{-1/2} G G0^{-1/2})|_F <= rho (G, G0 of equal det).
Exact arithmetic (QQ) for the polytope; RIF for the final scalar inequalities.

Notation.  G0 = Gram of pi_v D5 (det 2, basis of pivd5_data.json), G = G0^{1/2} e^Lam G0^{1/2}, Lam traceless,
t = |Lam|_F, Y = G0^{1/2} Lam G0^{1/2} (so linear functionals of Lam are rational functionals of Y).
Phi~(G) = kappa_raw(G)^2 det(G)^{-1/2} eta~(G); log Phi~(G) - log 6 is bounded as follows (all bounds global):
  q_s(G) = s^T G s:     log q_s(G) <= log q_s(G0) + l_s(Y) + t^2/4      (Hoeffding: spread of spec(Lam) <= sqrt2 t)
  flags f (36):         log h_f(G) <= log h_f(G0) + L_f(Y) + 3 tr E     (log det concave; E = e^Lam - I - Lam >= 0,
                                                                         tr E <= t^2/2 exp(sqrt(3/4) t))
  kappa_raw(G) = sup_phi min_s k_s, k_s = q_s/(1 - ||s.phi||^2); 6 'half' and 7 'integer' active s at phi0.
  U = c t.  Exterior (some half h with dist(h.phi, 1/2 + Z) >= U):  k_h <= q_h / (3/4 + U - U^2), giving
     2 log(kappa/2) + eta part <= B t + t^2/2 + 3/2 t^2 e^{..} - 2 log(1 + 4/3 (U - U^2)),
     B >= max_{|Lam|=1} [2 l_h + min_f L_f]  (any mu in the simplex: B_h <= |2 l_h + sum mu_f L_f|_*).
  Interior (all |h.Delta| < U, phi = phi0 + Delta, a single component since the half vectors generate Z^4):
     half: psi(1/2 - u) <= log(4/3) - theta (4/3) u,  theta = 1 - (5/3) c rho
     int s = +-h_a +-h_b:  psi(|s.Delta|) <= eps (u_a + u_b),  eps = 2 c rho / (1 - 4 c^2 rho^2)
     so 2 log(kappa/2) + eta part <= D(Y) + t^2/2 + 3/2 t^2 e^{..},  D(Y) = min_{g in Gamma} g(Y) <= -r t,
     Gamma = Gamma_kappa(theta, eps) + conv(L_f) (LP duality), r = inradius of Gamma in the dual norm (exact facets).
Data: pivd5_data.json, written by pivd5_sharp.py --dump (every property used is re-verified here exactly).
Usage: sage pivd5_lemma.sage rho c      (rationals, e.g. 1/20 11/20; about 30 min)
"""
import json, sys, time, os
D = json.load(open(os.path.join(os.path.dirname(os.path.abspath(sys.argv[0])), "pivd5_data.json")))
rho = QQ(sys.argv[1]); c = QQ(sys.argv[2])
B5 = matrix(ZZ, D["G5"]); G5 = B5 * B5.T
# B5 is a basis of D5 = {x in Z^5 : sum x even}: integral rows with even coordinate sums span a sublattice of D5,
# and [Z^5 : D5] = 2 = |det B5| makes it all of D5.  G0 below is pi_v D5 for v = the first row (a root).
assert B5.nrows() == 5 and B5.ncols() == 5 and all(sum(row) % 2 == 0 for row in B5.rows()) and abs(B5.det()) == 2
assert B5[0] * B5[0] == 2
G0 = G5[1:, 1:] - matrix(QQ, 4, 1, G5[0, 1:].list()) * matrix(QQ, 1, 4, G5[0, 1:].list()) / G5[0, 0]
G0i = G0.inverse()
H = [vector(ZZ, s) for s in D["half"]]; I = [vector(ZZ, s) for s in D["int"]]; F = [[vector(ZZ, v) for v in f] for f in D["flags"]]
print("G0 =", G0.list(), " det", G0.det())
q = lambda s: s * G0 * s
assert all(q(h) == 3/2 for h in H) and all(q(s) == 2 for s in I)
phi0 = vector(QQ, [0, 1/2, 0, 0])
assert all(h * phi0 - 1/2 in ZZ for h in H) and all(s * phi0 in ZZ for s in I)
print("k_s(G0, phi0) = 2 for the 6 half and 7 integer vectors")
# The half vectors generate Z^4: integer left inverse K (4x6), K H = 1, where H is the 6x4 matrix of rows h_a.
# Single interior component.  Let phi be interior: h_a.phi in 1/2 + Z + (-U, U) for all a.  With Delta = phi - phi0,
# H Delta = n + e, n in Z^6, |e|_inf < U.  Then Delta = K H Delta = K n + K e with K n in Z^4, so phi = phi0 + Delta'
# mod Z^4, Delta' = K e.  Now n + e = H Delta = H K n + H K e, i.e. n - H K n = H K e - e: the left side is
# integral and the right side has |.|_inf <= (|HK|_inf + 1) U < 1, so both vanish and H Delta' = H K e = e.
# Hence every interior phi is phi0 + Delta' mod Z^4 with |h_a.Delta'| < U for all a (the case analysed below),
# and Delta' is unique (H has full column rank).  The script asserts (|HK|_inf + 1) U < 1 at the end.
HM = matrix(ZZ, [list(h) for h in H])
Sm, Lm, Rm = HM.smith_form(); assert all(Sm[i, i] == 1 for i in range(4))
K = Rm * Sm.T * Lm; assert K * HM == 1
nHK = max(sum(abs(x) for x in row) for row in (HM * K).rows())
print("half vectors generate Z^4; |HK|_inf =", nHK)
pairs = []
for s in I:
    pr = [(a, b) for a in range(6) for b in range(a + 1, 6) for sa in (1, -1) for sb in (1, -1) if sa * H[a] + sb * H[b] == s]
    assert pr; pairs.append(pr[0])
print("integer vectors as +-h_a +-h_b:", pairs)
# flags: primitive prefixes (3 vectors extend to a basis), equal base value
for f in F:
    assert len(f) == 3, "flags must have 3 vectors (the minors(3) primitivity test assumes it)"
    assert gcd((matrix(ZZ, f)).minors(3)) == 1
hv = set(prod((matrix(QQ, f[:k]) * G0 * matrix(QQ, f[:k]).T).det() for k in range(1, 4)) for f in F)
assert len(hv) == 1; h0 = hv.pop()
base2 = 2^4 / G0.det() * h0^2 / G0.det()^3
print("%d flags, prod_k det G_k = %s; Phi~(G0)^2 = %s" % (len(F), h0, base2)); assert base2 == 36
# traceless Y basis, norm Gram N (|Lam|_F^2 = y^T N y), functionals
Ys = []
for i in range(4):
    for j in range(i, 4):
        E = matrix(QQ, 4, 4); E[i, j] = 1; E[j, i] = 1
        Ys.append(E - (G0i * E).trace() / 4 * G0)
Ys = [Ys[p] for p in matrix(QQ, [Y.list() for Y in Ys]).pivot_rows()]; assert len(Ys) == 9
N = matrix(QQ, 9, 9, lambda a, b: (G0i * Ys[a] * G0i * Ys[b]).trace()); Ni = N.inverse()
fun = lambda W: vector(QQ, [(W * Y).trace() for Y in Ys])
ell = lambda s: fun(matrix(QQ, 4, 1, list(s)) * matrix(QQ, 1, 4, list(s)) / q(s))
LH = [ell(h) for h in H]; LI = [ell(s) for s in I]
def Lflag(f):
    W = matrix(QQ, 4, 4)
    for k in range(1, 4):
        U = matrix(QQ, f[:k]); W += U.T * (U * G0 * U.T).inverse() * U
    return fun(W)
LF = [Lflag(f) for f in F]
dn2 = lambda g: g * Ni * g

theta = 1 - 5 / 3 * c * rho; eps = 2 * c * rho / (1 - 4 * c^2 * rho^2)
print("rho = %s, c = %s, theta = %s, eps = %s" % (rho, c, theta, eps))
# exterior constant B: B_h <= |2 l_h + sum mu_f L_f|_* for a rational mu (simplex), from a float QP.  The QP only
# picks mu (any mu in the simplex is valid); its starts are seeded (default_rng(k), k < 4), so B is reproducible.
import numpy as np
from scipy.optimize import minimize
LFn = np.array([[float(x) for x in g] for g in LF]); Nin = np.array(Ni, dtype=float)
Bs = []
for h in range(6):
    g2 = np.array([float(x) for x in 2 * LH[h]])
    obj = lambda m: (lambda v: v @ Nin @ v)(g2 + (np.abs(m) / np.abs(m).sum()) @ LFn)
    m = min((minimize(obj, np.random.default_rng(k).random(36), method="Powell", options={"maxiter": 40000}) for k in range(4)), key=lambda r: r.fun).x
    from fractions import Fraction
    mu = [QQ(Fraction(abs(float(x))).limit_denominator(int(1000000))) for x in m]; sm = sum(mu); mu = [x / sm for x in mu]
    Bh = sqrt(RIF(dn2(2 * LH[h] + sum(mu[f] * LF[f] for f in range(len(F))))))
    Bs.append(Bh.upper())
Bmax = RIF(max(Bs))
print("B <=", Bmax.upper())

t0 = time.time()
nvar = 19; ieqs = []; eqns = []
for t_ in range(13):
    r_ = [0] * (nvar + 1); r_[1 + t_] = 1; ieqs.append(r_)
for h in range(6):
    p = [0] * (nvar + 1); p[1 + 7 + h] = theta * 4 / 3
    for si, (a, b) in enumerate(pairs):
        if h in (a, b): p[1 + si] -= eps
    for sg in (1, -1):
        r_ = list(p); r_[1 + 13 + h] = -sg; ieqs.append(r_)
eqns.append([-2] + [1] * 13 + [0] * 6)
for col in range(4):
    eqns.append([0] * 14 + [H[h][col] for h in range(6)])
Pl = Polyhedron(ieqs=ieqs, eqns=eqns, base_ring=QQ)
Lall = LI + LH
Gk = Polyhedron(vertices=[sum(v[t_] * Lall[t_] for t_ in range(13)) for v in Pl.vertices_list()], base_ring=QQ)
Ge = Polyhedron(vertices=LF, base_ring=QQ)
Gm = Gk + Ge
assert Gm.dim() == 9
IE = Gm.inequalities_list()      # b + a.y >= 0; 0 strictly interior <=> b > 0 for every facet
assert all(ie[0] > 0 for ie in IE), "0 is not in the interior of Gamma: no first-order decrease"
r2 = min(ie[0]^2 / (vector(QQ, ie[1:]) * N * vector(QQ, ie[1:])) for ie in IE)
r = sqrt(RIF(r2))
print("Gamma: %d vertices, %d facets (%.0fs); inradius^2 = %s, r = %s" % (Gm.n_vertices(), Gm.n_inequalities(), time.time() - t0, r2, r))

R = RIF; rh = R(rho); cc = R(c)
rem = rh / 2 + R(3) / 2 * rh * exp(sqrt(R(3) / 4) * rh)            # (t/2 + 3/2 t e^{sqrt(3/4) t}) at t = rho
interior = r - rem
exterior = R(8) / 3 * cc * (1 - cc * rh) - R(16) / 9 * cc^2 * rh - (Bmax + rem)
uniq = (nHK + 1) * cc * rh
print("interior margin r - rem = %s   exterior margin = %s   (|HK|+1) U = %s" % (interior, exterior, uniq))
assert interior > 0 and exterior > 0 and uniq < 1 and 2 * cc * rh <= R(1) / 2
print("LEMMA HOLDS: Phi_5(G) <= 6 for d(G, G0) <= %s (equality only at G0)" % rho)
