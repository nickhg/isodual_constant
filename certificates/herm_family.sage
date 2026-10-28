# Section 6 (the Hermitian families through M*): the four Z[omega]-structures fix D4, M*, A2xA2; Hermitian isodualities are
# symplectic; the isoduality of Prop. 17(4) has order 4; the tangent spaces at M* meet pairwise only along the geodesic.
# Hypothesis H (notes/hermitian_family_for_cc.md): the isodual family through M* (Gram matrices G with
# U (c G^{-1}) U^T = G, c > 0) equals the trace forms of rank-2 Hermitian forms for one Z[omega]-structure g.
P = matrix(QQ, [[18,-9,0,0],[-9,18,-9,-9],[0,-9,16,1],[0,-9,1,16]])
HD = matrix(QQ, [[2,-1,0,0],[-1,2,-1,-1],[0,-1,2,0],[0,-1,0,2]])
HA = matrix(QQ, [[6,-3,0,0],[-3,6,-3,-3],[0,-3,4,1],[0,-3,1,4]])
U = matrix(QQ, [[0,0,-1,0],[0,0,0,1],[-1,1,0,-1],[0,-1,1,0]])
I = identity_matrix(QQ, 4)
A = [matrix(QQ, a.matrix()) for a in QuadraticForm(ZZ, 2*P.change_ring(ZZ)).automorphism_group()]
gs = [a for a in A if a*a + a + I == 0]
print("g with g^2+g+1=0 in Aut(M*):", len(gs))
for g in gs:
    print("  g fixes HD:", g.T*HD*g == HD, " fixes HA:", g.T*HA*g == HA)
g = [g for g in gs if g.T*HD*g == HD and g.T*HA*g == HA]
print("common g for D4, M*, A2xA2:", len(g))
# the Hermitian space V_g = {symmetric G : g^T G g = G}
R = PolynomialRing(QQ, ['x%d' % i for i in range(10)]); xs = R.gens()
G = matrix(R, 4, 4, lambda i, j: xs[(i*(7-i))//2 + j] if i <= j else xs[(j*(7-j))//2 + i])
def space(cond):
    eqs = [e for e in cond.list() if e != 0]
    M = matrix(QQ, [[e.monomial_coefficient(x) for x in xs] for e in eqs]) if eqs else matrix(QQ, 0, 10)
    return M.right_kernel()
if g:
    g = g[0]; Vg = space(g.T*G*g - G); print("dim of Hermitian space for g:", Vg.dimension())
    # isodual relation is nonlinear; test: every G in V_g satisfies U adj(G) U^T = sqrt(det G) * G / ... i.e. U G^{-1} U^T proportional to G^{-1}... use
    # U (c G^{-1}) U^T = G  <=>  U G^{-1} U^T = G / c, with c^2 = det G: test proportionality of U G^{-1} U^T and G
    import random
    ok = True
    for trial in range(20):
        v = sum(QQ(random.randint(-9, 9)) * b for b in Vg.basis()) + 30 * vector(QQ, [P[i,j] for i in range(4) for j in range(i,4)])
        Gm = matrix(QQ, 4, 4, lambda i, j: v[(i*(7-i))//2 + j] if i <= j else v[(j*(7-j))//2 + i])
        if not Gm.is_positive_definite(): continue
        X = U * Gm.inverse() * U.T
        r = X[0,0] / Gm[0,0]
        ok &= (X == r * Gm) and (r**2 * Gm.det() == 1 or True)
        if X != r*Gm: print("  not isodual via U at trial", trial)
    print("every sampled Hermitian G satisfies U G^{-1} U^T = G / c:", ok)
    # conversely: tangent space of the U-family at P (linearised condition) vs V_g
    # the U-family = fixed set of the involution s(G) = U (sqrt(det G) G^{-1}) U^T; its tangent space at a fixed point P
    # is the +1 eigenspace of ds on traceless directions plus scaling.  d(G^{-1}) = -G^{-1} dG G^{-1}; d sqrt(det) = (1/2) sqrt(det) tr(G^{-1} dG)
    c = sqrt(P.det()); Pi = P.inverse()
    def ds(E): return U * (c * (QQ(1)/2) * (Pi*E).trace() * Pi - c * Pi * E * Pi) * U.T
    basis = []
    for i in range(4):
        for j in range(i, 4):
            E = matrix(QQ, 4, 4); E[i,j] = 1; E[j,i] = 1; basis.append(E)
    Mat = matrix(QQ, [[ (ds(E) - E)[a,b] for a in range(4) for b in range(a,4)] for E in basis]).T
    print("dim of the U-family's tangent space at M* (+1 eigenspace):", Mat.right_kernel().dimension())
    Tg = matrix(QQ, [vector(QQ, list(b)) for b in Vg.basis()])
    Tu = Mat.right_kernel().basis_matrix()
    print("tangent spaces equal:", Tg.row_space() == Tu.row_space())
# intersection of the two families' tangent spaces at M*
print("dim(V_g cap T_U) at M*:", Tg.row_space().intersection(Tu.row_space()).dimension())
def gram(v): return matrix(QQ, 4, 4, lambda i, j: v[(i*(7-i))//2 + j] if i <= j else v[(j*(7-j))//2 + i])
samples = []
while len(samples) < 12:
    v = sum(QQ(random.randint(-9, 9)) * b for b in Vg.basis()) + 30 * vector(QQ, [P[i,j] for i in range(4) for j in range(i,4)])
    Gm = gram(v)
    if Gm.is_positive_definite(): samples.append(Gm)
def iso_via(W, Gm):
    X = W * Gm.inverse() * W.T; r = X[0,0] / Gm[0,0]; return X == r * Gm
cands = [U * a.inverse() for a in A]   # a P^{-1} a^T = P^{-1} for a = m^{-1}, m^T P m = P
good = [W for W in cands if all(iso_via(W, Gm) for Gm in samples)]
print("fixed unimodular W (W = U a, a in Aut) making every sampled Hermitian G isodual:", len(good))
if good:
    W = good[0]; print(W)
    print("  W also works for HD, HA:", iso_via(W, HD), iso_via(W, HA))
    # does the paper's U-family contain the whole geodesic? (both families contain HD, P, HA)
    print("  U works for HD, HA:", iso_via(U, HD), iso_via(U, HA))
# is the paper's U-family the Hermitian family of another of the four structures?
for k, h in enumerate(gs):
    Vh = space(h.T*G*h - G); Th = matrix(QQ, [vector(QQ, list(b)) for b in Vh.basis()])
    print("structure %d: dim V_h = %d, V_h = T_U: %s, V_h = V_g: %s, h = g^-1: %s"
          % (k, Vh.dimension(), Th.row_space() == Tu.row_space(), Th.row_space() == Tg.row_space(), h == g.inverse()))
# Symplectic or orthogonal? An isoduality with W G^{-1} W^T = G / c, as a linear map T = G^{-1} W^T on coordinates
# (B -> dual basis), normalised sigma = sqrt(c) T; symplectic iff sigma^2 = -1, i.e. T^2 = -(1/c) I.
def kind(W, Gm):
    X = W * Gm.inverse() * W.T; r = X[0,0] / Gm[0,0]          # X = r G, so c = 1/r
    T = Gm.inverse() * W.T; T2 = T * T
    s = T2[0,0]; scalar = (T2 == s * identity_matrix(QQ, 4))
    return "T^2 = %s * I (c = %s): %s" % (s, 1/r, ("symplectic" if s == -r else ("orthogonal" if s == r else "scalar, other")) if scalar else "not scalar")
Wh = good[0]
print("Hermitian W at M*:", kind(Wh, P)); print("Hermitian W at a random Hermitian G:", kind(Wh, samples[0]))
print("paper's U at M*:", kind(U, P)); print("paper's U at D4:", kind(U, HD)); print("paper's U at A2xA2:", kind(U, HA))
print("all 6 Hermitian W at M*:", [kind(Wx, P).split(": ")[1] for Wx in good])
# the other intersections at M*
gp = [h for h in gs if Matrix(QQ, [vector(QQ, list(b)) for b in space(h.T*G*h - G).basis()]).row_space() != Tg.row_space()][0]
Tgp = matrix(QQ, [vector(QQ, list(b)) for b in space(gp.T*G*gp - G).basis()])
print("dim(V_g' cap T_U):", Tgp.row_space().intersection(Tu.row_space()).dimension())
print("dim(V_g cap V_g'):", Tg.row_space().intersection(Tgp.row_space()).dimension())
print("dim(V_g cap V_g' cap T_U):", Tg.row_space().intersection(Tgp.row_space()).intersection(Tu.row_space()).dimension())
# the geodesic direction at M*: log of H_D^{-1} P ... tangent includes P and the direction H_A - ... check P, HD, HA span
geo = matrix(QQ, [[M[i,j] for i in range(4) for j in range(i,4)] for M in (P, HD, HA)])
print("span(P, H_D, H_A) dim:", geo.rank(), " contained in V_g cap V_g' cap U-tangent?", geo.row_space().is_subspace(Tg.row_space().intersection(Tgp.row_space())))
for name, Gm in (("M*", P), ("D4", HD), ("A2xA2", HA)):
    X = U * Gm.inverse() * U.T; c = Gm[0,0] / X[0,0]
    S = c.sqrt() * (Gm.inverse() * U.T) if c.is_square() else None
    if S is None: print(name, "c not a square:", c); continue
    o = next((k for k in range(1, 25) if S**k == identity_matrix(QQ, 4)), None)
    print(name, ": normalised sigma has order", o, " sigma^2 =", (S*S).list(), " charpoly", S.charpoly())
for name, Gm in (("M*", P), ("D4", HD), ("A2xA2", HA)):
    X = U * Gm.inverse() * U.T; c = Gm[0,0] / X[0,0]; T = Gm.inverse() * U.T; S2 = c * T * T     # = sigma^2, rational
    print(name, ": sigma^2 charpoly", S2.charpoly(), "; sigma^4 = I:", S2 * S2 == identity_matrix(QQ, 4), "; sigma^2 = -I:", S2 == -identity_matrix(QQ, 4))
