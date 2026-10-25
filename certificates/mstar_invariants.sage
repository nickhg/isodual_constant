# M_* as an object (WC's M4 and the reader's report): theta series, |Aut|, strong modularity, the isodual family through
# M_* and whether A2 (x) A2 can lie on the same family.
import itertools
P0 = matrix(ZZ, [[18,-9,0,0],[-9,18,-9,-9],[0,-9,16,1],[0,-9,1,16]])
K = matrix(ZZ, [[2,-1],[-1,2]]).tensor_product(matrix(ZZ, [[2,-1],[-1,2]]))      # A2 (x) A2, det 81
for nm, A in (("M*", P0), ("A2xA2", K)):
    Q = QuadraticForm(ZZ, 2 * A)                       # Q(x) = x^T A x
    th = Q.theta_series(40).list()
    print(nm, "det", A.det(), "theta (norm: count) ", [(i, th[i]) for i in range(len(th)) if th[i]][:12])
    G = pari(A).qfauto(); print("  |Aut| =", G[0], " (pari qfauto)")
# strong modularity (Quebbemann): for every exact divisor d of N = 135, sqrt(d) (L^# cap (1/d) L) is isometric to L
N = 135; Ls = P0.inverse()                             # Gram of L^# in the dual basis
def sublattice_gram(G, M):                             # Gram of the lattice spanned by rows of M (coefficients w.r.t. G's basis)
    return M * G * M.transpose()
V = QQ^4
Lmod = V.span(identity_matrix(QQ, 4).rows(), ZZ)                       # L, coordinates in L's basis
Ldual = V.span((P0.inverse()).rows(), ZZ)                              # L^#: rows of P0^{-1} in L's coordinates
for d in [d for d in divisors(N) if gcd(d, N // d) == 1]:
    I = Ldual.intersection(V.span((identity_matrix(QQ, 4) / d).rows(), ZZ))
    Bm = I.basis_matrix()
    S = d * Bm * P0.change_ring(QQ) * Bm.transpose()                   # sqrt(d) (L^# cap L/d)
    integral = all(x in ZZ for x in S.list())
    same = integral and QuadraticForm(ZZ, 2 * S.change_ring(ZZ)).is_globally_equivalent_to(QuadraticForm(ZZ, 2 * P0))
    print("  d =", d, ": det", S.det(), " isometric to M*:", same)
# isodual family through M*: tau(G) = U (135 G^{-1}) U^T with the U of Prop. mstarprops(4); tangent = 1-eigenspace of D tau
U = matrix(ZZ, [[0,0,-1,0],[0,0,0,1],[-1,1,0,-1],[0,-1,1,0]])
assert U * (135 * P0.inverse()) * U.transpose() == P0
Pi = P0.inverse()
basis = []
for i in range(4):
    for j in range(i, 4):
        E = matrix(QQ, 4, 4); E[i, j] = E[j, i] = 1; basis.append(E)
def vec(E): return vector(QQ, [E[i, j] for i in range(4) for j in range(i, 4)])
D = matrix(QQ, [vec(-135 * U * Pi * E * Pi * U.transpose()) for E in basis]).transpose()
k1 = (D - 1).right_kernel().dimension(); km1 = (D + 1).right_kernel().dimension()
print("isodual family through M*: dim of the 1-eigenspace of D tau =", k1, "(including scaling? scaling E = P0 has D tau(P0) =", "-P0" if D * vec(P0) == -vec(P0) else "?", ")")
print("  eigenvalues of D tau:", sorted(D.charpoly().roots(QQbar), key=lambda r: r[0].real()))
# congruence invariant of the isoduality: char poly of U^{-1} U^T
print("M*: charpoly(U^{-1} U^T) =", (U.inverse() * U.transpose()).charpoly())
# isodualities of A2 x A2: W with W (9 K^{-1}) W^T = K (isoduality scale det^{2/n} = 9)
W = QuadraticForm(ZZ, 2 * (9 * K.inverse()).change_ring(ZZ)).is_globally_equivalent_to(QuadraticForm(ZZ, 2 * K), return_matrix=True); ok = W is not False
print("A2xA2 isodual:", ok)
if ok:
    W = matrix(ZZ, W); W = W.transpose() if W * (9 * K.inverse()) * W.transpose() != K else W
    assert W * (9 * K.inverse()) * W.transpose() == K
    # all isodualities = W composed with Aut(K); collect char polys of W'^{-1} W'^T
    AutK = [matrix(ZZ, g) for g in pari(K).qfauto()[1]]
    from itertools import product
    gens = [matrix(ZZ, 4, 4, [int(x) for x in g.list()]) for g in AutK]
    # enumerate the group generated (|Aut| small)
    grp = {tuple(identity_matrix(ZZ, 4).list())}; frontier = [identity_matrix(ZZ, 4)]
    while frontier:
        new = []
        for a in frontier:
            for g in gens:
                b = a * g; t = tuple(b.list())
                if t not in grp: grp.add(t); new.append(b)
        frontier = new
    cps = set()
    for t in grp:
        a = matrix(ZZ, 4, 4, list(t)); V = a * W
        if V * (9 * K.inverse()) * V.transpose() != K: continue
        Vt = V.transpose(); cps.add(str((V.inverse() * Vt).charpoly()))
    print("A2xA2: |Aut| =", len(grp), "; charpolys of V^{-1} V^T over its isodualities:", cps)
