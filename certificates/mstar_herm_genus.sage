# Hermitian structures in the Z-genus of M*: (1) Aut-classes of g with g^2+g+1=0, with and without identifying g ~ g^2;
# (2) Hermitian Gram matrix over Z[w] and its Z-trace form, to compare local invariants at 3 (ramified) and 5 (inert).
P = matrix(ZZ, [[18,-9,0,0],[-9,18,-9,-9],[0,-9,16,1],[0,-9,1,16]])
reps = Genus(P).representatives(); I = identity_matrix(ZZ, 4)
K.<w> = CyclotomicField(3)
for idx, R in enumerate(reps):
    Q = matrix(ZZ, R)
    A = [matrix(ZZ, a.matrix()) for a in QuadraticForm(ZZ, 2*Q).automorphism_group()]
    gs = [a for a in A if a*a + a + I == 0]
    if not gs: continue
    cls = []
    for a in gs:
        if not any(any(c == u.inverse()*a*u for u in A) for c in cls): cls.append(a)
    g = cls[0]; conj_same = any(g*g == u.inverse()*g*u for u in A)
    # h(x,y) = (1/3) sum_k b(x, g^k y) w^{-k}... choose h with Tr_{E/Q} h = b  (so h(x,x) = b(x,x)/2)
    b = lambda x, y: x * Q * y
    h = lambda x, y: sum(b(g^k * x, y) * w^k for k in range(3)) / 3
    # Z[w]-basis e, f with e, ge, f, gf a Z-basis
    basis = None
    for e in I.columns():
        for f in I.columns():
            if abs(matrix(ZZ, [e, g*e, f, g*f]).det()) == 1: basis = (e, f); break
        if basis: break
    e, f = basis
    H = matrix(K, 2, 2, [h(e, e), h(e, f), h(f, e), h(f, f)])
    # check Tr h = b on the Z-basis
    Zb = [e, g*e, f, g*f]
    T = matrix(QQ, 4, 4, lambda i, j: (lambda z: z.trace())(sum(b(g^k*Zb[i], Zb[j]) * w^k for k in range(3)) / 3))
    print("rep %d: min %d, #g = %d, Aut-classes = %d, g ~ g^2 in Aut: %s" % (idx, QuadraticForm(ZZ, 2*Q).theta_series(40).exponents()[1], len(gs), len(cls), conj_same))
    print("   H =", H.list(), " det H =", H.det(), " Tr-form det =", T.det(), " Tr-form = Gram in basis? ", T == matrix(ZZ, Zb).transpose().transpose()*0 + matrix(QQ, 4, 4, lambda i, j: b(Zb[i], Zb[j])))
    print("   Z-genus of Tr-form at 3, 5:", Genus(T.change_ring(ZZ)).local_symbol(3), "|", Genus(T.change_ring(ZZ)).local_symbol(5))
