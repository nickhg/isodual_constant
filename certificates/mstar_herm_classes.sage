# For each class in the Z-genus of M*, count Aut-conjugacy classes of g with g^2+g+1 = 0 (Z[omega]-structures),
# and record their Hermitian determinants, to see how many Hermitian lattices share M*'s Z-genus.
P = matrix(ZZ, [[18,-9,0,0],[-9,18,-9,-9],[0,-9,16,1],[0,-9,1,16]])
reps = Genus(P).representatives()
I = identity_matrix(ZZ, 4)
K.<w> = CyclotomicField(3)
def herm_det(Q, g):
    # Z[omega]-module via omega -> g (acting on columns); pick x with x, gx a Z-basis of a rank-2 Z[omega]-submodule
    # Hermitian form h with h(x,x) = Q(x)/2-type normalisation: h(x,y) = (2 b(x,y) - b(x,gy)(1-w^2)... use generic
    # construction: h(x,y) = sum over the group {1,g,g^2} b(x, g^k y) w^k  (sesquilinear, h(x,x) = 3/2 b(x,x) - ... )
    b = lambda x, y: x * Q * y
    h = lambda x, y: sum(b(x, g^k * y) * w^k for k in range(3))
    # Z[omega]-basis: try standard vectors
    for e1 in I.columns():
        for e2 in I.columns():
            M = matrix(ZZ, [e1, g*e1, e2, g*e2]).transpose()
            if abs(M.det()) == 1:
                H = matrix(K, 2, 2, [h(e1,e1), h(e1,e2), h(e2,e1), h(e2,e2)])
                return H.det()
    return None
for idx, R in enumerate(reps):
    Q = matrix(ZZ, R)
    A = [matrix(ZZ, a.matrix()) for a in QuadraticForm(ZZ, 2*Q).automorphism_group()]
    gs = [a for a in A if a*a + a + I == 0]
    # conjugacy classes under Aut
    classes = []
    for a in gs:
        if not any(any(c == u.inverse()*a*u for u in A) for c in classes): classes.append(a)
    isM = Q.is_similar(P) if hasattr(Q, 'is_similar') else None
    print("rep %d: |Aut| = %d, theta min = %s, #g = %d, Aut-classes of g = %d, dets = %s"
          % (idx, len(A), QuadraticForm(ZZ, 2*Q).theta_series(20).exponents()[1:3], len(gs), len(classes),
             [herm_det(Q, c) for c in classes]))
