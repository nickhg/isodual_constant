# The kappa-maximiser of M* is unique up to sign: exact check (Web Claude's argument, phi65_final_reply_for_cc.md).
# Unscaled Gram G of M*: minimal norm 8/9 (9 pairs), next norm 1, and kappa_raw(G) = 1 at phi0 = (0,0,1/3,1/3).
# Step 1: min_s k_s(phi) >= 1 iff ||phi(s)|| >= 1/3 for the 9 minimal pairs (every other s has q_s >= 1, so k_s >= 1).
# Step 2 (here): the set {phi in R^4/Z^4 : ||s.phi|| >= 1/3 for all 9 minimal s} is exactly {phi0, -phi0}.
# Method: pick 4 minimal vectors forming a Z-basis, so phi <-> a = (s_i.phi mod 1) is a bijection of tori; then
# a_i in [1/3, 2/3], and each remaining minimal s = sum c_i s_i needs sum c_i a_i in [n + 1/3, n + 2/3] for an integer
# n.  Enumerate the n's over exact rational polyhedra (QQ) and list every feasible region.
import itertools
G = matrix(QQ, [[1, -1/2, 0, 0], [-1/2, 1, -1/2, -1/2], [0, -1/2, 8/9, 1/18], [0, -1/2, 1/18, 8/9]])
# short vectors: |v_i|^2 <= q (G^-1)_ii < 68/15 for q < 1, so |v_i| <= 2 and the box R = 3 contains them all
R = 3
vecs = {}
for v in itertools.product(range(-R, R + 1), repeat=4):
    v = vector(QQ, v)
    if v == 0 or next(x for x in v if x != 0) < 0: continue
    q = v * G * v
    if q < 1: vecs.setdefault(q, []).append(v)
mins = vecs[min(vecs)]
print("norms below 1:", {k: len(w) for k, w in vecs.items()}, " minimal pairs:", len(mins))
assert min(vecs) == 8/9 and len(mins) == 9 and len(vecs) == 1
# a Z-basis of minimal vectors
B = next(matrix(QQ, c) for c in itertools.combinations(mins, 4) if abs(matrix(QQ, c).det()) == 1)
others = [s for s in mins if s not in B.rows()]
C = [B.solve_left(s) for s in others]                 # s = c * B
assert all(c in ZZ^4 for c in C)
print("basis:", [list(r) for r in B.rows()]); print("others in that basis:", [list(c) for c in C])
# a in [1/3, 2/3]^4; for each other s, c.a ranges over an interval: enumerate integer shifts n with c.a in [n+1/3, n+2/3]
third = QQ(1) / 3
cube = [[-third] + [1 if j == i else 0 for j in range(4)] for i in range(4)] + \
       [[2 * third] + [-1 if j == i else 0 for j in range(4)] for i in range(4)]
regions = [Polyhedron(ieqs=cube, base_ring=QQ)]
for c in C:
    lo = sum(min(x * third, x * 2 * third) for x in c); hi = sum(max(x * third, x * 2 * third) for x in c)
    new = []
    for P in regions:
        for n in range(floor(lo) - 1, ceil(hi) + 1):
            Q = P.intersection(Polyhedron(ieqs=[[-(n + third)] + list(c), [n + 2 * third] + [-x for x in c]],
                                          base_ring=QQ))
            if not Q.is_empty(): new.append(Q)
    regions = new
print("feasible regions:", len(regions))
sols = set()
for P in regions:
    assert P.dim() == 0, P                                  # each region is a single point
    a = vector(QQ, P.vertices()[0])
    phi = B.solve_right(a)                                  # B * phi = a
    sols.add(tuple(x - floor(x) for x in phi))
print("phi mod Z^4:", sols)
phi0 = (0, 0, third, third); mphi0 = (0, 0, 2 * third, 2 * third)
assert sols == {phi0, mphi0}
print("UNIQUE: the kappa-maximisers of M* are exactly +-phi0 = +-(0,0,1/3,1/3) mod Z^4")
