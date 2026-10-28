# W(E_m) acts irreducibly on traceless symmetric forms (commutant of dimension 1, exact over QQ), m = 6, 7, 8;
# and for m = 8 the interval check of the radius 1/20 in Proposition e6e7.  sage schur_check.sage <m>

import sys
m = int(sys.argv[1]); C = CartanMatrix(['E', m]); G = matrix(QQ, C)
n = m
# simple reflections on coordinates in the root basis: s_i(x) = x - <x,a_i> a_i, as matrices acting on forms B -> S^T B S
refl = []
for i in range(n):
    S = identity_matrix(QQ, n)
    for j in range(n): S[i, j] -= G[i, j]
    refl.append(S)
for S in refl: assert S.transpose()*G*S == G
# basis of symmetric forms B with tr(G^{-1}B)=0
Gi = G.inverse()
sym = []
for i in range(n):
    for j in range(i, n):
        E = matrix(QQ, n, n); E[i, j] = 1; E[j, i] = 1; sym.append(E)
V = matrix(QQ, [list(E.list()) for E in sym])  # 36 x 64
tr = vector(QQ, [ (Gi*E).trace() for E in sym ])
K = matrix(QQ, [tr]).right_kernel_matrix()      # 35 x 36 coefficient vectors
basis = [sum(c*E for c, E in zip(row, sym)) for row in K.rows()]
Bm = matrix(QQ, [b.list() for b in basis]).transpose()   # 64 x 35
def act(S):
    cols = [ (S.transpose()*b*S).list() for b in basis ]
    M = matrix(QQ, cols).transpose()
    return Bm.solve_right(M)          # 35 x 35
R = [act(S) for S in refl]
D = n*(n+1)//2 - 1
# commutant: X with R_i X = X R_i
eqs = []
for Ri in R:
    eqs.append(Ri.tensor_product(identity_matrix(QQ, D)) - identity_matrix(QQ, D).tensor_product(Ri.transpose()))
A = block_matrix([[e] for e in eqs])
print("commutant dim:", A.right_kernel().dimension())
if n == 8:
  RIF = RealIntervalField(200); rho = RIF(1)/20
  lhs = RIF(7)/2*rho*exp(sqrt(RIF(7)/8)*rho); r = sqrt(RIF(42))/35
  print("lhs", lhs.upper(), "r", r.lower(), "ok", lhs < r)
