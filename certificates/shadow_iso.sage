# Exact identifications for the landscape table: pi_v A4 = (Z^3, I + J/2) and pi_u(pi_v E6) = A2 (x) A2 (u shortest),
# as isometries of the primitive integral forms (up to scaling).
def proj(G, u):
    M = matrix(ZZ, [u]); D, U, V = M.smith_form(); B = V.inverse()
    H = B * G * B.transpose(); n = G.nrows()
    return matrix(QQ, n-1, n-1, lambda i, j: H[i+1, j+1] - H[0, i+1]*H[0, j+1]/H[0, 0])
def prim(G):
    d = lcm([x.denominator() for x in G.list()]); H = (d*G).change_ring(ZZ); return H / gcd(H.list())
def iso(A, B): return QuadraticForm(ZZ, 2*prim(A)).is_globally_equivalent_to(QuadraticForm(ZZ, 2*prim(B)))
A4 = matrix(QQ, 4, 4, lambda i, j: 2 if i == j else (-1 if abs(i-j) == 1 else 0))
print("pi_v A4 ~ I + J/2:", iso(proj(A4, vector(ZZ, [1,0,0,0])), identity_matrix(QQ, 3) + matrix(QQ, 3, 3, [1/2]*9)))
E6 = matrix(QQ, [[2,-1,0,0,0,0],[-1,2,-1,0,0,0],[0,-1,2,-1,0,-1],[0,0,-1,2,-1,0],[0,0,0,-1,2,0],[0,0,-1,0,0,2]])
P = proj(E6, vector(ZZ, [1,0,0,0,0,0]))
import itertools
vs = [vector(ZZ, x) for x in itertools.product(range(-2, 3), repeat=5) if any(x)]
m = min(v*P*v for v in vs); us = [v for v in vs if v*P*v == m]
A2 = matrix(QQ, [[2,-1],[-1,2]])
print("pi_u(pi_v E6) ~ A2 x A2 for all", len(us), "shortest u:", all(iso(proj(P, u), A2.tensor_product(A2)) for u in us))
