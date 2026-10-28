# Is M* a Hermitian lattice over Z[omega]?  Look for g in Aut(P0) with g^2 + g + 1 = 0 (order 3, no fixed vector).
P = matrix(ZZ, [[18,-9,0,0],[-9,18,-9,-9],[0,-9,16,1],[0,-9,1,16]])
A = QuadraticForm(ZZ, 2*P).automorphism_group()
G = [matrix(ZZ, g.matrix()) for g in A]
print("|Aut| =", len(G))
I = identity_matrix(ZZ, 4)
E = [g for g in G if g*g + g + I == 0]
print("elements with g^2+g+1=0:", len(E))
for g in G: assert g.transpose()*P*g == P or g*P*g.transpose() == P
if E:
    g = E[0]; print(g)
    # Hermitian form: h(x,y) = <x,y> - omega-twisted part; rank-2 over Z[omega] since g has no fixed vector
    print("fixed space of g:", (g - I).kernel().dimension())
print("theta series:", QuadraticForm(ZZ, 2*P).theta_series(40))
