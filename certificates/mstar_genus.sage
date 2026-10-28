# Genus of M* (P0) and its class number; the Z[omega]-Hermitian form attached to an order-3 automorphism.
P = matrix(ZZ, [[18,-9,0,0],[-9,18,-9,-9],[0,-9,16,1],[0,-9,1,16]])
# scale: P is 9*? check content and the scaled-down even form
print("gcd of entries:", gcd(P.list()), " det:", P.det().factor())
g = matrix(ZZ, [[0,-1,1,0],[1,-1,1,-1],[0,0,0,-1],[0,0,1,-1]])
assert g*g + g + 1 == 0
# which side does g act on? check g^T P g == P
print("g^T P g == P:", g.transpose()*P*g == P, "  g P g^T == P:", g*P*g.transpose() == P)
G = Genus(P)
print("genus symbol:", G)
reps = G.representatives()
print("class number of the genus:", len(reps))
