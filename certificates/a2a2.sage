# Exact values at A2 (x) A2 (the shadow of pi_v E6 for Phi_5): eta by exact Psi-minimisation; kappa = 3/2 at det 1
# (proof in the paper: the 9 minimal pairs alpha (x) beta force phi = +-(1/3,1/3,1/3,1/3) on the basis alpha_i (x) beta_j),
# checked here by evaluating min_s k_s at that phi exactly and on a grid of torsion points.
load("psi.sage")
import itertools
A2 = matrix(QQ, [[2, -1], [-1, 2]]); K = A2.tensor_product(A2)
d = K.det(); p2, _ = psi2_min(K); eta = p2 / d^(5/2)
print("det", d, "eta(A2 x A2) =", eta)
vs = [vector(ZZ, x) for x in itertools.product(range(-3, 4), repeat=4) if any(x)]
m = min(v * K * v for v in vs); mins = [v for v in vs if v * K * v == m]
print("min", m, "minimal vectors", len(mins), "next", min(v * K * v for v in vs if v * K * v > m))
def kval(phi):
    return min((v * K * v) / d^(1/4) / (1 - min(frac(v * phi), 1 - frac(v * phi))^2) for v in vs if v * K * v <= 3 * m)
phi0 = vector(QQ, [1/3] * 4)
k0 = kval(phi0); print("min_s k_s(phi0) =", k0, "= 3/2 * ... check:", bool(k0 == 3/2))
print("Phi_5 =", k0^2 * eta)
best = max(kval(vector(QQ, a) / 6) for a in itertools.product(range(6), repeat=4))
print("max over 6-torsion points:", best)
