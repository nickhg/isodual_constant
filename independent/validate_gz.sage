# Independent validation of the parametrisation used by stage 1 (x -> S = mu R mu^T) and stage 2 (z -> G(z)), from first
# principles: build an explicit basis and compare Gram matrices, at random points, in 200-bit arithmetic.
#   stage 1: a basis b_1..b_4 with Gram-Schmidt norms r_k = e^{y_{k+1} - y_k} (y_0 = y_4 = 0) and coefficients mu_ij
#            has Gram matrix S; check S == mu R mu^T and that G_k = e^{y_k} are its prefix determinants.
#   stage 2: v = sqrt(r) e_0 and b_i' = r^{-1/8} (0, b_i) + phi_i v in R^5; check Gram(v, b_1', .., b_4') == G(z) and det = 1.
import random
R200 = RealField(200); random.seed(int(11)); worst1 = worst2 = 0
for trial in range(int(100000)):
    y = [R200(random.uniform(-0.6, 0.6)) for _ in range(3)]; Y = [R200(0)] + y + [R200(0)]
    mu = matrix(R200, 4, 4, 0)
    for i in range(4):
        mu[i, i] = 1
        for j in range(i): mu[i, j] = R200(random.uniform(-0.5, 0.5))
    rk = [exp(Y[k + 1] - Y[k]) for k in range(4)]
    # explicit basis: b_i = sum_j mu_ij sqrt(r_j) e_j  (orthonormal e_j), so its Gram-Schmidt norms are r_j
    B = matrix(R200, 4, 4, lambda i, j: mu[i, j] * sqrt(rk[j]))
    S_explicit = B * B.transpose()
    S_formula = mu * diagonal_matrix(R200, [exp(y[0]), exp(y[1] - y[0]), exp(y[2] - y[1]), exp(-y[2])]) * mu.transpose()
    worst1 = max(worst1, max(abs(x) for x in (S_explicit - S_formula).list()))
    for k in range(1, 4): worst1 = max(worst1, abs(S_explicit[:k, :k].det() - exp(Y[k])))
    r = R200(random.uniform(1.3, 1.6)); phi = vector(R200, [random.uniform(-0.5, 0.5) for _ in range(4)])
    v = vector(R200, [sqrt(r), 0, 0, 0, 0])
    rows = [v] + [r ** R200(-1/8) * vector(R200, [0] + list(B.row(i))) + phi[i] * v for i in range(4)]
    G_explicit = matrix(R200, rows) * matrix(R200, rows).transpose()
    G_formula = matrix(R200, 5, 5, 0); G_formula[0, 0] = r
    for i in range(4):
        G_formula[0, i + 1] = G_formula[i + 1, 0] = r * phi[i]
        for j in range(4): G_formula[i + 1, j + 1] = r ** R200(-1/4) * S_formula[i, j] + r * phi[i] * phi[j]
    worst2 = max(worst2, max(abs(x) for x in (G_explicit - G_formula).list()), abs(G_formula.det() - 1))
print("100000 random points: max |S_explicit - mu R mu^T| (and prefix dets) = %.3e;  max |Gram(v, b') - G(z)| (and |det - 1|) = %.3e" % (worst1, worst2))
