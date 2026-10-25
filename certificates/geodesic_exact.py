"""M_* lies exactly on the geodesic from D4 to A2 (x) A2 in the isodual family through it: with H_D, H_A (the copies found by
family_exact.py, in the coordinates of P0), X = H_D^{-1} P0 and Y = H_D^{-1} H_A commute, and their det-normalised eigenvalues
are lambda = sqrt2, 1/sqrt2 (each twice) and mu = sqrt(6/5), sqrt(5/6) (each twice), so mu = lambda^t with t = log_2(6/5):
d(D4, M_*) = ln(6/5), d(M_*, A2 x A2) = ln(5/3), d(D4, A2 x A2) = ln 2."""
import sympy as sp
HD = sp.Matrix([[2,-1,0,0],[-1,2,-1,-1],[0,-1,2,0],[0,-1,0,2]])
HA = sp.Matrix([[6,-3,0,0],[-3,6,-3,-3],[0,-3,4,1],[0,-3,1,4]])
P0 = sp.Matrix([[18,-9,0,0],[-9,18,-9,-9],[0,-9,16,1],[0,-9,1,16]])
X = HD.inv() * P0; Y = HD.inv() * HA
print("det HD, HA, P0:", HD.det(), HA.det(), P0.det())
# simultaneous eigenvectors: eigenvectors of Y + pi X (generic combination) diagonalise both
Z = Y + sp.Rational(7, 3) * X
ev = Z.eigenvects()
print("eigenvalue pairs (lambda of Y, mu of X) with multiplicity:")
pairs = []
for val, mult, vecs in ev:
    for v in vecs:
        lam = sp.nsimplify((Y * v)[next(i for i in range(4) if v[i] != 0)] / v[next(i for i in range(4) if v[i] != 0)])
        mu = sp.nsimplify((X * v)[next(i for i in range(4) if v[i] != 0)] / v[next(i for i in range(4) if v[i] != 0)])
        pairs.append((sp.simplify(lam), sp.simplify(mu))); print("  ", sp.simplify(lam), sp.simplify(mu))
# normalise to det 1: lambda_i / det(Y)^{1/4}, mu_i / det(X)^{1/4}; t_i = log(mu_i n) / log(lambda_i n)
dY, dX = Y.det(), X.det()
for lam, mu in pairs:
    ln, mn = lam / dY ** sp.Rational(1, 4), mu / dX ** sp.Rational(1, 4)
    print("   log mu / log lambda (normalised) = %.15f" % (sp.log(mn) / sp.log(ln)).evalf(30) if abs(sp.log(ln).evalf()) > 1e-12 else "   (lambda normalised = 1)", "  lambda_n %.12f mu_n %.12f" % (ln.evalf(20), mn.evalf(20)))
