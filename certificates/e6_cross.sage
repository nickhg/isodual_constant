# Exact lower bound on the E6 radius by a cross-polytope: for an H^{-1}-orthogonal rational basis e_j of the
# functional space, t_j^{+-} e_j in Gamma = conv{L_f} is certified by an explicit convex combination (primal only;
# HiGHS proposes the support, exact rational solve confirms lambda >= 0).  Then r >= (sum_j 1/(t_j^2 |e_j|^2))^{-1/2}.
import numpy as np
from scipy.optimize import linprog
A, H, Ys = load("e6_A.sobj")
Hi = H.inverse(); N = A.nrows()
# Gram-Schmidt of the standard basis in the dual inner product H^{-1}
E = []
for i in range(20):
    v = vector(QQ, [int(i == j) for j in range(20)])
    for e in E: v -= (v * Hi * e) / (e * Hi * e) * e
    E.append(v)
Af = np.array(A.change_ring(RDF)).T          # 20 x N
def best_t(e):
    # max t s.t. sum lam_f a_f = t e, sum lam = 1, lam >= 0
    ef = np.array([float(x) for x in e])
    c = np.zeros(N + 1); c[-1] = -1
    Aeq = np.zeros((21, N + 1)); Aeq[:20, :N] = Af; Aeq[:20, N] = -ef; Aeq[20, :N] = 1
    beq = np.zeros(21); beq[20] = 1
    res = linprog(c, A_eq=Aeq, b_eq=beq, bounds=[(0, None)] * N + [(None, None)], method="highs")
    supp = [f for f in range(N) if res.x[f] > 1e-9]
    # exact: solve sum_{f in supp} lam_f a_f - t e = 0, sum lam = 1 (least-squares-free: pick independent columns)
    M = matrix(QQ, [list(A.row(f)) + [1] for f in supp] + [list(-e) + [0]]).T     # 21 x (|supp|+1)
    rhs = vector(QQ, [0] * 20 + [1])
    sol = M.solve_right(rhs)
    lam, t = sol[:-1], sol[-1]
    assert all(l >= 0 for l in lam) and sum(lam) == 1
    assert sum(l * A.row(f) for l, f in zip(lam, supp)) == t * e
    return t, res.x[-1]
inv2 = 0
for j, e in enumerate(E):
    tp, fp = best_t(e); tm, fm = best_t(-e)
    t = min(tp, tm); inv2 += 1 / (t**2 * (e * Hi * e))
    print(j, "t+ = %.6f (lp %.6f)  t- = %.6f" % (tp, fp, tm), flush=True)
r2 = 1 / inv2
print("certified r^2 >=", r2.n(digits=10), " r >=", sqrt(r2).n(digits=8))
save(r2, "e6_cross_r2.sobj")
