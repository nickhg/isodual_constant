# Numerical check for the reviewers' Q4: for m = 6, 7, 8, is there a Psi-minimal basis whose first vector is longer
# than gamma_m Delta^(1/m)?  That would show mu_m > gamma_m.  For each sample lattice (exact rational Gram matrix):
# all Psi-minimal bases (psi2_all at cap Psi_min^2), the largest first norm r_1 among them, the ratio
# q = (r_1^m / Delta) / gamma_m^m (q > 1 iff the basis beats Hermite), and whether r_1 > lambda_1.
# Samples: random perturbations of A_m, D_m, E_m and their duals at several scales (q > 1 needs a nearly critical
# lattice), then a hill climb on q from the two best samples.  Exact arithmetic for Psi; floats only to choose perturbations.
# Usage: sage mu_check.sage M NSAMP NCLIMB SEED
import sys, os, random
load(os.environ.get('PSI', '../certificates/psi.sage'))   # run from numerics/
m, NS, NC, seed = (int(a) for a in sys.argv[1:5]); random.seed(seed); set_random_seed(seed)
GM = {6: QQ(64)/3, 7: QQ(64), 8: QQ(256)}[m]          # gamma_m^m

def q_of(G):
    p2, _ = psi2_min(G)
    bases = psi2_all(G, p2)
    r1 = max(U[0] * G * U[0] for _, U in bases)
    return r1 ** m / G.det() / GM, r1 > lambda1(G), len(bases)

def perturb(G, eps):
    E = matrix(QQ, m, m, lambda i, j: 0)
    for i in range(m):
        for j in range(i, m):
            e = QQ(random.randint(-1000, 1000)) / 1000 * eps
            E[i, j] += e; E[j, i] += e if i != j else 0
    H = G + E
    return H if H.is_positive_definite() else perturb(G, eps / 2)     # retry smaller: never return G unchanged

bases0 = []
for t in [('A', m), ('D', m)] + ([('E', m)] if m in (6, 7, 8) else []):
    R = root_gram(t[0], m); bases0 += [(t, R), (t + ('*',), R.inverse() * R.det() ** (QQ(1) / 1 if False else 1))]
stats = dict(n=0, notshortest=0); best = []
for s in range(NS):
    name, R = random.choice(bases0); eps = random.choice([QQ(1)/200, QQ(1)/50, QQ(1)/20, QQ(1)/5])
    G = perturb(R, eps * R[0, 0])
    t0 = walltime(); q, ns, nb = q_of(G); stats['n'] += 1; stats['notshortest'] += ns
    best.append((q, name, G))
    print("sample %d %s eps %s: q %.6f, %d Psi-minimal bases, first vector not shortest: %s (%.0f s)" % (s, "".join(map(str, name)), eps, float(q), nb, ns, walltime(t0))); sys.stdout.flush()
best.sort(key=lambda x: -x[0])
print("m=%d samples %d, first vector of some Psi-minimal basis not shortest: %d" % (m, stats['n'], stats['notshortest']))
print("top q:", ["%.6f %s" % (float(q), "".join(map(str, n))) for q, n, _ in best[:5]])
# hill climb on q from the five best
for q0, name, G in best[:2]:
    q = q0; eps = G[0, 0] / 20
    for it in range(NC):
        H = perturb(G, eps)
        qh = q_of(H)[0]
        if qh > q: G, q = H, qh; print("  climb step %d: q %.6f" % (it, float(q))); sys.stdout.flush()
        elif it % 20 == 19: eps /= 2
    print("climb from %s: q %.6f -> %.6f" % ("".join(map(str, name)), float(q0), float(q)))
    sys.stdout.flush()
