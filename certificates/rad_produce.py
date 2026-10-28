# Producer (untrusted) for the certificate  r(D5) >= 4/(3 sqrt5):  for every traceless symmetric Lam,
#   max_f <C_f, Lam> >= r0 |Lam|_F,   r0 = 4/(3 sqrt5),   r0^2 = 16/45.
# Lam = sum_i x_i E_i (E_i orthogonal, rational; m_i = |E_i|^2).  Every Lam != 0 is a positive multiple of a point of
# the boundary of [-1,1]^14; signed permutations (Aut D5) permute the flags, so only cubes meeting the fundamental
# region R = {lam_1 >= ... >= lam_5, Lam_1j >= 0} are needed.  Dyadic cubes, bisected; each leaf is one of
#   O j            : inequality j of R is < 0 on the whole cube
#   K sigma        : cap lemma: sum_i min_{x in {lo,hi}} (a0_i x - m_i x^2/(2 sigma)) >= 10 c0^2 sigma
#   F f sigma      : same with a_i = <C_f, E_i>, right side 8 sigma / 45
#   L sigma n f_1 l_1 ... f_n l_n : same with a = sum l_k A[f_k], l_k >= 0 rational, sum l_k = 1
# Internal nodes: "S k" (split coordinate k at the midpoint; lower half first).  Every leaf is pre-checked exactly here.
# Usage: python3 rad_produce.py OUTDIR [workers]
import sys, os, json, time, itertools
from fractions import Fraction as Fr
from multiprocessing import Pool
import numpy as np
from scipy.optimize import linprog

HERE = os.path.dirname(os.path.abspath(__file__))
CF = [[Fr(x) for x in c] for c in json.load(open(os.path.join(HERE, "cf_orth.json")))]
def mat(v): return [v[5 * i:5 * i + 5] for i in range(5)]
E = []
for k in range(1, 5):
    M = [[Fr(0)] * 5 for _ in range(5)]
    for t in range(k): M[t][t] = Fr(1)
    M[k][k] = Fr(-k); E.append(M)
OFF = [(i, j) for i in range(5) for j in range(i + 1, 5)]
for i, j in OFF:
    M = [[Fr(0)] * 5 for _ in range(5)]; M[i][j] = M[j][i] = Fr(1); E.append(M)
def ip(X, Y): return sum(X[i][j] * Y[i][j] for i in range(5) for j in range(5))
mQ = [ip(Ek, Ek) for Ek in E]
AQ = [[ip(mat(c), Ek) for Ek in E] for c in CF]                   # 720 x 14 exact
L0 = [[Fr(0 if i == j else 1) for j in range(5)] for i in range(5)]
a0Q = [ip(L0, Ek) for Ek in E]
C0 = Fr(9271, 10729)                                              # cos(theta0), tan(theta0/2) = 27/100
INEQ = [[E[k][t][t] - E[k][t + 1][t + 1] for k in range(14)] for t in range(4)]
for j in range(1, 5):
    INEQ.append([Fr(1) if k == 4 + OFF.index((0, j)) else Fr(0) for k in range(14)])
m = np.array([float(x) for x in mQ]); A = np.array([[float(x) for x in r] for r in AQ]); a0 = np.array([float(x) for x in a0Q])
IN = np.array([[float(x) for x in r] for r in INEQ])
r0 = 4 / (3 * np.sqrt(5)); c0 = float(C0); nL0 = np.sqrt(20)

def exact_sep(a, lo, hi, sigma):
    s = Fr(0)
    for i in range(14):
        s += min(a[i] * x - mQ[i] * x * x / (2 * sigma) for x in (lo[i], hi[i]))
    return s

def sepmin(a, c, h, tau, coef):
    lo = c - h; hi = c + h
    return np.minimum(a * lo - coef * m * lo ** 2 / (2 * tau), a * hi - coef * m * hi ** 2 / (2 * tau)).sum(axis=-1)

def lpz(c, h, tau):
    lo = c - h; hi = c + h; F = A.shape[0]; nv = F + 14
    cost = np.zeros(nv); cost[F:] = -1; rows = []; rhs = []
    for i in range(14):
        for x in (lo[i], hi[i]):
            r = np.zeros(nv); r[:F] = -A[:, i] * x; r[F + i] = 1; rows.append(r); rhs.append(-r0 * m[i] * x * x / (2 * tau))
    Aeq = np.zeros((1, nv)); Aeq[0, :F] = 1
    res = linprog(cost, A_ub=np.array(rows), b_ub=np.array(rhs), A_eq=Aeq, b_eq=[1],
                  bounds=[(0, None)] * F + [(None, None)] * 14, method="highs")
    if res.status != 0 or -res.fun - r0 * tau / 2 < 0: return None
    return res.x[:F]

def frac_of(x):   # exact rational of a dyadic double
    return Fr(x)

def leaf(c, h):
    lo = c - h; hi = c + h
    loQ = [frac_of(x) for x in lo]; hiQ = [frac_of(x) for x in hi]
    for j, w in enumerate(IN):
        if w @ c + np.abs(w) @ h < 0:
            if sum(max(INEQ[j][i] * loQ[i], INEQ[j][i] * hiQ[i]) for i in range(14)) < 0: return "O %d" % j
    qc = np.sum(m * c ** 2); qa = qc + np.sum(m * h ** 2) / 3
    taus = (np.sqrt(qa), np.sqrt(qc), np.sqrt(np.sum(m * (np.abs(c) + h) ** 2)))
    for tau in taus:
        # cap: kappa = c0 sqrt20; l0 - q/(2 sigma) >= kappa^2 sigma / 2 = 10 c0^2 sigma  implies  l0 >= kappa sqrt(q)
        sig = Fr(tau / (c0 * nL0)).limit_denominator(4096)
        if (sig > 0 and sepmin(a0, c, h, float(sig), 1.0) >= 10 * c0 * c0 * float(sig) - 1e-9
                and exact_sep(a0Q, loQ, hiQ, sig) >= 10 * C0 * C0 * sig):
            return "K %s" % sig
        sig = Fr(tau / r0).limit_denominator(4096)
        v = sepmin(A, c, h, float(sig), 1.0) - 8 * float(sig) / 45
        f = int(np.argmax(v))
        if v[f] >= -1e-9 and exact_sep(AQ[f], loQ, hiQ, sig) >= 8 * sig / 45:
            return "F %d %s" % (f, sig)
    tau = taus[0]
    lam = lpz(c, h, tau)
    if lam is not None:
        sig = Fr(tau / r0).limit_denominator(4096)
        idx = [f for f in range(len(lam)) if lam[f] > 1e-10]
        lq = [Fr(lam[f]).limit_denominator(10 ** 6) for f in idx]
        keep = [(f, l) for f, l in zip(idx, lq) if l > 0]
        tot = sum(l for _, l in keep); keep = [(f, l / tot) for f, l in keep]
        a = [sum(l * AQ[f][i] for f, l in keep) for i in range(14)]
        if exact_sep(a, loQ, hiQ, sig) >= 8 * sig / 45:
            return "L %s %d %s" % (sig, len(keep), " ".join("%d %s" % (f, l) for f, l in keep))
    return None

def split_coord(c, h):
    return int(np.argmax(h * np.sqrt(m)))

def run_job(job):
    """job = (face k, sign, list of pre-splits) -> (log lines, stats)"""
    k, sg, path, outdir = job
    c = np.zeros(14); h = np.ones(14); c[k] = sg; h[k] = 0
    for (kk, side) in path:
        h[kk] /= 2; c[kk] += side * h[kk]
    name = os.path.join(outdir, "f%02d%s_%s.log" % (k, "p" if sg > 0 else "m", "".join("%d%s" % (kk, "u" if s > 0 else "l") for kk, s in path) or "root"))
    if os.path.exists(name) and open(name).read().endswith("END\n"): return name, None
    t0 = time.time(); lines = []; stats = {}
    stack = [(c, h, 0)]
    maxd = 0
    while stack:
        c, h, d = stack.pop()
        if d > 120: raise RuntimeError("depth")
        L = leaf(c, h)
        if L:
            lines.append(L); stats[L[0]] = stats.get(L[0], 0) + 1; maxd = max(maxd, d); continue
        kk = split_coord(c, h); lines.append("S %d" % kk)
        h2 = h.copy(); h2[kk] /= 2
        c1 = c.copy(); c1[kk] -= h2[kk]; c2 = c.copy(); c2[kk] += h2[kk]
        stack.append((c2, h2, d + 1)); stack.append((c1, h2, d + 1))     # lower half popped first
    with open(name + ".tmp", "w") as fh:
        fh.write("FACE %d %d PATH %s\n" % (k, sg, " ".join("%d %d" % (kk, s) for kk, s in path)))
        fh.write("\n".join(lines) + "\nEND\n")
    os.replace(name + ".tmp", name)
    stats["secs"] = round(time.time() - t0, 1); stats["maxdepth"] = maxd
    return name, stats

if __name__ == "__main__":
    outdir = sys.argv[1]; workers = int(sys.argv[2]) if len(sys.argv) > 2 else 8
    os.makedirs(outdir, exist_ok=True)
    jobs = []
    for k in range(14):
        for sg in (1, -1):
            free = [i for i in range(14) if i != k][:4]          # pre-split 4 coordinates -> 16 jobs per face
            for sides in itertools.product((-1, 1), repeat=4):
                jobs.append((k, sg, list(zip(free, sides)), outdir))
    t0 = time.time(); tot = {}
    with Pool(workers) as P:
        for name, st in P.imap_unordered(run_job, jobs):
            if st:
                for key, v in st.items():
                    if key in "OKFL": tot[key] = tot.get(key, 0) + v
                print(os.path.basename(name), st, flush=True)
    print("TOTAL", tot, "%.0fs" % (time.time() - t0))
