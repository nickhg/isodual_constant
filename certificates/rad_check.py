# Exact checker for the certificate  r(D5) = 4/(3 sqrt5)  (lower bound part).  Trusted; all arithmetic in Fractions.
#
# Claim.  For every traceless symmetric 5x5 Lam:  max_f <C_f, Lam> >= r0 |Lam|_F,  r0 = 4/(3 sqrt5),
# where f runs over Psi-minimal flags of D5 (in Z^5, orthonormal coordinates), C_f = sum_{k<=4} P_{V_k} - 2I.
# Equivalently (App. A), the gradient polytope Gamma at D5 contains the dual-norm ball of radius r0.
#
# Steps.
#  1. Flags: from d5_bases.json (root coordinates), check each U is unimodular and its prefix Gram determinants are
#     2, 3, 4, 4 (det 4), so Psi^2 = 384 = eta(D5) det^3: each flag is Psi-minimal.  Build C_f exactly.
#  2. Symmetry: the set {C_f} is invariant under conjugation by signed permutations (generators checked).
#     Hence it suffices to treat Lam in R = {lam_1 >= ... >= lam_5, Lam_1j >= 0 (j >= 2)}: any Lam is moved into R by
#     a permutation (sorting the diagonal) and then a sign change diag(1, e_2, ..., e_5) with e_j = sign Lam_1j.
#  3. Cap lemma at Lam0 = J - I: the 120 flags with <C_f, Lam0> = 8/3 (= r0 |Lam0|) have tangential parts whose hull
#     contains a ball of radius s (inscribed cross-polytope, witnesses in cap_witness.json), and tan(theta0/2) = 27/100
#     <= s/r0.  Then every Lam at angle <= theta0 from Lam0 satisfies the claim, and cos(theta0) = 9271/10729.
#  4. Cover: Lam = sum x_i E_i in an orthogonal rational basis; every Lam != 0 is a positive multiple of a point on the
#     boundary of [-1,1]^14, i.e. on one of 28 faces.  Each face is split into 16 boxes (first 4 free coordinates) and
#     each box by the logged tree of midpoint bisections.  Every leaf box Q is shown to satisfy one of
#       O: some inequality of R is negative on all of Q (Q misses R);
#       K: l0(x) - q(x)/(2 sigma) >= 10 c0^2 sigma on Q, so <Lam0, Lam> >= c0 sqrt20 |Lam| (inside the cap);
#       F/L: l(x) - q(x)/(2 sigma) >= 8 sigma/45 on Q, with l = <z, Lam>, z = C_f or a convex combination of C_f's,
#            so <z, Lam> >= q/(2 sigma) + r0^2 sigma/2 >= r0 |Lam| (AM-GM), and max_f <C_f, Lam> >= <z, Lam>.
#     Here q = |Lam|^2 = sum m_i x_i^2.  l - q/(2 sigma) is separable and concave, so its minimum over Q is the sum of
#     coordinate-wise minima over the two endpoints: exact.
# Usage: python3 rad_check.py CERTDIR
import sys, os, json, itertools
from fractions import Fraction as Fr

HERE = os.path.dirname(os.path.abspath(__file__))
def die(msg): print("FAIL", msg); sys.exit(1)

def matmul(X, Y): return [[sum(X[i][t] * Y[t][j] for t in range(len(Y))) for j in range(len(Y[0]))] for i in range(len(X))]
def transpose(X): return [list(r) for r in zip(*X)]
def det(M):
    M = [[Fr(x) for x in r] for r in M]; n = len(M); d = Fr(1)
    for c in range(n):
        p = next((r for r in range(c, n) if M[r][c] != 0), None)
        if p is None: return Fr(0)
        if p != c: M[c], M[p] = M[p], M[c]; d = -d
        d *= M[c][c]
        for r in range(c + 1, n):
            f = M[r][c] / M[c][c]; M[r] = [a - f * b for a, b in zip(M[r], M[c])]
    return d
def inv(M):
    n = len(M); A = [[Fr(x) for x in r] + [Fr(int(i == j)) for j in range(n)] for i, r in enumerate(M)]
    for c in range(n):
        p = next(r for r in range(c, n) if A[r][c] != 0); A[c], A[p] = A[p], A[c]
        pv = A[c][c]; A[c] = [x / pv for x in A[c]]
        for r in range(n):
            if r != c and A[r][c] != 0:
                f = A[r][c]; A[r] = [x - f * y for x, y in zip(A[r], A[c])]
    return [r[n:] for r in A]
def ip(X, Y): return sum(X[i][j] * Y[i][j] for i in range(5) for j in range(5))

# ---- 1. flags
B0 = [[1, -1, 0, 0, 0], [0, 1, -1, 0, 0], [0, 0, 1, -1, 0], [0, 0, 0, 1, -1], [0, 0, 0, 1, 1]]   # D5 simple roots
G0 = matmul(B0, transpose(B0))
if [G0[i][i] for i in range(5)] != [2] * 5 or det(G0) != 4: die("D5 basis")
bases = json.load(open(os.path.join(HERE, "d5_bases.json")))
Cset = set()
for u in bases:
    U = [u[5 * i:5 * i + 5] for i in range(5)]
    if abs(det(U)) != 1: die("not unimodular")
    V = matmul(U, B0)
    if [det(matmul(V[:k], transpose(V[:k]))) for k in range(1, 6)] != [2, 3, 4, 4, 4]: die("not Psi-minimal")
    C = [[Fr(-2) if i == j else Fr(0) for j in range(5)] for i in range(5)]
    for k in range(1, 5):
        Vk = V[:k]; P = matmul(matmul(transpose(Vk), inv(matmul(Vk, transpose(Vk)))), Vk)
        C = [[C[i][j] + P[i][j] for j in range(5)] for i in range(5)]
    Cset.add(tuple(x for r in C for x in r))
CF = sorted(Cset)
if len(CF) != 720: die("expected 720 distinct C_f, got %d" % len(CF))
Cm = [[list(c[5 * i:5 * i + 5]) for i in range(5)] for c in CF]
print("1. 720 Psi-minimal flags (prefix determinants 2,3,4,4,4), 720 distinct C_f")

# ---- 2. symmetry: generators of the signed permutation group
def conj_perm(c, p, s):     # (D C D^T) with D = signed permutation matrix e_i -> s_i e_{p(i)}
    M = [list(c[5 * i:5 * i + 5]) for i in range(5)]
    N = [[None] * 5 for _ in range(5)]
    for i in range(5):
        for j in range(5): N[p[i]][p[j]] = s[i] * s[j] * M[i][j]
    return tuple(x for r in N for x in r)
for p, s in [((1, 0, 2, 3, 4), (1,) * 5), ((1, 2, 3, 4, 0), (1,) * 5), ((0, 1, 2, 3, 4), (-1, 1, 1, 1, 1))]:
    if set(conj_perm(c, p, s) for c in CF) != Cset: die("flag set not invariant")
print("2. {C_f} invariant under signed permutations")

# ---- 3. cap lemma
J = [[Fr(0 if i == j else 1) for j in range(5)] for i in range(5)]
act = [Cm[f] for f in range(720) if ip(Cm[f], J) == Fr(8, 3)]
if len(act) != 120 or max(ip(C, J) for C in Cm) != Fr(8, 3): die("active set")
gens = []
for i in range(4):
    X = [[Fr(0)] * 5 for _ in range(5)]; X[i][i] = Fr(1); X[4][4] = Fr(-1); gens.append(X)
for i in range(5):
    for j in range(i + 1, 5):
        X = [[Fr(0)] * 5 for _ in range(5)]; X[i][j] = X[j][i] = Fr(1); gens.append(X)
Bs = []
for X in [J] + gens:
    Y = X
    for Z in Bs:
        f = ip(Y, Z) / ip(Z, Z); Y = [[Y[i][j] - f * Z[i][j] for j in range(5)] for i in range(5)]
    if any(Y[i][j] != 0 for i in range(5) for j in range(5)): Bs.append(Y)
if len(Bs) != 14: die("tangent basis")
T = Bs[1:]                                                     # orthogonal, all orthogonal to J and to I
for Tk in T:
    if sum(Tk[i][i] for i in range(5)) != 0 or ip(Tk, J) != 0: die("tangent basis not traceless / not tangent")
Gt = [[ip(C, Tk) for Tk in T] for C in act]
W = json.load(open(os.path.join(HERE, "cap_witness.json")))
tmin = {}
for w in W["witnesses"]:
    k, sg, t = w["k"], w["sign"], Fr(w["t"])
    lam = {int(f): Fr(l) for f, l in w["lam"].items()}
    if any(l < 0 for l in lam.values()) or sum(lam.values()) != 1: die("cap witness weights")
    pt = [sum(l * Gt[f][kk] for f, l in lam.items()) for kk in range(13)]
    if pt != [t * sg if kk == k else 0 for kk in range(13)] or t <= 0: die("cap witness point")
    tmin[k] = min(tmin.get(k, t), t)
if sorted(tmin) != list(range(13)): die("cap axes")
# half-axes a_k = t_k / |T_k| along orthonormal directions; ball radius s with s^-2 = sum a_k^-2
s2 = 1 / sum(ip(T[k], T[k]) / tmin[k] ** 2 for k in range(13))
th = Fr(W["tan_half_theta0"])
if th != Fr(27, 100) or not th * th * Fr(16, 45) <= s2: die("cap angle")     # tan(theta0/2)^2 <= s^2 / r0^2
c0 = (1 - th * th) / (1 + th * th)
if c0 != Fr(9271, 10729): die("cos theta0")
print("3. cap lemma: s^2 >= %s (s >= %.5f), tan(theta0/2) = 27/100 <= s/r0, cos(theta0) = %s" % (s2, float(s2) ** 0.5, c0))

# ---- 4. the cover
E = []
for k in range(1, 5):
    M = [[Fr(0)] * 5 for _ in range(5)]
    for t in range(k): M[t][t] = Fr(1)
    M[k][k] = Fr(-k); E.append(M)
OFF = [(i, j) for i in range(5) for j in range(i + 1, 5)]
for i, j in OFF:
    M = [[Fr(0)] * 5 for _ in range(5)]; M[i][j] = M[j][i] = Fr(1); E.append(M)
for a, b in itertools.combinations(range(14), 2):
    if ip(E[a], E[b]) != 0: die("basis not orthogonal")
if any(sum(Ek[i][i] for i in range(5)) != 0 for Ek in E): die("basis not traceless")
mQ = [ip(Ek, Ek) for Ek in E]
AQ = [[ip(C, Ek) for Ek in E] for C in Cm]
a0 = [ip(J, Ek) for Ek in E]
INEQ = [[E[k][t][t] - E[k][t + 1][t + 1] for k in range(14)] for t in range(4)]       # lam_t - lam_{t+1}
for j in range(1, 5):
    INEQ.append([E[k][0][j] for k in range(14)])                                         # Lam_1j
R2 = Fr(16, 45)

def sep(a, lo, hi, sigma):
    return sum(min(a[i] * x - mQ[i] * x * x / (2 * sigma) for x in (lo[i], hi[i])) for i in range(14))

def check_leaf(tok, lo, hi):
    t = tok[0]
    if t == "O":
        w = INEQ[int(tok[1])]
        return sum(max(w[i] * lo[i], w[i] * hi[i]) for i in range(14)) < 0
    if t == "K":
        sg = Fr(tok[1]); return sg > 0 and sep(a0, lo, hi, sg) >= 10 * c0 * c0 * sg
    if t == "F":
        f, sg = int(tok[1]), Fr(tok[2]); return sg > 0 and sep(AQ[f], lo, hi, sg) >= R2 * sg / 2
    if t == "L":
        sg, n = Fr(tok[1]), int(tok[2])
        pairs = [(int(tok[3 + 2 * q]), Fr(tok[4 + 2 * q])) for q in range(n)]
        if sg <= 0 or any(l < 0 for _, l in pairs) or sum(l for _, l in pairs) != 1: return False
        a = [sum(l * AQ[f][i] for f, l in pairs) for i in range(14)]
        return sep(a, lo, hi, sg) >= R2 * sg / 2
    return False

certdir = sys.argv[1]; counts = {}; files = 0
for k in range(14):
    for sg in (1, -1):
        free = [i for i in range(14) if i != k][:4]
        for sides in itertools.product((-1, 1), repeat=4):
            path = list(zip(free, sides))
            name = os.path.join(certdir, "f%02d%s_%s.log" % (k, "p" if sg > 0 else "m", "".join("%d%s" % (kk, "u" if s > 0 else "l") for kk, s in path)))
            if not os.path.exists(name): die("missing " + name)
            lines = open(name).read().split("\n")
            if lines[0] != "FACE %d %d PATH %s" % (k, sg, " ".join("%d %d" % (kk, s) for kk, s in path)): die("header " + name)
            lo = [Fr(-1)] * 14; hi = [Fr(1)] * 14; lo[k] = hi[k] = Fr(sg)
            for kk, s in path:
                mid = (lo[kk] + hi[kk]) / 2
                if s > 0: lo[kk] = mid
                else: hi[kk] = mid
            stack = [(lo, hi)]; pos = 1
            while stack:
                lo, hi = stack.pop()
                tok = lines[pos].split(); pos += 1
                if tok[0] == "S":
                    kk = int(tok[1])
                    if kk == k or lo[kk] == hi[kk]: die("bad split " + name)
                    mid = (lo[kk] + hi[kk]) / 2
                    hi1 = list(hi); hi1[kk] = mid; lo2 = list(lo); lo2[kk] = mid
                    stack.append((lo2, hi)); stack.append((lo, hi1))       # lower half next
                else:
                    if not check_leaf(tok, lo, hi): die("leaf %s in %s line %d" % (" ".join(tok[:3]), name, pos))
                    counts[tok[0]] = counts.get(tok[0], 0) + 1
            if lines[pos] != "END" or any(l.strip() for l in lines[pos + 1:]): die("trailing " + name)
            files += 1
print("4. %d logs, every leaf proved: %s" % (files, counts))
print("THEOREM HOLDS: max_f <C_f, Lam> >= (4/(3 sqrt5)) |Lam| for all traceless Lam; r(D5) >= 4/(3 sqrt5)")
