"""check_lift2: an independent re-implementation of the stage-2 leaf checker, written from the log-format specification
in Appendix E of the paper (and the header of check_lift.c), in Python with Arb ball arithmetic (python-flint), so
that it shares no code, no arithmetic and no language with check_lift.c.  Only one table is taken from check_lift.c: the
integral targets K_0 = Cartan(D5) and K_1 = 4 Gram(D5^*) of the (B) leaves (T_t = K_t / det(K_t)^{1/5}), since the
specification fixes them only up to the choice of basis.

Coordinates z = (r, phi_1..phi_4, y_1, y_2, y_3, mu_10, mu_20, mu_21, mu_30, mu_31, mu_32);
  S = mu R mu^T, R = diag(e^{y1}, e^{y2-y1}, e^{y3-y2}, e^{-y3}), mu unit lower triangular;
  G(z) = [[r, r phi^T], [r phi, r^{-1/4} S + r phi phi^T]]   (det G = 1).
Log: "ROOT2 ver=..", 14 lines "lo hi" (hex), then the tree in pre-order: "S t" splits coordinate t at 0.5 (lo + hi)
(round to nearest), lower half first; leaves
  N k c1..c4          w = (-k, c) (c != 0): w^T G w < G_00 on the cell;
  F flag              some bound B_f(U) < log 6 on the cell, flag = "#j" (j-th distinct flag written in this log) or
                      25 integers u_11..u_55 (U's columns are a basis); here B_f is decided as the minimum over f of
                      the forms below (any one below log 6 proves the leaf);
  P m, m lines "f lambda flag": sum lambda_i B_{f_i}(U_i) / sum lambda_i < log 6;
  B t w_11..w_55      W (rows a basis): delta = |L^{-1} W G W^T L^{-T} - I|_F <= 1/11 with T_t = L L^T, so
                      d(W G W^T, T_t) <= 1/10 (Theorem D5 and its dual);
  "END complete".
Flag bounds (A_j = leading j x j block of U^T G U, det G = 1):
  f = 0:  log eta <= sum_{j=1}^4 log det A_j;      f = k >= 1:  log eta <= log eta_k + log eta_{5-k} + (5/2) log det A_k,
  eta_1 = 1, eta_2 = 2/sqrt3, eta_3 = 3/2, eta_4 = 3.
Method: every bound is reduced to the supremum over the cell of a linear functional <Q, G(z)> (for flags via the
tangent inequality of the concave log det at the cell centre, valid for any positive definite matrix), and that
supremum is bounded by the mean-value form <Q, G(c)> + sum_t |<Q, dG/dz_t(cell)>| hw_t, with dG/dz_t enclosed by ball
evaluation over the whole cell.  A leaf that fails as given is bisected (largest term of the mean-value sum; for an
(F) leaf, one of three heuristic rules for choosing among the five forms, tried in turn, see prove) up to MAXSPLIT times, each part with the same certificate; this is sound for the same reason as check_lift's EXTRA.
Usage: python3 check_lift2.py [--sample K] [--seed s] [--types NFPB] LOG...  (prints one line per checked leaf
failure, and a summary)."""
import sys, random
from fractions import Fraction as Fr
from flint import arb, arb_mat, fmpq, ctx
ctx.prec = 100
N, NZ, MAXSPLIT = 5, 14, 60
LOG6 = arb(6).log()
LOGETA = [None, arb(0), (arb(2) / arb(3).sqrt()).log(), (arb(3) / 2).log(), arb(3).log()]
K = [[[2, -1, 0, 0, 0], [-1, 2, -1, 0, 0], [0, -1, 2, -1, -1], [0, 0, -1, 2, 0], [0, 0, -1, 0, 2]],
     [[4, 0, 0, 0, 2], [0, 4, 0, 0, 2], [0, 0, 4, 0, 2], [0, 0, 0, 4, 2], [2, 2, 2, 2, 5]]]


def idet(M):
    """exact integer determinant (fraction-free Gaussian elimination)"""
    A = [[Fr(x) for x in r] for r in M]; n = len(A); d = Fr(1)
    for c in range(n):
        p = next((r for r in range(c, n) if A[r][c] != 0), None)
        if p is None: return 0
        if p != c: A[c], A[p] = A[p], A[c]; d = -d
        d *= A[c][c]
        for r in range(c + 1, n):
            f = A[r][c] / A[c][c]; A[r] = [a - f * b for a, b in zip(A[r], A[c])]
    return int(d)


def chol(T):
    """lower-triangular ball Cholesky factor of a (ball) symmetric positive definite matrix, or None"""
    L = [[arb(0)] * N for _ in range(N)]
    for i in range(N):
        for j in range(i + 1):
            s = T[i][j] - sum((L[i][k] * L[j][k] for k in range(j)), arb(0))
            if i == j:
                if not s > 0: return None
                L[i][i] = s.sqrt()
            else:
                L[i][j] = s / L[j][j]
    return L


TL = []
for t in range(2):
    dk = idet(K[t]); sc = arb(dk) ** (arb(-1) / 5)
    TL.append(chol([[arb(K[t][i][j]) * sc for j in range(N)] for i in range(N)]))
TLinv = [arb_mat([[TL[t][i][j] for j in range(N)] for i in range(N)]).inv() for t in range(2)]


def ball(lo, hi):
    """a ball containing the exact interval [lo, hi] of two doubles"""
    c = arb(fmpq(Fr(lo).numerator, Fr(lo).denominator) + fmpq(Fr(hi).numerator, Fr(hi).denominator)) / 2
    r = (Fr(hi) - Fr(lo)) / 2
    return c + arb(0, float(r) * (1 + 2.0 ** -40) + 1e-300)


def exact(x):
    """the exact midpoint of a cell coordinate, as a ball of radius 0 (rational)"""
    return arb(fmpq(x.numerator, x.denominator))


def gram_and_derivs(z, derivs):
    """G(z) for ball coordinates z, and (if derivs) the 14 partial derivatives dG/dz_t, all as 5x5 lists of balls"""
    r, ph, y1, y2, y3 = z[0], z[1:5], z[5], z[6], z[7]
    mu = [[arb(1), arb(0), arb(0), arb(0)], [z[8], arb(1), arb(0), arb(0)], [z[9], z[10], arb(1), arb(0)],
          [z[11], z[12], z[13], arb(1)]]
    R = [y1.exp(), (y2 - y1).exp(), (y3 - y2).exp(), (-y3).exp()]
    def sandwich(D, A=mu, B=mu):            # A diag(D) B^T
        return [[sum((A[i][k] * D[k] * B[j][k] for k in range(4)), arb(0)) for j in range(4)] for i in range(4)]
    S = sandwich(R)
    rq = (-(r.log()) / 4).exp()             # r^{-1/4}
    G = [[None] * N for _ in range(N)]
    G[0][0] = r
    for i in range(4):
        G[0][i + 1] = G[i + 1][0] = r * ph[i]
        for j in range(4): G[i + 1][j + 1] = rq * S[i][j] + r * ph[i] * ph[j]
    if not derivs: return G, None
    Z = arb(0); dG = []
    # r
    D = [[Z] * N for _ in range(N)]; D[0][0] = arb(1)
    for i in range(4):
        D[0][i + 1] = D[i + 1][0] = ph[i]
        for j in range(4): D[i + 1][j + 1] = -rq / (4 * r) * S[i][j] + ph[i] * ph[j]
    dG.append(D)
    # phi_k
    for k in range(4):
        D = [[Z] * N for _ in range(N)]; D[0][k + 1] = D[k + 1][0] = r
        for i in range(4):
            for j in range(4):
                v = Z
                if i == k: v += r * ph[j]
                if j == k: v += r * ph[i]
                D[i + 1][j + 1] = v
        dG.append(D)
    # y_1..y_3: dR/dy
    for dR in ([R[0], -R[1], Z, Z], [Z, R[1], -R[2], Z], [Z, Z, R[2], -R[3]]):
        dS = sandwich(dR); D = [[Z] * N for _ in range(N)]
        for i in range(4):
            for j in range(4): D[i + 1][j + 1] = rq * dS[i][j]
        dG.append(D)
    # mu_ab (a > b): dS = E_ab R mu^T + mu R E_ba
    for a, b in ((1, 0), (2, 0), (2, 1), (3, 0), (3, 1), (3, 2)):
        D = [[Z] * N for _ in range(N)]
        for i in range(4):
            for j in range(4):
                v = Z
                if i == a: v += R[b] * mu[j][b]
                if j == a: v += mu[i][b] * R[b]
                D[i + 1][j + 1] = rq * v
        dG.append(D)
    return G, dG


def pair(Q, G): return sum((Q[i][j] * G[i][j] for i in range(N) for j in range(N)), arb(0))


class Cell:
    def __init__(self, lo, hi):
        self.lo, self.hi = lo, hi
        self.c = [(Fr(a) + Fr(b)) / 2 for a, b in zip(lo, hi)]
        self.hw = [(Fr(b) - Fr(a)) / 2 for a, b in zip(lo, hi)]
        self.Gc, _ = gram_and_derivs([exact(x) for x in self.c], False)
        _, self.dG = gram_and_derivs([ball(a, b) for a, b in zip(lo, hi)], True)
    def sup_linear(self, Q):
        """upper bound of <Q, G(z)> over the cell, and the per-coordinate mean-value terms"""
        terms = [pair(Q, self.dG[t]).abs_upper() * exact(self.hw[t]) for t in range(NZ)]
        return (pair(Q, self.Gc) + sum(terms, arb(0))).upper(), terms
    def split(self, t):
        m = 0.5 * (self.lo[t] + self.hi[t])
        a = Cell(self.lo, self.hi[:t] + [m] + self.hi[t + 1:]); b = Cell(self.lo[:t] + [m] + self.lo[t + 1:], self.hi)
        return a, b


def flag_forms(cell, U):
    """[(const, Q)] for the forms f = 0..4 of a flag U (columns a basis): B_f(z) <= const + <Q, G(z)>"""
    cols = [[U[i][j] for i in range(N)] for j in range(N)]
    out = []
    pref = []
    for j in range(1, N):
        Uj = cols[:j]
        A = [[sum((Uj[a][p] * cell.Gc[p][q] * Uj[b][q] for p in range(N) for q in range(N)), arb(0)) for b in range(j)]
             for a in range(j)]
        Am = arb_mat(A); dA = Am.det()
        if not dA > 0: return None
        Y = Am.inv()
        Q = [[sum((Uj[a][p] * Y[a, b] * Uj[b][q] for a in range(j) for b in range(j)), arb(0)) for q in range(N)]
             for p in range(N)]
        pref.append((dA.log() - j, Q))                    # log det A_j <= log det A_j(c) - j + <Q_j, G>
    c0 = sum((p[0] for p in pref), arb(0)); Q0 = [[sum((p[1][a][b] for p in pref), arb(0)) for b in range(N)] for a in range(N)]
    out.append((c0, Q0))
    for k in range(1, N):
        ck, Qk = pref[k - 1]
        out.append((LOGETA[k] + LOGETA[N - k] + arb(5) / 2 * ck, [[arb(5) / 2 * Qk[a][b] for b in range(N)] for a in range(N)]))
    return out


def prove(cell, leaf):
    """True if the leaf's claim holds on the cell (bisecting up to MAXSPLIT times).  For an (F) leaf the choice of
    split coordinate is a heuristic, and up to three rules are tried in turn; any proof found is valid."""
    rules = ("alt", "norm", "score") if leaf[0] == "F" else ("alt",)
    return any(_prove(cell, leaf, 0, r) for r in rules)


def _prove(cell, leaf, depth, rule):
    typ = leaf[0]; worst_terms = None; ok = False
    if typ == "N":
        k, c = leaf[1], leaf[2]
        if not any(c): return False
        w = [-k] + c
        Q = [[arb(w[i] * w[j] - (1 if i == j == 0 else 0)) for j in range(N)] for i in range(N)]
        s, worst_terms = cell.sup_linear(Q); ok = s < 0
    elif typ == "F":
        forms = flag_forms(cell, leaf[1])
        if forms is not None:
            info = []
            for c0, Q in forms:
                s, terms = cell.sup_linear(Q); v = c0 + s
                if v < LOG6: ok = True; break
                slack = (LOG6 - c0 - pair(Q, cell.Gc)).mid()          # slack of the form at the cell centre
                if slack > 0: info.append(((v - LOG6).mid() / slack, slack, terms))
            if not ok and info:
                # score: the form closest to provable relative to its slack; alt: that form and the form with the
                # most slack in turn; norm: the largest mean-value term relative to its form's slack
                if rule == "norm":
                    worst_terms = [max(x[2][t].mid() / x[1] for x in info) for t in range(NZ)]
                elif rule == "alt" and depth % 2 == 1:
                    worst_terms = max(info, key=lambda x: x[1])[2]
                else:
                    worst_terms = min(info, key=lambda x: x[0])[2]
    elif typ == "P":
        num_c = arb(0); numQ = [[arb(0)] * N for _ in range(N)]; den = arb(0); bad = False
        for f, lam, U in leaf[1]:
            if lam < 0: return False
            forms = flag_forms(cell, U)
            if forms is None: bad = True; break
            c0, Q = forms[f]; L = arb(fmpq(Fr(lam).numerator, Fr(lam).denominator))
            num_c += L * c0; den += L
            for a in range(N):
                for b in range(N): numQ[a][b] += L * Q[a][b]
        if not bad and den > 0:
            s, worst_terms = cell.sup_linear(numQ); ok = (num_c + s) / den < LOG6
    elif typ == "B":
        t, W = leaf[1], leaf[2]
        if idet(W) not in (1, -1): return False
        Li = TLinv[t]
        M = [[sum((Li[i, k] * W[k][j] for k in range(N)), arb(0)) for j in range(N)] for i in range(N)]   # L^{-1} W
        tot = arb(0); worst_terms = [arb(0)] * NZ
        for c in range(N):
            for d in range(N):
                Q = [[M[c][p] * M[d][q] for q in range(N)] for p in range(N)]
                centre = pair(Q, cell.Gc) - (1 if c == d else 0)
                terms = [pair(Q, cell.dG[t2]).abs_upper() * exact(cell.hw[t2]) for t2 in range(NZ)]
                bound = centre.abs_upper() + sum(terms, arb(0))
                tot += bound * bound
                worst_terms = [a + b for a, b in zip(worst_terms, terms)]
        ok = tot < arb(1) / 121
    if ok: return True
    if depth >= MAXSPLIT or worst_terms is None: return False
    t = max(range(NZ), key=lambda i: float(worst_terms[i]) if isinstance(worst_terms[i], float) else worst_terms[i].mid())
    if not (cell.lo[t] < 0.5 * (cell.lo[t] + cell.hi[t]) < cell.hi[t]): return False
    a, b = cell.split(t)
    return _prove(a, leaf, depth + 1, rule) and _prove(b, leaf, depth + 1, rule)


def parse_flag(tok, table):
    if len(tok) == 1 and tok[0].startswith("#"): return table[int(tok[0][1:])]
    vals = [int(x) for x in tok]; U = [vals[5 * i:5 * i + 5] for i in range(N)]
    if idet(U) not in (1, -1): raise ValueError("flag not unimodular")
    table.append(U); return U


def leaves(lines):
    """yield (box_lo, box_hi, leaf) for every leaf of a log, replaying the tree exactly"""
    it = iter(lines); head = next(it)
    if not head.startswith("ROOT2 "): raise ValueError("bad header")
    lo, hi = [], []
    for _ in range(NZ):
        a, b = next(it).split(); lo.append(float.fromhex(a)); hi.append(float.fromhex(b))
    stack = [(lo, hi)]; table = []
    for l in it:
        tok = l.split()
        if not tok: continue
        if tok[0] == "END":
            if tok[1:] != ["complete"] or stack: raise ValueError("log incomplete")
            return
        if tok[0] == "S":
            t = int(tok[1]); lo, hi = stack.pop(); m = 0.5 * (lo[t] + hi[t])
            stack.append((lo[:t] + [m] + lo[t + 1:], hi)); stack.append((lo, hi[:t] + [m] + hi[t + 1:]))
            continue
        lo, hi = stack.pop()
        if tok[0] == "N": leaf = ("N", int(tok[1]), [int(x) for x in tok[2:6]])
        elif tok[0] == "F": leaf = ("F", parse_flag(tok[1:], table))
        elif tok[0] == "B": leaf = ("B", int(tok[1]), [[int(x) for x in tok[2 + 5 * i:7 + 5 * i]] for i in range(N)])
        elif tok[0] == "P":
            m = int(tok[1]); parts = []
            for _ in range(m):
                p = next(it).split(); parts.append((int(p[0]), float.fromhex(p[1]), parse_flag(p[2:], table)))
            leaf = ("P", parts)
        else: raise ValueError("unknown line " + l[:40])
        yield lo, hi, leaf
    raise ValueError("no END line")


def check_log(lines, select, name=""):
    """check the leaves chosen by select(index, type) -> bool; return counts {type: (checked, failed)}"""
    res = {}
    for i, (lo, hi, leaf) in enumerate(leaves(lines)):
        if not select(i, leaf[0]): continue
        ok = prove(Cell(lo, hi), leaf)
        c, f = res.get(leaf[0], (0, 0)); res[leaf[0]] = (c + 1, f + (not ok))
        if not ok: print("FAIL %s leaf %d (%s)" % (name, i, leaf[0]), flush=True)
    return res


if __name__ == "__main__":
    args = sys.argv[1:]; K_s = None; seed = 1; types = "NFPB"
    while args and args[0].startswith("--"):
        if args[0] == "--sample": K_s = float(args[1]); args = args[2:]
        elif args[0] == "--seed": seed = int(args[1]); args = args[2:]
        elif args[0] == "--types": types = args[1]; args = args[2:]
    rng = random.Random(seed); tot = {}
    for fn in args:
        sel = (lambda i, t: t in types and (K_s is None or rng.random() < K_s))
        r = check_log(open(fn).read().splitlines(), sel, fn)
        for k, (c, f) in r.items(): a, b = tot.get(k, (0, 0)); tot[k] = (a + c, b + f)
    print("checked", {k: v[0] for k, v in tot.items()}, "failed", {k: v[1] for k, v in tot.items()},
          "FAILURES", sum(v[1] for v in tot.values()))
