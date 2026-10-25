/* certify.c -- interval branch-and-bound certificate that eta_n <= c (C port of certify.sage).
 *
 *   gcc -O2 -o certify certify.c -lm
 *   ./certify n c [max_boxes] [balls_file]
 *
 * Coordinates of a Psi-minimal basis B of a det-1 lattice: y_k = log G_k (k = 1..n-1) and mu_ij (i > j),
 * with r_i = exp(y_i - y_{i-1}) and Gram = M diag(r) M^T.  A box is
 *   DISCARDED  if the proved window constraints exclude every point of it (no Psi-minimal basis lies there):
 *              r_i <= mu_L (G_j/G_{i-1})^(1/L), r_j >= (G_j/G_{i-1})^(1/L)/mu_L for every window i..j of length L,
 *              and Lagrange r_{i+1} >= (1 - mu_{i+1,i}^2) r_i;
 *   CERTIFIED  if some integer basis U has (interval) sup of Psi^2_U <= c over the box, via Cauchy-Binet:
 *              det(U_k G U_k^T) = sum_{|S|=k} det((U_k M)[:,S])^2 prod_{i in S} r_i.
 *              U need not be primitive: a non-primitive flag over-estimates its saturation's value;
 *   BALL       if the box lies in a certified ball around a known maximiser (balls_file, from curv.sage);
 *   otherwise it is split.  Interval arithmetic uses outward rounding (nextafter) after every operation.
 * Candidate flags are found at box centres by a floating-point Psi-minimisation (exactness not needed).
 */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <fenv.h>
#include <stdint.h>
#include <time.h>
#include <glpk.h>

/* source version, stamped into every certificate log: build with -DGITREV=\"$(git rev-parse HEAD)\" */
#ifndef GITREV
#define GITREV "unknown"
#endif
#define NMAX 6
#define DMAX (NMAX - 1 + NMAX * (NMAX - 1) / 2)

/* ---------------- interval arithmetic ---------------- */
typedef struct { double lo, hi; } I;
/* dn/up = nextafter(x, -/+INFINITY), bit for bit (incl. +-0, +-inf, NaN), inline: the libm call was most of the run time */
static inline double dn(double x) {
    if (!(x > -INFINITY)) return x;
    if (x == 0) return -0x1p-1074;
    int64_t b; memcpy(&b, &x, 8); b += x > 0 ? -1 : 1; memcpy(&x, &b, 8); return x;
}
static inline double up(double x) {
    if (!(x < INFINITY)) return x;
    if (x == 0) return 0x1p-1074;
    int64_t b; memcpy(&b, &x, 8); b += x > 0 ? 1 : -1; memcpy(&x, &b, 8); return x;
}
static inline I iv(double a, double b) { I r = {a, b}; return r; }
static inline I ipt(double a) { return iv(a, a); }
static inline I iadd(I a, I b) { return iv(dn(a.lo + b.lo), up(a.hi + b.hi)); }
static inline I isub(I a, I b) { return iv(dn(a.lo - b.hi), up(a.hi - b.lo)); }
static inline I imul(I a, I b) {
    double p[4] = {a.lo * b.lo, a.lo * b.hi, a.hi * b.lo, a.hi * b.hi};
    double lo = p[0], hi = p[0];
    for (int i = 1; i < 4; i++) { if (p[i] < lo) lo = p[i]; if (p[i] > hi) hi = p[i]; }
    return iv(dn(lo), up(hi));
}
static inline I isqr(I a) {
    double l = fabs(a.lo), h = fabs(a.hi);
    double mx = l > h ? l : h, mn = (a.lo <= 0 && a.hi >= 0) ? 0 : (l < h ? l : h);
    return iv(mn == 0 ? 0 : dn(mn * mn), up(mx * mx));
}
static inline I iscale(I a, double s) { /* s exact, s > 0 */ return iv(dn(a.lo * s), up(a.hi * s)); }
static inline I iexp(I a) { return iv(dn(dn(exp(a.lo))), up(up(exp(a.hi)))); }
static inline I ilog(I a) { return iv(dn(dn(log(a.lo))), up(up(log(a.hi)))); }

/* ---------------- problem setup ---------------- */
static int n, nmu, DIM;
static int use_sens = 0, use_comb = 1;   /* sensitivity splitting is unstable: off unless SENS is set */
static int MUI[NMAX * NMAX], MUJ[NMAX * NMAX];
static double LOGMU_HI[NMAX + 1];     /* upper bounds on log mu_L */
static double LOGETA_HI[NMAX + 1];    /* upper bounds on log eta_L (proved values), for window potentials */
static double logc_lo;
static int use_dualcut = 1;        /* NODUAL: off (needed when collecting windows, which need not satisfy it) */
/* Localisation (n = 5): a Psi-minimal basis with eta > c has both 4-windows Psi-minimal with eta > c / mu_5^(5/2),
 * since eta(L) = r_1^(5/2) eta(tail) = G_4^(5/2) eta(head) and r_1, G_4 <= mu_5.  S (SREGION file) is a set of
 * boxes in n = 4 coordinates containing every Psi-minimal 4-dim basis with eta > c4, produced by COLLECT. */
#define SMAX 2000000
static int nS = 0;
static float (*SLO)[9], (*SHI)[9];
static double collect_w = 0;       /* COLLECT: emit unresolved n = 4 boxes narrower than this instead of splitting */
static FILE *collect_f;
static double ETA_N_HI = INFINITY;
/* SHORTEST (n = 5): search bases with b_1 a shortest vector and pi_{b_1}(b_2..b_n) Psi-minimal for pi_{b_1} L.
 * Every lattice has one; such a basis is size-reduced (we choose) and LLL-reduced with delta = 1, and
 * eta(L) <= gamma(L)^{n/2} eta(pi L) with gamma(L) = r_1 <= gamma_n.  Rules that use Psi-minimality across b_1
 * are restricted to the tail window b_2..b_n or to moves fixing b_1. */
static int shortest_mode = 0;
static double LOGGAM[NMAX + 1];     /* upper bounds on log gamma_L (Hermite) */   /* log of a proved upper bound on eta_n (n = 5: 3 2^(5/4)) */

typedef struct { double lo[DMAX], hi[DMAX]; } Box;

static void unpack(const Box *b, I *y /* n+1 */, I mu[NMAX][NMAX]) {
    y[0] = ipt(0); y[n] = ipt(0);
    for (int k = 0; k < n - 1; k++) y[k + 1] = iv(b->lo[k], b->hi[k]);
    for (int i = 0; i < n; i++) for (int j = 0; j < n; j++) mu[i][j] = ipt(i == j ? 1 : 0);
    for (int t = 0; t < nmu; t++) mu[MUI[t]][MUJ[t]] = iv(b->lo[n - 1 + t], b->hi[n - 1 + t]);
}

/* Deep insertions (PotLLL moves) in closed form.  Moving b_k to position i < k multiplies Psi^2 by
 *   prod_{j=i}^{k-1} |pi_j(b_k)|^2 / r_j = prod_j ( mu_kj^2 + sum_{l=j+1}^{k} mu_kl^2 r_l / r_j ),   mu_kk = 1,
 * so a Psi-minimal basis has every such product >= 1.  Only differences of log r enter, which keeps the
 * interval evaluation tight (comparing two separately bounded potentials loses the correlation).
 * Returns 1 if some insertion lowers Psi^2 at every point of the box. */
static int use_pot = 1;
static double potbest;   /* diagnostics: smallest upper bound on a log insertion ratio seen in the last pot_beaten call */
static int pot_beaten(I mu[NMAX][NMAX], const I *logr) {
    for (int k = 1; k < n; k++)
        for (int i = 0; i < k; i++) {
            I tot = ipt(0);
            for (int j = i; j < k; j++) {
                I a = isqr(mu[k][j]);
                for (int l = j + 1; l <= k; l++) {
                    I m2 = (l == k) ? ipt(1) : isqr(mu[k][l]);
                    a = iadd(a, imul(m2, iexp(isub(logr[l], logr[j]))));
                }
                if (a.hi <= 0) return 1;
                tot = iadd(tot, a.lo > 0 ? ilog(a) : iv(-INFINITY, ilog(iv(a.hi, a.hi)).hi));
            }
            if (tot.hi < potbest) potbest = tot.hi;
            if (tot.hi < 0) return 1;
        }
    return 0;
}
/* the same test on the reversed dual basis: mu'_{ab} = (M^{-1})_{n-1-b, n-1-a}, log r'_a = -log r_{n-1-a} */
static int pot_beaten_dual(I mu[NMAX][NMAX], const I *logr) {
    I Mi[NMAX][NMAX];
    for (int i = 0; i < n; i++) {
        for (int j = 0; j < n; j++) Mi[i][j] = ipt(i == j ? 1 : 0);
        for (int j = i - 1; j >= 0; j--) {      /* (M^{-1})_{ij} = - sum_{t=j}^{i-1} mu_{it} (M^{-1})_{tj} */
            I acc = ipt(0);
            for (int t = j; t < i; t++) acc = iadd(acc, imul(mu[i][t], Mi[t][j]));
            Mi[i][j] = isub(ipt(0), acc);
        }
    }
    I mud[NMAX][NMAX], lrd[NMAX];
    for (int a = 0; a < n; a++) {
        lrd[a] = isub(ipt(0), logr[n - 1 - a]);
        for (int b = 0; b < n; b++) mud[a][b] = (b <= a) ? Mi[n - 1 - b][n - 1 - a] : ipt(0);
    }
    return pot_beaten(mud, lrd);
}

/* SHORTEST mode: discard the box if some integer v != +-e_1 has |sum v_i b_i|^2 < r_1 everywhere on it.
 * f_v = |v B|^2 / r_1 = sum_k exp(lr_k - lr_0) c_k^2, c_k = sum_i v_i mu_ik (mu_kk = 1), bounded by the
 * mean-value form f(centre) + sum_t sup|df/dx_t| halfwidth_t (exact gradients, evaluated in intervals). */
static int nshortv = 0;
static int SHORTV[256][NMAX];
static void gen_shortv(void) {
    int tot = 1; for (int i = 0; i < n; i++) tot *= 3;
    for (int code = 0; code < tot; code++) {
        int v[NMAX], c = code, first = -1, nz = 0;
        for (int i = 0; i < n; i++) { v[i] = c % 3 - 1; c /= 3; }
        for (int i = 0; i < n; i++) if (v[i]) { if (first < 0) first = i; nz++; }
        if (first < 0 || v[first] < 0) continue;
        if (nz == 1 && v[0]) continue;
        memcpy(SHORTV[nshortv++], v, sizeof(v));
    }
}
static void fv_eval(const int *v, const I *lr, I mu[NMAX][NMAX], I *val, I *grad /* DIM */) {
    I ck[NMAX], ek[NMAX], glr[NMAX];
    for (int k = 0; k < n; k++) {
        I s2 = ipt(0);
        for (int i = k; i < n; i++) if (v[i]) s2 = iadd(s2, imul(ipt(v[i]), mu[i][k]));
        ck[k] = s2;
        ek[k] = k == 0 ? ipt(1) : iexp(isub(lr[k], lr[0]));
    }
    I f = ipt(0);
    for (int k = 0; k < n; k++) f = iadd(f, imul(ek[k], isqr(ck[k])));
    *val = f;
    if (!grad) return;
    glr[0] = ipt(0);
    for (int k = 1; k < n; k++) { glr[k] = imul(ek[k], isqr(ck[k])); glr[0] = isub(glr[0], glr[k]); }
    /* y_m (coordinate m-1, m = 1..n-1): lr_{m-1} = y_m - y_{m-1}, lr_m = y_{m+1} - y_m */
    for (int m = 1; m < n; m++) grad[m - 1] = isub(glr[m - 1], glr[m]);
    for (int t = 0; t < nmu; t++) {
        int i = MUI[t], j = MUJ[t];
        grad[n - 1 + t] = v[i] ? imul(ipt(2 * v[i]), imul(ek[j], ck[j])) : ipt(0);
    }
}
static long sv_calls = 0, sv_centre = 0, sv_hits = 0; static double sv_mintot = 1e9;
static void sv_report(void) { if (sv_calls) fprintf(stderr, "shortv: calls %ld centre<1 %ld hits %ld min tot %.3f\n", sv_calls, sv_centre, sv_hits, sv_mintot); }
static int shorter_than_b1(const Box *b, const I *lr, I mu[NMAX][NMAX]) {
    sv_calls++;
    Box cb;
    for (int t = 0; t < DIM; t++) { double m = 0.5 * (b->lo[t] + b->hi[t]); cb.lo[t] = cb.hi[t] = m; }
    I yc[NMAX + 1], muc[NMAX][NMAX], lrc[NMAX];
    unpack(&cb, yc, muc);
    for (int i = 0; i < n; i++) lrc[i] = isub(yc[i + 1], yc[i]);
    for (int q = 0; q < nshortv; q++) {
        I fc, fb, g[DMAX];
        fv_eval(SHORTV[q], lrc, muc, &fc, NULL);
        if (fc.hi >= 1) continue;
        sv_centre++;
        fv_eval(SHORTV[q], lr, mu, &fb, g);
        if (fb.hi < 1) return 1;
        I tot = fc;
        for (int t = 0; t < DIM; t++) {
            double hw = up(0.5 * (b->hi[t] - b->lo[t]));
            tot = iadd(tot, ipt(up(hw * fmax(fabs(g[t].lo), fabs(g[t].hi)))));
        }
        if (tot.hi < sv_mintot) sv_mintot = tot.hi;
        if (tot.hi < 1) { sv_hits++; return 1; }
    }
    return 0;
}

static int infeasible(const I *y, I mu[NMAX][NMAX]) {
    /* duality symmetry: eta(L) = eta(L*), and the reversed dual of a Psi-minimal basis is Psi-minimal with
       y*_k = y_{n-k}; so we may assume y_1 >= y_{n-1} */
    if (use_dualcut && !shortest_mode && n >= 3 && y[1].hi < y[n - 1].lo) return 1;
    /* proved global bound on eta_n (Lemma split): a Psi-minimal basis has Psi^2 = eta(L) <= ETA_N_HI */
    if (ETA_N_HI < INFINITY && !shortest_mode) {
        I sy = ipt(0);
        for (int k = 1; k < n; k++) sy = iadd(sy, y[k]);
        if (sy.lo > ETA_N_HI) return 1;
    }
    I logr[NMAX];
    for (int i = 0; i < n; i++) logr[i] = isub(y[i + 1], y[i]);
    if (use_pot && !shortest_mode && (pot_beaten(mu, logr) || pot_beaten_dual(mu, logr))) return 1;
    if (use_pot && shortest_mode) {
        /* the tail is a Psi-minimal basis of pi_{b_1} L: same tests on its GS data */
        I mu4[NMAX][NMAX], lr4[NMAX];
        for (int i = 0; i < n - 1; i++) { lr4[i] = logr[i + 1]; for (int j = 0; j < n - 1; j++) mu4[i][j] = mu[i + 1][j + 1]; }
        int n0 = n; n = n0 - 1;
        int hit = pot_beaten(mu4, lr4) || pot_beaten_dual(mu4, lr4);
        n = n0;
        if (hit) return 1;
    }
    if (shortest_mode) {
        /* b_1 shortest: |b_i|^2 = r_i + sum_{j<i} mu_ij^2 r_j >= r_1 */
        for (int i = 1; i < n; i++) {
            I s2 = iexp(isub(logr[i], logr[0]));
            for (int j = 1; j < i; j++) s2 = iadd(s2, imul(isqr(mu[i][j]), iexp(isub(logr[j], logr[0]))));
            s2 = iadd(s2, isqr(mu[i][0]));
            if (s2.hi < 1) return 1;
        }
    }
    for (int L = 2; L <= n; L++)
        for (int i = 0; i + L - 1 < n; i++) {
            int j = i + L - 1;
            I blk = isub(y[j + 1], y[i]);
            blk = iv(dn(blk.lo / L), up(blk.hi / L));
            /* window from b_1 in SHORTEST mode: r_1 = lambda_1(L)^2 <= gamma_L G_L^(1/L); no lower bound on r_j */
            I a = isub(isub(logr[i], ipt(shortest_mode && i == 0 ? LOGGAM[L] : LOGMU_HI[L])), blk);
            if (a.lo > 0) return 1;
            if (shortest_mode && i == 0) continue;
            I b2 = isub(isub(blk, ipt(LOGMU_HI[L])), logr[j]);
            if (b2.lo > 0) return 1;
        }
    /* every window i..j (length L < n) of a Psi-minimal basis is Psi-minimal, so its normalised potential
       sum_{k=i}^{j-1} [ (y_{k+1} - y_i) - (k-i+1)/L (y_{j+1} - y_i) ] is at most log eta_L */
    for (int L = 2; L < n; L++)
        for (int i = 0; i + L - 1 < n; i++) {
            int j = i + L - 1;
            if (shortest_mode && i == 0) continue;
            I tot = ipt(0);
            I span = isub(y[j + 1], y[i]);
            for (int k = i; k < j; k++) {
                I term = isub(isub(y[k + 1], y[i]), iv(dn(span.lo * (k - i + 1) / L), up(span.hi * (k - i + 1) / L)));
                tot = iadd(tot, term);
            }
            if (tot.lo > LOGETA_HI[L]) return 1;
        }
    for (int i = 0; i + 1 < n; i++) {
        I m2 = isqr(mu[i + 1][i]);
        I one_m = isub(ipt(1), m2);
        if (one_m.lo > 0) {
            I d = isub(logr[i + 1], logr[i]);
            I lg = ilog(one_m);
            if (d.hi < lg.lo) return 1;
        }
    }
    return 0;
}

/* determinant of k x k interval matrix by Laplace expansion (k <= 5) */
static I idet(I A[NMAX][NMAX], int k, const int *cols) {
    if (k == 1) return A[0][cols[0]];
    I s = ipt(0);
    int sub[NMAX];
    for (int c = 0; c < k; c++) {
        int t = 0;
        for (int d = 0; d < k; d++) if (d != c) sub[t++] = cols[d];
        /* minor on rows 1..k-1 */
        I B[NMAX][NMAX];
        for (int r = 1; r < k; r++) for (int q = 0; q < NMAX; q++) B[r - 1][q] = A[r][q];
        I m = idet(B, k - 1, sub);
        I term = imul(A[0][cols[c]], m);
        s = (c % 2 == 0) ? iadd(s, term) : isub(s, term);
    }
    return s;
}

static int NSUB[NMAX + 1];
static int SUBS[NMAX + 1][64][NMAX];
static void gen_subsets(void) {
    for (int k = 1; k < n; k++) {
        NSUB[k] = 0;
        for (int mask = 0; mask < (1 << n); mask++) {
            if (__builtin_popcount(mask) != k) continue;
            int t = 0;
            for (int i = 0; i < n; i++) if (mask & (1 << i)) SUBS[k][NSUB[k]][t++] = i;
            NSUB[k]++;
        }
    }
}

typedef struct { int u[NMAX][NMAX]; } Flag;

static double log_psi2_upper(const Flag *U, const I *y, I mu[NMAX][NMAX]) {
    I A[NMAX][NMAX], logr[NMAX];
    for (int i = 0; i < n; i++) logr[i] = isub(y[i + 1], y[i]);
    for (int a = 0; a < n; a++)
        for (int b = 0; b < n; b++) {
            I s = ipt(0);
            for (int t = b; t < n; t++) if (U->u[a][t]) s = iadd(s, imul(ipt(U->u[a][t]), mu[t][b]));
            A[a][b] = s;
        }
    I total = ipt(0);
    for (int k = 1; k < n; k++) {
        I s = ipt(0);
        for (int q = 0; q < NSUB[k]; q++) {
            const int *S = SUBS[k][q];
            I d = idet(A, k, S);
            I lr = ipt(0);
            for (int t = 0; t < k; t++) lr = iadd(lr, logr[S[t]]);
            s = iadd(s, imul(isqr(d), iexp(lr)));
        }
        if (s.lo <= 0) return INFINITY;
        total = iadd(total, ilog(s));
    }
    return total.hi;
}

/* Split bound (Lemma split, any sublattice): eta(L) <= eta_k eta_{n-k} det(M)^{n/2} at det L = 1, for every
 * rank-k sublattice M, in particular M = span of the first k rows of U.  Valid at every point of the box,
 * Psi-minimal or not.  Returns min over k of the log upper bound. */
static int use_split = 1;
static double split_bound(const Flag *U, const I *y, I mu[NMAX][NMAX]) {
    I A[NMAX][NMAX], logr[NMAX];
    for (int i = 0; i < n; i++) logr[i] = isub(y[i + 1], y[i]);
    for (int a = 0; a < n; a++)
        for (int b = 0; b < n; b++) {
            I s = ipt(0);
            for (int t = b; t < n; t++) if (U->u[a][t]) s = iadd(s, imul(ipt(U->u[a][t]), mu[t][b]));
            A[a][b] = s;
        }
    double best = INFINITY;
    for (int k = 1; k < n; k++) {
        I s = ipt(0);
        for (int q = 0; q < NSUB[k]; q++) {
            const int *S = SUBS[k][q];
            I d = idet(A, k, S);
            I lr = ipt(0);
            for (int t = 0; t < k; t++) lr = iadd(lr, logr[S[t]]);
            s = iadd(s, imul(isqr(d), iexp(lr)));
        }
        if (s.lo <= 0) continue;
        double v = up(up(LOGETA_HI[k] + LOGETA_HI[n - k]) + up(0.5 * n * ilog(s).hi));
        if (v < best) best = v;
    }
    return best;
}


/* interval inverse of a k x k matrix (k <= 4) via adjugate / determinant */
static int iinv(I A[NMAX][NMAX], int k, I B[NMAX][NMAX]) {
    int all[NMAX];
    for (int i = 0; i < k; i++) all[i] = i;
    I d = idet(A, k, all);
    if (d.lo <= 0 && d.hi >= 0) return 0;
    for (int i = 0; i < k; i++) for (int j = 0; j < k; j++) {
        /* cofactor C_ji: delete row j, column i */
        I Mn[NMAX][NMAX];
        int rr = 0;
        for (int r = 0; r < k; r++) {
            if (r == j) continue;
            int cc = 0;
            for (int c2 = 0; c2 < k; c2++) { if (c2 == i) continue; Mn[rr][cc++] = A[r][c2]; }
            rr++;
        }
        I cof = ipt(1);
        if (k > 1) { int cols[NMAX]; for (int t = 0; t < k - 1; t++) cols[t] = t; cof = idet(Mn, k - 1, cols); }
        if ((i + j) % 2) cof = isub(ipt(0), cof);
        /* divide by d */
        I q;
        double a = cof.lo / d.lo, b = cof.lo / d.hi, c3 = cof.hi / d.lo, e = cof.hi / d.hi;
        double lo = fmin(fmin(a, b), fmin(c3, e)), hi = fmax(fmax(a, b), fmax(c3, e));
        q = iv(dn(dn(lo)), up(up(hi)));
        B[i][j] = q;
    }
    return 1;
}

/* mean-value bound: f(center) + sum_t sup|df/dx_t| * halfwidth_t, f = log Psi^2_U, or
   f = log Psi^2_U - log Psi^2_id = log Psi^2_U - sum_k y_k when sub_id is set */
static double lastcontrib[DMAX];
static I lastgrad[DMAX];     /* interval gradient from the most recent centred evaluation */
static double combcontrib[DMAX];   /* per-coordinate terms of the last LP combination bound */
static int combvalid = 0;
static double lastfc;        /* centre value (upper bound) from the most recent centred evaluation */
static int use_tangent = 1;   /* NOTANGENT: off */
static double log_psi2_centered_g(const Flag *U, const Box *b, int sub_id) {
    /* f at the centre (tight: point intervals) */
    Box cb;
    for (int t = 0; t < DIM; t++) { double m = 0.5 * (b->lo[t] + b->hi[t]); cb.lo[t] = m; cb.hi[t] = m; }
    I yc[NMAX + 1], muc[NMAX][NMAX];
    unpack(&cb, yc, muc);
    double fc = log_psi2_upper(U, yc, muc);
    if (!isfinite(fc)) return INFINITY;
    if (sub_id) {
        I sy = ipt(0);
        for (int k = 1; k < n; k++) sy = iadd(sy, yc[k]);
        fc = up(fc - sy.lo);
    }
    lastfc = fc;
    /* gradient over the whole box */
    I y[NMAX + 1], mu[NMAX][NMAX], r[NMAX], G[NMAX][NMAX];
    unpack(b, y, mu);
    for (int i = 0; i < n; i++) r[i] = iexp(isub(y[i + 1], y[i]));
    for (int i = 0; i < n; i++) for (int j = 0; j <= i; j++) {
        I acc = ipt(0);
        for (int k = 0; k <= j; k++) acc = iadd(acc, imul(imul(mu[i][k], mu[j][k]), r[k]));
        G[i][j] = G[j][i] = acc;
    }
    /* X = U M (interval); A_k = X_k D X_k^T.  det and cofactors of A_k by (generalised) Cauchy-Binet:
         det A_k        = sum_{|S|=k}   det(X_k[:,S])^2                 prod_{i in S} r_i      (sum of squares: stays > 0)
         cof_{ij}(A_k)  = (-1)^{i+j} sum_{|S|=k-1} det(X_k^{-i}[:,S]) det(X_k^{-j}[:,S]) prod_{i in S} r_i
       and A_k^{-1} = cof^T / det. */
    I X[NMAX][NMAX];
    for (int a = 0; a < n; a++) for (int c2 = 0; c2 < n; c2++) {
        I acc = ipt(0);
        for (int t = c2; t < n; t++) if (U->u[a][t]) acc = iadd(acc, imul(ipt(U->u[a][t]), mu[t][c2]));
        X[a][c2] = acc;
    }
    I S[NMAX][NMAX];
    for (int i = 0; i < n; i++) for (int j = 0; j < n; j++) S[i][j] = ipt(0);
    for (int k = 1; k < n; k++) {
        I det = ipt(0);
        for (int q = 0; q < NSUB[k]; q++) {
            const int *SS = SUBS[k][q];
            I d = idet(X, k, SS);
            I pr = ipt(1);
            for (int t = 0; t < k; t++) pr = imul(pr, r[SS[t]]);
            det = iadd(det, imul(isqr(d), pr));
        }
        if (det.lo <= 0) return INFINITY;
        I Ai[NMAX][NMAX];
        if (k == 1) Ai[0][0] = iv(dn(1 / det.hi), up(1 / det.lo));
        else {
            /* minors of X_k with row i removed, over column subsets of size k-1 */
            static I mnr[NMAX][64];
            for (int i = 0; i < k; i++) {
                I Xi[NMAX][NMAX];
                int rr = 0;
                for (int a = 0; a < k; a++) if (a != i) { for (int c2 = 0; c2 < n; c2++) Xi[rr][c2] = X[a][c2]; rr++; }
                for (int q = 0; q < NSUB[k - 1]; q++) mnr[i][q] = idet(Xi, k - 1, SUBS[k - 1][q]);
            }
            for (int i = 0; i < k; i++) for (int j = 0; j < k; j++) {
                I cof = ipt(0);
                for (int q = 0; q < NSUB[k - 1]; q++) {
                    const int *SS = SUBS[k - 1][q];
                    I pr = ipt(1);
                    for (int t = 0; t < k - 1; t++) pr = imul(pr, r[SS[t]]);
                    cof = iadd(cof, imul(imul(mnr[i][q], mnr[j][q]), pr));
                }
                if ((i + j) % 2) cof = isub(ipt(0), cof);
                /* (A^{-1})_{ji} = cof_{ij} / det; A symmetric so cof symmetric */
                double a1 = cof.lo / det.lo, a2 = cof.lo / det.hi, a3 = cof.hi / det.lo, a4 = cof.hi / det.hi;
                Ai[j][i] = iv(dn(dn(fmin(fmin(a1, a2), fmin(a3, a4)))), up(up(fmax(fmax(a1, a2), fmax(a3, a4)))));
            }
        }
        /* S += U_k^T Ai U_k */
        for (int p = 0; p < n; p++) for (int q = 0; q < n; q++) {
            I acc = ipt(0);
            for (int a = 0; a < k; a++) if (U->u[a][p]) for (int c2 = 0; c2 < k; c2++) if (U->u[c2][q])
                acc = iadd(acc, imul(ipt(U->u[a][p] * U->u[c2][q]), Ai[a][c2]));
            S[p][q] = iadd(S[p][q], acc);
        }
    }
    /* SM = S M, T = M^T S M */
    I SM[NMAX][NMAX], T[NMAX][NMAX];
    for (int p = 0; p < n; p++) for (int j = 0; j < n; j++) {
        I acc = ipt(0);
        for (int t = j; t < n; t++) acc = iadd(acc, imul(S[p][t], mu[t][j]));
        SM[p][j] = acc;
    }
    for (int i = 0; i < n; i++) {
        I acc = ipt(0);
        for (int p = i; p < n; p++) acc = iadd(acc, imul(mu[p][i], SM[p][i]));
        T[i][i] = acc;
    }
    double bound = fc;
    for (int k = 1; k < n; k++) {      /* y_k, box coordinate k-1 */
        I g = isub(imul(T[k - 1][k - 1], r[k - 1]), imul(T[k][k], r[k]));
        if (sub_id) g = isub(g, ipt(1));
        lastgrad[k - 1] = g;
        double gm = fmax(fabs(g.lo), fabs(g.hi));
        double cb_ = up(gm * 0.5 * (b->hi[k - 1] - b->lo[k - 1]));
        lastcontrib[k - 1] = cb_;
        bound = up(bound + cb_);
    }
    for (int t = 0; t < nmu; t++) {   /* d f / d mu_ij = 2 (S M D)_ij = 2 SM_ij r_j */
        int i = MUI[t], j = MUJ[t];
        I g = iscale(imul(SM[i][j], r[j]), 2.0);
        lastgrad[n - 1 + t] = g;
        double gm = fmax(fabs(g.lo), fabs(g.hi));
        double cb_ = up(gm * 0.5 * (b->hi[n - 1 + t] - b->lo[n - 1 + t]));
        lastcontrib[n - 1 + t] = cb_;
        bound = up(bound + cb_);
    }
    return bound;
}
static double log_psi2_centered(const Flag *U, const Box *b) { return log_psi2_centered_g(U, b, 0); }

/* tangent-plane bound (as in check_phi): f(G) = sum_{k<n} log det(U_k G U_k^T) is concave in G, so
   f(G(x)) <= f(G*) + F(x) - F(x*) with F(x) = tr(P M R M^T), P = sum_k U_k^T (U_k G* U_k^T)^{-1} U_k at the centre;
   F(x) - F(x*) by the mean value theorem with dF/dmu_ab = 2 (P M R)_ab, dF/dy_s = r_{s-1}(M^T P M)_{s-1,s-1}
   - r_s (M^T P M)_{ss} enclosed over the box.  Only the outer function's gradient is enclosed, so the dependency
   loss of log_psi2_centered_g (whose gradient of f itself is enclosed) is avoided.  sub_id: bound f - sum_k y_k. */
static double log_psi2_tangent(const Flag *U, const Box *b, int sub_id) {
    Box cb;
    for (int t = 0; t < DIM; t++) { double m = 0.5 * (b->lo[t] + b->hi[t]); cb.lo[t] = m; cb.hi[t] = m; }
    I yc[NMAX + 1], muc[NMAX][NMAX], rc[NMAX], Gc[NMAX][NMAX];
    unpack(&cb, yc, muc);
    for (int i = 0; i < n; i++) rc[i] = iexp(isub(yc[i + 1], yc[i]));
    for (int i = 0; i < n; i++) for (int j = 0; j <= i; j++) {
        I acc = ipt(0);
        for (int k = 0; k <= j; k++) acc = iadd(acc, imul(imul(muc[i][k], muc[j][k]), rc[k]));
        Gc[i][j] = Gc[j][i] = acc;
    }
    I A[NMAX][NMAX];
    for (int a = 0; a < n; a++) for (int c2 = 0; c2 < n; c2++) {
        I acc = ipt(0);
        for (int i = 0; i < n; i++) if (U->u[a][i]) for (int j = 0; j < n; j++) if (U->u[c2][j])
            acc = iadd(acc, imul(ipt((double)U->u[a][i] * U->u[c2][j]), Gc[i][j]));
        A[a][c2] = acc;
    }
    /* Cholesky of the leading (n-1) x (n-1) block once: prefixes share it */
    I L[NMAX][NMAX], val = ipt(0), ld = ipt(0), P[NMAX][NMAX], Wt[NMAX][NMAX];
    for (int i = 0; i < n; i++) for (int j = 0; j < n; j++) P[i][j] = ipt(0);
    for (int c = 0; c < n - 1; c++) {
        I sq = A[c][c];
        for (int q = 0; q < c; q++) sq = isub(sq, isqr(L[c][q]));
        if (!(sq.lo > 0)) return INFINITY;
        L[c][c] = iv(dn(sqrt(sq.lo)), up(sqrt(sq.hi)));
        for (int i = c + 1; i < n - 1; i++) {
            I t = A[i][c];
            for (int q = 0; q < c; q++) t = isub(t, imul(L[i][q], L[c][q]));
            double a1 = t.lo / L[c][c].lo, a2 = t.lo / L[c][c].hi, a3 = t.hi / L[c][c].lo, a4 = t.hi / L[c][c].hi;
            L[i][c] = iv(dn(fmin(fmin(a1, a2), fmin(a3, a4))), up(fmax(fmax(a1, a2), fmax(a3, a4))));
        }
        ld = iadd(ld, iscale(ilog(L[c][c]), 2));     /* log det A_{c+1} */
        val = iadd(val, ld);
        /* row c of L^{-1} U: Wt[c][i] = (U[c][i] - sum_{q<c} L[c][q] Wt[q][i]) / L[c][c]; prefix k = c+1 adds
           sum_{c' <= c} Wt[c'] Wt[c']^T, so P = sum_c (n-1-c) Wt[c] Wt[c]^T */
        for (int i = 0; i < n; i++) {
            I t = ipt(U->u[c][i]);
            for (int q = 0; q < c; q++) t = isub(t, imul(L[c][q], Wt[q][i]));
            double a1 = t.lo / L[c][c].lo, a2 = t.lo / L[c][c].hi, a3 = t.hi / L[c][c].lo, a4 = t.hi / L[c][c].hi;
            Wt[c][i] = iv(dn(fmin(fmin(a1, a2), fmin(a3, a4))), up(fmax(fmax(a1, a2), fmax(a3, a4))));
        }
        for (int i = 0; i < n; i++) for (int j = 0; j < n; j++)
            P[i][j] = iadd(P[i][j], iscale(imul(Wt[c][i], Wt[c][j]), (double)(n - 1 - c)));
    }
    if (sub_id) for (int k = 1; k < n; k++) val = isub(val, yc[k]);
    /* gradient of F over the box */
    I y[NMAX + 1], mu[NMAX][NMAX], r[NMAX];
    unpack(b, y, mu);
    for (int i = 0; i < n; i++) r[i] = iexp(isub(y[i + 1], y[i]));
    I PM[NMAX][NMAX], Q[NMAX];
    for (int i = 0; i < n; i++) for (int k = 0; k < n; k++) {
        I acc = ipt(0);
        for (int j = k; j < n; j++) acc = iadd(acc, imul(P[i][j], mu[j][k]));
        PM[i][k] = acc;
    }
    for (int k = 0; k < n; k++) {
        I acc = ipt(0);
        for (int i = k; i < n; i++) acc = iadd(acc, imul(mu[i][k], PM[i][k]));
        Q[k] = acc;
    }
    I tot = val;
    for (int s1 = 1; s1 < n; s1++) {
        I g = isub(imul(r[s1 - 1], Q[s1 - 1]), imul(r[s1], Q[s1]));
        if (sub_id) g = isub(g, ipt(1));
        I dx = iv(dn(b->lo[s1 - 1] - cb.lo[s1 - 1]), up(b->hi[s1 - 1] - cb.lo[s1 - 1]));
        tot = iadd(tot, imul(g, dx));
    }
    for (int t = 0; t < nmu; t++) {
        I g = iscale(imul(PM[MUI[t]][MUJ[t]], r[MUJ[t]]), 2);
        I dx = iv(dn(b->lo[n - 1 + t] - cb.lo[n - 1 + t]), up(b->hi[n - 1 + t] - cb.lo[n - 1 + t]));
        tot = iadd(tot, imul(g, dx));
    }
    return tot.hi;
}



/* sup over the box of min_U l_U <= min over w in the simplex of
     sum_U w_U l_U(c) + sum_t hw_t * sup |sum_U w_U g_{U,t}|     (mean value theorem per flag)
   We search w on pairs of flags (ternary search; the bound is convex in w). */
#define CMAX 48
static double comb_bound(const Flag *F, int nf, const Box *b, int *nused) {
    static double fcs[CMAX];
    static I gr[CMAX][DMAX];
    double hw[DMAX];
    int m = 0;
    for (int t = 0; t < DIM; t++) hw[t] = 0.5 * (b->hi[t] - b->lo[t]);
    for (int q = 0; q < nf && m < CMAX; q++) {
        double bq = log_psi2_centered(&F[q], b);
        if (!isfinite(bq)) continue;
        fcs[m] = lastfc;
        memcpy(gr[m], lastgrad, sizeof(I) * DIM);
        m++;
    }
    *nused = m;
    double best = INFINITY;
    for (int a = 0; a < m; a++) {
        double v = fcs[a];
        for (int t = 0; t < DIM; t++) v = up(v + up(hw[t] * fmax(fabs(gr[a][t].lo), fabs(gr[a][t].hi))));
        if (v < best) best = v;
    }
    for (int a = 0; a < m; a++) for (int c2 = a + 1; c2 < m; c2++) {
        double lo = 0, hi = 1;
        for (int it = 0; it < 40; it++) {
            double l1 = lo + (hi - lo) / 3, l2 = hi - (hi - lo) / 3, v1 = 0, v2 = 0;
            for (int pass = 0; pass < 2; pass++) {
                double lam = pass ? l2 : l1;
                double v = (1 - lam) * fcs[a] + lam * fcs[c2];
                for (int t = 0; t < DIM; t++) {
                    double glo = (1 - lam) * gr[a][t].lo + lam * gr[c2][t].lo;
                    double ghi = (1 - lam) * gr[a][t].hi + lam * gr[c2][t].hi;
                    v += hw[t] * fmax(fabs(glo), fabs(ghi));
                }
                if (pass) v2 = v; else v1 = v;
            }
            if (v1 < v2) hi = l2; else lo = l1;
        }
        /* rigorous evaluation at the chosen lambda */
        double lam = 0.5 * (lo + hi);
        I L = ipt(lam), Lc = isub(ipt(1), L);
        I v = iadd(imul(Lc, ipt(fcs[a])), imul(L, ipt(fcs[c2])));
        for (int t = 0; t < DIM; t++) {
            I g = iadd(imul(Lc, gr[a][t]), imul(L, gr[c2][t]));
            v = iadd(v, imul(ipt(hw[t]), ipt(fmax(fabs(g.lo), fabs(g.hi)))));
        }
        if (v.hi < best) best = v.hi;
    }
    return best;
}


/* LP version: minimise over w in the simplex  sum_U w_U fc_U + sum_t hw_t u_t,
   u_t >= sum_U w_U ghi_{U,t},  u_t >= -sum_U w_U glo_{U,t}.  The bound is then re-evaluated rigorously at the
   returned w (interval arithmetic), so it is valid whatever the LP solver's accuracy. */
static double lp_comb_bound(const Flag *F, int nf, const Box *b) {
    static double fcs[CMAX];
    static I gr[CMAX][DMAX];
    double hw[DMAX];
    int m = 0;
    for (int t = 0; t < DIM; t++) hw[t] = 0.5 * (b->hi[t] - b->lo[t]);
    for (int q = 0; q < nf && m < CMAX; q++) {
        double bq = log_psi2_centered(&F[q], b);
        if (!isfinite(bq)) continue;
        fcs[m] = lastfc;
        memcpy(gr[m], lastgrad, sizeof(I) * DIM);
        m++;
    }
    if (m < 2) return INFINITY;
    glp_prob *lp = glp_create_prob();
    glp_set_obj_dir(lp, GLP_MIN);
    int nv = m + DIM;
    glp_add_cols(lp, nv);
    for (int q = 0; q < m; q++) { glp_set_col_bnds(lp, q + 1, GLP_LO, 0, 0); glp_set_obj_coef(lp, q + 1, fcs[q]); }
    for (int t = 0; t < DIM; t++) { glp_set_col_bnds(lp, m + t + 1, GLP_LO, 0, 0); glp_set_obj_coef(lp, m + t + 1, hw[t]); }
    int nr = 2 * DIM + 1;
    glp_add_rows(lp, nr);
    static int ia[1 + (2 * DMAX + 1) * (CMAX + 1)], ja[1 + (2 * DMAX + 1) * (CMAX + 1)];
    static double ar[1 + (2 * DMAX + 1) * (CMAX + 1)];
    int ne = 0;
    for (int t = 0; t < DIM; t++) {
        /* u_t - sum w ghi >= 0 */
        glp_set_row_bnds(lp, 2 * t + 1, GLP_LO, 0, 0);
        for (int q = 0; q < m; q++) { ne++; ia[ne] = 2 * t + 1; ja[ne] = q + 1; ar[ne] = -gr[q][t].hi; }
        ne++; ia[ne] = 2 * t + 1; ja[ne] = m + t + 1; ar[ne] = 1;
        /* u_t + sum w glo >= 0 */
        glp_set_row_bnds(lp, 2 * t + 2, GLP_LO, 0, 0);
        for (int q = 0; q < m; q++) { ne++; ia[ne] = 2 * t + 2; ja[ne] = q + 1; ar[ne] = gr[q][t].lo; }
        ne++; ia[ne] = 2 * t + 2; ja[ne] = m + t + 1; ar[ne] = 1;
    }
    glp_set_row_bnds(lp, nr, GLP_FX, 1, 1);
    for (int q = 0; q < m; q++) { ne++; ia[ne] = nr; ja[ne] = q + 1; ar[ne] = 1; }
    glp_load_matrix(lp, ne, ia, ja, ar);
    glp_smcp parm;
    glp_init_smcp(&parm);
    parm.msg_lev = GLP_MSG_OFF;
    parm.it_lim = 5000;      /* the simplex can cycle on degenerate LPs (seen in slice 469); a capped LP gives no bound */
    double result = INFINITY;
    if (glp_simplex(lp, &parm) == 0 && glp_get_status(lp) == GLP_OPT) {
        double w[CMAX], sw = 0;
        for (int q = 0; q < m; q++) { w[q] = glp_get_col_prim(lp, q + 1); if (w[q] < 0) w[q] = 0; sw += w[q]; }
        if (sw > 0) {
            /* rigorous evaluation with normalised weights (normalisation done in intervals) */
            I S = ipt(0);
            for (int q = 0; q < m; q++) S = iadd(S, ipt(w[q]));
            I v = ipt(0);
            I wq[CMAX];
            for (int q = 0; q < m; q++) {
                wq[q] = iv(dn(w[q] / S.hi), up(w[q] / S.lo));
                if (w[q] == 0) wq[q] = ipt(0);
                v = iadd(v, imul(wq[q], ipt(fcs[q])));
            }
            /* the normalised weights may sum to slightly more than 1; bound sum_U w_U l_U >= min_U l_U needs
               sum w = 1 exactly, so we use w_q / sum and the interval enclosure above, then require the
               resulting bound to hold for all weight vectors in the enclosure: v is an enclosure of
               sum_q w_q fc_q, and each gradient term is bounded below. */
            for (int t = 0; t < DIM; t++) {
                I g = ipt(0);
                for (int q = 0; q < m; q++) if (w[q] > 0) g = iadd(g, imul(wq[q], gr[q][t]));
                double term = hw[t] * fmax(fabs(g.lo), fabs(g.hi));
                combcontrib[t] = term;
                v = iadd(v, imul(ipt(hw[t]), ipt(fmax(fabs(g.lo), fabs(g.hi)))));
            }
            combvalid = 1;
            result = v.hi;
        }
    }
    glp_delete_prob(lp);
    return result;
}

static int is_identity(const Flag *U) {
    for (int i = 0; i < n; i++) for (int j = 0; j < n; j++) if (U->u[i][j] != (i == j)) return 0;
    return 1;
}

/* (A) no point of the box is a Psi-minimal basis if some U beats the identity everywhere on it */
static int use_elem_cg = 1;   /* NOELEMCG: skip the centred test below (closed no box in any PHI-mode sample, ~40% of the time) */
static int beats_identity(const Flag *U, const Box *b) {
    if (is_identity(U)) return 0;
    if (shortest_mode) {   /* only moves that keep b_1 */
        if (abs(U->u[0][0]) != 1) return 0;
        for (int c = 1; c < n; c++) if (U->u[0][c]) return 0;
    }
    I y[NMAX + 1], mu[NMAX][NMAX];
    unpack(b, y, mu);
    I sy = ipt(0);
    for (int k = 1; k < n; k++) sy = iadd(sy, y[k]);
    if (log_psi2_upper(U, y, mu) < sy.lo) return 1;     /* naive interval form: robust on wide boxes */
    if (use_tangent && log_psi2_tangent(U, b, 1) < 0) return 1;
    return use_elem_cg && log_psi2_centered_g(U, b, 1) < 0;
}

static double flag_bound(const Flag *U, const Box *b, const I *y, I mu[NMAX][NMAX]) {
    double a = log_psi2_upper(U, y, mu);
    if (a <= logc_lo) return a;
    double c2 = log_psi2_centered(U, b);
    if (use_tangent && c2 > logc_lo) { double c3 = log_psi2_tangent(U, b, 0); if (c3 < c2) c2 = c3; }
    return a < c2 ? a : c2;
}

/* ---------------- floating-point Psi minimisation (candidate flags) ---------------- */
static double bestval;
static int bestflag[NMAX][NMAX], curflag[NMAX][NMAX];
#define KBEST 128
static int kbest_cap = 32;     /* flags kept by the centre enumeration (COUPLED raises it for its working set) */
static int nkbest;
static double kbestval[KBEST];
static int kbestflag[KBEST][NMAX][NMAX];
static void kbest_insert(double val) {
    /* keep the kbest_cap smallest distinct values (bestval = worst kept, for pruning) */
    for (int q = 0; q < nkbest; q++) if (fabs(kbestval[q] - val) < 1e-12 * (1 + fabs(val))) {
        int same = 1;
        for (int i = 0; i < n && same; i++) for (int j = 0; j < n; j++) if (kbestflag[q][i][j] != curflag[i][j]) { same = 0; break; }
        if (same) return;
    }
    if (nkbest < kbest_cap) { kbestval[nkbest] = val; memcpy(kbestflag[nkbest], curflag, sizeof(curflag)); nkbest++; }
    else {
        int w = 0;
        for (int q = 1; q < kbest_cap; q++) if (kbestval[q] > kbestval[w]) w = q;
        if (val >= kbestval[w]) return;
        kbestval[w] = val; memcpy(kbestflag[w], curflag, sizeof(curflag));
    }
}
static double kbest_worst(void) {
    if (nkbest < kbest_cap) return INFINITY;
    double w = kbestval[0];
    for (int q = 1; q < kbest_cap; q++) if (kbestval[q] > w) w = kbestval[q];
    return w;
}

static long gcdl(long a, long b) { a = labs(a); b = labs(b); while (b) { long t = a % b; a = b; b = t; } return a; }

/* enumerate x in Z^m, x != 0, first nonzero > 0, with x^T G x <= cap; store up to maxout */
static int fp_enum(int m, double G[NMAX][NMAX], double cap, int out[][NMAX], double *nrm, int maxout) {
    double R[NMAX][NMAX] = {{0}};           /* Cholesky G = R^T R, R upper */
    for (int i = 0; i < m; i++) {
        double s = G[i][i];
        for (int k = 0; k < i; k++) s -= R[k][i] * R[k][i];
        if (s <= 0) return 0;
        R[i][i] = sqrt(s);
        for (int j = i + 1; j < m; j++) {
            double t = G[i][j];
            for (int k = 0; k < i; k++) t -= R[k][i] * R[k][j];
            R[i][j] = t / R[i][i];
        }
    }
    int cnt = 0, x[NMAX];
    double part[NMAX + 1];
    part[m] = 0;
    /* recursive enumeration from last coordinate */
    int lvl = m - 1;
    double c[NMAX];
    int lo_[NMAX], hi_[NMAX];
    /* compute range at level */
    #define RANGE(L) do { double cc = 0; for (int j = (L) + 1; j < m; j++) cc += R[(L)][j] * x[j]; cc /= R[(L)][(L)]; \
        c[(L)] = cc; double rem = cap - part[(L) + 1]; if (rem < 0) rem = 0; double w = sqrt(rem) / R[(L)][(L)]; \
        lo_[(L)] = (int)ceil(-cc - w - 1e-9); hi_[(L)] = (int)floor(-cc + w + 1e-9); x[(L)] = lo_[(L)]; } while (0)
    RANGE(lvl);
    while (1) {
        if (x[lvl] > hi_[lvl]) {
            lvl++;
            if (lvl >= m) break;
            x[lvl]++;
            continue;
        }
        double v = R[lvl][lvl] * (x[lvl] + c[lvl]);
        part[lvl] = part[lvl + 1] + v * v;
        if (part[lvl] > cap * (1 + 1e-12)) { x[lvl]++; continue; }
        if (lvl == 0) {
            int nz = -1;
            for (int i = 0; i < m; i++) if (x[i]) { nz = i; break; }
            if (nz >= 0 && x[nz] > 0 && cnt < maxout) {
                long g = 0;
                for (int i = 0; i < m; i++) g = gcdl(g, x[i]);
                if (g == 1) { for (int i = 0; i < m; i++) out[cnt][i] = x[i]; nrm[cnt] = part[0]; cnt++; }
            }
            x[0]++;
            continue;
        }
        lvl--;
        RANGE(lvl);
    }
    #undef RANGE
    return cnt;
}

static double rest_lb(double first_sq, int m) {   /* lower bound on log Psi^2 of m-dim LLL(delta=1) tail */
    double s = 0;
    for (int j = 0; j < m; j++) s += j * (m - j) * log(0.75) + (m - j) * log(first_sq);
    return s;
}

/* W: m x n integer rows (lifts in original coordinates), G: m x m Gram of the projected lattice */
static void fp_psimin(int m, double G[NMAX][NMAX], int W[NMAX][NMAX], double prev_sq, double acc, int depth) {
    if (m == 1) {
        double val = acc + log(G[0][0]);
        for (int i = 0; i < n; i++) curflag[depth][i] = W[0][i];
        kbest_insert(val);
        if (val < bestval) { bestval = val; memcpy(bestflag, curflag, sizeof(curflag)); }
        return;
    }
    double det = 1;   /* det via Cholesky */
    {
        double T[NMAX][NMAX]; memcpy(T, G, sizeof(T));
        for (int i = 0; i < m; i++) {
            for (int k = 0; k < i; k++) for (int j = i; j < m; j++) T[i][j] -= T[k][i] * T[k][j] / T[k][k];
            det *= T[i][i];
        }
    }
    double logcap = ((m - 1) / 2.0) * log(4.0 / 3.0) + log(det) / m + 1e-9;
    static int cand[NMAX][4096][NMAX];
    static double cn[NMAX][4096];
    int nc = fp_enum(m, G, exp(logcap), cand[m], cn[m], 4096);
    /* sort by norm */
    for (int a = 1; a < nc; a++) {
        double kn = cn[m][a]; int kx[NMAX]; memcpy(kx, cand[m][a], sizeof(kx));
        int b = a - 1;
        while (b >= 0 && cn[m][b] > kn) { cn[m][b + 1] = cn[m][b]; memcpy(cand[m][b + 1], cand[m][b], sizeof(kx)); b--; }
        cn[m][b + 1] = kn; memcpy(cand[m][b + 1], kx, sizeof(kx));
    }
    for (int a = 0; a < nc; a++) {
        double n2 = cn[m][a];
        if (n2 < 0.75 * prev_sq * (1 - 1e-12)) continue;
        if (acc + m * log(n2) + rest_lb(0.75 * n2, m - 1) >= kbest_worst()) break;
        /* complete x to unimodular T with first row x: column Euclid on x, tracking C^{-1} */
        long xs[NMAX], Ci[NMAX][NMAX];
        for (int i = 0; i < m; i++) { xs[i] = cand[m][a][i]; for (int j = 0; j < m; j++) Ci[i][j] = (i == j); }
        while (1) {
            int jmin = -1;
            for (int i = 0; i < m; i++) if (xs[i] && (jmin < 0 || labs(xs[i]) < labs(xs[jmin]))) jmin = i;
            int others = 0;
            for (int i = 0; i < m; i++) if (i != jmin && xs[i]) {
                long q = xs[i] / xs[jmin];
                xs[i] -= q * xs[jmin];
                for (int t = 0; t < m; t++) Ci[jmin][t] += q * Ci[i][t];   /* row j of C^{-1} += q row i */
                if (xs[i]) others = 1;
            }
            if (!others) {
                if (jmin != 0) {
                    long tmp = xs[0]; xs[0] = xs[jmin]; xs[jmin] = tmp;
                    for (int t = 0; t < m; t++) { long tt = Ci[0][t]; Ci[0][t] = Ci[jmin][t]; Ci[jmin][t] = tt; }
                }
                if (xs[0] < 0) { xs[0] = -xs[0]; for (int t = 0; t < m; t++) Ci[0][t] = -Ci[0][t]; }
                break;
            }
        }
        /* new rows W' = Ci * W ; Gram H = Ci G Ci^T */
        int W2[NMAX][NMAX];
        for (int i = 0; i < m; i++) for (int t = 0; t < n; t++) {
            long s = 0;
            for (int j = 0; j < m; j++) s += Ci[i][j] * W[j][t];
            W2[i][t] = (int)s;
        }
        double H[NMAX][NMAX];
        for (int i = 0; i < m; i++) for (int j = 0; j < m; j++) {
            double s = 0;
            for (int p = 0; p < m; p++) for (int q = 0; q < m; q++) s += Ci[i][p] * G[p][q] * Ci[j][q];
            H[i][j] = s;
        }
        double P[NMAX][NMAX];
        for (int i = 1; i < m; i++) for (int j = 1; j < m; j++) P[i - 1][j - 1] = H[i][j] - H[i][0] * H[0][j] / H[0][0];
        int Wr[NMAX][NMAX];
        for (int i = 1; i < m; i++) for (int t = 0; t < n; t++) Wr[i - 1][t] = W2[i][t];
        for (int t = 0; t < n; t++) curflag[depth][t] = W2[0][t];
        fp_psimin(m - 1, P, Wr, H[0][0], acc + m * log(H[0][0]), depth + 1);
    }
}

static void center_flag(const Box *b, Flag *out) {
    double y[NMAX + 1], r[NMAX], M[NMAX][NMAX] = {{0}}, G[NMAX][NMAX];
    y[0] = 0; y[n] = 0;
    for (int k = 0; k < n - 1; k++) y[k + 1] = 0.5 * (b->lo[k] + b->hi[k]);
    for (int i = 0; i < n; i++) { r[i] = exp(y[i + 1] - y[i]); M[i][i] = 1; }
    for (int t = 0; t < nmu; t++) M[MUI[t]][MUJ[t]] = 0.5 * (b->lo[n - 1 + t] + b->hi[n - 1 + t]);
    for (int i = 0; i < n; i++) for (int j = 0; j < n; j++) {
        double s = 0;
        for (int k = 0; k < n; k++) s += M[i][k] * r[k] * M[j][k];
        G[i][j] = s;
    }
    int W[NMAX][NMAX];
    for (int i = 0; i < n; i++) for (int j = 0; j < n; j++) W[i][j] = (i == j);
    bestval = INFINITY;
    nkbest = 0;
    for (int i = 0; i < n; i++) for (int j = 0; j < n; j++) bestflag[i][j] = (i == j);
    fp_psimin(n, G, W, 0, 0, 0);
    memcpy(out->u, bestflag, sizeof(bestflag));
}

/* ---------------- certified balls around maximisers ---------------- */
#define MAXBALL 4096
static int nball = 0;
static double BALLG[MAXBALL][NMAX][NMAX], BALLGI[MAXBALL][NMAX][NMAX], BALLR[MAXBALL];

static void load_balls(const char *fn) {
    FILE *f = fopen(fn, "r");
    if (!f) { fprintf(stderr, "cannot open %s\n", fn); exit(1); }
    while (nball < MAXBALL) {
        double rho;
        if (fscanf(f, "%lf", &rho) != 1) break;
        BALLR[nball] = rho;
        for (int i = 0; i < n; i++) for (int j = 0; j < n; j++) if (fscanf(f, "%lf", &BALLG[nball][i][j]) != 1) exit(2);
        for (int i = 0; i < n; i++) for (int j = 0; j < n; j++) if (fscanf(f, "%lf", &BALLGI[nball][i][j]) != 1) exit(2);
        nball++;
    }
    fclose(f);
    fprintf(stderr, "loaded %d balls\n", nball);
}

/* Is every Gram in the box (scaled into the slice tr(G0^{-1} G) = n) within distance rho of G0?
 * Distance = Euclidean norm of the upper-triangular entries of the difference (as in curv.sage). */
/* INVBALL: balls in the affine-invariant metric d(G, G') = |log(G^{-1/2} G' G^{-1/2})|_F.  The file gives, per
 * ball, the radius r, G_j and W_j = G_j^{-1/2}.  With Y = W_j (G - G_j) W_j and delta = |Y|_F < 1,
 * d(G, G_j) = |log(I + Y)|_F <= delta / (1 - delta), so sup delta <= r / (1 + r) puts the box inside the ball. */
static int use_invball = 0;
static int in_invball(I G[NMAX][NMAX], int s) {
    double dmax = BALLR[s] / (1 + BALLR[s]);
    I E[NMAX][NMAX], T[NMAX][NMAX];
    for (int i = 0; i < n; i++) for (int j = 0; j < n; j++) E[i][j] = isub(G[i][j], ipt(BALLG[s][i][j]));
    for (int i = 0; i < n; i++) for (int j = 0; j < n; j++) {
        I acc = ipt(0);
        for (int k = 0; k < n; k++) acc = iadd(acc, imul(ipt(BALLGI[s][i][k]), E[k][j]));
        T[i][j] = acc;
    }
    I d2 = ipt(0);
    for (int i = 0; i < n; i++) for (int j = 0; j < n; j++) {
        I acc = ipt(0);
        for (int k = 0; k < n; k++) acc = iadd(acc, imul(T[i][k], ipt(BALLGI[s][k][j])));
        d2 = iadd(d2, isqr(acc));
        if (d2.lo > dmax * dmax) return 0;
    }
    return d2.hi <= dmax * dmax;
}

/* mean-value form of delta = |W (G - G_j) W|_F over the box: Y = W (G - G_j) W at the centre plus, per entry,
 * sum_t sup|dY/dx_t| hw_t, with dG/dx_t evaluated in intervals over the box (x = (y_1..y_{n-1}, mu)). */
static int in_invball_mv(const Box *b, int s) {
    double dmax = BALLR[s] / (1 + BALLR[s]);
    I y[NMAX + 1], mu[NMAX][NMAX], yc[NMAX + 1], muc[NMAX][NMAX];
    unpack(b, y, mu);
    Box cb;
    for (int t = 0; t < DIM; t++) { double m = 0.5 * (b->lo[t] + b->hi[t]); cb.lo[t] = cb.hi[t] = m; }
    unpack(&cb, yc, muc);
    I r[NMAX], rc[NMAX];
    for (int i = 0; i < n; i++) { r[i] = iexp(isub(y[i + 1], y[i])); rc[i] = iexp(isub(yc[i + 1], yc[i])); }
    /* centre G and Y */
    I Gc[NMAX][NMAX];
    for (int i = 0; i < n; i++) for (int j = 0; j < n; j++) {
        I acc = ipt(0);
        for (int k = 0; k <= (i < j ? i : j); k++) acc = iadd(acc, imul(imul(muc[i][k], muc[j][k]), rc[k]));
        Gc[i][j] = isub(acc, ipt(BALLG[s][i][j]));
    }
    double absY[NMAX][NMAX];
    I Yc[NMAX][NMAX];
    for (int a = 0; a < n; a++) for (int c = 0; c < n; c++) {
        I acc = ipt(0);
        for (int i = 0; i < n; i++) for (int j = 0; j < n; j++)
            acc = iadd(acc, imul(ipt(BALLGI[s][a][i] * BALLGI[s][j][c]), Gc[i][j]));
        Yc[a][c] = acc; absY[a][c] = 0;
    }
    {   /* cheap rejection: the centre itself must be well inside the ball */
        double c2 = 0;
        for (int a = 0; a < n; a++) for (int c = 0; c < n; c++) { double m = fmax(fabs(Yc[a][c].lo), fabs(Yc[a][c].hi)); c2 += m * m; }
        if (c2 > dmax * dmax) return 0;
    }
    /* derivative terms */
    for (int t = 0; t < DIM; t++) {
        double hw = 0.5 * (b->hi[t] - b->lo[t]);
        if (hw == 0) continue;
        I dG[NMAX][NMAX];
        for (int i = 0; i < n; i++) for (int j = 0; j < n; j++) {
            I acc = ipt(0);
            int lim = i < j ? i : j;
            if (t < n - 1) {               /* y_m, m = t+1: d/dy_m = d/dlr_{m-1} - d/dlr_m, d g_ij/d lr_k = mu_ik mu_jk r_k */
                int m = t + 1;
                if (m - 1 <= lim) acc = iadd(acc, imul(imul(mu[i][m - 1], mu[j][m - 1]), r[m - 1]));
                if (m <= lim) acc = isub(acc, imul(imul(mu[i][m], mu[j][m]), r[m]));
            } else {                       /* mu_pq: d g_ij = [i==p] mu_jq r_q + [j==p] mu_iq r_q */
                int p = MUI[t - (n - 1)], q = MUJ[t - (n - 1)];
                if (i == p && q <= j) acc = iadd(acc, imul(mu[j][q], r[q]));
                if (j == p && q <= i) acc = iadd(acc, imul(mu[i][q], r[q]));
            }
            dG[i][j] = acc;
        }
        for (int a = 0; a < n; a++) for (int c = 0; c < n; c++) {
            I acc = ipt(0);
            for (int i = 0; i < n; i++) for (int j = 0; j < n; j++)
                acc = iadd(acc, imul(ipt(BALLGI[s][a][i] * BALLGI[s][j][c]), dG[i][j]));
            absY[a][c] = up(absY[a][c] + up(hw * fmax(fabs(acc.lo), fabs(acc.hi))));
        }
    }
    double d2 = 0;
    for (int a = 0; a < n; a++) for (int c = 0; c < n; c++) {
        double m = up(fmax(fabs(Yc[a][c].lo), fabs(Yc[a][c].hi)) + absY[a][c]);
        d2 = up(d2 + up(m * m));
        if (d2 > dmax * dmax) return 0;
    }
    return d2 <= dmax * dmax;
}

static int in_ball(const I *y, I mu[NMAX][NMAX], int s) {
    I r[NMAX], G[NMAX][NMAX];
    for (int i = 0; i < n; i++) r[i] = iexp(isub(y[i + 1], y[i]));
    for (int i = 0; i < n; i++) for (int j = 0; j <= i; j++) {
        I acc = ipt(0);
        for (int k = 0; k <= j; k++) acc = iadd(acc, imul(imul(mu[i][k], mu[j][k]), r[k]));
        G[i][j] = G[j][i] = acc;
    }
    if (use_invball) return in_invball(G, s);
    I tr = ipt(0);
    for (int i = 0; i < n; i++) for (int j = 0; j < n; j++) tr = iadd(tr, imul(ipt(BALLGI[s][i][j]), G[j][i]));
    if (tr.lo <= 0) return 0;
    I scale = iv(dn(n / tr.hi), up(n / tr.lo));
    I d2 = ipt(0);
    for (int i = 0; i < n; i++) for (int j = i; j < n; j++) {
        I d = isub(imul(scale, G[i][j]), ipt(BALLG[s][i][j]));
        d2 = iadd(d2, isqr(d));
        if (d2.lo > BALLR[s] * BALLR[s]) return 0;
    }
    return d2.hi <= BALLR[s] * BALLR[s];
}

/* ---------------- main loop ---------------- */
#define LIB 64
static Flag lib[LIB];
static int nlib = 0, libpos = 0;

typedef struct { Box b; Flag inh; int has; } Item;

static double vol0;
static double box_vol(const Box *b) { double v = 1; for (int t = 0; t < DIM; t++) v *= (b->hi[t] - b->lo[t]); return v; }

/* elementary moves: deep insertions (b_k moved to position i) and b_i += a b_j, a = +-1.  A Psi-minimal basis
   is not beaten by any of them, so a box on which one of them beats the identity everywhere is discarded. */
#define ELMAX 128
static Flag elem[ELMAX];
static int nelem = 0, use_elem = 1;
static void gen_elem(void) {
    for (int i = 0; i < n; i++) for (int k = 0; k < n; k++) if (i != k) {
        int p[NMAX], t = 0;
        for (int q = 0; q < n; q++) if (q != k) p[t++] = q;
        for (int q = n - 1; q > i; q--) p[q] = p[q - 1];
        p[i] = k;
        Flag *F = &elem[nelem++];
        for (int r = 0; r < n; r++) for (int c = 0; c < n; c++) F->u[r][c] = (c == p[r]);
    }
    for (int i = 0; i < n; i++) for (int j = 0; j < n; j++) if (i != j)
        for (int a = -1; a <= 1; a += 2) {
            Flag *F = &elem[nelem++];
            for (int r = 0; r < n; r++) for (int c = 0; c < n; c++) F->u[r][c] = (r == c);
            F->u[i][j] = a;
        }
}

static int window_in_S(const I *w /* 9 intervals: y'_1..y'_3, mu'(1,0),(2,0),(2,1),(3,0),(3,1),(3,2) */) {
    for (int q = 0; q < nS; q++) {
        int ok = 1;
        for (int t = 0; t < 9 && ok; t++) if (w[t].hi < SLO[q][t] || w[t].lo > SHI[q][t]) ok = 0;
        if (ok) return 1;
    }
    return 0;
}
static int localised_out(const I *y, I mu[NMAX][NMAX]) {
    /* head window b1..b4: y'_k = y_k - (k/4) y_4 ; tail window pi(b2..b5): y''_k = y_{k+1} - y_1 + (k/4) y_1 */
    I w[9];
    for (int k = 1; k <= 3; k++) w[k - 1] = isub(y[k], iv(dn(y[4].lo * k / 4.0), up(y[4].hi * k / 4.0)));
    int t = 3;
    for (int i = 1; i < 4; i++) for (int j = 0; j < i; j++) w[t++] = mu[i][j];
    if (!shortest_mode && !window_in_S(w)) return 1;
    for (int k = 1; k <= 3; k++) w[k - 1] = iadd(isub(y[k + 1], y[1]), iv(dn(y[1].lo * k / 4.0), up(y[1].hi * k / 4.0)));
    t = 3;
    for (int i = 1; i < 4; i++) for (int j = 0; j < i; j++) w[t++] = mu[i + 1][j + 1];
    if (!window_in_S(w)) return 1;
    return 0;
}
static void load_S(const char *fn) {
    FILE *f = fopen(fn, "r");
    if (!f) { fprintf(stderr, "cannot open %s\n", fn); exit(1); }
    SLO = malloc(sizeof(*SLO) * SMAX); SHI = malloc(sizeof(*SHI) * SMAX);
    double a, b2;
    while (nS < SMAX) {
        int t;
        for (t = 0; t < 9; t++) { if (fscanf(f, "%lf %lf", &a, &b2) != 2) break; SLO[nS][t] = nextafterf((float)a, -INFINITY); SHI[nS][t] = nextafterf((float)b2, INFINITY); }
        if (t < 9) break;
        nS++;
    }
    fclose(f);
    fprintf(stderr, "loaded %d localisation boxes\n", nS);
}

/* PHI mode (n = 4): certify Phi(M) = kappa(M)^2 eta(M) <= PHI_C, which gives eta_5 <= PHI_C (notes/eta5_reply6).
 *   kappa(M) = max_phi min_{s != 0} |s|^2 / (1 - ||phi(s)||^2)   (det M = 1).
 * Per box: kappa_bar from finitely many short s = sum k_i b_i, k in [-2,2]^4, with interval upper bounds n_s on |s|^2
 * (fewer s can only raise the max, so this is sound), first the cheap kappa <= (4/3) min_s n_s, then a bisection
 * whose infeasibility tests are a branch-and-bound over the phi-torus.  The box is then certified at the per-box
 * threshold log PHI_C - 2 log kappa_bar.  Run with c = PHI_C / (16/9 gamma_4^2) so the initial box is valid. */
static int phi_mode = 0;
static double PHI_C = 7;
static int NSHV = 0;
static int SHK[400][NMAX];
static long st_phi_refine = 0, st_phi_bb = 0;
static void gen_phi_vectors(void) {
    for (int code = 0; code < 625; code++) {
        int k[4], c = code, first = -1;
        for (int i = 0; i < 4; i++) { k[i] = c % 5 - 2; c /= 5; }
        for (int i = 0; i < 4; i++) if (k[i]) { first = i; break; }
        if (first < 0 || k[first] < 0) continue;
        memcpy(SHK[NSHV++], k, sizeof(k));
    }
}
static double PHN[400];
static I PN_GC[NMAX][NMAX], PN_DG[DMAX][NMAX][NMAX];   /* from the last phi_norms: G at the box centre, dG/dx_t over the box */
static double PN_HW[DMAX];
/* upper bounds on |s|^2 = k^T G k over the box: min of the direct interval form and the mean-value form
 * (centre value + sum_t hw_t sup |k^T dG/dx_t k|, with dG/dx_t in intervals over the box) */
static void phi_norms(const Box *b, const I *y, I mu[NMAX][NMAX]) {
    I r[NMAX];
    for (int i = 0; i < n; i++) r[i] = iexp(isub(y[i + 1], y[i]));
    Box cb;
    for (int t = 0; t < DIM; t++) { double m = 0.5 * (b->lo[t] + b->hi[t]); cb.lo[t] = cb.hi[t] = m; }
    I yc[NMAX + 1], muc[NMAX][NMAX], rc[NMAX];
    unpack(&cb, yc, muc);
    for (int i = 0; i < n; i++) rc[i] = iexp(isub(yc[i + 1], yc[i]));
    I Gc[NMAX][NMAX], dG[DMAX][NMAX][NMAX];
    double hw[DMAX];
    for (int i = 0; i < n; i++) for (int j = 0; j < n; j++) {
        I acc = ipt(0);
        for (int k = 0; k <= (i < j ? i : j); k++) acc = iadd(acc, imul(imul(muc[i][k], muc[j][k]), rc[k]));
        Gc[i][j] = acc;
    }
    for (int t = 0; t < DIM; t++) {
        hw[t] = up(0.5 * (b->hi[t] - b->lo[t]));
        for (int i = 0; i < n; i++) for (int j = 0; j < n; j++) {
            I acc = ipt(0);
            int lim = i < j ? i : j;
            if (t < n - 1) {
                int m = t + 1;
                if (m - 1 <= lim) acc = iadd(acc, imul(imul(mu[i][m - 1], mu[j][m - 1]), r[m - 1]));
                if (m <= lim) acc = isub(acc, imul(imul(mu[i][m], mu[j][m]), r[m]));
            } else {
                int p = MUI[t - (n - 1)], q = MUJ[t - (n - 1)];
                if (i == p && q <= j) acc = iadd(acc, imul(mu[j][q], r[q]));
                if (j == p && q <= i) acc = iadd(acc, imul(mu[i][q], r[q]));
            }
            dG[t][i][j] = acc;
        }
    }
    memcpy(PN_GC, Gc, sizeof(Gc)); memcpy(PN_DG, dG, sizeof(dG)); memcpy(PN_HW, hw, sizeof(hw));
    for (int q = 0; q < NSHV; q++) {
        const int *k = SHK[q];
        I tot = ipt(0);
        for (int j = 0; j < n; j++) {
            I c = ipt(0);
            for (int i = j; i < n; i++) if (k[i]) c = iadd(c, imul(ipt(k[i]), mu[i][j]));
            tot = iadd(tot, imul(r[j], isqr(c)));
        }
        I cen = ipt(0);
        for (int i = 0; i < n; i++) for (int j = 0; j < n; j++) if (k[i] && k[j]) cen = iadd(cen, imul(ipt(k[i] * k[j]), Gc[i][j]));
        double mv = cen.hi;
        for (int t = 0; t < DIM; t++) {
            if (hw[t] == 0) continue;
            I d = ipt(0);
            for (int i = 0; i < n; i++) for (int j = 0; j < n; j++) if (k[i] && k[j]) d = iadd(d, imul(ipt(k[i] * k[j]), dG[t][i][j]));
            mv = up(mv + up(hw[t] * fmax(fabs(d.lo), fabs(d.hi))));
        }
        PHN[q] = fmin(tot.hi, mv);
    }
}
/* is there phi with ||phi(s)||^2 >= 1 - n_s / c for all s?  0 = proved infeasible, 1 = possibly feasible */
static int phi_feasible(double c) {
    int act[400], na = 0; double tau[400];
    for (int q = 0; q < NSHV; q++) if (PHN[q] < c) { act[na] = q; tau[na] = sqrt(1 - PHN[q] / c) * (1 - 1e-12) - 1e-12; if (tau[na] > 0) na++; }
    if (na == 0) return 1;
    static double stk[200000][8];
    int sp = 0;
    /* phi -> -phi symmetry: phi_1 in [0, 1/2] */
    double b0[8] = {0, 0.5, 0, 1, 0, 1, 0, 1};
    memcpy(stk[sp++], b0, sizeof(b0));
    long iters = 0;
    while (sp > 0) {
        double *bx = stk[--sp], lo4[4], hi4[4];
        for (int i = 0; i < 4; i++) { lo4[i] = bx[2 * i]; hi4[i] = bx[2 * i + 1]; }
        iters++;
        int excluded = 0, allok = 1;
        for (int a = 0; a < na && !excluded; a++) {
            const int *k = SHK[act[a]];
            double lo = 0, hi = 0;
            for (int i = 0; i < 4; i++) {
                if (k[i] > 0) { lo += k[i] * lo4[i]; hi += k[i] * hi4[i]; }
                else if (k[i] < 0) { lo += k[i] * hi4[i]; hi += k[i] * lo4[i]; }
            }
            double t = tau[a];
            lo -= 1e-12; hi += 1e-12;            /* conservative against rounding in the phi intervals */
            double m = floor(0.5 * (lo + hi) + 0.5);
            if (lo > m - t && hi < m + t) { excluded = 1; break; }
            double f = floor(lo);
            if (!(lo >= f + t && hi <= f + 1 - t)) allok = 0;
        }
        if (excluded) continue;
        if (allok) { st_phi_bb += iters; return 1; }
        int d = 0; double w = -1;
        for (int i = 0; i < 4; i++) if (hi4[i] - lo4[i] > w) { w = hi4[i] - lo4[i]; d = i; }
        if (w < 1.0 / 1024 || sp + 2 >= 200000) { st_phi_bb += iters; return 1; }
        double mid = 0.5 * (lo4[d] + hi4[d]);
        double c1[8], c2[8];
        memcpy(c1, bx, sizeof(c1)); memcpy(c2, bx, sizeof(c2));
        c1[2 * d + 1] = mid; c2[2 * d] = mid;
        memcpy(stk[sp++], c1, sizeof(c1)); memcpy(stk[sp++], c2, sizeof(c2));
    }
    st_phi_bb += iters;
    return 0;
}
static double phi_kappa_cheap(void) {
    double m = INFINITY; for (int q = 0; q < NSHV; q++) if (PHN[q] < m) m = PHN[q];
    return up(m * 4.0 / 3.0);
}
/* best-first branch-and-bound for kappa_bar = max_phi min_s n_s / (1 - ||phi(s)||^2) over the torus (phi_1 in [0,1/2]).
 * On a phi-box, U = min_s n_s / (1 - dmax_s^2), where dmax_s >= ||phi(s)|| on the box, is an upper bound; the global
 * upper bound is the max of U over open boxes.  Stops when it is within rel_tol of the best centre value. */
#define PHQ 40000
static double PQ_U[PHQ], PQ_B[PHQ][8];
static int pq_n;
static void pq_push(double U, const double *bx) {
    int i = pq_n++;
    PQ_U[i] = U; memcpy(PQ_B[i], bx, sizeof(PQ_B[i]));
    while (i > 0) { int p = (i - 1) / 2; if (PQ_U[p] >= PQ_U[i]) break;
        double tu = PQ_U[p]; PQ_U[p] = PQ_U[i]; PQ_U[i] = tu;
        double tb[8]; memcpy(tb, PQ_B[p], sizeof(tb)); memcpy(PQ_B[p], PQ_B[i], sizeof(tb)); memcpy(PQ_B[i], tb, sizeof(tb));
        i = p; }
}
static void pq_pop(double *U, double *bx) {
    *U = PQ_U[0]; memcpy(bx, PQ_B[0], sizeof(PQ_B[0]));
    pq_n--; PQ_U[0] = PQ_U[pq_n]; memcpy(PQ_B[0], PQ_B[pq_n], sizeof(PQ_B[0]));
    int i = 0;
    while (1) { int l = 2 * i + 1, r = l + 1, m = i;
        if (l < pq_n && PQ_U[l] > PQ_U[m]) m = l;
        if (r < pq_n && PQ_U[r] > PQ_U[m]) m = r;
        if (m == i) break;
        double tu = PQ_U[m]; PQ_U[m] = PQ_U[i]; PQ_U[i] = tu;
        double tb[8]; memcpy(tb, PQ_B[m], sizeof(tb)); memcpy(PQ_B[m], PQ_B[i], sizeof(tb)); memcpy(PQ_B[i], tb, sizeof(tb));
        i = m; }
}
static int PH_ACT[400], ph_na;
static double phi_box_upper(const double *bx, double *centre_val) {
    double U = INFINITY, C = INFINITY;
    for (int a = 0; a < ph_na; a++) {
        const int *k = SHK[PH_ACT[a]];
        double lo = 0, hi = 0, mid = 0;
        for (int i = 0; i < 4; i++) {
            double l = bx[2 * i], h = bx[2 * i + 1];
            if (k[i] > 0) { lo += k[i] * l; hi += k[i] * h; } else if (k[i] < 0) { lo += k[i] * h; hi += k[i] * l; }
            mid += k[i] * 0.5 * (l + h);
        }
        lo -= 1e-12; hi += 1e-12;
        /* max over [lo,hi] of the distance to Z: 1/2 if a half-integer lies inside, else the larger endpoint value
           (between consecutive half-integers the distance is V-shaped) */
        double dmax;
        double kh = ceil(lo - 0.5);
        if (kh + 0.5 <= hi) dmax = 0.5;
        else dmax = fmax(fabs(lo - floor(lo + 0.5)), fabs(hi - floor(hi + 0.5)));
        double u = up(PHN[PH_ACT[a]] / dn(1 - up(dmax * dmax)));
        if (u < U) U = u;
        double dm = fabs(mid - floor(mid + 0.5));
        double cv = PHN[PH_ACT[a]] / (1 - dm * dm);
        if (cv < C) C = cv;
    }
    *centre_val = C;
    return U;
}
static double phi_kappa_bb(double cheap, double rel_tol) {
    ph_na = 0;
    for (int q = 0; q < NSHV; q++) if (PHN[q] < cheap) PH_ACT[ph_na++] = q;
    if (ph_na == 0) return cheap;
    pq_n = 0;
    double b0[8] = {0, 0.5, 0, 1, 0, 1, 0, 1}, cv;
    double U0 = phi_box_upper(b0, &cv);
    double best = cv;
    pq_push(fmin(U0, cheap), b0);
    long it = 0;
    while (pq_n > 0) {
        double U, bx[8];
        pq_pop(&U, bx);
        it++;
        if (U <= best * (1 + rel_tol) || pq_n + 2 >= PHQ || it > 20000) { st_phi_bb += it; return U; }
        int d = 0; double w = -1;
        for (int i = 0; i < 4; i++) if (bx[2 * i + 1] - bx[2 * i] > w) { w = bx[2 * i + 1] - bx[2 * i]; d = i; }
        double mid = 0.5 * (bx[2 * d] + bx[2 * d + 1]);
        for (int h = 0; h < 2; h++) {
            double c[8]; memcpy(c, bx, sizeof(c));
            if (h) c[2 * d] = mid; else c[2 * d + 1] = mid;
            double u = phi_box_upper(c, &cv);
            if (cv > best) best = cv;
            pq_push(fmin(u, U), c);
        }
    }
    st_phi_bb += it;
    return cheap;            /* unreachable (each pop pushes two boxes); cheap is a valid upper bound */
}
static double phi_kappa_refined(double cheap) {
    if (!getenv("PHIBISECT")) { st_phi_refine++; return phi_kappa_bb(cheap, 2e-3); }
    double m = INFINITY; for (int q = 0; q < NSHV; q++) if (PHN[q] < m) m = PHN[q];
    double a = m, b = cheap;            /* kappa >= ... a is only a bracket guess; b is a proved upper bound */
    if (phi_feasible(b)) return b;
    for (int it = 0; it < 10; it++) { double mid = 0.5 * (a + b); if (phi_feasible(mid)) a = mid; else b = mid; }
    st_phi_refine++;
    return b;
}
/* ---------------- COUPLED: per-phi-cell coupled bound on 2 log kappa + log eta ----------------
 * For lam_s >= 0 (sum 2) over a set of vectors s, mu_f >= 0 (sum 1) over a set of flags f, and every G in the box
 * and phi in a cell C of the torus:
 *   2 log min_s q_s/(1 - ||s.phi||^2) + log min_f h_f  <=  sum_s lam_s (log q_s + psi_s(phi)) + sum_f mu_f log h_f,
 * q_s = s^T G s, psi_s(phi) = -log(1 - ||s.phi||^2), h_f = prod_{k<=3} det(U_k G U_k^T).  Each term is bounded by a
 * tangent of a concave function, at an arbitrary point (so floating-point tangent points cost nothing in rigour):
 *   log q <= log a + q/a - 1                    (any a > 0)
 *   log det A <= -log det T + tr(T A) - k       (any symmetric positive definite k x k matrix T)
 *   psi_s(phi) <= c0_s + sl_s s.(phi - phi_C)   (the chord of psi(|t - m|), convex where s.C has no half-integer;
 *                                                else sl_s = 0 and c0_s = psi(dmax_s))
 * The G-dependence collects into tr(W G) with W = sum lam_s s s^T / a_s + sum mu_f sum_k U_k^T T_k U_k, bounded on the
 * box by its centre value plus sum_t hw_t sup|tr(W dG/dx_t)| (mean value, dG in intervals over the box).  lam, mu,
 * a_s, T_k come from an LP linearised at the box centre (floating point); everything after that is in intervals.
 * The cell's 2 log kappa + log eta, hence the box's Phi, is <= the result.  Cells whose kappa bound alone suffices
 * (with the best single-flag eta bound) are closed without the LP. */
static int use_coupled = 0;
static long st_cpl_boxes = 0, st_cpl_ok = 0, st_cpl_lp = 0, st_cpl_cells = 0;
#define CPF 128
static int cp_nf;
static Flag cp_flag[CPF];
static double cp_lh[CPF], cp_g[CPF][DMAX];           /* float: log h_f and its gradient at the centre */
static I cp_P[CPF][NMAX][NMAX], cp_F[CPF];          /* sum_k U_k^T T_k U_k (intervals) and sum_k (-log det T_k - k) */
static double cp_qc[400], cp_gq[400][DMAX];          /* float: q_s and its gradient at the centre */
static double cp_Gc[NMAX][NMAX], cp_dGc[DMAX][NMAX][NMAX];
static void fgram(const double *x, double G[NMAX][NMAX], double dG[DMAX][NMAX][NMAX]) {
    double y[NMAX + 1], r[NMAX], M[NMAX][NMAX] = {{0}};
    y[0] = 0; y[n] = 0;
    for (int k = 0; k < n - 1; k++) y[k + 1] = x[k];
    for (int i = 0; i < n; i++) { r[i] = exp(y[i + 1] - y[i]); M[i][i] = 1; }
    for (int t = 0; t < nmu; t++) M[MUI[t]][MUJ[t]] = x[n - 1 + t];
    for (int i = 0; i < n; i++) for (int j = 0; j < n; j++) {
        double s = 0; for (int k = 0; k < n; k++) s += M[i][k] * r[k] * M[j][k]; G[i][j] = s;
        for (int t = 0; t < DIM; t++) {
            double d = 0;
            if (t < n - 1) { int m = t + 1; d = M[i][m - 1] * M[j][m - 1] * r[m - 1] - M[i][m] * M[j][m] * r[m]; }
            else { int p = MUI[t - (n - 1)], q = MUJ[t - (n - 1)]; if (i == p) d += M[j][q] * r[q]; if (j == p) d += M[i][q] * r[q]; }
            dG[t][i][j] = d;
        }
    }
}
/* float inverse of a k x k SPD matrix (Gauss-Jordan), symmetrised */
static int finv(int k, double A[NMAX][NMAX], double B[NMAX][NMAX]) {
    double T[NMAX][2 * NMAX];
    for (int i = 0; i < k; i++) for (int j = 0; j < k; j++) { T[i][j] = A[i][j]; T[i][k + j] = (i == j); }
    for (int c = 0; c < k; c++) {
        if (!(T[c][c] > 0)) return 0;
        double p = T[c][c];
        for (int j = 0; j < 2 * k; j++) T[c][j] /= p;
        for (int i = 0; i < k; i++) if (i != c) { double f = T[i][c]; for (int j = 0; j < 2 * k; j++) T[i][j] -= f * T[c][j]; }
    }
    for (int i = 0; i < k; i++) for (int j = 0; j < k; j++) B[i][j] = 0.5 * (T[i][k + j] + T[j][k + i]);
    return 1;
}
/* add flag U to the working set: float data at the centre and the rigorous tangent data */
static int cp_add_flag(const Flag *U) {
    if (cp_nf >= CPF) return 0;
    int f = cp_nf;
    double lh = 0, g[DMAX] = {0};
    I P[NMAX][NMAX], F = ipt(0);
    for (int i = 0; i < n; i++) for (int j = 0; j < n; j++) P[i][j] = ipt(0);
    for (int k = 1; k < n; k++) {
        double A[NMAX][NMAX], T[NMAX][NMAX];
        for (int a = 0; a < k; a++) for (int b = 0; b < k; b++) {
            double s = 0; for (int i = 0; i < n; i++) for (int j = 0; j < n; j++) s += U->u[a][i] * cp_Gc[i][j] * U->u[b][j]; A[a][b] = s;
        }
        if (!finv(k, A, T)) return 0;
        /* T must be positive definite: leading principal minors > 0 (in intervals) */
        I TI[NMAX][NMAX]; int cols[NMAX];
        for (int a = 0; a < k; a++) for (int b = 0; b < k; b++) TI[a][b] = ipt(T[a][b]);
        I det = ipt(1);
        for (int m = 1; m <= k; m++) { for (int c = 0; c < m; c++) cols[c] = c; det = idet(TI, m, cols); if (!(det.lo > 0)) return 0; }
        F = isub(isub(F, ilog(det)), ipt(k));
        double ld = 0; { double D[NMAX][NMAX]; memcpy(D, A, sizeof(D)); for (int c = 0; c < k; c++) { ld += log(D[c][c]); for (int i = c + 1; i < k; i++) { double q = D[i][c] / D[c][c]; for (int j = c; j < k; j++) D[i][j] -= q * D[c][j]; } } }
        lh += ld;
        for (int i = 0; i < n; i++) for (int j = 0; j < n; j++) {
            I acc = ipt(0); double fa = 0;
            for (int a = 0; a < k; a++) for (int b = 0; b < k; b++) if (U->u[a][i] && U->u[b][j]) {
                acc = iadd(acc, imul(ipt(U->u[a][i] * U->u[b][j]), TI[a][b])); fa += U->u[a][i] * U->u[b][j] * T[a][b];
            }
            P[i][j] = iadd(P[i][j], acc);
            for (int t = 0; t < DIM; t++) g[t] += fa * cp_dGc[t][i][j];
        }
    }
    cp_flag[f] = *U; cp_lh[f] = lh; memcpy(cp_g[f], g, sizeof(g)); memcpy(cp_P[f], P, sizeof(P)); cp_F[f] = F;
    cp_nf++;
    return 1;
}
/* psi(d) = -log(1 - d^2), upper bound, 0 <= d <= 1/2 */
static double psi_up(double d) { I one_m = iv(dn(1 - up(d * d)), up(1 - dn(d * d))); return -ilog(one_m).lo; }
static int cp_act[400], cp_na;
static double cp_lastv;
/* rigorous value of the coupled bound on one cell for given multipliers */
static double cp_eval(const double *lo4, const double *hi4, int na, const int *act, const double *lam, const double *c0,
                      const double *sl, int nf, const int *fl, const double *mu) {
    I L = ipt(0), Mu = ipt(0);
    for (int a = 0; a < na; a++) L = iadd(L, ipt(lam[a]));
    for (int q = 0; q < nf; q++) Mu = iadd(Mu, ipt(mu[q]));
    if (!(L.lo > 0) || !(Mu.lo > 0)) return INFINITY;
    I W[NMAX][NMAX], val = ipt(0), ph[4] = {ipt(0), ipt(0), ipt(0), ipt(0)};
    for (int i = 0; i < n; i++) for (int j = 0; j < n; j++) W[i][j] = ipt(0);
    for (int a = 0; a < na; a++) {
        I w = iv(dn(2 * lam[a] / L.hi), up(2 * lam[a] / L.lo));       /* lam_a / sum lam * 2, enclosed */
        const int *k = SHK[act[a]]; double as = cp_qc[act[a]];
        I term = iadd(isub(ilog(ipt(as)), ipt(1)), ipt(c0[a]));
        val = iadd(val, imul(w, term));
        I wa = iv(dn(w.lo / as), up(w.hi / as));
        for (int i = 0; i < n; i++) for (int j = 0; j < n; j++) if (k[i] && k[j]) W[i][j] = iadd(W[i][j], imul(wa, ipt(k[i] * k[j])));
        if (sl[a] != 0) for (int i = 0; i < 4; i++) if (k[i]) ph[i] = iadd(ph[i], imul(imul(w, ipt(sl[a])), ipt(k[i])));
    }
    for (int q = 0; q < nf; q++) {
        I w = iv(dn(mu[q] / Mu.hi), up(mu[q] / Mu.lo));
        val = iadd(val, imul(w, cp_F[fl[q]]));
        for (int i = 0; i < n; i++) for (int j = 0; j < n; j++) W[i][j] = iadd(W[i][j], imul(w, cp_P[fl[q]][i][j]));
    }
    I tc = ipt(0);
    for (int i = 0; i < n; i++) for (int j = 0; j < n; j++) tc = iadd(tc, imul(W[i][j], PN_GC[i][j]));
    double v = up(val.hi + tc.hi);
    for (int t = 0; t < DIM; t++) {
        if (PN_HW[t] == 0) continue;
        I d = ipt(0);
        for (int i = 0; i < n; i++) for (int j = 0; j < n; j++) d = iadd(d, imul(W[i][j], PN_DG[t][i][j]));
        v = up(v + up(PN_HW[t] * fmax(fabs(d.lo), fabs(d.hi))));
    }
    for (int i = 0; i < 4; i++) v = up(v + up(0.5 * (hi4[i] - lo4[i]) * fmax(fabs(ph[i].lo), fabs(ph[i].hi))));
    return v;
}
/* one cell: 1 = closed */
static int cp_cell(const double *bx, double logphi, double etab_log) {
    double lo4[4], hi4[4], mid4[4];
    for (int i = 0; i < 4; i++) { lo4[i] = bx[2 * i]; hi4[i] = bx[2 * i + 1]; mid4[i] = 0.5 * (lo4[i] + hi4[i]); }
    /* per-vector cell data: [a, b] = s.cell, dmax, chord */
    double U[400], Umin = INFINITY;
    double c0[400], sl[400];
    for (int a = 0; a < cp_na; a++) {
        const int *k = SHK[cp_act[a]];
        double l = 0, h = 0, mid = 0;
        for (int i = 0; i < 4; i++) {
            if (k[i] > 0) { l += k[i] * lo4[i]; h += k[i] * hi4[i]; } else if (k[i] < 0) { l += k[i] * hi4[i]; h += k[i] * lo4[i]; }
            mid += k[i] * mid4[i];
        }
        /* cell corners are dyadic with few bits and |k_i| <= 2: l, h, mid are exact */
        double dmax, m = floor(mid + 0.5);
        int halfin = (ceil(l - 0.5) + 0.5 <= h);
        if (halfin) dmax = 0.5; else dmax = fmax(fabs(l - m), fabs(h - m));
        U[a] = up(PHN[cp_act[a]] / dn(1 - up(dmax * dmax)));
        if (U[a] < Umin) Umin = U[a];
        if (halfin) { c0[a] = psi_up(dmax); sl[a] = 0; }
        else {
            double pa = psi_up(fabs(l - m)), pb = psi_up(fabs(h - m));
            double s_ = (h > l) ? (pb - pa) / (h - l) : 0;
            I ea = imul(ipt(s_), ipt(l - mid)), eb = imul(ipt(s_), ipt(h - mid));
            c0[a] = fmax(up(pa - ea.lo), up(pb - eb.lo));
            sl[a] = s_;
        }
    }
    /* kappa alone (with the box's single-flag eta bound) */
    if (up(2 * up(log(Umin)) + 1e-12) + etab_log <= logphi) return 1;
    /* LP at the centre: max 2z + w, z <= log q_s + c0_s + (gq_s/q_s).dx + sl_s s.dphi, w <= lh_f + g_f.dx */
    int act[400], na = 0;
    for (int a = 0; a < cp_na; a++) if (U[a] <= 1.25 * Umin) act[na++] = a;
    double lh0 = INFINITY; for (int f = 0; f < cp_nf; f++) if (cp_lh[f] < lh0) lh0 = cp_lh[f];
    int fl[CPF], nf = 0;
    for (int f = 0; f < cp_nf; f++) if (cp_lh[f] <= lh0 + 0.2) fl[nf++] = f;
    int ncol = DIM + 4 + 2, nrow = na + nf;
    glp_prob *lp = glp_create_prob();
    glp_set_obj_dir(lp, GLP_MAX);
    glp_add_cols(lp, ncol);
    for (int t = 0; t < DIM; t++) glp_set_col_bnds(lp, t + 1, PN_HW[t] > 0 ? GLP_DB : GLP_FX, -PN_HW[t], PN_HW[t]);
    for (int i = 0; i < 4; i++) glp_set_col_bnds(lp, DIM + i + 1, GLP_DB, -0.5 * (hi4[i] - lo4[i]), 0.5 * (hi4[i] - lo4[i]));
    glp_set_col_bnds(lp, DIM + 5, GLP_FR, 0, 0); glp_set_obj_coef(lp, DIM + 5, 2);
    glp_set_col_bnds(lp, DIM + 6, GLP_FR, 0, 0); glp_set_obj_coef(lp, DIM + 6, 1);
    glp_add_rows(lp, nrow);
    int *ia = malloc(sizeof(int) * (nrow * ncol + 1)), *ja = malloc(sizeof(int) * (nrow * ncol + 1));
    double *ar = malloc(sizeof(double) * (nrow * ncol + 1));
    int ne = 0;
    for (int r = 0; r < na; r++) {
        int s = cp_act[act[r]]; const int *k = SHK[s];
        glp_set_row_bnds(lp, r + 1, GLP_UP, 0, log(cp_qc[s]) + c0[act[r]]);
        for (int t = 0; t < DIM; t++) { ne++; ia[ne] = r + 1; ja[ne] = t + 1; ar[ne] = -cp_gq[s][t] / cp_qc[s]; }
        for (int i = 0; i < 4; i++) if (k[i] && sl[act[r]] != 0) { ne++; ia[ne] = r + 1; ja[ne] = DIM + i + 1; ar[ne] = -sl[act[r]] * k[i]; }
        ne++; ia[ne] = r + 1; ja[ne] = DIM + 5; ar[ne] = 1;
    }
    for (int q = 0; q < nf; q++) {
        int r = na + q;
        glp_set_row_bnds(lp, r + 1, GLP_UP, 0, cp_lh[fl[q]]);
        for (int t = 0; t < DIM; t++) { ne++; ia[ne] = r + 1; ja[ne] = t + 1; ar[ne] = -cp_g[fl[q]][t]; }
        ne++; ia[ne] = r + 1; ja[ne] = DIM + 6; ar[ne] = 1;
    }
    glp_load_matrix(lp, ne, ia, ja, ar);
    glp_smcp parm; glp_init_smcp(&parm); parm.msg_lev = GLP_MSG_OFF;
    double v = INFINITY;
    st_cpl_lp++;
    if (glp_simplex(lp, &parm) == 0 && glp_get_status(lp) == GLP_OPT) {
        double lam[400], mu[CPF], c0a[400], sla[400]; int acts[400], nl = 0, fls[CPF], nm = 0;
        for (int r = 0; r < na; r++) { double d = glp_get_row_dual(lp, r + 1); if (d > 1e-12) { lam[nl] = d; acts[nl] = cp_act[act[r]]; c0a[nl] = c0[act[r]]; sla[nl] = sl[act[r]]; nl++; } }
        for (int q = 0; q < nf; q++) { double d = glp_get_row_dual(lp, na + q + 1); if (d > 1e-12) { mu[nm] = d; fls[nm] = fl[q]; nm++; } }
        if (nl > 0 && nm > 0) v = cp_eval(lo4, hi4, nl, acts, lam, c0a, sla, nm, fls, mu);
    }
    cp_lastv = v;
    glp_delete_prob(lp); free(ia); free(ja); free(ar);
    return v <= logphi;
}
static double cp_cell0 = 1.0 / 16, cp_cellmin = 1.0 / 64;
static int cp_maxcells = 200000;
/* 1 = Phi <= PHI_C on the box (PHN, PN_* must be those of this box) */
static int coupled_bound(const Box *b, const I *y, I mu[NMAX][NMAX]) {
    st_cpl_boxes++;
    double xc[DMAX];
    for (int t = 0; t < DIM; t++) xc[t] = 0.5 * (b->lo[t] + b->hi[t]);
    fgram(xc, cp_Gc, cp_dGc);
    double logphi = dn(log(PHI_C) - 1e-12);
    /* vectors: n_s below the cheap kappa bound (the others cannot lower any cell's min) */
    double cheap = phi_kappa_cheap();
    cp_na = 0;
    for (int q = 0; q < NSHV; q++) if (PHN[q] < cheap) {
        const int *k = SHK[q]; double s = 0;
        for (int i = 0; i < n; i++) for (int j = 0; j < n; j++) s += k[i] * cp_Gc[i][j] * k[j];
        cp_qc[q] = s;
        for (int t = 0; t < DIM; t++) { double d = 0; for (int i = 0; i < n; i++) for (int j = 0; j < n; j++) d += k[i] * cp_dGc[t][i][j] * k[j]; cp_gq[q][t] = d; }
        if (s > 0) cp_act[cp_na++] = q;
    }
    /* flags: the centre enumeration's best ones, and the identity */
    int save = kbest_cap; kbest_cap = CPF - 1;
    Flag f; center_flag(b, &f);
    kbest_cap = save;
    cp_nf = 0;
    Flag id; for (int i = 0; i < n; i++) for (int j = 0; j < n; j++) id.u[i][j] = (i == j);
    cp_add_flag(&id);
    for (int q = 0; q < nkbest; q++) { Flag g; memcpy(g.u, kbestflag[q], sizeof(g.u)); cp_add_flag(&g); }
    /* rigorous single-flag eta bound for the kappa-only cell test: the best few flags */
    double etab = INFINITY;
    int ord[CPF]; for (int q = 0; q < cp_nf; q++) ord[q] = q;
    for (int a = 1; a < cp_nf; a++) { int t = ord[a], c = a - 1; while (c >= 0 && cp_lh[ord[c]] > cp_lh[t]) { ord[c + 1] = ord[c]; c--; } ord[c + 1] = t; }
    for (int q = 0; q < cp_nf && q < 4; q++) { double e = flag_bound(&cp_flag[ord[q]], b, y, mu); if (e < etab) etab = e; }
    /* cells of the torus (phi_1 in [0, 1/2] by phi -> -phi) */
    static double stk[4096][8];
    int sp = 0; long cells = 0;
    double b0[8] = {0, 0.5, 0, 1, 0, 1, 0, 1};
    memcpy(stk[sp++], b0, sizeof(b0));
    while (sp > 0) {
        double bx[8]; memcpy(bx, stk[--sp], sizeof(bx));
        if (++cells > cp_maxcells) { st_cpl_cells += cells; return 0; }
        double w = 0; int d = 0;
        for (int i = 0; i < 4; i++) if (bx[2 * i + 1] - bx[2 * i] > w) { w = bx[2 * i + 1] - bx[2 * i]; d = i; }
        /* kappa-only test on every cell; the LP once the cell is small enough */
        if (w > cp_cell0) {
            /* cheap test: reuse cp_cell's first stage by asking with etab only (no LP) via a large-cell shortcut */
            double U = INFINITY;
            for (int a = 0; a < cp_na; a++) {
                const int *k = SHK[cp_act[a]]; double l = 0, h = 0;
                for (int i = 0; i < 4; i++) { if (k[i] > 0) { l += k[i] * bx[2 * i]; h += k[i] * bx[2 * i + 1]; } else if (k[i] < 0) { l += k[i] * bx[2 * i + 1]; h += k[i] * bx[2 * i]; } }
                double dm = (ceil(l - 0.5) + 0.5 <= h) ? 0.5 : fmax(fabs(l - floor(l + 0.5)), fabs(h - floor(h + 0.5)));
                double u = up(PHN[cp_act[a]] / dn(1 - up(dm * dm))); if (u < U) U = u;
            }
            if (up(2 * up(log(U)) + 1e-12) + etab <= logphi) continue;
        } else {
            if (cp_cell(bx, logphi, etab)) continue;
            if (w <= cp_cellmin) {
                st_cpl_cells += cells;
                if (getenv("CPDBG")) {      /* failing cell: its excess, the box width, and the centre's distance to the nearest ball */
                    double wm = 0; for (int t = 0; t < DIM; t++) wm = fmax(wm, b->hi[t] - b->lo[t]);
                    double dbest = INFINITY;
                    for (int q = 0; q < nball; q++) {     /* |W_q (G_c - G_q) W_q|_F, first-order invariant distance */
                        double d2 = 0;
                        for (int a = 0; a < n; a++) for (int c = 0; c < n; c++) {
                            double y_ = 0; for (int i = 0; i < n; i++) for (int j = 0; j < n; j++) y_ += BALLGI[q][a][i] * (cp_Gc[i][j] - BALLG[q][i][j]) * BALLGI[q][j][c];
                            d2 += y_ * y_;
                        }
                        if (d2 < dbest) dbest = d2;
                    }
                    fprintf(stderr, "CPFAIL width %.4f dist %.4f excess %.5f nflags %d nact %d phi %.4f %.4f %.4f %.4f\n", wm, sqrt(dbest), cp_lastv - logphi, cp_nf, cp_na, bx[0], bx[2], bx[4], bx[6]);
                }
                return 0;
            }
        }
        if (sp + 2 >= 4096) { st_cpl_cells += cells; return 0; }
        double mid = 0.5 * (bx[2 * d] + bx[2 * d + 1]);
        memcpy(stk[sp], bx, sizeof(bx)); stk[sp][2 * d + 1] = mid; sp++;
        memcpy(stk[sp], bx, sizeof(bx)); stk[sp][2 * d] = mid; sp++;
    }
    st_cpl_cells += cells; st_cpl_ok++;
    return 1;
}
static void cpl_report(void) { if (st_cpl_boxes) fprintf(stderr, "coupled: boxes %ld closed %ld cells %ld LPs %ld\n", st_cpl_boxes, st_cpl_ok, st_cpl_cells, st_cpl_lp); }
static double last_kappa = 0;
static double hope_h = 0, last_phi_centre = 0;
static int hope_on = 0, last_hopeless = 0;

static Flag hint_flag; static int hint_has = 0;   /* flag that closed the last box (written to CERTLOG as a hint) */
static int classify_inner(const Box *b, const Flag *inh, int has, Flag *f);
static int classify(const Box *b, const Flag *inh, int has, Flag *f) {
    hint_has = 0;
    if (!phi_mode) return classify_inner(b, inh, has, f);
    I y[NMAX + 1], mu[NMAX][NMAX];
    unpack(b, y, mu);
    phi_norms(b, y, mu);
    double save = logc_lo;
    double kc = phi_kappa_cheap();
    /* eta(M) <= eta_4 = 3 for every 4-lattice (Theorem eta4), so kappa_bar^2 * 3 <= PHI_C certifies the box outright */
    if (up(up(kc * kc) * 3) <= PHI_C) { last_kappa = kc; return 1; }
    last_kappa = kc;
    logc_lo = dn(log(PHI_C) - 2 * log(kc) - 1e-12);   /* libm log is not directed: explicit margin */
    int r = classify_inner(b, inh, has, f);
    if (!r) {
        double kr = phi_kappa_refined(kc);
        last_kappa = kr;
        if (up(up(kr * kr) * 3) <= PHI_C) { logc_lo = save; return 1; }
        if (kr < kc) { logc_lo = dn(log(PHI_C) - 2 * log(kr) - 1e-12); r = classify_inner(b, inh, has, f); }
        if (!r && use_coupled) { phi_norms(b, y, mu); r = coupled_bound(b, y, mu); }
    }
    /* HOPELESS=h (with COLLECT): an unresolved box whose centre has Phi > PHI_C (1 + h), numerically, can never
       close, so it is collected at once instead of being bisected down to the COLLECT width.  Collecting is always
       sound; this is only a heuristic for when to stop. */
    last_hopeless = 0;
    if (!r && hope_on) {
        Box cb; for (int t = 0; t < DIM; t++) cb.lo[t] = cb.hi[t] = 0.5 * (b->lo[t] + b->hi[t]);
        I yc[NMAX + 1], mc[NMAX][NMAX]; unpack(&cb, yc, mc);
        phi_norms(&cb, yc, mc);
        double kcen = phi_kappa_bb(phi_kappa_cheap(), 1e-3);
        Flag fc; center_flag(b, &fc);
        last_phi_centre = kcen * kcen * exp(bestval);
        if (last_phi_centre > PHI_C * (1 + hope_h)) last_hopeless = 1;
    }
    if (getenv("PHIDBG")) {
        double wm = 0; for (int t = 0; t < DIM; t++) wm = fmax(wm, b->hi[t] - b->lo[t]);
        if (hope_on && !r) fprintf(stderr, "  centre Phi %.5f%s\n", last_phi_centre, last_hopeless ? " HOPELESS" : "");
        fprintf(stderr, "  width %.4f kappa_bar %.5f (cheap %.5f) -> eta threshold %.4f  %s\n", wm, last_kappa, kc, exp(logc_lo), r ? "RESOLVED" : "");
        if (!r && getenv("PHIDBG")[0] == '2') {
            Box cb; for (int t = 0; t < DIM; t++) cb.lo[t] = cb.hi[t] = 0.5 * (b->lo[t] + b->hi[t]);
            I yc[NMAX + 1], mc[NMAX][NMAX]; unpack(&cb, yc, mc);
            double rr[NMAX]; for (int i = 0; i < n; i++) rr[i] = exp(yc[i + 1].lo - yc[i].lo);
            fprintf(stderr, "  GRAM");
            for (int i = 0; i < n; i++) for (int j = i; j < n; j++) { double g = 0; for (int k = 0; k <= i; k++) g += mc[i][k].lo * mc[j][k].lo * rr[k]; fprintf(stderr, " %.10g", g); }
            fprintf(stderr, "\n");
        }
    }
    logc_lo = save;
    return r;
}

static int classify_inner(const Box *b, const Flag *inh, int has, Flag *f) {
    combvalid = 0;
    I y[NMAX + 1], mu[NMAX][NMAX];
    unpack(b, y, mu);
    if (infeasible(y, mu)) return 1;
    if (shortest_mode) {
        I lr[NMAX];
        for (int i = 0; i < n; i++) lr[i] = isub(y[i + 1], y[i]);
        if (shorter_than_b1(b, lr, mu)) return 1;
    }
    if (nS && n == 5 && localised_out(y, mu)) return 1;
    I sy = ipt(0);
    for (int k = 1; k < n; k++) sy = iadd(sy, y[k]);
    if (sy.hi <= logc_lo) return 1;
    if (use_elem) for (int q = 0; q < nelem; q++) if (beats_identity(&elem[q], b)) return 1;
    if (use_split) {
        Flag id;
        for (int i = 0; i < n; i++) for (int j = 0; j < n; j++) id.u[i][j] = (i == j);
        if (split_bound(&id, y, mu) <= logc_lo) return 1;
        if (has && split_bound(inh, y, mu) <= logc_lo) return 1;
        for (int q = 0; q < nlib; q++) if (split_bound(&lib[q], y, mu) <= logc_lo) return 1;
    }
    if (has && (flag_bound(inh, b, y, mu) <= logc_lo || beats_identity(inh, b))) { hint_flag = *inh; hint_has = 1; return 1; }
    for (int q = 0; q < nlib; q++) if (flag_bound(&lib[q], b, y, mu) <= logc_lo || beats_identity(&lib[q], b)) { hint_flag = lib[q]; hint_has = 1; return 1; }
    center_flag(b, f);
    if (flag_bound(f, b, y, mu) <= logc_lo || beats_identity(f, b)) {
        hint_flag = *f; hint_has = 1;
        lib[libpos] = *f; libpos = (libpos + 1) % LIB; if (nlib < LIB) nlib++; return 1;
    }
    if (use_split) {
        Flag g;
        for (int q = 0; q < nkbest; q++) {
            memcpy(g.u, kbestflag[q], sizeof(g.u));
            if (split_bound(&g, y, mu) <= logc_lo) { hint_flag = g; hint_has = 1; lib[libpos] = g; libpos = (libpos + 1) % LIB; if (nlib < LIB) nlib++; return 1; }
        }
    }
    if (use_comb) {
        Flag cand[CMAX];
        int nc = 0;
        for (int q = 0; q < nkbest && nc < CMAX - 2; q++) memcpy(cand[nc++].u, kbestflag[q], sizeof(kbestflag[q]));
        for (int i = 0; i < n; i++) for (int j = 0; j < n; j++) cand[nc].u[i][j] = (i == j);
        nc++;
        if (has) cand[nc++] = *inh;
        for (int q = 0; q < nlib && nc < CMAX; q++) cand[nc++] = lib[q];
        if (lp_comb_bound(cand, nc, b) <= logc_lo) return 1;
    }
    for (int s_ = 0; s_ < nball; s_++) if (in_ball(y, mu, s_) || (use_invball && !getenv("NOMVBALL") && in_invball_mv(b, s_))) return 1;
    log_psi2_centered(f, b);          /* fills lastcontrib for sensitivity splitting */
    return 0;
}

static int split_dim(const Box *b) {
    int ts = 0; double wmax = -1;
    for (int t = 0; t < DIM; t++) {
        double w = (b->hi[t] - b->lo[t]) * (t < n - 1 ? 0.5 : 1.0);
        if (use_sens && combvalid && isfinite(combcontrib[t])) w = combcontrib[t] + 1e-3 * w;
        else if (use_sens && isfinite(lastcontrib[t])) w = lastcontrib[t] + 1e-3 * w;   /* (B) sensitivity */
        if (w > wmax) { wmax = w; ts = t; }
    }
    return ts;
}

/* TARGET=j (with INVBALL): dive towards ball centre j instead of randomly, and report the depth */
static int target_on = 0;
static double TARGET[DMAX];
static void set_target(int j) {
    /* Cholesky of BALLG[j] -> y_k = log G_k, mu */
    double L[NMAX][NMAX] = {{0}}, A[NMAX][NMAX];
    memcpy(A, BALLG[j], sizeof(A));
    for (int c = 0; c < n; c++) {
        double s_ = A[c][c]; for (int k = 0; k < c; k++) s_ -= L[c][k] * L[c][k];
        L[c][c] = sqrt(s_);
        for (int i = c + 1; i < n; i++) { double t = A[i][c]; for (int k = 0; k < c; k++) t -= L[i][k] * L[c][k]; L[i][c] = t / L[c][c]; }
    }
    double y = 0;
    for (int k = 0; k < n - 1; k++) { y += 2 * log(L[k][k]); TARGET[k] = y; }
    for (int t = 0; t < nmu; t++) TARGET[n - 1 + t] = L[MUI[t]][MUJ[t]] / L[MUJ[t]][MUJ[t]];
    target_on = 1;
    fprintf(stderr, "target ball %d:", j); for (int t = 0; t < DIM; t++) fprintf(stderr, " %.4f", TARGET[t]); fprintf(stderr, "\n");
}

/* Knuth's estimator of the branch-and-bound tree size: random root-to-leaf dives */
static long dumpdeep = 0;
static void estimate(Box b0, int dives) {
    double sum = 0, sum2 = 0, csum = 0;
    long maxdepth = 0;
    srand(12345);
    for (int d = 0; d < dives; d++) {
        Box b = b0; Flag inh, f; int has = 0;
        double weight = 1, est = 0, cest = 0;
        long depth = 0;
        while (1) {
            est += weight;
            if (classify(&b, &inh, has, &f)) {
                if (dumpdeep && depth >= dumpdeep) {   /* print the leaf's centre Gram for offline study */
                    double y[NMAX + 1], r[NMAX], M[NMAX][NMAX] = {{0}};
                    y[0] = 0; y[n] = 0;
                    for (int k = 0; k < n - 1; k++) y[k + 1] = 0.5 * (b.lo[k] + b.hi[k]);
                    for (int i = 0; i < n; i++) { r[i] = exp(y[i + 1] - y[i]); M[i][i] = 1; }
                    for (int t = 0; t < nmu; t++) M[MUI[t]][MUJ[t]] = 0.5 * (b.lo[n - 1 + t] + b.hi[n - 1 + t]);
                    printf("DEEP %ld", depth);
                    for (int i = 0; i < n; i++) for (int j = 0; j < n; j++) {
                        double g = 0;
                        for (int k = 0; k < n; k++) g += M[i][k] * r[k] * M[j][k];
                        printf(" %.10g", g);
                    }
                    printf("\n");
                }
                break;
            }
            if (collect_w > 0) {      /* a collected box is a leaf of the tree (COLLECT width or HOPELESS) */
                double wm = 0; for (int t = 0; t < DIM; t++) wm = fmax(wm, b.hi[t] - b.lo[t]);
                if (wm < collect_w || last_hopeless) { cest = weight; break; }
            }
            int ts = split_dim(&b);
            double mid = 0.5 * (b.lo[ts] + b.hi[ts]);
            if (target_on) { if (TARGET[ts] < mid) b.hi[ts] = mid; else b.lo[ts] = mid; }
            else if (rand() & 1) b.hi[ts] = mid; else b.lo[ts] = mid;
            inh = f; has = 1;
            weight *= 2;
            if (++depth > 400) { est = INFINITY; break; }
        }
        if (depth > maxdepth) maxdepth = depth;
        sum += est; sum2 += est * est; csum += cest;
        if ((d + 1) % 100 == 0 || d + 1 == dives) {
            double m = sum / (d + 1), sd = sqrt(fabs(sum2 / (d + 1) - m * m) / (d + 1));
            printf("  dives=%d  estimated tree size %.3g  (+- %.2g)  collected %.3g  max depth %ld\n", d + 1, m, sd, csum / (d + 1), maxdepth);
            fflush(stdout);
        }
    }
}


/* diagnostic: one random dive, printing why each box fails */
static void diagnose(Box b0, int seed) {
    srand(seed);
    Box b = b0; Flag inh, f; int has = 0;
    for (int depth = 0; depth < 80; depth++) {
        if (classify(&b, &inh, has, &f)) { printf("depth %d: resolved\n", depth); return; }
        I y[NMAX + 1], mu[NMAX][NMAX];
        unpack(&b, y, mu);
        double naive = log_psi2_upper(&f, y, mu);
        double cen = log_psi2_centered(&f, &b);
        double gsum = cen - lastfc;
        double wmin = 1e9, wmax = 0;
        for (int t = 0; t < DIM; t++) { double w = b.hi[t] - b.lo[t]; if (w < wmin) wmin = w; if (w > wmax) wmax = w; }
        if (depth % 4 == 0 || depth > 40) {
            printf("depth %2d: logc=%.4f  fc(center, best flag)=%.4f  naive=%.4f  centred=%.4f (grad terms %.4f)  width[min,max]=[%.2e,%.2e]\n",
                   depth, logc_lo, lastfc, naive, cen, gsum, wmin, wmax);
            printf("          grad-term per coord:");
            for (int t = 0; t < DIM; t++) printf(" %.3f", lastcontrib[t]);
            printf("\n          grad interval widths: ");
            for (int t = 0; t < DIM; t++) printf(" %.2f", lastgrad[t].hi - lastgrad[t].lo);
            printf("\n");
        }
        {
            I lr[NMAX]; for (int i = 0; i < n; i++) lr[i] = isub(y[i + 1], y[i]);
            potbest = INFINITY; pot_beaten(mu, lr); double pb = potbest;
            potbest = INFINITY; pot_beaten_dual(mu, lr); double pd = potbest;
            Box pt; for (int t = 0; t < DIM; t++) pt.lo[t] = pt.hi[t] = 0.5 * (b.lo[t] + b.hi[t]);
            I py[NMAX + 1], pmu[NMAX][NMAX]; unpack(&pt, py, pmu); I plr[NMAX];
            for (int i = 0; i < n; i++) plr[i] = isub(py[i + 1], py[i]);
            potbest = INFINITY; pot_beaten(pmu, plr); double qb = potbest;
            potbest = INFINITY; pot_beaten_dual(pmu, plr); double qd = potbest;
            printf("   pot box primal %.3f dual %.3f | centre primal %.3f dual %.3f\n", pb, pd, qb, qd);
        }
        int ts = split_dim(&b);
        double mid = 0.5 * (b.lo[ts] + b.hi[ts]);
        if (rand() & 1) b.hi[ts] = mid; else b.lo[ts] = mid;
        inh = f; has = 1;
    }
}

int main(int argc, char **argv) {
    if (argc < 3) { fprintf(stderr, "usage: %s n c [max_boxes] [balls_file]\n", argv[0]); return 1; }
    n = atoi(argv[1]);
    double c = atof(argv[2]);
    long max_boxes = argc > 3 ? atol(argv[3]) : 100000000L;
    if (argc > 4 && strcmp(argv[4], "none") != 0) load_balls(argv[4]);
    if (getenv("SENS")) use_sens = 1;
    if (getenv("NOCOMB")) use_comb = 0;
    if (getenv("NOSPLIT")) use_split = 0;
    if (getenv("NOELEM")) use_elem = 0;
    if (getenv("NOPOT")) use_pot = 0;
    if (getenv("NODUAL")) use_dualcut = 0;
    if (getenv("HOPELESS")) { hope_h = atof(getenv("HOPELESS")); hope_on = 1; }
    if (getenv("NOTANGENT")) use_tangent = 0;
    if (getenv("NOELEMCG")) use_elem_cg = 0;
    if (getenv("COLLECT")) { collect_w = atof(getenv("COLLECT")); collect_f = fopen(getenv("COLLECT_FILE") ? getenv("COLLECT_FILE") : "S.txt", "w"); }
    if (getenv("SREGION")) load_S(getenv("SREGION"));
    if (getenv("SHORTEST")) shortest_mode = 1;
    if (getenv("INVBALL")) use_invball = 1;
    if (getenv("PHI") && *getenv("PHI")) {
        /* PHI=a/b: the double just below or at a/b (rounded down), so every leaf proves Phi < a/b strictly */
        phi_mode = 1;
        const char *ps = getenv("PHI"), *sl = strchr(ps, '/');
        if (sl) { volatile double a = atof(ps), d = atof(sl + 1); fesetround(FE_DOWNWARD); PHI_C = a / d; fesetround(FE_TONEAREST); }
        else PHI_C = atof(ps);
    }
    if (getenv("DUMPDEEP")) dumpdeep = atol(getenv("DUMPDEEP"));
    if (getenv("COUPLED")) { use_coupled = 1; atexit(cpl_report); }
    if (getenv("CPCELL0")) cp_cell0 = atof(getenv("CPCELL0"));
    if (getenv("CPCELLMIN")) cp_cellmin = atof(getenv("CPCELLMIN"));
    logc_lo = dn(log(c));
    nmu = 0;
    for (int i = 0; i < n; i++) for (int j = 0; j < i; j++) { MUI[nmu] = i; MUJ[nmu] = j; nmu++; }
    DIM = n - 1 + nmu;
    gen_subsets();
    gen_elem();
    gen_shortv(); atexit(sv_report);
    if (phi_mode) { if (n != 4) { fprintf(stderr, "PHI: n = 4 only\n"); return 1; } gen_phi_vectors(); use_dualcut = 0; }
    /* log mu_L upper bounds: mu_1 = 1, mu_2..mu_4 = gamma (Theorem K), mu_5 from (K) with eta_5 <= 3 2^(5/4) */
    LOGMU_HI[1] = 0;
    LOGMU_HI[2] = up(up(0.5 * log(4.0 / 3.0)));
    LOGMU_HI[3] = up(up(log(2.0) / 3));
    LOGMU_HI[4] = up(up(0.5 * log(2.0)));
    LOGMU_HI[5] = up(up((log(3.0) + 1.25 * log(2.0) + log(4.0 / 3.0) + log(8.0)) / 10));
    if (n > 5) { fprintf(stderr, "n <= 5 only\n"); return 1; }
    if (n == 5 && !getenv("NOETACUT")) ETA_N_HI = up(up(log(3.0) + 1.25 * log(2.0)));
    LOGGAM[1] = 0; LOGGAM[2] = LOGMU_HI[2]; LOGGAM[3] = LOGMU_HI[3]; LOGGAM[4] = LOGMU_HI[4];
    LOGGAM[5] = up(up(log(8.0) / 5));
    if (shortest_mode && n != 5) { fprintf(stderr, "SHORTEST: n = 5 only\n"); return 1; }
    LOGETA_HI[1] = 0;
    LOGETA_HI[2] = up(up(log(2.0 / sqrt(3.0)) + 1e-15));
    LOGETA_HI[3] = up(up(log(1.5)));
    LOGETA_HI[4] = up(up(log(3.0)));
    /* initial box */
    Box b0;
    double hi[NMAX];
    for (int k = 1; k < n; k++) {
        double lll = up(up(k * (n - k) / 2.0 * log(4.0 / 3.0)));
        hi[k - 1] = (k == 1 || k == n - 1) ? (LOGMU_HI[n] < lll ? LOGMU_HI[n] : lll) : lll;
        if (shortest_mode) hi[k - 1] = (k == 1) ? LOGGAM[n] : lll;
    }
    double sumhi = 0;
    for (int k = 0; k < n - 1; k++) sumhi = up(sumhi + hi[k]);
    for (int k = 0; k < n - 1; k++) { b0.lo[k] = dn(logc_lo - dn(sumhi - hi[k]) - 1e-12); b0.hi[k] = hi[k]; }
    /* eta <= eta_{n-1} gamma(L)^{n/2} (split, k = 1): r_1 > (c / eta_{n-1})^(2/n) */
    if (shortest_mode) { double l1 = dn(dn(2.0 / n * (logc_lo - LOGETA_HI[n - 1]))); if (l1 > b0.lo[0]) b0.lo[0] = l1; }
    for (int t = 0; t < nmu; t++) {
        b0.lo[n - 1 + t] = (MUI[t] == MUJ[t] + 1) ? 0.0 : -0.5;
        b0.hi[n - 1 + t] = 0.5;
    }
    vol0 = box_vol(&b0);
    double resolved = 0;
    /* optional parallel slicing: argv[5] = "k/K" processes slice k of K (split the first y-coordinate) */
    if (argc > 5) {
        int k, K;
        if (sscanf(argv[5], "%d/%d", &k, &K) == 2 && K > 1) {
            double w = (b0.hi[0] - b0.lo[0]) / K;
            double lo = b0.lo[0] + k * w, hi = (k == K - 1) ? b0.hi[0] : b0.lo[0] + (k + 1) * w;
            b0.lo[0] = lo; b0.hi[0] = hi;
            /* "k/K/j/J": sub-slice j of J of slice k; adjacent sub-slices share the same computed endpoint */
            int j, J;
            if (sscanf(argv[5], "%*d/%*d/%d/%d", &j, &J) == 2 && J > 1) {
                double w2 = (hi - lo) / J;
                if (j > 0) b0.lo[0] = lo + j * w2;
                if (j < J - 1) b0.hi[0] = lo + (j + 1) * w2;
            }
        }
    }
    if (getenv("ROOTBOX")) {        /* "lo_0 hi_0 ... lo_8 hi_8" (hex floats, as COLLECT writes them): start from this box */
        char *p = getenv("ROOTBOX");
        for (int t = 0; t < DIM; t++) { b0.lo[t] = strtod(p, &p); b0.hi[t] = strtod(p, &p); }
    }
    double myvol = box_vol(&b0);
    if (getenv("TARGET")) { set_target(atoi(getenv("TARGET"))); if (phi_mode && !getenv("KEEPBALLS")) nball = 0; }   /* balls file only supplies the target (KEEPBALLS: also excluded) */
    if (target_on && getenv("NEAR")) {      /* restrict the initial box to a neighbourhood of the target */
        double w = atof(getenv("NEAR"));
        for (int t = 0; t < DIM; t++) {
            double lo = TARGET[t] - w, hi = TARGET[t] + w;
            if (lo > b0.lo[t]) b0.lo[t] = lo;
            if (hi < b0.hi[t]) b0.hi[t] = hi;
        }
        fprintf(stderr, "NEAR %.3f: box restricted around target\n", w);
        target_on = 0;          /* random dives inside the restricted box */
    }
    if (argc > 6 && strncmp(argv[6], "estimate:", 9) == 0) { estimate(b0, atoi(argv[6] + 9)); return 0; }
    if (argc > 6 && strncmp(argv[6], "diag:", 5) == 0) { diagnose(b0, atoi(argv[6] + 5)); return 0; }
    long cap = 1 << 20;
    Item *stack = malloc(sizeof(Item) * cap);
    long sp = 0;
    stack[sp].b = b0; stack[sp].has = 0; sp++;
    long collected = 0;
    long processed = 0, certified = 0, discarded = 0, balled = 0, newflags = 0;
    int gave_up = 0;
    clock_t t0 = clock();
    /* CERTLOG=dir: write the search tree for an independent checker.  Pre-order: the root box (hex floats), then per
       node "S t" (split coordinate t at mid = 0.5*(lo+hi), round-to-nearest; the upper half comes first), "L k"
       (leaf, with the kappa bound used, 0 if not in PHI mode, then optionally the 16 entries of the flag that closed it,
       a hint only) or "C" (COLLECT: unresolved box handed to stage 2). */
    FILE *clog = NULL;
    if (getenv("CERTLOG")) {
        char fn[512], tag[64] = "all";
        if (argc > 5) { snprintf(tag, sizeof(tag), "%s", argv[5]); for (char *q = tag; *q; q++) if (*q == '/') *q = '_'; }
        snprintf(fn, sizeof(fn), "%s/slice_%s.log", getenv("CERTLOG"), tag);
        clog = fopen(fn, "w");
        if (!clog) { fprintf(stderr, "cannot open %s\n", fn); return 1; }
        fprintf(clog, "ROOT n=%d dim=%d phi=%a c=%.17g ver=%s\n", n, DIM, phi_mode ? PHI_C : 0.0, c, GITREV);
        for (int t = 0; t < DIM; t++) fprintf(clog, "%a %a\n", b0.lo[t], b0.hi[t]);
    }
    while (sp > 0) {
        Item it = stack[--sp];
        if (++processed > max_boxes) { gave_up = 1; break; }
        Flag f;
        if (classify(&it.b, &it.inh, it.has, &f)) {
            certified++; resolved += box_vol(&it.b);
            if (clog) {
                if (phi_mode) fprintf(clog, "L %a", last_kappa); else fputs("L", clog);
                /* closed without a single flag (LP combination, split bound, ball): hint the best flag at the centre,
                   which check_phi's own descent may not find near M*, where many flags nearly tie */
                if (!hint_has && phi_mode && n == 4) { center_flag(&it.b, &hint_flag); hint_has = 1; }
                if (hint_has && n == 4) for (int i = 0; i < n; i++) for (int j = 0; j < n; j++) fprintf(clog, " %d", hint_flag.u[i][j]);
                fputs("\n", clog);
            }
            continue;
        }
        if (collect_w > 0) {
            double wm = 0;
            for (int t = 0; t < DIM; t++) wm = fmax(wm, it.b.hi[t] - it.b.lo[t]);
            if (wm < collect_w || last_hopeless) {      /* handed to stage 2: exact (hex) box, and a "C" leaf in the certificate log */
                for (int t = 0; t < DIM; t++) fprintf(collect_f, "%a %a ", it.b.lo[t], it.b.hi[t]);
                fprintf(collect_f, "\n"); collected++; resolved += box_vol(&it.b);
                if (clog) fputs("C\n", clog);
                continue;
            }
        }
        newflags++;
        int ts = split_dim(&it.b);
        double mid = 0.5 * (it.b.lo[ts] + it.b.hi[ts]);
        if (clog) fprintf(clog, "S %d\n", ts);
        if (sp + 2 >= cap) { cap *= 2; stack = realloc(stack, sizeof(Item) * cap); }
        Item a = it, bb = it;
        a.b.hi[ts] = mid; bb.b.lo[ts] = mid;
        a.inh = f; bb.inh = f; a.has = bb.has = 1;
        stack[sp++] = a; stack[sp++] = bb;
        if (processed % (getenv("PROGRESS") ? atol(getenv("PROGRESS")) : 200000) == 0)
            fprintf(stderr, "  processed=%ld certified=%ld discarded=%ld balled=%ld newflags=%ld stack=%ld "
                    "resolved=%.6f%% %.0fs\n", processed, certified, discarded, balled, newflags, sp,
                    100 * resolved / myvol, (double)(clock() - t0) / CLOCKS_PER_SEC);
    }
    if (collect_w > 0) { fclose(collect_f); fprintf(stderr, "collected %ld boxes\n", collected); }
    if (clog) { fprintf(clog, "END %s\n", (sp == 0 && !gave_up) ? "complete" : "incomplete"); fclose(clog); }
    int ok = (sp == 0) && !gave_up;
    if (ok && collected > 0) printf("COMPLETE: %ld boxes handed to stage 2\n", collected);
    printf("n=%d c=%g: %s  processed=%ld certified=%ld discarded=%ld balled=%ld newflags=%ld resolved=%.6f%%  %.1fs\n",
           n, c, ok ? (collected ? "COMPLETE (collected boxes, see above)" : "CERTIFIED eta_n <= c") : "NOT certified", processed, certified, discarded, balled, newflags,
           100 * resolved / myvol, (double)(clock() - t0) / CLOCKS_PER_SEC);
    return ok ? 0 : 3;
}
