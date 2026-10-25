/* check_phi.c -- independent checker for the Phi <= PHI_C certificate of certify.c (PHI mode, n = 4).
 *
 *   gcc -O2 -frounding-math -o check_phi check_phi.c -lm
 *   ./check_phi slice_k_K.log [more logs...]
 *   MSTAR=mstar_U.txt ./check_phi ...   (also accept boxes inside the balls of Proposition [M*], see load_family)
 *   PVD5=pivd5_U.txt ./check_phi ...    (the pi_v D5 balls, radius 1/20, for the PHI = 6 run)
 *
 * A log may also contain "C" leaves (certify COLLECT mode): boxes handed to a later stage.  They are counted, and with
 * STAGE2=file written there exactly; together with the proved leaves they tile the root box by construction.
 *
 * Written separately from certify.c.  It shares only the coordinates and the search tree (read from the log), and
 * re-proves every leaf with its own arithmetic and its own bounds, subdividing a leaf further when its weaker bounds
 * need it.  Nothing about how certify.c decided a leaf is trusted.
 *
 * Coordinates: y_1..y_3 (y_k = log G_k, y_0 = y_4 = 0, so det = 1) and mu_10, mu_20, mu_21, mu_30, mu_31, mu_32 of a
 * Psi-minimal, size-reduced, sign-normalised basis of a 4-dimensional lattice M (every lattice has one).
 * A leaf box is accepted if
 *   (D) no Psi-minimal basis lies in it: a window bound, the window potential bound, the Lagrange bound, a deep
 *       insertion, or an elementary basis change U beating the identity everywhere; or
 *   (C) eta(M) <= E on the box (tangent bound of a concave function, see below) and kappa(M) < K := sqrt(PHI_C / E)
 *       on the box (branch-and-bound over the phi-torus).  Then Phi(M) = kappa^2 eta < PHI_C.
 * The proved inputs: mu_2 = gamma_2, mu_3 = gamma_3, mu_4 = gamma_4 (window bounds, Theorem K and eta_4 = 3),
 * eta_2 = 2/sqrt 3, eta_3 = 3/2 (window potentials, split lemma), and the LLL bound G_k <= (4/3)^{k(4-k)/2}.
 *
 * Arithmetic: the FPU stays in round-upward mode; lower ends are computed as -(upward(-x)).  exp/log are CORE-MATH's
 * correctly rounded cr_exp/cr_log (third_party/core-math), widened by one ulp; sqrt is IEEE correctly rounded (pow is not
 * used).  Splits are replayed in round-to-nearest, exactly as certify.c computed them (fesetround around the midpoint
 * only).  Build: gcc -O2 -frounding-math -ffp-contract=off (baseline x86-64, no FMA).
 *
 * eta bound: for a flag U (rows = new basis vectors), log det(U_k G U_k^T) is concave in G, so with G* the box
 * centre, f(G) <= f(G*) + <grad f(G*), G - G*>; f(G*) and grad f(G*) are enclosed in intervals, G - G* by the box.
 */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <fenv.h>
#ifndef GITREV
#define GITREV "unknown"
#endif

#define N 4
#define D 9

/* ---------------- interval arithmetic (FE_UPWARD throughout) ---------------- */
typedef struct { double lo, hi; } J;
static J jv(double a, double b) { J r = {a, b}; return r; }
static J jp(double a) { return jv(a, a); }
static J jadd(J a, J b) { return jv(-((-a.lo) - b.lo), a.hi + b.hi); }
static J jsub(J a, J b) { return jv(-(b.hi - a.lo), a.hi - b.lo); }
static J jmul(J a, J b) {
    double h1 = a.lo * b.lo, h2 = a.lo * b.hi, h3 = a.hi * b.lo, h4 = a.hi * b.hi;
    double l1 = -((-a.lo) * b.lo), l2 = -((-a.lo) * b.hi), l3 = -((-a.hi) * b.lo), l4 = -((-a.hi) * b.hi);
    return jv(fmin(fmin(l1, l2), fmin(l3, l4)), fmax(fmax(h1, h2), fmax(h3, h4)));
}
static J jsqr(J a) {
    if (a.lo >= 0) return jv(-((-a.lo) * a.lo), a.hi * a.hi);
    if (a.hi <= 0) return jv(-((-a.hi) * a.hi), a.lo * a.lo);
    double m = fmax(-a.lo, a.hi);
    return jv(0, m * m);
}
static J jdiv(J a, J b) {      /* b > 0 */
    if (!(b.lo > 0)) { fprintf(stderr, "jdiv: divisor not positive\n"); exit(3); }
    double h1 = a.lo / b.lo, h2 = a.lo / b.hi, h3 = a.hi / b.lo, h4 = a.hi / b.hi;
    double l1 = -((-a.lo) / b.lo), l2 = -((-a.lo) / b.hi), l3 = -((-a.hi) / b.lo), l4 = -((-a.hi) / b.hi);
    return jv(fmin(fmin(l1, l2), fmin(l3, l4)), fmax(fmax(h1, h2), fmax(h3, h4)));
}
/* exp and log: CORE-MATH's correctly rounded cr_exp, cr_log (third_party/core-math), correct in every rounding mode, so
   each result is within one ulp of the true value and one nextafter in each direction encloses it */
double cr_exp(double), cr_log(double);
static double up1(double x) { return nextafter(x, INFINITY); }
static double dn1(double x) { return nextafter(x, -INFINITY); }
static J jexp(J a) { return jv(fmax(0, dn1(cr_exp(a.lo))), up1(cr_exp(a.hi))); }
static J jlog(J a) { if (!(a.lo > 0)) return jv(-INFINITY, up1(cr_log(a.hi))); return jv(dn1(cr_log(a.lo)), up1(cr_log(a.hi))); }
static J jsqrt(J a) { return jv(a.lo > 0 ? nextafter(sqrt(a.lo), -INFINITY) : 0, sqrt(a.hi)); }
static J jscale(J a, double s) { return jmul(a, jp(s)); }

/* ---------------- proved constants (as intervals) ---------------- */
static J LOGMU[N + 1], LOGETA[N + 1], LOGPHI;

static void constants(double phi) {
    /* log gamma_2 = 1/2 log(4/3), log gamma_3 = 1/3 log 2, log gamma_4 = 1/2 log 2 */
    LOGMU[1] = jp(0);
    LOGMU[2] = jscale(jlog(jdiv(jp(4), jp(3))), 0.5);
    LOGMU[3] = jdiv(jlog(jp(2)), jp(3));
    LOGMU[4] = jscale(jlog(jp(2)), 0.5);
    LOGETA[1] = jp(0);
    LOGETA[2] = jsub(jlog(jp(2)), jscale(jlog(jp(3)), 0.5));   /* log(2/sqrt 3) */
    LOGETA[3] = jsub(jlog(jp(3)), jlog(jp(2)));                /* log(3/2) */
    LOGETA[4] = jlog(jp(3));
    LOGPHI = jlog(jp(phi));
}

/* ---------------- box -> intervals ---------------- */
static const int MI[6] = {1, 2, 2, 3, 3, 3}, MJ[6] = {0, 0, 1, 0, 1, 2};
typedef struct { double lo[D], hi[D]; } Box;
typedef struct { J y[N + 1], lr[N], r[N], mu[N][N], G[N][N]; } Data;

static void fill(const Box *b, Data *d) {
    d->y[0] = jp(0); d->y[N] = jp(0);
    for (int k = 1; k < N; k++) d->y[k] = jv(b->lo[k - 1], b->hi[k - 1]);
    for (int i = 0; i < N; i++) for (int j = 0; j < N; j++) d->mu[i][j] = jp(i == j ? 1 : 0);
    for (int t = 0; t < 6; t++) d->mu[MI[t]][MJ[t]] = jv(b->lo[3 + t], b->hi[3 + t]);
    for (int i = 0; i < N; i++) { d->lr[i] = jsub(d->y[i + 1], d->y[i]); d->r[i] = jexp(d->lr[i]); }
    for (int i = 0; i < N; i++) for (int j = 0; j <= i; j++) {
        J acc = jp(0);
        for (int k = 0; k <= j; k++) acc = jadd(acc, jmul(jmul(d->mu[i][k], d->mu[j][k]), d->r[k]));
        d->G[i][j] = d->G[j][i] = acc;
    }
}

/* ---------------- discards ---------------- */
static long n_disc[6];
/* deep insertion test on Gram-Schmidt data (mu unit lower triangular, lr = log r): 1 if some insertion of b_k at
   position i < k lowers Psi^2 everywhere on the box */
static int insertion_beats(J mu[N][N], const J *lr) {
    for (int k = 1; k < N; k++) for (int i = 0; i < k; i++) {
        J tot = jp(0);
        int ok = 1;
        for (int j = i; j < k; j++) {
            J a = jsqr(mu[k][j]);
            for (int l = j + 1; l <= k; l++) a = jadd(a, jmul(l == k ? jp(1) : jsqr(mu[k][l]), jexp(jsub(lr[l], lr[j]))));
            if (!(a.lo > 0)) { ok = 0; break; }
            tot = jadd(tot, jlog(a));
        }
        if (ok && tot.hi < 0) return 1;
    }
    return 0;
}
static int discard(const Data *d) {
    /* window bounds: in a Psi-minimal basis each window i..j (length L) satisfies
       r_i <= mu_L (G_j/G_{i-1})^{1/L} and r_j >= (G_j/G_{i-1})^{1/L} / mu_L */
    for (int L = 2; L <= N; L++) for (int i = 0; i + L - 1 < N; i++) {
        int j = i + L - 1;
        J blk = jdiv(jsub(d->y[j + 1], d->y[i]), jp(L));
        if (jsub(jsub(d->lr[i], LOGMU[L]), blk).lo > 0) { n_disc[0]++; return 1; }
        if (jsub(jsub(blk, LOGMU[L]), d->lr[j]).lo > 0) { n_disc[0]++; return 1; }
    }
    /* window potentials: a window of length L < N is Psi-minimal, so its normalised potential is <= log eta_L */
    for (int L = 2; L < N; L++) for (int i = 0; i + L - 1 < N; i++) {
        int j = i + L - 1;
        J span = jsub(d->y[j + 1], d->y[i]), tot = jp(0);
        for (int k = i; k < j; k++)
            tot = jadd(tot, jsub(jsub(d->y[k + 1], d->y[i]), jdiv(jscale(span, k - i + 1), jp(L))));
        if (tot.lo > LOGETA[L].hi) { n_disc[1]++; return 1; }
    }
    /* Lagrange: r_{i+1} >= (1 - mu_{i+1,i}^2) r_i */
    for (int i = 0; i + 1 < N; i++) {
        J om = jsub(jp(1), jsqr(d->mu[i + 1][i]));
        if (om.lo > 0 && jsub(d->lr[i + 1], d->lr[i]).hi < jlog(om).lo) { n_disc[2]++; return 1; }
    }
    /* deep insertion of b_k at position i multiplies Psi^2 by prod_{j=i}^{k-1} (mu_kj^2 + sum_{l=j+1}^{k} mu_kl^2 r_l/r_j) */
    J mu[N][N];
    memcpy(mu, d->mu, sizeof(mu));
    if (insertion_beats(mu, d->lr)) { n_disc[3]++; return 1; }
    /* the same on the reversed dual basis (Psi^2 of the reversed dual basis is Psi^2 / det^{...}, so a dual insertion
       that lowers it lowers Psi^2 of the corresponding primal basis): with Minv = M^{-1} (unit lower triangular),
       mu'_{ab} = Minv[n-1-b][n-1-a], log r'_a = -log r_{n-1-a} */
    J Mi[N][N];
    for (int i = 0; i < N; i++) {
        for (int j = 0; j < N; j++) Mi[i][j] = jp(i == j ? 1 : 0);
        for (int j = i - 1; j >= 0; j--) {
            J acc = jp(0);
            for (int t = j; t < i; t++) acc = jadd(acc, jmul(d->mu[i][t], Mi[t][j]));
            Mi[i][j] = jsub(jp(0), acc);
        }
    }
    J mud[N][N], lrd[N];
    for (int a = 0; a < N; a++) {
        lrd[a] = jsub(jp(0), d->lr[N - 1 - a]);
        for (int b = 0; b < N; b++) mud[a][b] = (b <= a) ? Mi[N - 1 - b][N - 1 - a] : jp(0);
    }
    if (insertion_beats(mud, lrd)) { n_disc[5]++; return 1; }
    return 0;
}

/* ---------------- eta bound via concavity ---------------- */
typedef struct { int u[N][N]; } Flag;        /* rows = new basis vectors in old coordinates */

/* interval Cholesky of a k x k matrix; returns 0 if positivity cannot be shown */
static int jchol(int k, J A[N][N], J L[N][N]) {
    for (int c = 0; c < k; c++) {
        J s = A[c][c];
        for (int q = 0; q < c; q++) s = jsub(s, jsqr(L[c][q]));
        if (!(s.lo > 0)) return 0;
        L[c][c] = jsqrt(s);
        for (int i = c + 1; i < k; i++) {
            J t = A[i][c];
            for (int q = 0; q < c; q++) t = jsub(t, jmul(L[i][q], L[c][q]));
            L[i][c] = jdiv(t, L[c][c]);
        }
    }
    return 1;
}

static J GC[N][N];      /* G at the box centre (tight intervals) */
static J DG[N][N];      /* G - G* over the box */
static const Data *BX;  /* the box's interval data (mean-value bound) */
static J DX[D];         /* x - x* over the box, per coordinate */
static int YSUB = 0;    /* 1: bound f(G(x)) - (y_1 + y_2 + y_3)(x) instead of f(G(x)) (comparisons with the identity) */
static J SY_BOX, SY_C;  /* y_1 + y_2 + y_3 over the box and at the centre */

/* <P, G(x) - G(x*)> over the box by the mean value theorem: F(x) = tr(P M R M^T) has
   dF/dmu_ab = 2 (P M R)_ab (a > b) and dF/dy_s = r_{s-1} (M^T P M)_{s-1,s-1} - r_s (M^T P M)_{ss};
   the gradient is enclosed over the whole box, so the bound loses only second-order terms */
static J mean_value_term(J P[N][N]) {
    J PM[N][N], Q[N];
    for (int i = 0; i < N; i++) for (int k = 0; k < N; k++) {
        J s = jp(0);
        for (int j = k; j < N; j++) s = jadd(s, jmul(P[i][j], BX->mu[j][k]));
        PM[i][k] = s;
    }
    for (int k = 0; k < N; k++) {
        J s = jp(0);
        for (int i = k; i < N; i++) s = jadd(s, jmul(BX->mu[i][k], PM[i][k]));
        Q[k] = s;
    }
    J tot = jp(0);
    for (int s1 = 1; s1 < N; s1++) {
        J g = jsub(jmul(BX->r[s1 - 1], Q[s1 - 1]), jmul(BX->r[s1], Q[s1]));
        if (YSUB) g = jsub(g, jp(1));          /* d/dy_s of -(y_1 + y_2 + y_3) */
        tot = jadd(tot, jmul(g, DX[s1 - 1]));
    }
    for (int t = 0; t < 6; t++) {
        J g = jscale(jmul(PM[MI[t]][MJ[t]], BX->r[MJ[t]]), 2);
        tot = jadd(tot, jmul(g, DX[3 + t]));
    }
    return tot;
}

/* upper bound over the box of  const + sum_k w_k log det(U_k G U_k^T); returns +inf on failure */
static double tangent_bound(const Flag *U, const double *w, J cst) {
    J A[N][N];
    for (int a = 0; a < N; a++) for (int b = 0; b < N; b++) {
        J s = jp(0);
        for (int i = 0; i < N; i++) for (int j = 0; j < N; j++)
            if (U->u[a][i] && U->u[b][j]) s = jadd(s, jscale(GC[i][j], U->u[a][i] * U->u[b][j]));
        A[a][b] = s;
    }
    J val = cst, P[N][N];
    for (int i = 0; i < N; i++) for (int j = 0; j < N; j++) P[i][j] = jp(0);
    for (int k = 1; k < N; k++) {
        if (w[k] == 0) continue;
        J L[N][N];
        if (!jchol(k, A, L)) return INFINITY;
        J ld = jp(0);
        for (int c = 0; c < k; c++) ld = jadd(ld, jscale(jlog(L[c][c]), 2));
        val = jadd(val, jscale(ld, w[k]));
        /* A_k^{-1} = L^{-T} L^{-1}; W = U_k^T L^{-T}, so U_k^T A_k^{-1} U_k = W W^T, W[i][c] = (L^{-1} U_k)[c][i] */
        J Wt[N][N];      /* Wt[c][i] = (L^{-1} U_k)[c][i] */
        for (int i = 0; i < N; i++) for (int c = 0; c < k; c++) {
            J s = jp(U->u[c][i]);
            for (int q = 0; q < c; q++) s = jsub(s, jmul(L[c][q], Wt[q][i]));
            Wt[c][i] = jdiv(s, L[c][c]);
        }
        for (int i = 0; i < N; i++) for (int j = 0; j < N; j++) {
            J s = jp(0);
            for (int c = 0; c < k; c++) s = jadd(s, jmul(Wt[c][i], Wt[c][j]));
            P[i][j] = jadd(P[i][j], jscale(s, w[k]));
        }
    }
    J lin = jp(0);
    for (int i = 0; i < N; i++) for (int j = 0; j < N; j++) lin = jadd(lin, jmul(P[i][j], DG[i][j]));
    if (YSUB) lin = jsub(lin, jsub(SY_BOX, SY_C));
    J mv = mean_value_term(P);
    return jadd(val, mv.hi < lin.hi ? mv : lin).hi;
}

/* best upper bound on log eta from flag U (Psi form and split forms) */
static double eta_bound(const Flag *U) {
    double w[N] = {0, 1, 1, 1};
    double best = tangent_bound(U, w, jp(0));          /* det = 1: log eta = sum_{k<4} log G_k(U) */
    for (int k = 1; k < N; k++) {
        double ws[N] = {0, 0, 0, 0}; ws[k] = 2;           /* eta <= eta_k eta_{4-k} G_k^{2} */
        double v = tangent_bound(U, ws, jadd(LOGETA[k], LOGETA[N - k]));
        if (v < best) best = v;
    }
    return best;
}

/* floating-point flag search at the centre (only proposes flags; their bounds are rigorous) */
static double Gf[N][N];
static double flog_eta(const Flag *U) {
    double A[N][N], L[N][N] = {{0}};
    for (int a = 0; a < N; a++) for (int b = 0; b < N; b++) {
        double s = 0;
        for (int i = 0; i < N; i++) for (int j = 0; j < N; j++) s += U->u[a][i] * Gf[i][j] * U->u[b][j];
        A[a][b] = s;
    }
    double tot = 0, ld = 0;
    for (int c = 0; c < N; c++) {
        double s = A[c][c];
        for (int q = 0; q < c; q++) s -= L[c][q] * L[c][q];
        if (!(s > 0)) return INFINITY;
        L[c][c] = sqrt(s);
        for (int i = c + 1; i < N; i++) { double t = A[i][c]; for (int q = 0; q < c; q++) t -= L[i][q] * L[c][q]; L[i][c] = t / L[c][c]; }
        ld += 2 * log(L[c][c]);
        if (c < N - 1) tot += ld;
    }
    return tot - 1.5 * ld;
}
static void descend(Flag *U) {
    double v = flog_eta(U);
    for (int it = 0; it < 40; it++) {
        Flag best = *U; double bv = v;
        for (int i = 0; i < N; i++) for (int j = 0; j < N; j++) if (i != j) for (int s = -1; s <= 1; s += 2) {
            Flag V = *U;
            for (int c = 0; c < N; c++) V.u[i][c] += s * U->u[j][c];
            double x = flog_eta(&V);
            if (x < bv - 1e-12) { bv = x; best = V; }
        }
        for (int i = 0; i < N; i++) for (int j = i + 1; j < N; j++) {
            Flag V = *U;
            for (int c = 0; c < N; c++) { int t = V.u[i][c]; V.u[i][c] = V.u[j][c]; V.u[j][c] = t; }
            double x = flog_eta(&V);
            if (x < bv - 1e-12) { bv = x; best = V; }
        }
        if (bv >= v - 1e-12) break;
        v = bv; *U = best;
    }
}

/* elementary changes: b_i +- b_j and transpositions; U beats the identity on the box => not Psi-minimal */
static int beaten(const Data *d) {
    J sy = jadd(jadd(d->y[1], d->y[2]), d->y[3]);
    for (int i = 0; i < N; i++) for (int j = 0; j < N; j++) if (i != j) for (int s = -1; s <= 1; s += 2) {
        Flag U; memset(&U, 0, sizeof(U)); for (int q = 0; q < N; q++) U.u[q][q] = 1;
        U.u[i][j] = s;
        double w[N] = {0, 1, 1, 1};
        if (tangent_bound(&U, w, jp(0)) < sy.lo) { n_disc[4]++; return 1; }
        /* the same comparison as one function, log Psi^2_U - log Psi^2_id, so that the first-order terms cancel
           (needed near ties, where both are equal at a nearby point) */
        YSUB = 1;
        double dv = tangent_bound(&U, w, jsub(jp(0), SY_C));
        YSUB = 0;
        if (dv < 0) { n_disc[4]++; return 1; }
    }
    /* moves of b_k to position i (any i != k), in the difference form */
    for (int i = 0; i < N; i++) for (int k = 0; k < N; k++) if (i != k) {
        int p[N], t = 0;
        for (int q = 0; q < N; q++) if (q != k) p[t++] = q;
        for (int q = N - 1; q > i; q--) p[q] = p[q - 1];
        p[i] = k;
        Flag U; memset(&U, 0, sizeof(U));
        for (int r = 0; r < N; r++) U.u[r][p[r]] = 1;
        double w[N] = {0, 1, 1, 1};
        YSUB = 1;
        double dv = tangent_bound(&U, w, jsub(jp(0), SY_C));
        YSUB = 0;
        if (dv < 0) { n_disc[4]++; return 1; }
    }
    return 0;
}

/* ---------------- kappa < K via branch-and-bound over the phi-torus ---------------- */
static int SK[400][N], nsk = 0;
static void gen_sk(void) {
    for (int code = 0; code < 625; code++) {
        int k[N], c = code, f = -1;
        for (int i = 0; i < N; i++) { k[i] = c % 5 - 2; c /= 5; }
        for (int i = 0; i < N; i++) if (k[i]) { f = i; break; }
        if (f < 0 || k[f] < 0) continue;
        memcpy(SK[nsk++], k, sizeof(k));
    }
}
/* is every phi covered, i.e. for every phi some s has ||phi(s)|| < tau_s? (proves kappa < K) */
static int torus_covered(int na, int act[], double tau[]) {
    static double st[100000][8];
    int sp = 0;
    double b0[8] = {0, 0.5, 0, 1, 0, 1, 0, 1};
    memcpy(st[sp++], b0, sizeof(b0));
    while (sp > 0) {
        double *b = st[--sp], bx[8];
        memcpy(bx, b, sizeof(bx));
        int covered = 0;
        for (int a = 0; a < na && !covered; a++) {
            const int *k = SK[act[a]];
            J v = jp(0);
            for (int i = 0; i < N; i++) v = jadd(v, jscale(jv(bx[2 * i], bx[2 * i + 1]), k[i]));
            /* nearest integer m to the midpoint (any integer works; the test is exact interval logic) */
            fesetround(FE_TONEAREST); double m = nearbyint(0.5 * (v.lo + v.hi)); fesetround(FE_UPWARD);
            J dv = jsub(v, jp(m));
            if (dv.lo > -tau[a] && dv.hi < tau[a]) covered = 1;
        }
        if (covered) continue;
        int dd = 0; double w = -1;
        for (int i = 0; i < N; i++) if (bx[2 * i + 1] - bx[2 * i] > w) { w = bx[2 * i + 1] - bx[2 * i]; dd = i; }
        if (w < 1.0 / 2048 || sp + 2 >= 100000) return 0;
        fesetround(FE_TONEAREST); double mid = 0.5 * (bx[2 * dd] + bx[2 * dd + 1]); fesetround(FE_UPWARD);
        double c1[8], c2[8];
        memcpy(c1, bx, sizeof(c1)); memcpy(c2, bx, sizeof(c2));
        c1[2 * dd + 1] = mid; c2[2 * dd] = mid;
        memcpy(st[sp++], c1, sizeof(c1)); memcpy(st[sp++], c2, sizeof(c2));
    }
    return 1;
}
static int kappa_below(const Data *d, J K) {
    /* |s|^2 = sum_j r_j (sum_{i >= j} k_i mu_ij)^2 */
    int act[400], na = 0; double tau[400];
    for (int q = 0; q < nsk; q++) {
        J tot = jp(0);
        for (int j = 0; j < N; j++) {
            J c = jp(0);
            for (int i = j; i < N; i++) if (SK[q][i]) c = jadd(c, jscale(d->mu[i][j], SK[q][i]));
            tot = jadd(tot, jmul(d->r[j], jsqr(c)));
        }
        /* need ||phi(s)||^2 < 1 - n_s/K for coverage: tau = lower bound of sqrt(1 - n_s.hi / K.lo) */
        J t2 = jsub(jp(1), jdiv(jp(tot.hi), jp(K.lo)));
        if (!(t2.lo > 0)) continue;
        double t = jsqrt(jp(t2.lo)).lo;
        if (t > 0.5) return 1;      /* ||phi(s)|| <= 1/2 < tau for every phi: the cheap bound kappa <= n_s/(3/4) < K */
        act[na] = q; tau[na] = t; na++;
    }
    if (na == 0) return 0;
    return torus_covered(na, act, tau);
}

/* ---------------- leaf verification ---------------- */
static Flag HINT; static int hint_ok = 0;   /* flag proposed by the log (any integer matrix: checked to be unimodular) */
static int unimodular(const Flag *U) {
    long a[N][N]; for (int i = 0; i < N; i++) for (int j = 0; j < N; j++) a[i][j] = U->u[i][j];
    /* integer determinant by fraction-free elimination (Bareiss) */
    long prev = 1; int sign = 1;
    for (int k = 0; k < N - 1; k++) {
        if (a[k][k] == 0) {
            int r = k + 1; while (r < N && a[r][k] == 0) r++;
            if (r == N) return 0;
            for (int j = 0; j < N; j++) { long t = a[k][j]; a[k][j] = a[r][j]; a[r][j] = t; }
            sign = -sign;
        }
        for (int i = k + 1; i < N; i++) for (int j = k + 1; j < N; j++) a[i][j] = (a[i][j] * a[k][k] - a[i][k] * a[k][j]) / prev;
        prev = a[k][k];
    }
    long det = sign * a[N - 1][N - 1];
    return det == 1 || det == -1;
}
/* ---------------- M* balls (MSTAR=file) ----------------
 * Proposition [M*]: Phi_5(G) <= 2^16/(3^4 5^3) whenever d(G, G0) <= 1/40, for G0 a det-1 Gram matrix of M*.  The file
 * lists bases of M* as unimodular integer matrices U (16 entries per line); then G0 = U P U^T / sqrt(135) with P the
 * integral Gram matrix of the proposition (det 135^2).  A box lies in the ball if delta^2 := tr((G0^{-1} E)^2) <= 1/41^2
 * on it, E = G - G0: delta = |W E W|_F (W = G0^{-1/2}) and d(G, G0) = |log(I + W E W)|_F <= delta / (1 - delta) <= 1/40.
 * G0^{-1} = sqrt(135) adj(Q) / 135^2, Q = U P U^T, with adj(Q) computed exactly in integers. */
#define MBMAX 64
static int nmb = 0;
static J MB_G[MBMAX][N][N], MB_GI[MBMAX][N][N];
static long n_mball = 0;
static long long det3(long long a[3][3]) {
    return a[0][0] * (a[1][1] * a[2][2] - a[1][2] * a[2][1]) - a[0][1] * (a[1][0] * a[2][2] - a[1][2] * a[2][0])
         + a[0][2] * (a[1][0] * a[2][1] - a[1][1] * a[2][0]);
}
static long long det4(long long A[N][N]) {
    long long s = 0;
    for (int c = 0; c < N; c++) {
        long long m[3][3];
        for (int i = 1; i < N; i++) for (int j = 0, k = 0; j < N; j++) if (j != c) m[i - 1][k++] = A[i][j];
        s += (c % 2 ? -1 : 1) * A[0][c] * det3(m);
    }
    return s;
}
/* Ball families: each lists bases U of one lattice (16 entries per line), with an integral Gram matrix P; the det-1
 * Gram matrix is G0 = U P U^T / s, s = (det P)^{1/4}, and the ball is d(G, G0) <= rho, accepted when
 * delta^2 <= (rho / (1 + rho))^2 (see above).
 *   MSTAR=file: M*, P of Proposition [M*] (det 135^2), s = sqrt(135), rho = 1/40 (Phi <= 65536/10125 inside);
 *   PVD5=file:  pi_v D5, P = 2 G0(pi_v D5) (det 32), s = 2^{5/4}, rho = 1/20 (Prop. [pi_v D5]: Phi <= 6 inside, equality
 *               only at G0). */
static J MB_LIM[MBMAX];
static int MB_FAM[MBMAX];
static long n_fball[2];
static void load_family(const char *fn, const long long P[N][N], long long detP, J s, J lim, int fam, const char *name) {
    FILE *f = fopen(fn, "r");
    if (!f) { fprintf(stderr, "cannot open %s\n", fn); exit(1); }
    int n0 = nmb;
    if (det4((long long (*)[N])P) != detP) { fprintf(stderr, "%s: det P\n", name); exit(2); }
    while (nmb < MBMAX) {
        long long U[N][N], Q[N][N], A[N][N];
        int ok = 1;
        for (int i = 0; i < N * N && ok; i++) if (fscanf(f, "%lld", &U[i / N][i % N]) != 1) ok = 0;
        if (!ok) break;
        long long du = det4(U);
        if (du != 1 && du != -1) { fprintf(stderr, "%s: matrix %d not unimodular\n", name, nmb - n0); exit(2); }
        for (int i = 0; i < N; i++) for (int j = 0; j < N; j++) {
            long long t = 0;
            for (int a = 0; a < N; a++) for (int b = 0; b < N; b++) t += U[i][a] * P[a][b] * U[j][b];
            Q[i][j] = t;
        }
        if (det4(Q) != detP) { fprintf(stderr, "%s: det\n", name); exit(2); }
        for (int i = 0; i < N; i++) for (int j = 0; j < N; j++) {     /* adj(Q)_{ij} = (-1)^{i+j} minor_{ji} */
            long long m[3][3];
            for (int a = 0, ra = 0; a < N; a++) if (a != j) { for (int b = 0, cb = 0; b < N; b++) if (b != i) m[ra][cb++] = Q[a][b]; ra++; }
            A[i][j] = ((i + j) % 2 ? -1 : 1) * det3(m);
        }
        for (int i = 0; i < N; i++) for (int j = 0; j < N; j++) {      /* check Q adj(Q) = det Q I */
            long long t = 0; for (int k = 0; k < N; k++) t += Q[i][k] * A[k][j];
            if (t != (i == j ? detP : 0)) { fprintf(stderr, "%s: adjugate\n", name); exit(2); }
            MB_G[nmb][i][j] = jdiv(jp((double)Q[i][j]), s);                       /* G0 = Q / s */
            MB_GI[nmb][i][j] = jdiv(jmul(jp((double)A[i][j]), s), jp((double)detP)); /* G0^{-1} = s adj(Q) / det Q */
        }
        MB_LIM[nmb] = lim; MB_FAM[nmb] = fam;
        nmb++;
    }
    fclose(f);
    fprintf(stderr, "%s: %d copies\n", name, nmb - n0);
}
static void load_mstar(const char *fn) {
    static const long long P[N][N] = {{18, -9, 0, 0}, {-9, 18, -9, -9}, {0, -9, 16, 1}, {0, -9, 1, 16}};
    load_family(fn, P, 135LL * 135LL, jsqrt(jp(135)), jdiv(jp(1), jp(41.0 * 41.0)), 0, "MSTAR");
}
static void load_pvd5(const char *fn) {
    static const long long P[N][N] = {{4, -2, 0, 0}, {-2, 3, -2, 0}, {0, -2, 4, -2}, {0, 0, -2, 4}};
    load_family(fn, P, 32, jsqrt(jsqrt(jp(32))), jdiv(jp(1), jp(21.0 * 21.0)), 1, "PVD5");
}
static int in_mstar_ball(const Data *d) {
    for (int q = 0; q < nmb; q++) {
        J T[N][N];
        for (int i = 0; i < N; i++) for (int j = 0; j < N; j++) {
            J acc = jp(0);
            for (int k = 0; k < N; k++) acc = jadd(acc, jmul(MB_GI[q][i][k], jsub(d->G[k][j], MB_G[q][k][j])));
            T[i][j] = acc;
        }
        J tr = jp(0);
        for (int i = 0; i < N; i++) for (int j = 0; j < N; j++) tr = jadd(tr, jmul(T[i][j], T[j][i]));
        if (tr.hi <= MB_LIM[q].lo) { n_mball++; n_fball[MB_FAM[q]]++; return 1; }
    }
    return 0;
}
static long n_leaf = 0, n_cert = 0, n_sub = 0, n_fail = 0, n_missing = 0;
static double worst_phi_hint = 0;

static int verify_box(const Box *b, int extra) {
    Data d;
    fill(b, &d);
    if (discard(&d)) return 1;
    if (nmb && in_mstar_ball(&d)) return 1;
    /* centre and G - G* */
    Box c;
    fesetround(FE_TONEAREST);
    for (int t = 0; t < D; t++) c.lo[t] = c.hi[t] = 0.5 * (b->lo[t] + b->hi[t]);
    fesetround(FE_UPWARD);
    Data dc;
    fill(&c, &dc);
    for (int i = 0; i < N; i++) for (int j = 0; j < N; j++) {
        GC[i][j] = dc.G[i][j];
        DG[i][j] = jsub(d.G[i][j], dc.G[i][j]);
        Gf[i][j] = 0.5 * (dc.G[i][j].lo + dc.G[i][j].hi);
    }
    BX = &d;
    SY_BOX = jadd(jadd(d.y[1], d.y[2]), d.y[3]);
    SY_C = jadd(jadd(dc.y[1], dc.y[2]), dc.y[3]);
    for (int t = 0; t < D; t++) DX[t] = jv(-(c.lo[t] - b->lo[t]), b->hi[t] - c.lo[t]);
    if (beaten(&d)) return 1;
    /* eta bound: identity, then descended flags from a few starts */
    J sy = jadd(jadd(d.y[1], d.y[2]), d.y[3]);
    double E = sy.hi;
    Flag U; memset(&U, 0, sizeof(U)); for (int q = 0; q < N; q++) U.u[q][q] = 1;
    double v = eta_bound(&U); if (v < E) E = v;
    descend(&U);
    /* the descended basis U (unimodular: built from row additions and swaps) beating the identity on the whole
       box, in the difference form, also shows that the box contains no Psi-minimal basis */
    {
        int id = 1;
        for (int i = 0; i < N; i++) for (int j = 0; j < N; j++) if (U.u[i][j] != (i == j)) id = 0;
        if (!id) {
            double w[N] = {0, 1, 1, 1};
            YSUB = 1;
            double dv = tangent_bound(&U, w, jsub(jp(0), SY_C));
            YSUB = 0;
            if (dv < 0) { n_disc[4]++; return 1; }
        }
    }
    v = eta_bound(&U); if (v < E) E = v;
    if (hint_ok && unimodular(&HINT)) {
        double w[N] = {0, 1, 1, 1};
        YSUB = 1;
        double dv = tangent_bound(&HINT, w, jsub(jp(0), SY_C));
        YSUB = 0;
        if (dv < 0) { n_disc[4]++; return 1; }
        v = eta_bound(&HINT); if (v < E) E = v;
    }
    /* eta(M) <= eta_4 = 3 for every 4-dimensional lattice (Theorem eta4) */
    if (LOGETA[4].hi < E) E = LOGETA[4].hi;
    /* K = sqrt(PHI / exp(E)), as a rigorous lower bound */
    J K = jexp(jscale(jsub(LOGPHI, jp(E)), 0.5));
    K.hi = K.lo;
    if (kappa_below(&d, K)) { n_cert++; return 1; }
    if (getenv("DBG") && extra == 0) {
        double fl = flog_eta(&U);
        double nmin = 1e9;
        for (int q = 0; q < nsk; q++) { double t = 0; for (int i = 0; i < N; i++) for (int j = 0; j < N; j++) t += SK[q][i] * Gf[i][j] * SK[q][j]; if (t < nmin) nmin = t; }
        fprintf(stderr, "  dbg: centre eta %.4f  E bound %.4f  K %.4f  lambda1^2 %.4f  (4/3)l1 %.4f  widths", exp(fl), exp(E), K.lo, nmin, nmin * 4 / 3);
        for (int t = 0; t < D; t++) fprintf(stderr, " %.3g", b->hi[t] - b->lo[t]);
        fprintf(stderr, "\n");
    }
    if (extra <= 0) return 0;
    /* subdivide the widest coordinate (relative to its range) */
    int ts = 0; double wm = -1;
    for (int t = 0; t < D; t++) { double w = (b->hi[t] - b->lo[t]) * (t < 3 ? 1.0 : 1.0); if (w > wm) { wm = w; ts = t; } }
    fesetround(FE_TONEAREST); double mid = 0.5 * (b->lo[ts] + b->hi[ts]); fesetround(FE_UPWARD);
    Box b1 = *b, b2 = *b;
    b1.hi[ts] = mid; b2.lo[ts] = mid;
    n_sub++;
    return verify_box(&b1, extra - 1) && verify_box(&b2, extra - 1);
}

/* ---------------- replay the tree ---------------- */
static FILE *lf;
static char line[256];
static int EXTRA = 12;
static long n_coll = 0;
static FILE *stage2 = NULL;   /* STAGE2=file: the collected boxes, exactly (hex) */
static int replay(Box *b) {
    if (!fgets(line, sizeof(line), lf) || line[0] == 'E') { n_missing++; return 0; }   /* truncated log: subtree unverified */
    if (line[0] == 'C') {          /* handed to stage 2: not proved here; recorded exactly */
        n_coll++;
        if (stage2) { for (int t = 0; t < D; t++) fprintf(stage2, "%a %a ", b->lo[t], b->hi[t]); fprintf(stage2, "\n"); }
        return 1;
    }
    if (line[0] == 'L') {
        n_leaf++;
        /* optional hint: the flag certify used (16 integers after the kappa value); only proposed, re-proved here */
        hint_ok = 0;
        {
            char *p = line + 1; strtod(p, &p);
            int v[16], k = 0;
            while (k < 16) { char *q; long x = strtol(p, &q, 10); if (q == p) break; v[k++] = (int)x; p = q; }
            if (k == 16) { for (int i = 0; i < 16; i++) HINT.u[i / 4][i % 4] = v[i]; hint_ok = 1; }
        }
        if (!verify_box(b, EXTRA)) {
            n_fail++;
            fprintf(stdout, "FAIL leaf:");
            for (int t = 0; t < D; t++) fprintf(stdout, " [%.6g,%.6g]", b->lo[t], b->hi[t]);
            fprintf(stdout, "\n");
        }
        if (n_leaf % 100000 == 0) { fprintf(stderr, "  %ld leaves (%ld fail)\n", n_leaf, n_fail); }
        return 1;
    }
    if (line[0] != 'S') { fprintf(stderr, "bad line: %s", line); exit(2); }
    int t = atoi(line + 2);
    if (t < 0 || t >= D) { fprintf(stderr, "bad split coordinate\n"); exit(2); }
    fesetround(FE_TONEAREST); double mid = 0.5 * (b->lo[t] + b->hi[t]); fesetround(FE_UPWARD);
    Box hi = *b, lo = *b;
    hi.lo[t] = mid; lo.hi[t] = mid;
    replay(&hi);          /* certify.c pops the upper half first */
    replay(&lo);
    return 1;
}

int main(int argc, char **argv) {
    if (argc < 2) { fprintf(stderr, "usage: %s logs...\n", argv[0]); return 1; }
    if (getenv("EXTRA")) EXTRA = atoi(getenv("EXTRA"));
    if (getenv("STAGE2")) { stage2 = fopen(getenv("STAGE2"), "w"); if (!stage2) { fprintf(stderr, "cannot open STAGE2 file\n"); return 1; } }
    fesetround(FE_UPWARD);
    gen_sk();
    double phi = 0;
    for (int a = 1; a < argc; a++) {
        lf = fopen(argv[a], "r");
        if (!lf) { fprintf(stderr, "cannot open %s\n", argv[a]); return 1; }
        int n, dim; double c;
        /* parse with round-to-nearest: a decimal %.17g header read upward can come out one ulp high */
        fesetround(FE_TONEAREST);
        int hdr = fgets(line, sizeof(line), lf) && sscanf(line, "ROOT n=%d dim=%d phi=%lf c=%lf", &n, &dim, &phi, &c) == 4;
        fesetround(FE_UPWARD);
        char line0[sizeof(line)]; snprintf(line0, sizeof(line0), "%s", line); line0[strcspn(line0, "\n")] = 0;
        if (!hdr || n != 4 || dim != D) { fprintf(stderr, "%s: bad header\n", argv[a]); return 1; }
        if (getenv("MSTAR") && !(jp(phi).hi <= jdiv(jp(65536), jp(10125)).lo)) {
            fprintf(stderr, "%s: MSTAR needs PHI <= 65536/10125 (got %a)\n", argv[a], phi); return 1;
        }
        if (getenv("PVD5") && !(phi <= 6.0)) { fprintf(stderr, "%s: PVD5 needs PHI <= 6 (got %a)\n", argv[a], phi); return 1; }
        constants(phi);
        if (!nmb && getenv("MSTAR")) load_mstar(getenv("MSTAR"));
        if (!nmb && getenv("PVD5")) load_pvd5(getenv("PVD5"));
        Box b;
        for (int t = 0; t < D; t++) {
            if (!fgets(line, sizeof(line), lf)) return 1;
            b.lo[t] = strtod(line, NULL);
            char *sp = strchr(line, ' ');
            b.hi[t] = strtod(sp + 1, NULL);
        }
        const char *lv = strstr(line0, " ver=");
        printf("%s: root y1 in [%a, %a]  certify %s  check_phi ver=%s\n", argv[a], b.lo[0], b.hi[0], lv ? lv + 1 : "ver=unknown", GITREV);
        replay(&b);
        if (n_missing) { printf("%s: log incomplete (%ld subtrees missing)\n", argv[a], n_missing); n_fail++; }
        else if (!fgets(line, sizeof(line), lf) || strncmp(line, "END complete", 12) != 0) { printf("%s: log not complete\n", argv[a]); n_fail++; }
        n_missing = 0;
        fclose(lf);
    }
    printf("leaves %ld, certified by kappa/eta %ld, extra subdivisions %ld, discards [window %ld, potential %ld, "
           "Lagrange %ld, insertion %ld, dual insertion %ld, elementary %ld], FAILURES %ld\n",
           n_leaf, n_cert, n_sub, n_disc[0], n_disc[1], n_disc[2], n_disc[3], n_disc[5], n_disc[4], n_fail);
    if (getenv("MSTAR")) printf("boxes inside an M* ball (Proposition [M*]): %ld\n", n_fball[0]);
    if (getenv("PVD5")) printf("boxes inside a pi_v D5 ball (Proposition [pi_v D5]): %ld\n", n_fball[1]);
    if (getenv("MSTAR") && !n_fail && !n_coll) printf("so Phi <= 65536/10125 on the root boxes, with equality only at M* (PHI = %a)\n", phi);
    if (n_coll) printf("collected boxes (handed to stage 2, not proved here): %ld\n", n_coll);
    if (stage2) fclose(stage2);
    (void)worst_phi_hint;
    return n_fail ? 4 : 0;
}
