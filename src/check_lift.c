/* check_lift.c -- independent checker for the stage-2 certificates of stage2.c (eta_5 = 6).
 *
 *   gcc -O2 -frounding-math -o check_lift check_lift.c -lm
 *   ./check_lift root.log [more...]        (EXTRA=k: up to k extra bisections of a leaf that fails as given, default 12)
 *
 * Written separately from stage2.c: it shares only the coordinates and the log format.  G(z) and its derivatives are
 * recomputed here from the formula, in interval arithmetic with outward rounding (the primitives of check_phi.c);
 * nothing about how stage2.c enclosed G or chose a leaf is trusted.
 *
 * Setting (notes/stage2_design_for_web_claude.md).  det L = 1, v shortest, r = |v|^2, M = pi_v L rescaled to det 1 with
 * a Gram matrix S in a stage-1 collected box X, lifts b_i' = b_i + phi_i v.  In the basis (v, b_1', ..., b_4')
 *      G(z) = [[r, r phi^T], [r phi, r^{-1/4} S + r phi phi^T]],   det G = 1,
 * z = (r, phi_1..phi_4, y_1..y_3, mu_10, mu_20, mu_21, mu_30, mu_31, mu_32), S = mu R mu^T,
 * R = diag(e^{y_1}, e^{y_2 - y_1}, e^{y_3 - y_2}, e^{-y_3}).  A root is X x [0,1/2] x [-1/2,1/2]^3 x [r_lo, r_hi] with
 * r_lo <= 2^{2/5} (below it eta <= 3 r^{5/2} < 6) and r_hi >= 8^{1/5} = gamma_5; phi modulo Z^4 and +-.
 *
 * Log: "ROOT2 ver=..", 14 lines "lo hi" (hex; coordinate order as z), then the tree in pre-order: "S t" splits
 * coordinate t at mid = 0.5 (lo + hi) (round to nearest), LOWER half first; leaves:
 *   N k c1..c4      not shortest: with w = (-k, c) (c != 0), w^T G w < G_00 = r on the cell, so v is not shortest
 *                   there (strictly), and every lattice also has a parametrisation with v shortest;
 *   F u_11..u_55    a unimodular U (columns: a basis); one of the forms of the flag bound below is < log 6;
 *   P m, then m lines "f lambda u_11..u_55": lambda_i >= 0; sum_i lambda_i B_{f_i}(U_i) / sum lambda < log 6, since
 *                   log eta <= min_i B_i <= that weighted mean;
 *   B t w_11..w_55  a unimodular W (rows: a basis); d(W G W^T, T_t) <= 1/10 for T_0 = Cartan(D5) / 4^{1/5},
 *                   T_1 = Gram(e_1..e_4, h/2 1) 4^{1/5} (det 1): eta <= 6 by Theorem D5 (and its dual), equality only at
 *                   the centre.  The only leaf type that may contain eta = 6.
 *   "END complete".
 * Flag bound forms for U, A_j = leading j x j block of U^T G U:
 *   f = 0:      log eta <= sum_{j=1}^{4} log det A_j                                (det G = 1)
 *   f = k >= 1: log eta <= log eta_k + log eta_{5-k} + (5/2) log det A_k            (split lemma, det G = 1)
 * with eta_1 = 1, eta_2 = 2/sqrt 3, eta_3 = 3/2, eta_4 = 3.  Each log det A_j is bounded by the tangent inequality of
 * the concave log det, valid for ANY positive definite Y_j:   log det A <= -log det Y_j - j + tr(Y_j A),
 * with Y_j = (A_j at the cell centre)^{-1} computed in floating point (any Y_j works; its log det is enclosed by an
 * interval Cholesky).  So every bound is c + l(G) with l linear in the entries of G.
 * Enclosure of l(G(z)) over the cell: the mean-value form l(G(z_c)) + sum_t l(dG/dz_t over the cell) (Z_t - z_c),
 * intersected with the naive interval l(G(cell)).  For B, delta^2 = sum_cd Y_cd^2 with
 * Y = L^{-1} (W G W^T) L^{-T} - I (T = L L^T; orthogonally similar to T^{-1/2}(.)T^{-1/2}), delta <= 1/11 gives
 * d <= delta / (1 - delta) <= 1/10.
 * Arithmetic: FE_UPWARD throughout, lower ends as -(upward(-x)); exp/log are CORE-MATH's correctly rounded cr_exp/cr_log
 * (third_party/core-math), widened by one ulp; sqrt is IEEE correctly rounded; splits replayed in round-to-nearest
 * (fesetround around the midpoint only).  Build: gcc -O2 -frounding-math -ffp-contract=off (baseline x86-64, no FMA).
 */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <fenv.h>
#ifndef GITREV
#define GITREV "unknown"
#endif

#define N 5
#define NZ 14

/* ---------------- interval arithmetic (FE_UPWARD throughout), as in check_phi.c ---------------- */
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
static J jhull(J a, J b) { return jv(fmin(a.lo, b.lo), fmax(a.hi, b.hi)); }
static J jmeet(J a, J b) { return jv(fmax(a.lo, b.lo), fmin(a.hi, b.hi)); }

/* ---------------- constants ---------------- */
static J LOGETA[N], LOG6, DLIM;
static J TSCALE[2];                 /* sigma_t: T_t = sigma_t K_t with K_t integral */
static const long long K_T[2][N][N] = {
    {{2, -1, 0, 0, 0}, {-1, 2, -1, 0, 0}, {0, -1, 2, -1, -1}, {0, 0, -1, 2, 0}, {0, 0, -1, 0, 2}},     /* Cartan(D5), det 4 */
    {{4, 0, 0, 0, 2}, {0, 4, 0, 0, 2}, {0, 0, 4, 0, 2}, {0, 0, 0, 4, 2}, {2, 2, 2, 2, 5}}};            /* 4 Gram(D5^*), det 256 */
static J LK[2][N][N];               /* interval Cholesky factors of K_t */
static long long idet(long long A0[N][N]);
static void constants(void) {
    LOGETA[1] = jp(0);
    LOGETA[2] = jsub(jlog(jp(2)), jscale(jlog(jp(3)), 0.5));   /* log(2/sqrt 3) */
    LOGETA[3] = jsub(jlog(jp(3)), jlog(jp(2)));                /* log(3/2) */
    LOGETA[4] = jlog(jp(3));
    LOG6 = jlog(jp(6));
    DLIM = jdiv(jp(1), jp(121));                                /* delta^2 <= 1/11^2 */
    TSCALE[0] = jexp(jscale(jlog(jp(4)), -0.2));               /* 4^{-1/5}: det(sigma K_0) = sigma^5 4 = 1 */
    TSCALE[1] = jexp(jscale(jlog(jp(4)), -0.8));               /* 4^{-4/5}: sigma^5 256 = 1 */
    for (int t = 0; t < 2; t++) {
        long long K[N][N]; memcpy(K, K_T[t], sizeof(K));
        if (idet(K) != (t ? 256 : 4)) { fprintf(stderr, "target det\n"); exit(2); }
        for (int j = 0; j < N; j++) {
            J s = jp((double)K_T[t][j][j]);
            for (int k = 0; k < j; k++) s = jsub(s, jsqr(LK[t][j][k]));
            if (!(s.lo > 0)) { fprintf(stderr, "target Cholesky\n"); exit(2); }
            LK[t][j][j] = jsqrt(s);
            for (int i = j + 1; i < N; i++) {
                J u = jp((double)K_T[t][i][j]);
                for (int k = 0; k < j; k++) u = jsub(u, jmul(LK[t][i][k], LK[t][j][k]));
                LK[t][i][j] = jdiv(u, LK[t][j][j]);
            }
            for (int i = 0; i < j; i++) LK[t][i][j] = jp(0);
        }
    }
}

/* exact integer determinant (Bareiss) */
static long long idet(long long A0[N][N]) {
    long long a[N][N]; memcpy(a, A0, sizeof(a));
    long long prev = 1; int sign = 1;
    for (int k = 0; k < N - 1; k++) {
        if (a[k][k] == 0) {
            int p = -1; for (int i = k + 1; i < N; i++) if (a[i][k]) { p = i; break; }
            if (p < 0) return 0;
            for (int j = 0; j < N; j++) { long long t = a[k][j]; a[k][j] = a[p][j]; a[p][j] = t; }
            sign = -sign;
        }
        for (int i = k + 1; i < N; i++) for (int j = k + 1; j < N; j++) a[i][j] = (a[i][j] * a[k][k] - a[i][k] * a[k][j]) / prev;
        prev = a[k][k];
    }
    return sign * a[N - 1][N - 1];
}

/* ---------------- G on a cell: second-order Taylor jets in interval arithmetic ----------------
 * A jet carries a value, the gradient in z and the Hessian.  Evaluated at the centre (point intervals) it gives G and
 * its exact gradient there; evaluated on the cell (z_t = [lo_t, hi_t]) its Hessian part encloses every second
 * derivative over the cell.  For a linear functional l of G, Taylor's theorem with the Lagrange remainder gives
 *   l(G(z)) in l(G(z_c)) + grad l(z_c) . (z - z_c) + 1/2 (z - z_c)^T H(cell) (z - z_c). */
static const int MI[6] = {1, 2, 2, 3, 3, 3}, MJ[6] = {0, 0, 1, 0, 1, 2};
static const int RC[4][3] = {{1, 0, 0}, {-1, 1, 0}, {0, -1, 1}, {0, 0, -1}};
typedef struct { double lo[NZ], hi[NZ]; } Box;
/* m: the coordinates the jet depends on (g and h are zero outside m, and are read only inside it); h is stored for
   s <= t only; with jt_hess = 0 the Hessian is not computed (centre jets need only value and gradient) */
typedef struct { J v, g[NZ], h[NZ][NZ]; unsigned m; } Jet;
static int jt_hess = 1;
#define GA(a, s) (((a)->m >> (s)) & 1 ? (a)->g[s] : jp(0))
#define HA(a, s, t) ((((a)->m >> (s)) & ((a)->m >> (t)) & 1) ? (a)->h[s][t] : jp(0))
static void jt_const(Jet *r, J c) { r->v = c; r->m = 0; }
static void jt_var(Jet *r, int t, J x) { r->v = x; r->m = 1u << t; r->g[t] = jp(1); r->h[t][t] = jp(0); }
static void jt_add(Jet *r, const Jet *a, const Jet *b) {
    Jet o; o.v = jadd(a->v, b->v); o.m = a->m | b->m;
    for (int s = 0; s < NZ; s++) if (o.m >> s & 1) {
        o.g[s] = jadd(GA(a, s), GA(b, s));
        if (jt_hess) for (int t = s; t < NZ; t++) if (o.m >> t & 1) o.h[s][t] = jadd(HA(a, s, t), HA(b, s, t));
    }
    *r = o;
}
static void jt_neg(Jet *r, const Jet *a) {
    Jet o = *a; o.v = jsub(jp(0), a->v);
    for (int s = 0; s < NZ; s++) if (o.m >> s & 1) {
        o.g[s] = jsub(jp(0), a->g[s]);
        if (jt_hess) for (int t = s; t < NZ; t++) if (o.m >> t & 1) o.h[s][t] = jsub(jp(0), a->h[s][t]);
    }
    *r = o;
}
static void jt_mul(Jet *r, const Jet *a, const Jet *b) {
    Jet o; o.v = jmul(a->v, b->v); o.m = a->m | b->m;
    for (int s = 0; s < NZ; s++) if (o.m >> s & 1) {
        J as = GA(a, s), bs = GA(b, s);
        o.g[s] = jadd(jmul(a->v, bs), jmul(b->v, as));
        if (jt_hess) for (int t = s; t < NZ; t++) if (o.m >> t & 1)
            o.h[s][t] = jadd(jadd(jmul(a->v, HA(b, s, t)), jmul(b->v, HA(a, s, t))), jadd(jmul(as, GA(b, t)), jmul(bs, GA(a, t))));
    }
    *r = o;
}
/* f(a) for a scalar function with f, f', f'' enclosed on the range of a */
static void jt_fun(Jet *r, const Jet *a, J f, J f1, J f2) {
    Jet o; o.v = f; o.m = a->m;
    for (int s = 0; s < NZ; s++) if (o.m >> s & 1) {
        o.g[s] = jmul(f1, a->g[s]);
        if (jt_hess) for (int t = s; t < NZ; t++) if (o.m >> t & 1) o.h[s][t] = jadd(jmul(f2, jmul(a->g[s], a->g[t])), jmul(f1, a->h[s][t]));
    }
    *r = o;
}
static void jt_exp(Jet *r, const Jet *a) { J e = jexp(a->v); jt_fun(r, a, e, e, e); }
static void jt_rq(Jet *r, const Jet *a) {       /* a^{-1/4}: f' = -f/(4a), f'' = 5 f/(16 a^2) */
    if (!(a->v.lo > 0)) { fprintf(stderr, "r not positive\n"); exit(3); }
    J f = jexp(jscale(jlog(a->v), -0.25));
    jt_fun(r, a, f, jdiv(jscale(f, -0.25), a->v), jdiv(jscale(f, 5.0 / 16), jsqr(a->v)));
}

/* G as jets, with every z_t = Z[t] (intervals) */
typedef struct { Jet G[N][N], S[4][4]; } JGram;     /* S: the Gram of M, for the N leaves */
static void gram(const J *Z, JGram *g) {
    Jet z[NZ], mu[4][4], R[4], S[4][4], s, t1, t2;
    for (int t = 0; t < NZ; t++) jt_var(&z[t], t, Z[t]);
    for (int i = 0; i < 4; i++) for (int j = 0; j < 4; j++) jt_const(&mu[i][j], jp(i == j ? 1 : 0));
    for (int q = 0; q < 6; q++) mu[MI[q]][MJ[q]] = z[8 + q];
    for (int k = 0; k < 4; k++) {
        Jet l; jt_const(&l, jp(0));
        for (int a = 0; a < 3; a++) if (RC[k][a]) { Jet y = z[5 + a]; if (RC[k][a] < 0) jt_neg(&y, &y); jt_add(&l, &l, &y); }
        jt_exp(&R[k], &l);
    }
    jt_rq(&s, &z[0]);
    for (int i = 0; i < 4; i++) for (int j = 0; j <= i; j++) {
        Jet acc; jt_const(&acc, jp(0));
        for (int k = 0; k <= j; k++) { jt_mul(&t1, &mu[i][k], &mu[j][k]); jt_mul(&t2, &t1, &R[k]); jt_add(&acc, &acc, &t2); }
        S[i][j] = S[j][i] = acc;
    }
    memcpy(g->S, S, sizeof(S));
    g->G[0][0] = z[0];
    for (int i = 0; i < 4; i++) { jt_mul(&g->G[0][1 + i], &z[0], &z[1 + i]); g->G[1 + i][0] = g->G[0][1 + i]; }
    for (int i = 0; i < 4; i++) for (int j = 0; j <= i; j++) {
        jt_mul(&t1, &s, &S[i][j]); jt_mul(&t2, &z[1 + i], &z[1 + j]); jt_mul(&t2, &z[0], &t2);
        jt_add(&g->G[1 + i][1 + j], &t1, &t2); g->G[1 + j][1 + i] = g->G[1 + i][1 + j];
    }
}

/* the cell: jets at the centre (value, gradient) and on the cell (Hessian, naive value), and Z_t - z_c */
typedef struct { JGram c, cell; J dz[NZ]; double Gmid[N][N]; } Cell;
static void cell_fill(const Box *b, Cell *c) {
    J Z[NZ], zc[NZ];
    fesetround(FE_TONEAREST);
    double m[NZ]; for (int t = 0; t < NZ; t++) m[t] = 0.5 * (b->lo[t] + b->hi[t]);
    fesetround(FE_UPWARD);
    for (int t = 0; t < NZ; t++) { Z[t] = jv(b->lo[t], b->hi[t]); zc[t] = jp(m[t]); c->dz[t] = jsub(Z[t], zc[t]); }
    gram(Z, &c->cell);
    jt_hess = 0; gram(zc, &c->c); jt_hess = 1;
    for (int i = 0; i < N; i++) for (int j = 0; j < N; j++) c->Gmid[i][j] = 0.5 * (c->c.G[i][j].v.lo + c->c.G[i][j].v.hi);
}
/* enclosure of sum_ij Q_ij G_ij(z) over the cell: second-order Taylor, intersected with the naive interval */
static J lin(const Cell *c, J Q[N][N]) {
    J v = jp(0), nv = jp(0), g[NZ], H[NZ][NZ];
    for (int t = 0; t < NZ; t++) { g[t] = jp(0); for (int u = 0; u < NZ; u++) H[t][u] = jp(0); }
    for (int i = 0; i < N; i++) for (int j = 0; j < N; j++) {
        J q = Q[i][j];
        if (q.lo == 0 && q.hi == 0) continue;
        v = jadd(v, jmul(q, c->c.G[i][j].v));
        nv = jadd(nv, jmul(q, c->cell.G[i][j].v));
        const Jet *jc = &c->c.G[i][j], *jz = &c->cell.G[i][j];
        for (int t = 0; t < NZ; t++) {
            if (jc->m >> t & 1) g[t] = jadd(g[t], jmul(q, jc->g[t]));
            if (jz->m >> t & 1) for (int u = t; u < NZ; u++) if (jz->m >> u & 1) H[t][u] = jadd(H[t][u], jmul(q, jz->h[t][u]));
        }
    }
    for (int t = 0; t < NZ; t++) {
        v = jadd(v, jmul(g[t], c->dz[t]));
        v = jadd(v, jmul(jscale(H[t][t], 0.5), jsqr(c->dz[t])));
        for (int u = t + 1; u < NZ; u++) v = jadd(v, jmul(H[t][u], jmul(c->dz[t], c->dz[u])));
    }
    J r = jmeet(v, nv);
    return r.lo <= r.hi ? r : v;
}

/* enclosure of c^T S c over the cell (second-order Taylor in the jets of S) */
static J quadS(const Cell *c, const int *cv) {
    J v = jp(0), nv = jp(0), g[NZ], H[NZ][NZ];
    for (int t = 0; t < NZ; t++) { g[t] = jp(0); for (int u = 0; u < NZ; u++) H[t][u] = jp(0); }
    for (int i = 0; i < 4; i++) for (int j = 0; j < 4; j++) {
        if (!cv[i] || !cv[j]) continue;
        J q = jp((double)(cv[i] * cv[j]));
        v = jadd(v, jmul(q, c->c.S[i][j].v)); nv = jadd(nv, jmul(q, c->cell.S[i][j].v));
        const Jet *jc = &c->c.S[i][j], *jz = &c->cell.S[i][j];
        for (int t = 0; t < NZ; t++) {
            if (jc->m >> t & 1) g[t] = jadd(g[t], jmul(q, jc->g[t]));
            if (jz->m >> t & 1) for (int u = t; u < NZ; u++) if (jz->m >> u & 1) H[t][u] = jadd(H[t][u], jmul(q, jz->h[t][u]));
        }
    }
    for (int t = 0; t < NZ; t++) {
        v = jadd(v, jmul(g[t], c->dz[t]));
        v = jadd(v, jmul(jscale(H[t][t], 0.5), jsqr(c->dz[t])));
        for (int u = t + 1; u < NZ; u++) v = jadd(v, jmul(H[t][u], jmul(c->dz[t], c->dz[u])));
    }
    J r = jmeet(v, nv);
    return r.lo <= r.hi ? r : v;
}

/* ---------------- leaf certificates ---------------- */
typedef struct { long long u[N][N]; } Mat;
static int unimodular(const Mat *U) { long long A[N][N]; memcpy(A, U->u, sizeof(A)); long long d = idet(A); return d == 1 || d == -1; }

/* flag U (columns = basis), form f: the bound as const + l_Q(G); 0 if a Y_j is not positive definite */
static int flag_form(const Cell *c, const Mat *U, int f, J *cst, J Q[N][N]) {
    double A[N][N];
    for (int a = 0; a < N; a++) for (int b = 0; b < N; b++) {      /* A = U^T Gmid U, floating point (only to choose Y) */
        double s = 0; for (int i = 0; i < N; i++) for (int j = 0; j < N; j++) s += U->u[i][a] * c->Gmid[i][j] * U->u[j][b];
        A[a][b] = s;
    }
    *cst = jp(0);
    for (int i = 0; i < N; i++) for (int j = 0; j < N; j++) Q[i][j] = jp(0);
    int j0 = f == 0 ? 1 : f, j1 = f == 0 ? 4 : f;
    double w = f == 0 ? 1 : 2.5;
    if (f) *cst = jadd(LOGETA[f], LOGETA[N - f]);
    for (int j = j0; j <= j1; j++) {
        /* Y = (A_j)^{-1} by Gauss-Jordan in floating point, symmetrised */
        double M[N][2 * N];
        for (int a = 0; a < j; a++) for (int b = 0; b < j; b++) { M[a][b] = A[a][b]; M[a][j + b] = (a == b); }
        for (int k = 0; k < j; k++) {
            int p = k; for (int i = k + 1; i < j; i++) if (fabs(M[i][k]) > fabs(M[p][k])) p = i;
            if (fabs(M[p][k]) < 1e-300) return 0;
            for (int b = 0; b < 2 * j; b++) { double t = M[k][b]; M[k][b] = M[p][b]; M[p][b] = t; }
            double piv = M[k][k];
            for (int b = 0; b < 2 * j; b++) M[k][b] /= piv;
            for (int i = 0; i < j; i++) if (i != k) { double fct = M[i][k]; for (int b = 0; b < 2 * j; b++) M[i][b] -= fct * M[k][b]; }
        }
        double Y[N][N];
        for (int a = 0; a < j; a++) for (int b = 0; b < j; b++) Y[a][b] = 0.5 * (M[a][j + b] + M[b][j + a]);
        /* log det Y by an interval Cholesky (fails unless Y is positive definite) */
        J L[N][N], ld = jp(0);
        for (int b = 0; b < j; b++) {
            J s = jp(Y[b][b]);
            for (int k = 0; k < b; k++) s = jsub(s, jsqr(L[b][k]));
            if (!(s.lo > 0)) return 0;
            L[b][b] = jsqrt(s);
            ld = jadd(ld, jlog(s));
            for (int i = b + 1; i < j; i++) {
                J u = jp(Y[i][b]);
                for (int k = 0; k < b; k++) u = jsub(u, jmul(L[i][k], L[b][k]));
                L[i][b] = jdiv(u, L[b][b]);
            }
        }
        /* -log det Y - j + tr(Y U_j^T G U_j): Q += U_j Y U_j^T */
        *cst = jadd(*cst, jscale(jsub(jsub(jp(0), ld), jp(j)), w));
        for (int p = 0; p < N; p++) for (int q = 0; q < N; q++) {
            J acc = jp(0);
            for (int a = 0; a < j; a++) for (int b = 0; b < j; b++)
                if (U->u[p][a] && U->u[q][b]) acc = jadd(acc, jscale(jp(Y[a][b]), (double)(U->u[p][a] * U->u[q][b])));
            Q[p][q] = jadd(Q[p][q], jscale(acc, w));
        }
    }
    return 1;
}

typedef struct { int type, k, c[4], tk, m; Mat U; int pf[64]; double pl[64]; Mat pU[64]; } Leaf;
static long n_type[5], n_sub = 0, n_fail = 0, n_leaf = 0;
static int EXTRA = 12;

static int verify_cell(const Box *b, const Leaf *lf) {
    static Cell c; cell_fill(b, &c);      /* not re-entrant: verify_cell does not recurse */
    if (lf->type == 1) {                                   /* N */
        /* |sum c_i b_i' - k v|^2 - |v|^2 = r^{-1/4} c^T S c - r (1 - u^2), u = c.phi - k: decreasing in r while
           u^2 < 1, so it is < 0 on the cell if it is < 0 at r = r_lo; u is linear in phi, so its interval is exact */
        J u = jp((double)-lf->k);
        for (int i = 0; i < 4; i++) u = jadd(u, jscale(jv(b->lo[1 + i], b->hi[1 + i]), lf->c[i]));
        J om = jsub(jp(1), jsqr(u));
        if (!(om.lo > 0)) return 0;
        J r0 = jp(b->lo[0]);
        J s0 = jexp(jscale(jlog(r0), -0.25));
        return jsub(jmul(s0, quadS(&c, lf->c)), jmul(r0, om)).hi < 0;
    }
    if (lf->type == 2) {                                   /* F: any of the five forms */
        double best = INFINITY;
        for (int f = 0; f < 5; f++) {
            J cst, Q[N][N];
            if (!flag_form(&c, &lf->U, f, &cst, Q)) continue;
            J v = jadd(cst, lin(&c, Q));
            if (v.hi < LOG6.lo) return 1;
            if (v.hi < best) best = v.hi;
        }
        if (getenv("DBG")) {
            double wm = 0; for (int t = 0; t < NZ; t++) wm = fmax(wm, b->hi[t] - b->lo[t]);
            /* sampled max over the cell of the best form evaluated exactly at points (a lower bound for any valid bound) */
            double smax = -INFINITY;
            for (int k = 0; k < 300; k++) {
                Box pt; for (int t = 0; t < NZ; t++) pt.lo[t] = pt.hi[t] = b->lo[t] + (b->hi[t] - b->lo[t]) * (rand() / (double)RAND_MAX);
                static Cell cp; cell_fill(&pt, &cp);
                double mn = INFINITY;
                for (int f = 0; f < 5; f++) { J cst, Q[N][N]; if (flag_form(&cp, &lf->U, f, &cst, Q)) mn = fmin(mn, jadd(cst, lin(&cp, Q)).hi); }
                smax = fmax(smax, mn);
            }
            fprintf(stderr, "F miss %.3e maxw %.4f sampled-max-minus-log6 %.3e\n", best - LOG6.lo, wm, smax - LOG6.lo);
        }
        return 0;
    }
    if (lf->type == 3) {                                   /* P: one weighted functional */
        J cst = jp(0), Q[N][N], sl = jp(0);
        for (int i = 0; i < N; i++) for (int j = 0; j < N; j++) Q[i][j] = jp(0);
        for (int q = 0; q < lf->m; q++) {
            J cq, Qq[N][N];
            if (!flag_form(&c, &lf->pU[q], lf->pf[q], &cq, Qq)) return 0;
            J l = jp(lf->pl[q]);
            cst = jadd(cst, jmul(l, cq)); sl = jadd(sl, l);
            for (int i = 0; i < N; i++) for (int j = 0; j < N; j++) Q[i][j] = jadd(Q[i][j], jmul(l, Qq[i][j]));
        }
        if (!(sl.lo > 0)) return 0;
        return jdiv(jadd(cst, lin(&c, Q)), sl).hi < LOG6.lo;
    }
    if (lf->type == 4) {                                   /* B */
        int t = lf->tk;
        J V[N][N];                                         /* V = L_K^{-1} W, forward substitution */
        for (int col = 0; col < N; col++) for (int i = 0; i < N; i++) {
            J s = jp((double)lf->U.u[i][col]);
            for (int k = 0; k < i; k++) s = jsub(s, jmul(LK[t][i][k], V[k][col]));
            V[i][col] = jdiv(s, LK[t][i][i]);
        }
        J inv = jdiv(jp(1), TSCALE[t]), d2 = jp(0);
        for (int a = 0; a < N; a++) for (int bb = a; bb < N; bb++) {
            J Q[N][N];
            for (int i = 0; i < N; i++) for (int j = 0; j < N; j++) Q[i][j] = jmul(inv, jmul(V[a][i], V[bb][j]));
            J y = jsub(lin(&c, Q), jp(a == bb ? 1 : 0));
            J y2 = jsqr(y);
            d2 = jadd(d2, a == bb ? y2 : jscale(y2, 2));
        }
        return d2.hi <= DLIM.lo;
    }
    return 0;
}

static int verify_leaf(const Box *b, const Leaf *lf, int depth) {
    if (verify_cell(b, lf)) return 1;
    if (depth >= EXTRA) return 0;
    /* bisect the widest coordinate (y_1..y_3 counted double) and prove both halves with the same certificate; an N leaf
       only y and mu (r enters at r_lo and phi exactly through u, so splitting them cannot tighten c^T S c) */
    int ts = 0; double wm = -1;
    for (int t = lf->type == 1 ? 5 : 0; t < NZ; t++) { double w = (b->hi[t] - b->lo[t]) * (t >= 5 && t < 8 ? 2 : 1); if (w > wm) { wm = w; ts = t; } }
    fesetround(FE_TONEAREST); double mid = 0.5 * (b->lo[ts] + b->hi[ts]); fesetround(FE_UPWARD);
    Box lo = *b, hi = *b; lo.hi[ts] = mid; hi.lo[ts] = mid;
    n_sub++;
    return verify_leaf(&lo, lf, depth + 1) && verify_leaf(&hi, lf, depth + 1);
}

/* ---------------- replay ---------------- */
static FILE *lf_in;
static char line[4096];
static long n_missing = 0;
/* a flag: 25 integers (then the next entry of this log's flag table), or "#k" = the k-th flag written in full */
#define FTAB 65536
static Mat ftab[FTAB]; static int nftab = 0;
static int read_mat(char *p, Mat *U) {
    while (*p == ' ') p++;
    if (*p == '#') { int k = atoi(p + 1); if (k < 0 || k >= nftab) return 0; *U = ftab[k]; return 1; }
    for (int i = 0; i < N; i++) for (int j = 0; j < N; j++) { char *q; U->u[i][j] = strtoll(p, &q, 10); if (q == p) return 0; p = q; }
    if (nftab < FTAB) ftab[nftab++] = *U;
    return 1;
}
static int lf_ended = 0;     /* the log's END line came where a node was expected (a producer timeout) */
static void replay(const Box *b) {
    if (lf_ended || !fgets(line, sizeof(line), lf_in)) { n_missing++; return; }
    if (line[0] == 'E') { lf_ended = 1; n_missing++; return; }   /* this and every later subtree are missing */
    if (line[0] == 'S') {
        int t = atoi(line + 2);
        if (t < 0 || t >= NZ) { fprintf(stderr, "bad split coordinate\n"); exit(2); }
        fesetround(FE_TONEAREST); double mid = 0.5 * (b->lo[t] + b->hi[t]); fesetround(FE_UPWARD);
        if (!(mid > b->lo[t] && mid < b->hi[t])) { fprintf(stderr, "degenerate split\n"); exit(2); }
        Box lo = *b, hi = *b; lo.hi[t] = mid; hi.lo[t] = mid;
        replay(&lo);           /* stage2.c recurses into the lower half first */
        replay(&hi);
        return;
    }
    Leaf L; memset(&L, 0, sizeof(L));
    int ok = 1;
    switch (line[0]) {
    case 'N':
        L.type = 1;
        ok = sscanf(line + 1, "%d %d %d %d %d", &L.k, &L.c[0], &L.c[1], &L.c[2], &L.c[3]) == 5 && (L.c[0] || L.c[1] || L.c[2] || L.c[3]);
        break;
    case 'F': L.type = 2; ok = read_mat(line + 1, &L.U) && unimodular(&L.U); break;
    case 'P':
        L.type = 3; L.m = atoi(line + 1);
        if (L.m < 1 || L.m > 64) { ok = 0; break; }
        for (int q = 0; q < L.m && ok; q++) {
            if (!fgets(line, sizeof(line), lf_in)) { fprintf(stderr, "truncated P leaf\n"); exit(2); }
            char *p = line, *e;
            L.pf[q] = (int)strtol(p, &e, 10); p = e;
            fesetround(FE_TONEAREST); L.pl[q] = strtod(p, &e); fesetround(FE_UPWARD); p = e;
            ok = L.pf[q] >= 0 && L.pf[q] <= 4 && L.pl[q] >= 0 && read_mat(p, &L.pU[q]) && unimodular(&L.pU[q]);
        }
        break;
    case 'B':
        L.type = 4; L.tk = atoi(line + 1);
        { int n0 = nftab; ok = (L.tk == 0 || L.tk == 1) && line[3] != '#' && read_mat(line + 3, &L.U) && unimodular(&L.U); nftab = n0; }
        break;
    default: fprintf(stderr, "bad line: %s", line); exit(2);
    }
    n_leaf++;
    if (ok) n_type[L.type]++;
    if (!ok || !verify_leaf(b, &L, 0)) {
        n_fail++;
        printf("FAIL %c leaf:", line[0]);
        for (int t = 0; t < NZ; t++) printf(" [%.6g,%.6g]", b->lo[t], b->hi[t]);
        printf("\n");
    }
}

int main(int argc, char **argv) {
    if (argc < 2) { fprintf(stderr, "usage: %s logs...\n", argv[0]); return 1; }
    if (getenv("EXTRA")) EXTRA = atoi(getenv("EXTRA"));
    fesetround(FE_UPWARD);
    constants();
    /* negative controls only (printed, and the run is marked as a test): a lower eta threshold or a smaller ball */
    if (getenv("TEST_ETA")) { LOG6 = jlog(jp(atof(getenv("TEST_ETA")))); printf("TEST: threshold eta < %s\n", getenv("TEST_ETA")); }
    if (getenv("TEST_BALLDEN")) { double d = atof(getenv("TEST_BALLDEN")); DLIM = jdiv(jp(1), jp(d * d)); printf("TEST: ball delta <= 1/%s\n", getenv("TEST_BALLDEN")); }
    if (getenv("TESTZ")) {          /* G at the point z (14 numbers); TESTD: gradient and Hessian vs differences */
        J z[NZ]; char *p = getenv("TESTZ");
        for (int t = 0; t < NZ; t++) z[t] = jp(strtod(p, &p));
        static JGram g; gram(z, &g);
        for (int i = 0; i < N; i++) { for (int j = 0; j < N; j++) printf(getenv("TESTHEX") ? " %a %a" : " %.12f", g.G[i][j].v.lo, g.G[i][j].v.hi); printf("\n"); }
        if (getenv("TESTD")) {
            double eg = 0, eh = 0, h = 1e-5;
            for (int t = 0; t < NZ; t++) {
                J zp[NZ], zm[NZ]; memcpy(zp, z, sizeof(zp)); memcpy(zm, z, sizeof(zm));
                zp[t] = jp(z[t].lo + h); zm[t] = jp(z[t].lo - h);
                static JGram gp, gm; gram(zp, &gp); gram(zm, &gm);
                for (int i = 0; i < N; i++) for (int j = 0; j < N; j++) {
                    eg = fmax(eg, fabs((gp.G[i][j].v.lo - gm.G[i][j].v.lo) / (2 * h) - GA(&g.G[i][j], t).lo));
                    for (int u = 0; u < NZ; u++) eh = fmax(eh, fabs((GA(&gp.G[i][j], u).lo - GA(&gm.G[i][j], u).lo) / (2 * h) - HA(&g.G[i][j], t < u ? t : u, t < u ? u : t).lo));
                }
            }
            printf("max |grad - FD| = %.2e, max |Hessian - FD| = %.2e\n", eg, eh);
        }
        return 0;
    }
    for (int a = 1; a < argc; a++) {
        lf_in = fopen(argv[a], "r");
        if (!lf_in) { fprintf(stderr, "cannot open %s\n", argv[a]); return 1; }
        if (!fgets(line, sizeof(line), lf_in) || strncmp(line, "ROOT2 ", 6) != 0) { fprintf(stderr, "%s: bad header\n", argv[a]); return 1; }
        line[strcspn(line, "\n")] = 0;
        nftab = 0;
        char hdr[256]; snprintf(hdr, sizeof(hdr), "%s", line + 6);
        Box b;
        fesetround(FE_TONEAREST);
        for (int t = 0; t < NZ; t++) {
            if (!fgets(line, sizeof(line), lf_in)) { fprintf(stderr, "%s: short header\n", argv[a]); return 1; }
            char *p; b.lo[t] = strtod(line, &p); b.hi[t] = strtod(p, NULL);
        }
        fesetround(FE_UPWARD);
        /* the fibre must be complete: r over [<= 2^{2/5}, >= 8^{1/5}], phi over [0,1/2] x [-1/2,1/2]^3 */
        /* r_lo <= 2^{2/5} iff r_lo^5 <= 4, r_hi >= 8^{1/5} iff r_hi^5 >= 8 (powers in outward rounding) */
        J l5 = jp(b.lo[0]), h5 = jp(b.hi[0]);
        for (int k = 0; k < 4; k++) { l5 = jmul(l5, jp(b.lo[0])); h5 = jmul(h5, jp(b.hi[0])); }
        int root_ok = b.lo[0] > 0 && l5.hi <= 4 && h5.lo >= 8 && b.lo[1] == 0 && b.hi[1] == 0.5;
        for (int i = 2; i <= 4; i++) root_ok = root_ok && b.lo[i] == -0.5 && b.hi[i] == 0.5;
        long f0 = n_fail;
        if (!root_ok) { printf("%s: root fibre incomplete\n", argv[a]); n_fail++; }
        lf_ended = 0; replay(&b);
        if (n_missing) { printf("%s: log incomplete (%ld subtrees missing)\n", argv[a], n_missing); n_fail++; n_missing = 0; }
        else if (!fgets(line, sizeof(line), lf_in) || strncmp(line, "END complete", 12) != 0) { printf("%s: log not complete\n", argv[a]); n_fail++; }
        else if (fgets(line, sizeof(line), lf_in)) { printf("%s: trailing data\n", argv[a]); n_fail++; }
        fclose(lf_in);
        /* the M-box X (hex), for the coverage check against the merged collected boxes */
        printf("%s: %s check_lift ver=%s X", argv[a], hdr, GITREV);
        for (int t = 5; t < NZ; t++) printf(" %a %a", b.lo[t], b.hi[t]);
        printf(" %s\n", n_fail == f0 ? "OK" : "FAILED");
    }
    if (getenv("TEST_ETA") || getenv("TEST_BALLDEN")) printf("TEST RUN (thresholds changed): not a certificate check\n");
    printf("leaves %ld [N %ld, F %ld, P %ld, B %ld], extra subdivisions %ld, FAILURES %ld\n",
           n_leaf, n_type[1], n_type[2], n_type[3], n_type[4], n_sub, n_fail);
    return n_fail ? 4 : 0;
}
