/* stage2.c -- stage-2 prototype for eta_5 = 6: eta(L) <= 6 for every lift L of a stage-1 collected box X.
 *
 *   gcc -O2 -o stage2 stage2.c -lglpk -lm
 *   ./stage2 COLFILE [first [count]]        one line of statistics per collected box (lines first..first+count-1)
 *   env: BOXCAP (max cells per box, default 2e6), TCAP (seconds per box, default 600), MERGE=k (treat the line as the
 *        union of 2^k boxes: not used; merging is done by runs/merge_col.py), NOBALL, NOLP, VERBOSE
 *
 * Reduction.  Let det L = 1, v a shortest vector, r = |v|^2, M = pi_v L rescaled to det 1, and b_i' = b_i + phi_i v
 * the lift of a basis b_i of M.  In the basis (v, b_1', ..., b_4') the Gram matrix is
 *      G = [[r, r phi^T], [r phi, r^{-1/4} S + r phi phi^T]],   det G = 1,
 * with S = mu R mu^T the Gram of M in stage-1 coordinates (y_1, y_2, y_3, mu_10, mu_20, mu_21, mu_30, mu_31, mu_32),
 * R = diag(e^{y_1}, e^{y_2 - y_1}, e^{y_3 - y_2}, e^{-y_3}).  phi is taken modulo Z^4 and +- (phi_1 in [0, 1/2]),
 * and r in [(6/eta_4)^{2/5}, gamma_5] = [2^{2/5}, 8^{1/5}].  A cell z = (r, phi, X-coordinates) closes when
 *   (a) v is not shortest on it (an explicit shorter vector, see not_shortest);
 *   (b) one flag of L bounds log eta < log 6 over it (concave tangent plane at the centre + affine enclosure of G(z));
 *   (c) an LP over several flags' tangent planes has no point above log 6;
 *   (d) (B) it lies within invariant distance 1/10 of a Gram matrix of D5 or D5^* (Theorem D5 and its dual).  The
 *       copy is found at the cell centre by matching short vectors to the target Gram (any unimodular U gives an
 *       exact copy U^{-1} T U^{-T}), and d <= delta / (1 - delta) with delta = |T^{-1/2}(U G U^T - T)T^{-1/2}|_F.
 *
 * NOTE: floating point, no outward rounding.  Exploratory: it measures cost per collected box, it proves nothing.
 */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <time.h>
#include <glpk.h>
#ifndef GITREV
#define GITREV "unknown"
#endif

#define N 5
#define NZ 14        /* z = (r, phi_1..phi_4, y_1..y_3, mu_10, mu_20, mu_21, mu_30, mu_31, mu_32) */
#define FMAX 48

static double logc;
static double LOGETA[N + 1];
static double R1LO, R1HI;

typedef struct { signed char u[N][N]; } Flag;   /* columns are the new basis vectors */
typedef struct { double lo[NZ], hi[NZ]; int nf; Flag f[FMAX]; int nU; signed char U[2][N][N]; signed char Ut[2]; } Box;

/* the certificate of the last cell closed, for the log (LOGDIR); see check_lift.c for what each type proves */
enum { LF_N = 1, LF_F, LF_P, LF_B };
#define PMAX 400
static int leaf_t, leaf_k, leaf_c[4], leaf_tk, leaf_m;
static Flag leaf_flag; static signed char leaf_W[N][N];
static Flag leaf_pf[PMAX]; static int leaf_pk[PMAX]; static double leaf_pl[PMAX];
static FILE *clog = NULL;

static Flag PERM[120];
static int nperm = 0;
static long st_cert = 0, st_split = 0, st_ball = 0;
static int verbose = 0;

/* ---------------- small linear algebra ---------------- */
static void jacobi(int m, double A0[4][4], double w[4], double V[4][4]) {
    double A[4][4];
    memcpy(A, A0, sizeof(A));
    for (int i = 0; i < 4; i++) for (int j = 0; j < 4; j++) V[i][j] = (i == j);
    for (int sweep = 0; sweep < 80; sweep++) {
        double off = 0;
        for (int p = 0; p < m; p++) for (int q = p + 1; q < m; q++) off += A[p][q] * A[p][q];
        if (off < 1e-32) break;
        for (int p = 0; p < m; p++) for (int q = p + 1; q < m; q++) {
            if (fabs(A[p][q]) < 1e-300) continue;
            double th = (A[q][q] - A[p][p]) / (2 * A[p][q]);
            double t = (th >= 0 ? 1 : -1) / (fabs(th) + sqrt(th * th + 1));
            double c = 1 / sqrt(t * t + 1), s = t * c;
            for (int k = 0; k < m; k++) { double a = A[k][p], b = A[k][q]; A[k][p] = c * a - s * b; A[k][q] = s * a + c * b; }
            for (int k = 0; k < m; k++) { double a = A[p][k], b = A[q][k]; A[p][k] = c * a - s * b; A[q][k] = s * a + c * b; }
            for (int k = 0; k < m; k++) { double a = V[k][p], b = V[k][q]; V[k][p] = c * a - s * b; V[k][q] = s * a + c * b; }
        }
    }
    for (int i = 0; i < m; i++) w[i] = A[i][i];
}

static int chol(int m, double A[N][N], double L[N][N]) {
    for (int i = 0; i < m; i++) for (int j = 0; j < m; j++) L[i][j] = 0;
    for (int j = 0; j < m; j++) {
        double s = A[j][j];
        for (int k = 0; k < j; k++) s -= L[j][k] * L[j][k];
        if (!(s > 0)) return 0;
        L[j][j] = sqrt(s);
        for (int i = j + 1; i < m; i++) {
            double t = A[i][j];
            for (int k = 0; k < j; k++) t -= L[i][k] * L[j][k];
            L[i][j] = t / L[j][j];
        }
    }
    return 1;
}

static void congr(const Flag *U, double G[N][N], double A[N][N]) {
    double T[N][N];
    for (int i = 0; i < N; i++) for (int j = 0; j < N; j++) {
        double s = 0; for (int k = 0; k < N; k++) s += G[i][k] * U->u[k][j]; T[i][j] = s;
    }
    for (int i = 0; i < N; i++) for (int j = 0; j < N; j++) {
        double s = 0; for (int k = 0; k < N; k++) s += U->u[k][i] * T[k][j]; A[i][j] = s;
    }
}

static int prefix(const Flag *U, double G[N][N], double ld[N + 1], double L[N][N]) {
    double A[N][N];
    congr(U, G, A);
    if (!chol(N, A, L)) return 0;
    ld[0] = 0;
    for (int k = 0; k < N; k++) ld[k + 1] = ld[k] + 2 * log(L[k][k]);
    return 1;
}

static double log_eta_U(const Flag *U, double G[N][N]) {
    double ld[N + 1], L[N][N];
    if (!prefix(U, G, ld, L)) return INFINITY;
    return ld[1] + ld[2] + ld[3] + ld[4] - 2 * ld[5];
}

/* P_k = U_k A_k^{-1} U_k^T, cumulative over k */
static void prefix_grads(const Flag *U, double L[N][N], double P[N + 1][N][N]) {
    double W[N][N];
    for (int r = 0; r < N; r++) for (int j = 0; j < N; j++) {
        double s = U->u[r][j];
        for (int k = 0; k < j; k++) s -= L[j][k] * W[r][k];
        W[r][j] = s / L[j][j];
    }
    memset(P[0], 0, sizeof(P[0]));
    for (int k = 1; k <= N; k++)
        for (int a = 0; a < N; a++) for (int b = 0; b < N; b++) P[k][a][b] = P[k - 1][a][b] + W[a][k - 1] * W[b][k - 1];
}

/* ---------------- local Psi search ---------------- */
static int flag_ok(const Flag *U) {
    for (int i = 0; i < N; i++) for (int j = 0; j < N; j++) if (abs(U->u[i][j]) > 60) return 0;
    return 1;
}
static double local_search(Flag *U, double G[N][N]) {
    double best = log_eta_U(U, G);
    for (int it = 0; it < 60; it++) {
        Flag bu = *U; double bv = best;
        for (int i = 0; i < N; i++) for (int j = 0; j < N; j++) if (i != j) for (int s = -1; s <= 1; s += 2) {
            Flag V = *U;
            for (int r = 0; r < N; r++) V.u[r][i] += s * U->u[r][j];
            if (!flag_ok(&V)) continue;
            double v = log_eta_U(&V, G);
            if (v < bv - 1e-12) { bv = v; bu = V; }
        }
        for (int k = 0; k < N; k++) for (int i = 0; i < N; i++) if (i != k) {
            int order[N], m = 0, ord2[N], p = 0;
            for (int q = 0; q < N; q++) if (q != k) order[m++] = q;
            for (int q = 0; q < N; q++) ord2[q] = (q == i) ? k : order[p++];
            Flag V;
            for (int r = 0; r < N; r++) for (int q = 0; q < N; q++) V.u[r][q] = U->u[r][ord2[q]];
            double v = log_eta_U(&V, G);
            if (v < bv - 1e-12) { bv = v; bu = V; }
        }
        if (bv >= best - 1e-12) break;
        best = bv; *U = bu;
    }
    return best;
}

/* beam search over flags built from short vectors (coefficients in [-2,2]^5): a much stronger upper bound for
 * eta than descent from permutations.  Returns the best log eta found and the flag in *U. */
#define NSV 3124
#define SVK 100
#define BEAM 30
static int SV[NSV][N], nsv = 0;
static void gen_sv(void) {
    int tot = 1; for (int i = 0; i < N; i++) tot *= 5;
    for (int code = 0; code < tot; code++) {
        int v[N], c = code, first = -1;
        for (int i = 0; i < N; i++) { v[i] = c % 5 - 2; c /= 5; }
        for (int i = 0; i < N; i++) if (v[i]) { first = i; break; }
        if (first < 0 || v[first] < 0) continue;
        memcpy(SV[nsv++], v, sizeof(v));
    }
}
static double idet5(int C[N][N]) {     /* integer determinant, columns C[.][j] */
    double A[N][N];
    for (int i = 0; i < N; i++) for (int j = 0; j < N; j++) A[i][j] = C[i][j];
    double d = 1;
    for (int k = 0; k < N; k++) {
        int p = k; for (int i = k + 1; i < N; i++) if (fabs(A[i][k]) > fabs(A[p][k])) p = i;
        if (fabs(A[p][k]) < 1e-9) return 0;
        if (p != k) { for (int j = 0; j < N; j++) { double t = A[k][j]; A[k][j] = A[p][j]; A[p][j] = t; } d = -d; }
        d *= A[k][k];
        for (int i = k + 1; i < N; i++) { double f = A[i][k] / A[k][k]; for (int j = k; j < N; j++) A[i][j] -= f * A[k][j]; }
    }
    return d;
}
typedef struct { int idx[N]; double s, L[N][N]; } BeamSt;
static int cmp_beam(const void *a, const void *b) { double x = ((const BeamSt *)a)->s, y = ((const BeamSt *)b)->s; return (x < y) ? -1 : (x > y); }
static double beam_search(double G[N][N], Flag *U) {
    static int ord[NSV]; static double nrm[NSV];
    int m = 0;
    for (int q = 0; q < nsv; q++) {
        double s = 0;
        for (int i = 0; i < N; i++) for (int j = 0; j < N; j++) s += SV[q][i] * G[i][j] * SV[q][j];
        nrm[q] = s;
    }
    /* partial selection of the SVK shortest */
    static int used[NSV]; memset(used, 0, sizeof(used));
    for (m = 0; m < SVK; m++) {
        int b = -1; for (int q = 0; q < nsv; q++) if (!used[q] && (b < 0 || nrm[q] < nrm[b])) b = q;
        used[b] = 1; ord[m] = b;
    }
    static BeamSt cur[BEAM * SVK], nxt[BEAM * SVK];
    int nc = 1; memset(&cur[0], 0, sizeof(BeamSt));
    for (int k = 0; k < N; k++) {
        int nn = 0;
        for (int b = 0; b < nc; b++) for (int q = 0; q < SVK; q++) {
            const int *v = SV[ord[q]];
            BeamSt *t = &nxt[nn];
            *t = cur[b]; t->idx[k] = ord[q];
            /* extend the Cholesky factor of the k+1 prefix Gram */
            double g[N];
            for (int j = 0; j <= k; j++) {
                const int *w = SV[t->idx[j]];
                double s = 0; for (int a = 0; a < N; a++) for (int c = 0; c < N; c++) s += v[a] * G[a][c] * w[c];
                g[j] = s;
            }
            for (int j = 0; j < k; j++) {
                double s = g[j]; for (int l = 0; l < j; l++) s -= t->L[k][l] * t->L[j][l];
                t->L[k][j] = s / t->L[j][j];
            }
            double d = g[k]; for (int l = 0; l < k; l++) d -= t->L[k][l] * t->L[k][l];
            if (d < 1e-9 * nrm[ord[q]]) continue;
            t->L[k][k] = sqrt(d);
            double lk = 0; for (int l = 0; l <= k; l++) lk += 2 * log(t->L[l][l]);
            t->s = cur[b].s + (k < N - 1 ? lk : 0);
            if (k == N - 1) {
                int C[N][N];
                for (int i = 0; i < N; i++) for (int j = 0; j < N; j++) C[i][j] = SV[t->idx[j]][i];
                if (fabs(fabs(idet5(C)) - 1) > 1e-6) continue;
            }
            nn++;
        }
        if (nn == 0) return INFINITY;
        qsort(nxt, nn, sizeof(BeamSt), cmp_beam);
        nc = nn < BEAM ? nn : BEAM;
        memcpy(cur, nxt, sizeof(BeamSt) * nc);
    }
    for (int i = 0; i < N; i++) for (int j = 0; j < N; j++) U->u[i][j] = SV[cur[0].idx[j]][i];
    return log_eta_U(U, G);
}
static int use_beam = 1;

/* ---------------- affine enclosure of G(z) on a box ---------------- */
typedef struct { double c, a[NZ], e, ea[NZ]; } AF;     /* c + a.(z - zc) + [-e, e]; ea: e attributed per coordinate (heuristic) */

static double af_rad(const AF *x, const double *hw) {
    double s = 0; for (int t = 0; t < NZ; t++) s += fabs(x->a[t]) * hw[t]; return s;
}
static AF af_mul(const AF *x, const AF *y, const double *hw) {
    AF r; r.c = x->c * y->c;
    for (int t = 0; t < NZ; t++) r.a[t] = x->c * y->a[t] + y->c * x->a[t];
    double rx = af_rad(x, hw), ry = af_rad(y, hw);
    r.e = rx * ry + fabs(x->c) * y->e + fabs(y->c) * x->e + (rx + x->e) * y->e + x->e * ry;
    double fx = fabs(y->c) + ry + y->e, fy = fabs(x->c) + rx;
    for (int t = 0; t < NZ; t++)
        r.ea[t] = 0.5 * (fabs(x->a[t]) * hw[t] * ry + rx * fabs(y->a[t]) * hw[t]) + fx * x->ea[t] + fy * y->ea[t];
    return r;
}

/* ---------------- G(z) on a cell ---------------- */
static const int MI[6] = {1, 2, 2, 3, 3, 3}, MJ[6] = {0, 0, 1, 0, 1, 2};   /* z[8 + q] = mu[MI[q]][MJ[q]] */
static const double RCOEF[4][3] = {{1, 0, 0}, {-1, 1, 0}, {0, -1, 1}, {0, 0, -1}};   /* log R_k in (y_1, y_2, y_3) */

static AF af_const(double c) { AF r; memset(&r, 0, sizeof(r)); r.c = c; return r; }
static AF af_var(int t, double c) { AF r = af_const(c); r.a[t] = 1; return r; }
static AF af_add(const AF *x, const AF *y) {
    AF r; r.c = x->c + y->c; r.e = x->e + y->e;
    for (int t = 0; t < NZ; t++) { r.a[t] = x->a[t] + y->a[t]; r.ea[t] = x->ea[t] + y->ea[t]; }
    return r;
}
/* exp(sum_t coef_t z_t): remainder e^D - 1 - D <= e^rho - 1 - rho for |D| <= rho */
static AF af_exp(const double *coef, const double *zc, const double *hw) {
    double l = 0, rho = 0, q = 0;
    for (int t = 0; t < NZ; t++) { l += coef[t] * zc[t]; rho += fabs(coef[t]) * hw[t]; q += coef[t] * coef[t] * hw[t] * hw[t]; }
    AF r = af_const(exp(l));
    for (int t = 0; t < NZ; t++) r.a[t] = r.c * coef[t];
    r.e = r.c * (exp(rho) - 1 - rho);
    for (int t = 0; t < NZ; t++) r.ea[t] = q > 0 ? r.e * coef[t] * coef[t] * hw[t] * hw[t] / q : 0;
    return r;
}

static void s_point(const double *z, double S[4][4]) {
    double mu[4][4] = {{0}}, R[4];
    for (int i = 0; i < 4; i++) mu[i][i] = 1;
    for (int q = 0; q < 6; q++) mu[MI[q]][MJ[q]] = z[8 + q];
    for (int k = 0; k < 4; k++) R[k] = exp(RCOEF[k][0] * z[5] + RCOEF[k][1] * z[6] + RCOEF[k][2] * z[7]);
    for (int i = 0; i < 4; i++) for (int j = 0; j < 4; j++) {
        double s = 0; for (int k = 0; k < 4; k++) s += mu[i][k] * mu[j][k] * R[k]; S[i][j] = s;
    }
}

static void gram_at(const double *z, double G[N][N]) {
    double S[4][4]; s_point(z, S);
    double s = pow(z[0], -0.25);
    G[0][0] = z[0];
    for (int i = 0; i < 4; i++) G[0][1 + i] = G[1 + i][0] = z[0] * z[1 + i];
    for (int i = 0; i < 4; i++) for (int j = 0; j < 4; j++) G[1 + i][1 + j] = s * S[i][j] + z[0] * z[1 + i] * z[1 + j];
}

static AF SA[4][4];   /* affine form of r^{-1/4} S on the last cell (set by gram_affine) */
static void gram_affine(const Box *b, AF GA[N][N], double Gc[N][N], double *hw) {
    double zc[NZ];
    for (int t = 0; t < NZ; t++) { zc[t] = 0.5 * (b->lo[t] + b->hi[t]); hw[t] = 0.5 * (b->hi[t] - b->lo[t]); }
    AF r1 = af_var(0, zc[0]), phi[4], mu[4][4], R[4];
    for (int i = 0; i < 4; i++) phi[i] = af_var(1 + i, zc[1 + i]);
    for (int i = 0; i < 4; i++) for (int j = 0; j < 4; j++) mu[i][j] = af_const(i == j);
    for (int q = 0; q < 6; q++) mu[MI[q]][MJ[q]] = af_var(8 + q, zc[8 + q]);
    for (int k = 0; k < 4; k++) {
        double coef[NZ] = {0};
        for (int a = 0; a < 3; a++) coef[5 + a] = RCOEF[k][a];
        R[k] = af_exp(coef, zc, hw);
    }
    /* s = r^{-1/4}, with the second-order remainder of the convex function on [lo, hi] */
    AF s = af_const(pow(zc[0], -0.25));
    s.a[0] = -0.25 * pow(zc[0], -1.25);
    s.e = 0.5 * (5.0 / 16) * pow(b->lo[0], -2.25) * hw[0] * hw[0];
    s.ea[0] = s.e;
    for (int i = 0; i < 4; i++) for (int j = i; j < 4; j++) {
        AF S = af_const(0);
        for (int k = 0; k <= i; k++) {      /* mu is unit lower triangular, i <= j */
            AF p = af_mul(&mu[i][k], &mu[j][k], hw);
            AF q = af_mul(&p, &R[k], hw);
            S = af_add(&S, &q);
        }
        SA[i][j] = SA[j][i] = af_mul(&s, &S, hw);
    }
    GA[0][0] = r1;
    for (int i = 0; i < 4; i++) GA[0][1 + i] = af_mul(&r1, &phi[i], hw);
    for (int i = 0; i < 4; i++) for (int j = i; j < 4; j++) {
        AF rm = af_mul(&r1, &phi[i], hw);
        AF rmm = af_mul(&rm, &phi[j], hw);
        GA[1 + i][1 + j] = af_add(&SA[i][j], &rmm);
    }
    for (int i = 0; i < N; i++) for (int j = i; j < N; j++) { Gc[i][j] = Gc[j][i] = GA[i][j].c; GA[j][i] = GA[i][j]; }
}

/* ---------------- (B): balls about D5 and D5^* ---------------- */
static double TG[2][N][N], TK[2][N][N];     /* targets at det 1 (0: D5 Cartan, 1: D5^*) and T^{-1/2} */
static const double BALLR = 0.1;            /* Theorem D5 (invariant metric), and its dual */

static void jacobi5(double A0[N][N], double w[N], double V[N][N]) {
    double A[N][N];
    memcpy(A, A0, sizeof(A));
    for (int i = 0; i < N; i++) for (int j = 0; j < N; j++) V[i][j] = (i == j);
    for (int sweep = 0; sweep < 100; sweep++) {
        double off = 0;
        for (int p = 0; p < N; p++) for (int q = p + 1; q < N; q++) off += A[p][q] * A[p][q];
        if (off < 1e-30) break;
        for (int p = 0; p < N; p++) for (int q = p + 1; q < N; q++) {
            if (fabs(A[p][q]) < 1e-300) continue;
            double th = (A[q][q] - A[p][p]) / (2 * A[p][q]);
            double t = (th >= 0 ? 1 : -1) / (fabs(th) + sqrt(th * th + 1));
            double c = 1 / sqrt(t * t + 1), s = t * c;
            for (int k = 0; k < N; k++) { double a = A[k][p], b = A[k][q]; A[k][p] = c * a - s * b; A[k][q] = s * a + c * b; }
            for (int k = 0; k < N; k++) { double a = A[p][k], b = A[q][k]; A[p][k] = c * a - s * b; A[q][k] = s * a + c * b; }
            for (int k = 0; k < N; k++) { double a = V[k][p], b = V[k][q]; V[k][p] = c * a - s * b; V[k][q] = s * a + c * b; }
        }
    }
    for (int i = 0; i < N; i++) w[i] = A[i][i];
}

static void init_targets(void) {
    double D5[N][N] = {{2, -1, 0, 0, 0}, {-1, 2, -1, 0, 0}, {0, -1, 2, -1, -1}, {0, 0, -1, 2, 0}, {0, 0, -1, 0, 2}};
    double sc = pow(4.0, -0.2);     /* det 4 -> 1 */
    for (int i = 0; i < N; i++) for (int j = 0; j < N; j++) TG[0][i][j] = D5[i][j] * sc;
    /* D5^* = Z^5 u (Z^5 + h), h = (1/2,...,1/2): basis e_1..e_4, h; det 1/4 -> 1 */
    for (int i = 0; i < N; i++) for (int j = 0; j < N; j++)
        TG[1][i][j] = ((i == j) ? (i < 4 ? 1.0 : 1.25) : ((i == 4 || j == 4) ? 0.5 : 0.0)) / sc;
    for (int k = 0; k < 2; k++) {
        double w[N], V[N][N];
        jacobi5(TG[k], w, V);
        for (int i = 0; i < N; i++) for (int j = 0; j < N; j++) {
            double s = 0; for (int l = 0; l < N; l++) s += V[i][l] * V[j][l] / sqrt(w[l]); TK[k][i][j] = s;
        }
    }
}

/* invariant distance from U G U^T to target k (rows of U are the new basis in coordinates of G) */
static double dist_target(int k, signed char U[N][N], double G[N][N]) {
    double M[N][N], Y[N][N], T[N][N], w[N], V[N][N];
    for (int a = 0; a < N; a++) for (int b = 0; b < N; b++) {
        double s = 0; for (int i = 0; i < N; i++) for (int j = 0; j < N; j++) s += U[a][i] * G[i][j] * U[b][j]; M[a][b] = s;
    }
    for (int a = 0; a < N; a++) for (int b = 0; b < N; b++) { double s = 0; for (int c = 0; c < N; c++) s += TK[k][a][c] * M[c][b]; T[a][b] = s; }
    for (int a = 0; a < N; a++) for (int b = 0; b < N; b++) { double s = 0; for (int c = 0; c < N; c++) s += T[a][c] * TK[k][c][b]; Y[a][b] = s; }
    jacobi5(Y, w, V);
    double d = 0;
    for (int i = 0; i < N; i++) { if (!(w[i] > 0)) return INFINITY; d += log(w[i]) * log(w[i]); }
    return sqrt(d);
}

/* short vectors x (one of each +-x) with x^T G x <= B, Fincke-Pohst */
#define MAXSV 600
static int svx[MAXSV][N], nsvx; static double svn[MAXSV];
static double FPq[N][N];
static void fp_rec(int i, int *x, double rem, double G[N][N]) {
    if (nsvx >= MAXSV) return;
    if (i < 0) {
        int nz = 0, first = 0;
        for (int k = 0; k < N; k++) if (x[k]) { nz = 1; first = x[k]; break; }
        if (!nz || first < 0) return;
        double q = 0; for (int a = 0; a < N; a++) for (int b = 0; b < N; b++) q += x[a] * G[a][b] * x[b];
        memcpy(svx[nsvx], x, sizeof(int) * N); svn[nsvx++] = q;
        return;
    }
    double c = 0; for (int j = i + 1; j < N; j++) c += FPq[i][j] * x[j];
    double rr = sqrt(fmax(rem, 0) / FPq[i][i]);
    int lo = (int)ceil(-c - rr - 1e-9), hi = (int)floor(-c + rr + 1e-9);
    for (int v = lo; v <= hi; v++) {
        x[i] = v; double t = v + c;
        fp_rec(i - 1, x, rem - FPq[i][i] * t * t, G);
    }
    x[i] = 0;
}
static void short_vectors(double G[N][N], double B) {
    double A[N][N]; memcpy(A, G, sizeof(A));
    memset(FPq, 0, sizeof(FPq));
    for (int i = 0; i < N; i++) for (int j = i; j < N; j++) FPq[i][j] = A[i][j];
    for (int i = 0; i < N; i++) {
        for (int j = i + 1; j < N; j++) { FPq[j][i] = FPq[i][j]; FPq[i][j] /= FPq[i][i]; }
        for (int k = i + 1; k < N; k++) for (int l = k; l < N; l++) FPq[k][l] -= FPq[k][i] * FPq[i][l];
    }
    int x[N] = {0}; nsvx = 0;
    fp_rec(N - 1, x, B, G);
}

/* match a basis of the lattice with Gram G to target k: backtrack over short vectors on the target's Gram pattern */
static double mbest; static signed char mU[N][N]; static int mq[N], msg[N]; static long mnodes;
static void match_rec(int k, int tk, double G[N][N], double tol) {
    if (++mnodes > 200000) return;
    if (k == N) {
        signed char U[N][N]; int C[N][N];
        for (int a = 0; a < N; a++) for (int i = 0; i < N; i++) { U[a][i] = msg[a] * svx[mq[a]][i]; C[i][a] = U[a][i]; }
        if (fabs(fabs(idet5(C)) - 1) > 1e-6) return;
        double d = dist_target(tk, U, G);
        if (d < mbest) { mbest = d; memcpy(mU, U, sizeof(U)); }
        return;
    }
    for (int q = 0; q < nsvx; q++) {
        if (fabs(svn[q] - TG[tk][k][k]) > tol * TG[tk][k][k]) continue;
        double ipj[N];
        for (int j = 0; j < k; j++) {
            double ip = 0;
            for (int a = 0; a < N; a++) for (int b = 0; b < N; b++) ip += svx[q][a] * G[a][b] * svx[mq[j]][b];
            ipj[j] = ip * msg[j];
        }
        for (int sg = 1; sg >= -1; sg -= 2) {
            if (sg < 0 && k == 0) continue;
            int ok = 1;
            for (int j = 0; j < k && ok; j++) if (fabs(sg * ipj[j] - TG[tk][k][j]) > tol * sqrt(TG[tk][k][k] * TG[tk][j][j])) ok = 0;
            if (!ok) continue;
            mq[k] = q; msg[k] = sg;
            match_rec(k + 1, tk, G, tol);
        }
        /* D5: every first-vector candidate is a root of any copy within 2 BALLR (its T-norm is in [1.27, 2.98], and
           the norms of D5 are even), and Aut(D5) is transitive on the roots: one candidate finds every such copy.
           (Not for D5^*: its norms 1 and 5/4 are not separated by the tolerance.) */
        if (k == 0 && tk == 0) break;
    }
}
/* best copy of target k near G: returns the distance and the basis in U (rows).  Only copies at distance < 2 BALLR are
   used, and there every eigenvalue of T^{-1/2} U G U^T T^{-1/2} is within a factor e^{2 BALLR} of 1, so
   |x^T (UGU^T - T) y| <= (e^{2 BALLR} - 1) sqrt(x^T T x y^T T y): the search prunes with exactly that tolerance */
static double find_copy_direct(int tk, double G[N][N], signed char U[N][N]) {
    double tmax = 0, tol = expm1(2 * BALLR); for (int i = 0; i < N; i++) tmax = fmax(tmax, TG[tk][i][i]);
    short_vectors(G, tmax * (1 + tol) * 1.001);
    mbest = INFINITY; mnodes = 0;
    match_rec(0, tk, G, tol);
    if (mbest < INFINITY) memcpy(U, mU, sizeof(mU));
    return mbest;
}
static int inv_d(double A[N][N], double X[N][N]) {      /* Gauss-Jordan with partial pivoting */
    double M[N][2 * N];
    for (int i = 0; i < N; i++) for (int j = 0; j < N; j++) { M[i][j] = A[i][j]; M[i][N + j] = (i == j); }
    for (int c = 0; c < N; c++) {
        int p = c; for (int i = c + 1; i < N; i++) if (fabs(M[i][c]) > fabs(M[p][c])) p = i;
        if (!(fabs(M[p][c]) > 1e-300)) return 0;
        for (int j = 0; j < 2 * N; j++) { double t = M[c][j]; M[c][j] = M[p][j]; M[p][j] = t; }
        double d = M[c][c]; for (int j = 0; j < 2 * N; j++) M[c][j] /= d;
        for (int i = 0; i < N; i++) if (i != c) { double f = M[i][c]; for (int j = 0; j < 2 * N; j++) M[i][j] -= f * M[c][j]; }
    }
    for (int i = 0; i < N; i++) for (int j = 0; j < N; j++) X[i][j] = M[i][N + j];
    return 1;
}
/* D5^* = (D5)^*: T1 = BD T0^{-1} BD^T (BD found once, at d ~ 0), so a basis V with V G^{-1} V^T ~ T0 gives
   U = BD V^{-T} with U G U^T ~ T1 at the same invariant distance (it is invariant under congruence and inversion):
   the D5^* search runs as a D5 search on G^{-1}, with its symmetry break */
static signed char BD[N][N]; static int have_bd = 0;
static double find_copy(int tk, double G[N][N], signed char U[N][N]) {
    if (tk == 0 || !have_bd) return find_copy_direct(tk, G, U);
    double Gi[N][N], Vd[N][N], Vi[N][N]; signed char V[N][N];
    if (!inv_d(G, Gi) || find_copy_direct(0, Gi, V) == INFINITY) return INFINITY;
    for (int i = 0; i < N; i++) for (int j = 0; j < N; j++) Vd[i][j] = V[i][j];
    if (!inv_d(Vd, Vi)) return INFINITY;
    int W[N][N];                                   /* W = V^{-T}, integral since V is unimodular */
    for (int i = 0; i < N; i++) for (int j = 0; j < N; j++) W[i][j] = (int)lround(Vi[j][i]);
    for (int i = 0; i < N; i++) for (int j = 0; j < N; j++) {    /* check W V^T = I exactly */
        int s = 0; for (int k = 0; k < N; k++) s += W[i][k] * V[j][k];
        if (s != (i == j)) return INFINITY;
    }
    for (int i = 0; i < N; i++) for (int j = 0; j < N; j++) {
        int s = 0; for (int k = 0; k < N; k++) s += BD[i][k] * W[k][j];
        if (s < -127 || s > 127) return INFINITY;
        U[i][j] = (signed char)s;
    }
    return dist_target(1, U, G);
}
static void init_bd(void) {
    double T0i[N][N];
    if (!inv_d(TG[0], T0i)) return;
    double d = find_copy_direct(1, T0i, BD);
    have_bd = d < 1e-9;
    if (!have_bd) fprintf(stderr, "init_bd: no D5^* basis of D5^{-1} (d = %g); D5^* search runs directly\n", d);
}

/* sup over the cell of delta = |T^{-1/2}(U G U^T - T)T^{-1/2}|_F (entrywise bounds), d <= delta/(1-delta) */
static int cell_in_ball(int tk, signed char U[N][N], AF GA[N][N], double Gc[N][N], const double *hw) {
    double V[N][N];      /* V = T^{-1/2} U */
    for (int c = 0; c < N; c++) for (int i = 0; i < N; i++) { double s = 0; for (int a = 0; a < N; a++) s += TK[tk][c][a] * U[a][i]; V[c][i] = s; }
    double d2 = 0;
    for (int c = 0; c < N; c++) for (int d = c; d < N; d++) {
        double cen = -(c == d), a[NZ] = {0}, e = 0;
        for (int i = 0; i < N; i++) for (int j = 0; j < N; j++) {
            double w = V[c][i] * V[d][j];
            if (w == 0) continue;
            cen += w * Gc[i][j]; e += fabs(w) * GA[i][j].e;
            for (int t = 0; t < NZ; t++) a[t] += w * GA[i][j].a[t];
        }
        double m = fabs(cen) + e; for (int t = 0; t < NZ; t++) m += fabs(a[t]) * hw[t];
        d2 += (c == d ? 1 : 2) * m * m;
    }
    return sqrt(d2) <= BALLR / (1 + BALLR);
}

/* upper bound of const + sum_ij Q_ij (G_ij - Gc_ij) over the box; also the per-coordinate sensitivity */
static double lin_upper(double Q[N][N], AF GA[N][N], const double *hw, double *sens) {
    double a[NZ] = {0}, e = 0, ea[NZ] = {0};
    for (int i = 0; i < N; i++) for (int j = 0; j < N; j++) {
        for (int t = 0; t < NZ; t++) { a[t] += Q[i][j] * GA[i][j].a[t]; ea[t] += fabs(Q[i][j]) * GA[i][j].ea[t]; }
        e += fabs(Q[i][j]) * GA[i][j].e;
    }
    double s = e;
    for (int t = 0; t < NZ; t++) { s += fabs(a[t]) * hw[t]; if (sens) sens[t] = fabs(a[t]) * hw[t] + ea[t]; }
    return s;
}

static double last_sens[NZ];

/* best upper bound for log eta over the box from flag U (Psi and split forms); updates last_sens */
static double flag_bound(const Flag *U, AF GA[N][N], double Gc[N][N], const double *hw, double *sens) {
    double ld[N + 1], L[N][N], P[N + 1][N][N], Q[N][N];
    if (!prefix(U, Gc, ld, L)) return INFINITY;
    prefix_grads(U, L, P);
    double best = INFINITY, s[NZ];
    for (int i = 0; i < N; i++) for (int j = 0; j < N; j++) Q[i][j] = P[1][i][j] + P[2][i][j] + P[3][i][j] + P[4][i][j];
    double v = ld[1] + ld[2] + ld[3] + ld[4] + lin_upper(Q, GA, hw, s);   /* det G = 1 */
    if (v < best) best = v;
    if (sens) memcpy(sens, s, sizeof(s));       /* split heuristic follows the Psi form */
    for (int k = 1; k < N; k++) {
        for (int i = 0; i < N; i++) for (int j = 0; j < N; j++) Q[i][j] = 2.5 * P[k][i][j];
        v = LOGETA[k] + LOGETA[N - k] + 2.5 * ld[k] + lin_upper(Q, GA, hw, s);
        if (v < best) best = v;
    }
    return best;
}

/* linear forms c + a.dz (+ e) of every bound attached to flag U (Psi and the four split forms) */
typedef struct { double c, a[NZ], e; } Row;
static int flag_rows(const Flag *U, AF GA[N][N], double Gc[N][N], Row *out) {
    double ld[N + 1], L[N][N], P[N + 1][N][N], Q[N][N];
    if (!prefix(U, Gc, ld, L)) return 0;
    prefix_grads(U, L, P);
    for (int f = 0; f < 5; f++) {
        Row *r = &out[f];
        double w = f == 0 ? 1 : 2.5;
        for (int i = 0; i < N; i++) for (int j = 0; j < N; j++)
            Q[i][j] = f == 0 ? P[1][i][j] + P[2][i][j] + P[3][i][j] + P[4][i][j] : w * P[f][i][j];
        r->c = f == 0 ? ld[1] + ld[2] + ld[3] + ld[4] : LOGETA[f] + LOGETA[N - f] + 2.5 * ld[f];
        memset(r->a, 0, sizeof(r->a)); r->e = 0;
        for (int i = 0; i < N; i++) for (int j = 0; j < N; j++) {
            for (int t = 0; t < NZ; t++) r->a[t] += Q[i][j] * GA[i][j].a[t];
            r->e += fabs(Q[i][j]) * GA[i][j].e;
        }
    }
    return 5;
}

/* LP: is there dz in the box with c_F + a_F.dz + e_F > log c for every row F?  1 = no (certified). */
static long st_lp = 0, st_lpcert = 0;
static int lp_joint(Row *R, int nr, const double *hw, double *wit) {
    glp_prob *lp = glp_create_prob();
    glp_add_cols(lp, NZ);
    for (int t = 0; t < NZ; t++) {
        if (hw[t] > 0) glp_set_col_bnds(lp, t + 1, GLP_DB, -hw[t], hw[t]);
        else glp_set_col_bnds(lp, t + 1, GLP_FX, 0, 0);
    }
    glp_add_rows(lp, nr);
    static int ia[1 + 400 * NZ], ja[1 + 400 * NZ]; static double ar[1 + 400 * NZ];
    int ne = 0;
    for (int q = 0; q < nr; q++) {
        glp_set_row_bnds(lp, q + 1, GLP_LO, logc - R[q].c - R[q].e, 0);
        for (int t = 0; t < NZ; t++) if (R[q].a[t] != 0) { ne++; ia[ne] = q + 1; ja[ne] = t + 1; ar[ne] = R[q].a[t]; }
    }
    glp_load_matrix(lp, ne, ia, ja, ar);
    glp_smcp parm; glp_init_smcp(&parm); parm.msg_lev = GLP_MSG_OFF;
    int ret = glp_simplex(lp, &parm), res = 0;
    st_lp++;
    if (ret == 0 && glp_get_prim_stat(lp) == GLP_NOFEAS) res = 1;
    else if (ret == 0) for (int t = 0; t < NZ; t++) wit[t] = glp_get_col_prim(lp, t + 1);
    glp_delete_prob(lp);
    return res;
}
static int use_lp = 1, ncuts = 4;
/* min over lambda in the simplex of max over the cell of sum_i lambda_i (c_i + a_i.dz + e_i): an LP in (lambda, u),
   minimise sum lambda_i (c_i + e_i) + sum_t hw_t u_t subject to |sum_i lambda_i a_it| <= u_t.  Its value is < log c
   exactly when lp_joint finds no point (LP duality); lambda is the certificate of a P leaf. */
static double lp_weights(Row *R, int nr, const double *hw, double *lam) {
    glp_prob *lp = glp_create_prob();
    glp_set_obj_dir(lp, GLP_MIN);
    glp_add_cols(lp, nr + NZ);
    for (int i = 0; i < nr; i++) { glp_set_col_bnds(lp, i + 1, GLP_LO, 0, 0); glp_set_obj_coef(lp, i + 1, R[i].c + R[i].e); }
    for (int t = 0; t < NZ; t++) { glp_set_col_bnds(lp, nr + t + 1, GLP_LO, 0, 0); glp_set_obj_coef(lp, nr + t + 1, hw[t]); }
    glp_add_rows(lp, 1 + 2 * NZ);
    static int ia[1 + 401 * (2 * NZ + 1)], ja[1 + 401 * (2 * NZ + 1)]; static double ar[1 + 401 * (2 * NZ + 1)];
    int ne = 0;
    glp_set_row_bnds(lp, 1, GLP_FX, 1, 1);
    for (int i = 0; i < nr; i++) { ne++; ia[ne] = 1; ja[ne] = i + 1; ar[ne] = 1; }
    for (int t = 0; t < NZ; t++) for (int sg = 0; sg < 2; sg++) {
        int row = 2 + 2 * t + sg;
        glp_set_row_bnds(lp, row, GLP_UP, 0, 0);
        for (int i = 0; i < nr; i++) if (R[i].a[t] != 0) { ne++; ia[ne] = row; ja[ne] = i + 1; ar[ne] = sg ? -R[i].a[t] : R[i].a[t]; }
        ne++; ia[ne] = row; ja[ne] = nr + t + 1; ar[ne] = -1;
    }
    glp_load_matrix(lp, ne, ia, ja, ar);
    glp_smcp parm; glp_init_smcp(&parm); parm.msg_lev = GLP_MSG_OFF;
    double v = INFINITY;
    if (glp_simplex(lp, &parm) == 0 && glp_get_status(lp) == GLP_OPT) {
        v = glp_get_obj_val(lp);
        for (int i = 0; i < nr; i++) lam[i] = glp_get_col_prim(lp, i + 1);
    }
    glp_delete_prob(lp);
    return v;
}

static void add_flag(Box *b, const Flag *U) {
    for (int q = 0; q < b->nf; q++) if (!memcmp(&b->f[q], U, sizeof(Flag))) return;
    if (b->nf < FMAX) b->f[b->nf++] = *U;
    else { memmove(&b->f[0], &b->f[1], sizeof(Flag) * (FMAX - 1)); b->f[FMAX - 1] = *U; }
}

/* b_1 shortest.  For a tail vector s = sum c_i s_i (c in Z^4), the lattice vectors s + (c.mu + j) b_1 have norm
 * |s|^2 + r1 (c.mu + j)^2 >= r1, so dist(c.mu, Z) >= sqrt(1 - |s|^2 / r1) whenever |s|^2 < r1.  c.mu is linear in
 * mu, so the test is exact in mu; |s|^2 is bounded above on the box via the affine form of S. */
#define NTV 320
static int TV[NTV][4], ntv = 0;
static long st_short = 0;
static int use_shortcut = 1;
static int TVR = 1;     /* tail vectors c in {-TVR..TVR}^4 (TV2: TVR = 2) */
static void gen_tv(void) {
    int b = 2 * TVR + 1, tot = b * b * b * b;
    for (int code = 0; code < tot; code++) {
        int c[4], x = code, first = -1;
        for (int i = 0; i < 4; i++) { c[i] = x % b - TVR; x /= b; }
        for (int i = 0; i < 4; i++) if (c[i]) { first = i; break; }
        if (first < 0 || c[first] < 0) continue;
        memcpy(TV[ntv++], c, sizeof(c));
    }
}
static int not_shortest(const Box *b, const double *hw) {
    double r1lo = b->lo[0];
    for (int q = 0; q < ntv; q++) {
        const int *c = TV[q];
        /* upper bound on c^T S c */
        double cen = 0, rad = 0, err = 0;
        for (int i = 0; i < 4; i++) for (int j = 0; j < 4; j++) {
            if (!c[i] || !c[j]) continue;
            double w = c[i] * c[j];
            cen += w * SA[i][j].c; err += fabs(w) * SA[i][j].e;
        }
        double a[NZ] = {0};
        for (int i = 0; i < 4; i++) for (int j = 0; j < 4; j++) if (c[i] && c[j])
            for (int t = 0; t < NZ; t++) a[t] += c[i] * c[j] * SA[i][j].a[t];
        for (int t = 0; t < NZ; t++) rad += fabs(a[t]) * hw[t];
        double s2 = cen + rad + err;
        if (s2 >= r1lo) continue;
        double tau = sqrt(1 - s2 / r1lo);
        /* range of c.mu */
        double lo = 0, hi = 0;
        for (int i = 0; i < 4; i++) {
            if (c[i] > 0) { lo += c[i] * b->lo[1 + i]; hi += c[i] * b->hi[1 + i]; }
            else if (c[i] < 0) { lo += c[i] * b->hi[1 + i]; hi += c[i] * b->lo[1 + i]; }
        }
        double k = floor(0.5 * (lo + hi) + 0.5);
        if (lo > k - tau && hi < k + tau) { st_short++; leaf_t = LF_N; leaf_k = (int)k; memcpy(leaf_c, c, sizeof(leaf_c)); return 1; }
    }
    return 0;
}

static int use_inherit = 0, use_cellrank = 0;
/* local descent (same moves as local_search) on the flag's upper bound over the cell instead of its centre value */
static double cell_search(Flag *U, AF GA[N][N], double Gc[N][N], const double *hw) {
    double best = flag_bound(U, GA, Gc, hw, NULL);
    for (int it = 0; it < 30; it++) {
        Flag bu = *U; double bv = best;
        for (int i = 0; i < N; i++) for (int j = 0; j < N; j++) if (i != j) for (int s = -1; s <= 1; s += 2) {
            Flag V = *U;
            for (int r = 0; r < N; r++) V.u[r][i] += s * U->u[r][j];
            if (!flag_ok(&V)) continue;
            double v = flag_bound(&V, GA, Gc, hw, NULL);
            if (v < bv - 1e-12) { bv = v; bu = V; }
        }
        if (bv >= best - 1e-12) break;
        best = bv; *U = bu;
    }
    return best;
}

static double last_centre;   /* best log eta found at the centre */
static int use_ball = 1;
static double BALLTRY;       /* look for a D5 / D5^* copy only when the centre's log eta exceeds this */
static long st_lpc = 0, st_balltry = 0;

static int classify(Box *b) {
    AF GA[N][N]; double Gc[N][N], hw[NZ];
    gram_affine(b, GA, Gc, hw);
    if (use_shortcut && not_shortest(b, hw)) return 1;
    double best = INFINITY, s[NZ];
    for (int q = 0; q < b->nf; q++) {
        double v = flag_bound(&b->f[q], GA, Gc, hw, s);
        if (v < best) best = v;
        if (v <= logc) { st_cert++; leaf_t = LF_F; leaf_flag = b->f[q]; return 1; }
    }
    /* an inherited ball copy */
    if (use_ball) for (int q = 0; q < b->nU; q++) if (cell_in_ball(b->Ut[q], b->U[q], GA, Gc, hw)) {
        st_ball++; leaf_t = LF_B; leaf_tk = b->Ut[q]; memcpy(leaf_W, b->U[q], sizeof(leaf_W)); return 1;
    }
    /* new flag from a local search at the centre, started from the best inherited flag (INHERIT; the permutations
       only when there is none) or from the best permutation */
    Flag bu = PERM[0]; double bv = INFINITY;
    if (!use_inherit || b->nf == 0)
        for (int q = 0; q < nperm; q++) { double v = log_eta_U(&PERM[q], Gc); if (v < bv) { bv = v; bu = PERM[q]; } }
    for (int q = 0; q < b->nf; q++) { double v = log_eta_U(&b->f[q], Gc); if (v < bv) { bv = v; bu = b->f[q]; } }
    if (use_beam) {
        Flag bb; double vb = beam_search(Gc, &bb);
        if (vb < log_eta_U(&bu, Gc)) bu = bb;
    }
    last_centre = local_search(&bu, Gc);
    if (use_cellrank) cell_search(&bu, GA, Gc, hw);   /* descend on the bound over the cell, from the centre's best */
    double v = flag_bound(&bu, GA, Gc, hw, s);
    add_flag(b, &bu);
    if (v < best) best = v;
    memcpy(last_sens, s, sizeof(s));
    if (v <= logc) { st_cert++; leaf_t = LF_F; leaf_flag = bu; return 1; }
    if (best == INFINITY) for (int t = 0; t < NZ; t++) last_sens[t] = hw[t];
    /* a new ball copy near the centre */
    if (use_ball && last_centre > BALLTRY && b->nU < 2) {
        st_balltry++;
        for (int tk = 0; tk < 2; tk++) {
            signed char U[N][N];
            double d = find_copy(tk, Gc, U);
            if (d < 2 * BALLR) {
                int q = b->nU < 2 ? b->nU++ : 1;
                memcpy(b->U[q], U, sizeof(U)); b->Ut[q] = tk;
                if (verbose > 1) fprintf(stderr, "  copy of %s at d = %.4f (centre log eta %.4f)\n", tk ? "D5*" : "D5", d, last_centre);
                if (cell_in_ball(tk, U, GA, Gc, hw)) { st_ball++; leaf_t = LF_B; leaf_tk = tk; memcpy(leaf_W, U, sizeof(leaf_W)); return 1; }
            }
        }
    }
    if (!use_lp) return 0;
    /* joint certificate over all flags of the cell, with cutting planes from local searches at LP witnesses */
    static Row R[400]; static Flag RF[400];
    int nr = 0;
    for (int q = 0; q < b->nf && nr + 5 <= 400; q++) { int k0 = nr; nr += flag_rows(&b->f[q], GA, Gc, &R[nr]); for (int i = k0; i < nr; i++) RF[i] = b->f[q]; }
    for (int round = 0; round <= ncuts; round++) {
        double wit[NZ];
        if (lp_joint(R, nr, hw, wit)) {
            /* the certificate: weights lambda on the rows with max over the cell of sum lambda_i row_i < log c */
            static double lam[400];
            if (lp_weights(R, nr, hw, lam) <= logc) {
                leaf_t = LF_P; leaf_m = 0;
                for (int i = 0; i < nr; i++) if (lam[i] > 1e-12 && leaf_m < PMAX) {
                    leaf_pf[leaf_m] = RF[i]; leaf_pk[leaf_m] = i % 5; leaf_pl[leaf_m] = lam[i]; leaf_m++;
                }
                st_cert++; st_lpc++; return 1;
            }
            return 0;
        }
        if (round == ncuts || nr + 5 > 400) break;
        double Gw[N][N], zw[NZ];
        for (int t = 0; t < NZ; t++) zw[t] = 0.5 * (b->lo[t] + b->hi[t]) + wit[t];
        gram_at(zw, Gw);
        Flag cu = bu; double cv = log_eta_U(&cu, Gw);
        for (int q = 0; q < b->nf; q++) { double v2 = log_eta_U(&b->f[q], Gw); if (v2 < cv) { cv = v2; cu = b->f[q]; } }
        if (use_beam) { Flag bb; double vb = beam_search(Gw, &bb); if (vb < cv) { cv = vb + 1e-9; cu = bb; } }
        double ls = local_search(&cu, Gw);
        if (ls >= cv - 1e-9 && use_beam) ls = cv - 1e-9;
        if (ls >= cv - 1e-9) break;
        add_flag(b, &cu);
        { int k0 = nr; nr += flag_rows(&cu, GA, Gc, &R[nr]); for (int i = k0; i < nr; i++) RF[i] = cu; }
    }
    return 0;
}

static double R1W = 1, MUW = 1;   /* optional pre-refinement of r and phi (off by default) */
static int split_dim(const Box *b) {
    int ts = -1; double bs = 0;
    if (b->hi[0] - b->lo[0] > R1W) { ts = 0; bs = (b->hi[0] - b->lo[0]) / R1W; }
    for (int i = 1; i <= 4; i++) { double w = (b->hi[i] - b->lo[i]) / MUW; if (w > 1 && w > bs) { bs = w; ts = i; } }
    if (ts >= 0) return ts;
    ts = 0; bs = -1;
    for (int t = 0; t < NZ; t++) {
        double w = b->hi[t] - b->lo[t];
        if (w <= 0) continue;
        double sc = last_sens[t] + 1e-6 * w;
        if (sc > bs) { bs = sc; ts = t; }
    }
    return ts;
}

/* flags are written in full once per log and then referred to as "#k" (k-th distinct flag written in this log) */
#define FTAB 65536
static Flag ftab[FTAB]; static int nftab = 0;
static void write_flag(const Flag *U) {
    for (int q = nftab - 1; q >= 0; q--) if (!memcmp(&ftab[q], U, sizeof(Flag))) { fprintf(clog, " #%d", q); return; }
    for (int i = 0; i < N; i++) for (int j = 0; j < N; j++) fprintf(clog, " %d", U->u[i][j]);
    if (nftab < FTAB) ftab[nftab++] = *U;
}
static void write_leaf(void) {
    switch (leaf_t) {
    case LF_N: fprintf(clog, "N %d %d %d %d %d\n", leaf_k, leaf_c[0], leaf_c[1], leaf_c[2], leaf_c[3]); break;
    case LF_F: fputs("F", clog); write_flag(&leaf_flag); fputs("\n", clog); break;
    case LF_P:
        fprintf(clog, "P %d\n", leaf_m);
        for (int i = 0; i < leaf_m; i++) { fprintf(clog, "%d %a", leaf_pk[i], leaf_pl[i]); write_flag(&leaf_pf[i]); fputs("\n", clog); }
        break;
    case LF_B:
        fprintf(clog, "B %d", leaf_tk);
        for (int i = 0; i < N; i++) for (int j = 0; j < N; j++) fprintf(clog, " %d", leaf_W[i][j]);
        fputs("\n", clog); break;
    default: fprintf(stderr, "leaf without certificate\n"); exit(2);
    }
    leaf_t = 0;
}
static long nboxes, maxboxes, maxdepth;
static long splits_at[NZ];
static double tcap; static clock_t tstart; static int timed_out;
static int run(Box *b, int depth) {
    if (++nboxes > maxboxes) return 0;
    if ((nboxes & 255) == 0 && (double)(clock() - tstart) / CLOCKS_PER_SEC > tcap) { timed_out = 1; return 0; }
    if (depth > maxdepth) maxdepth = depth;
    if (classify(b)) { if (clog) write_leaf(); return 1; }
    int ts = split_dim(b);
    splits_at[ts]++;
    if (clog) fprintf(clog, "S %d\n", ts);
    double mid = 0.5 * (b->lo[ts] + b->hi[ts]);
    Box *c = malloc(sizeof(Box));
    *c = *b; c->hi[ts] = mid; int ok = run(c, depth + 1);
    if (ok) { *c = *b; c->lo[ts] = mid; ok = run(c, depth + 1); }
    free(c);
    return ok;
}

static void perms(int k, int *p, int *used) {
    if (k == N) {
        Flag *U = &PERM[nperm++]; memset(U, 0, sizeof(*U));
        for (int j = 0; j < N; j++) U->u[p[j]][j] = 1;
        return;
    }
    for (int i = 0; i < N; i++) if (!used[i]) { used[i] = 1; p[k] = i; perms(k + 1, p, used); used[i] = 0; }
}

int main(int argc, char **argv) {
    if (argc < 2) { fprintf(stderr, "usage: %s COLFILE [first [count]]\n", argv[0]); return 1; }
    /* the producer closes a cell only with a margin below log 6; soundness rests on check_lift, not on this */
    logc = log(6.0) - (getenv("MARGIN") ? atof(getenv("MARGIN")) : 1e-9);
    if (getenv("VERBOSE")) verbose = atoi(getenv("VERBOSE"));
    use_beam = getenv("BEAM") != NULL;      /* descent flags (default) beat beam flags, see notes/stage2_prototype */
    R1W = 0.05;                             /* r pre-refined to width 0.05 (about 10% fewer cells) */
    if (getenv("NOLP")) use_lp = 0;
    if (getenv("NOBALL")) use_ball = 0;
    if (getenv("NOBEAM")) use_beam = 0;
    if (getenv("INHERIT")) use_inherit = 1;
    if (getenv("CELLRANK")) use_cellrank = 1;
    if (getenv("NOSHORT")) use_shortcut = 0;
    if (getenv("CUTS")) ncuts = atoi(getenv("CUTS"));
    if (getenv("MUW")) MUW = atof(getenv("MUW"));
    if (getenv("R1W")) R1W = atof(getenv("R1W"));
    BALLTRY = log(getenv("BALLTRY") ? atof(getenv("BALLTRY")) : 5.4);
    maxboxes = getenv("BOXCAP") ? atol(getenv("BOXCAP")) : 2000000;
    tcap = getenv("TCAP") ? atof(getenv("TCAP")) : 600;
    glp_term_out(GLP_OFF);
    if (getenv("TV2")) TVR = 2;
    gen_sv(); gen_tv(); init_targets(); init_bd();
    LOGETA[0] = 0; LOGETA[1] = 0; LOGETA[2] = log(2 / sqrt(3.0)); LOGETA[3] = log(1.5); LOGETA[4] = log(3.0);
    /* the root r-range must contain [2^{2/5}, 8^{1/5}] (check_lift verifies): rounded outward by two ulps */
    R1LO = nextafter(nextafter(pow(2.0, 0.4), 0), 0); R1HI = nextafter(nextafter(pow(8.0, 0.2), 2), 2);
    int p[N], used[N] = {0};
    perms(0, p, used);
    if (getenv("TESTZ")) {          /* print G at the point z (14 numbers) */
        double z[NZ], G[N][N]; char *p = getenv("TESTZ");
        for (int t = 0; t < NZ; t++) z[t] = strtod(p, &p);
        gram_at(z, G);
        for (int i = 0; i < N; i++) { for (int j = 0; j < N; j++) printf(" %.12f", G[i][j]); printf("\n"); }
        return 0;
    }
    if (!strcmp(argv[1], "test")) {      /* self-test: det 1, the affine enclosures, the copy finder */
        srand(1);
        double worst = 0, worstdet = 0;
        for (int it = 0; it < 2000; it++) {
            Box b; memset(&b, 0, sizeof(b));
            double cz[NZ] = {1.4, 0.2, -0.1, 0.3, 0.05, 0.3, 0.4, -0.15, 0.4, -0.5, -0.3, -0.5, 0.1, 0.4};
            for (int t = 0; t < NZ; t++) { double w = (t == 0 ? 0.05 : 0.12) * rand() / RAND_MAX; b.lo[t] = cz[t] - w; b.hi[t] = cz[t] + w; }
            AF GA[N][N]; double Gc[N][N], hw[NZ];
            gram_affine(&b, GA, Gc, hw);
            for (int k = 0; k < 20; k++) {
                double z[NZ], G[N][N];
                for (int t = 0; t < NZ; t++) z[t] = b.lo[t] + (b.hi[t] - b.lo[t]) * rand() / RAND_MAX;
                gram_at(z, G);
                double A[N][N], L[N][N]; memcpy(A, G, sizeof(A)); chol(N, A, L);
                double ld = 0; for (int i = 0; i < N; i++) ld += 2 * log(L[i][i]);
                worstdet = fmax(worstdet, fabs(ld));
                for (int i = 0; i < N; i++) for (int j = 0; j < N; j++) {
                    double lin = GA[i][j].c; for (int t = 0; t < NZ; t++) lin += GA[i][j].a[t] * (z[t] - 0.5 * (b.lo[t] + b.hi[t]));
                    double ex = fabs(G[i][j] - lin) - GA[i][j].e;
                    if (ex > worst) worst = ex;
                }
            }
        }
        printf("max |log det G| = %.2e; max enclosure violation = %.2e (should be <= ~1e-15)\n", worstdet, worst);
        for (int tk = 0; tk < 2; tk++) for (int it = 0; it < 5; it++) {
            int Ui[N][N]; signed char U[N][N]; double G[N][N];
            /* random unimodular: product of elementary operations */
            for (int i = 0; i < N; i++) for (int j = 0; j < N; j++) Ui[i][j] = (i == j);
            for (int k = 0; k < 12; k++) { int a = rand() % N, c = rand() % N; if (a == c) continue; int m = rand() % 3 - 1; for (int j = 0; j < N; j++) Ui[a][j] += m * Ui[c][j]; }
            for (int i = 0; i < N; i++) for (int j = 0; j < N; j++) {
                double s = 0; for (int a = 0; a < N; a++) for (int c = 0; c < N; c++) s += Ui[i][a] * TG[tk][a][c] * Ui[j][c]; G[i][j] = s;
            }
            double d = find_copy(tk, G, U);
            Flag fb; double le = beam_search(G, &fb); le = local_search(&fb, G);
            printf("%s copy %d: d = %.2e (svs %d, nodes %ld), eta = %.6f\n", tk ? "D5*" : "D5", it, d, nsvx, mnodes, exp(le));
        }
        return 0;
    }
    FILE *f = strcmp(argv[1], "-") ? fopen(argv[1], "r") : stdin;
    if (!f) { fprintf(stderr, "cannot open %s\n", argv[1]); return 1; }
    long first = argc > 2 ? atol(argv[2]) : 0, count = argc > 3 ? atol(argv[3]) : 1L << 40;
    char line[4096]; long ln = -1;
    printf("# line maxw_M vol_M result cells maxdepth flag lp short ball balltry secs\n");
    while (fgets(line, sizeof(line), f)) {
        ln++;
        if (ln < first) continue;
        if (ln >= first + count) break;
        Box *b = calloc(1, sizeof(Box));
        char *q = line; double wm = 0, vol = 1;
        for (int t = 0; t < 9; t++) {
            b->lo[5 + t] = strtod(q, &q); b->hi[5 + t] = strtod(q, &q);
            wm = fmax(wm, b->hi[5 + t] - b->lo[5 + t]); vol *= b->hi[5 + t] - b->lo[5 + t];
        }
        b->lo[0] = R1LO; b->hi[0] = R1HI;
        b->lo[1] = 0; b->hi[1] = 0.5;
        for (int i = 2; i <= 4; i++) { b->lo[i] = -0.5; b->hi[i] = 0.5; }
        nboxes = 0; maxdepth = 0; memset(splits_at, 0, sizeof(splits_at)); timed_out = 0; st_cert = st_lpc = st_short = st_ball = st_balltry = 0;
        tstart = clock();
        if (getenv("LOGDIR")) {     /* the certificate of this root: see check_lift.c for the format */
            char fn[1024]; snprintf(fn, sizeof(fn), "%s/%ld.log", getenv("LOGDIR"), ln);
            clog = fopen(fn, "w");
            if (!clog) { fprintf(stderr, "cannot open %s\n", fn); return 1; }
            fprintf(clog, "ROOT2 ver=%s\n", GITREV);
            nftab = 0;
            for (int t = 0; t < NZ; t++) fprintf(clog, "%a %a\n", b->lo[t], b->hi[t]);
        }
        int ok = run(b, 0);
        if (clog) { fprintf(clog, "END %s\n", ok ? "complete" : "incomplete"); fclose(clog); clog = NULL; }
        printf("%ld %.4f %.3e %s %ld %ld %ld %ld %ld %ld %ld %.1f\n", ln, wm, vol,
               ok ? "OK" : (timed_out ? "TIME" : "CAP"), nboxes, maxdepth, st_cert - st_lpc, st_lpc, st_short, st_ball,
               st_balltry, (double)(clock() - tstart) / CLOCKS_PER_SEC);
        if (verbose) { fprintf(stderr, "splits r phi1-4 y1-3 mu:"); for (int t = 0; t < NZ; t++) fprintf(stderr, " %ld", splits_at[t]); fprintf(stderr, "\n"); }
        fflush(stdout);
        free(b);
    }
    return 0;
}
