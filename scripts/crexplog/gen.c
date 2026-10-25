#include <stdio.h>
#include <stdlib.h>
#include <fenv.h>
#include <math.h>
double cr_exp(double), cr_log(double);
int main(int argc, char **argv) {
    long n = atol(argv[1]); srand(7);
    for (long i = 0; i < n; i++) {
        double u = rand() / (RAND_MAX + 1.0), v = rand() / (RAND_MAX + 1.0);
        double xe = (u - 0.5) * 60 * (i % 3 == 0 ? v : 1e-3 * v + (i % 3 == 1 ? 1 : 0));   /* exp args incl. small */
        double xl = exp((v - 0.5) * 40) * (1 + 1e-3 * u);                                 /* log args, positive */
        fesetround(FE_UPWARD); double eu = cr_exp(xe), lu = cr_log(xl);
        fesetround(FE_TONEAREST); double en = cr_exp(xe), ln = cr_log(xl);
        printf("%a %a %a %a %a %a\n", xe, eu, en, xl, lu, ln);
    }
    return 0;
}
