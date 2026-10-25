#!/bin/bash
# Tier 1 of VERIFY.md: every exact (non-search) certificate in the paper, with timings; outputs to expected/.
cd "$(dirname "$0")"; mkdir -p expected
run() { local name=$1; shift; local t0=$(date +%s); "$@" > expected/$name.txt 2>&1; echo "exit $? after $(( $(date +%s) - t0 )) s" >> expected/$name.txt; tail -n 3 expected/$name.txt | sed "s/^/[$name] /"; }
run d5_lemma          sage d5_lemma.sage 1/10 5            # Theorem D5 (radius 1/10; dual by isoduality)
run pivd5_lemma       sage pivd5_lemma.sage 1/20 11/20     # Proposition pivd5 (radius 1/20)
run mstar_lemma       sage mstar_lemma.sage 1/40 9/10      # Proposition mstar (radius 1/40)
run mstar_unique      sage mstar_unique.sage               # Prop. mstarprops (1): kappa(M_*) only at +-phi_0
run lstar_check       sage lstar_check.sage                # eta(L_*) = 486/125
run mstar_invariants  sage mstar_invariants.sage           # |Aut(M_*)|, theta series, modularity, isodual family
run a2a2              sage a2a2.sage                       # Phi_5(A2 x A2) = 16/3, kappa = 3/2
run shadow_iso        sage shadow_iso.sage                 # pi_v A4 = (Z^3, I + J/2), pi_u(pi_v E6) = A2 x A2
run family_exact      python3 family_exact.py              # D4 and A2 x A2 in M_*'s isodual family (exact)
run geodesic_exact    python3 geodesic_exact.py            # M_* on the D4 - A2 x A2 geodesic, parameter log_2(6/5)
