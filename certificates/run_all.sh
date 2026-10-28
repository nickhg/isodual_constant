#!/bin/bash
# Tier 1 of VERIFY.md: every exact (non-search) certificate in the paper, with timings; outputs to expected/.
cd "$(dirname "$0")"; mkdir -p expected
run() { local name=$1; shift; local t0=$(date +%s); "$@" > expected/$name.txt 2>&1; echo "exit $? after $(( $(date +%s) - t0 )) s" >> expected/$name.txt; tail -n 3 expected/$name.txt | sed "s/^/[$name] /"; }
run d5_lemma          sage d5_lemma.sage 1/10 5            # Theorem D5 (radius 1/10; dual by isoduality)
run radius_d4d5       sage radius_d4d5.sage                # App. A: D4 radius 2/3 (radius 33/100), D5 radius <= 4/(3 sqrt5)
run pivd5_lemma       sage pivd5_lemma.sage 1/20 11/20     # Proposition pivd5 (radius 1/20)
run mstar_lemma       sage mstar_lemma.sage 1/40 9/10      # Proposition mstar (radius 1/40)
run mstar_unique      sage mstar_unique.sage               # Prop. mstarprops (1): kappa(M_*) only at +-phi_0
run lstar_check       sage lstar_check.sage                # eta(L_*) = 486/125
run mstar_invariants  sage mstar_invariants.sage           # |Aut(M_*)|, theta series, modularity, isodual family
run a2a2              sage a2a2.sage                       # Phi_5(A2 x A2) = 16/3, kappa = 3/2
run shadow_iso        sage shadow_iso.sage                 # pi_v A4 = (Z^3, I + J/2), pi_u(pi_v E6) = A2 x A2
run family_exact      python3 family_exact.py              # D4 and A2 x A2 in M_*'s isodual family (exact)
run geodesic_exact    python3 geodesic_exact.py            # M_* on the D4 - A2 x A2 geodesic, parameter log_2(6/5)
run e6_flags          sage e6_flags.sage                   # Prop. e6e7: the 19,440 Psi-minimal flags of E6, design identities (writes e6_A.sobj)
run e6_cross          sage e6_cross.sage                   # Prop. e6e7: exact LP certificate r(E6) >= 3/10 (radius 1/10)
run e7_design         sage e7_design.sage                  # Prop. e6e7: E7 design identities over 544,320 flags; upper bound along omega_7
run schur_check_6     sage schur_check.sage 6              # Prop. e6e7: W(E6) irreducible on traceless symmetric forms
run schur_check_7     sage schur_check.sage 7              # Prop. e6e7: W(E7) irreducible
run schur_check_8     sage schur_check.sage 8              # Prop. e6e7: W(E8) irreducible; interval check of the radius 1/20
run d4_kappa          python3 d4_kappa_classes.py          # Prop. kapparoot: exact search, 48 classes at D4, 24 at A4
run mstar_genus       sage mstar_genus.sage                # Prop. mstarprops (7): genus of P_0, class number 6
run mstar_mass        sage mstar_mass.sage                 # Prop. mstarprops (7): mass check, minima of the six classes
run mstar_eisenstein  sage mstar_eisenstein.sage           # Prop. mstarprops (8): the four g with g^2 + g + 1 = 0
run mstar_herm_check  sage mstar_herm_check.sage           # Prop. mstarprops (8): h is Hermitian with trace P_0
run mstar_herm_genus  sage mstar_herm_genus.sage           # Section 6: Hermitian Gram matrices of the four structured classes
