# Numerical computations (not used in any proof)

These scripts produce the floating-point results that the paper labels as numerical. None of them is part of the
proof of any theorem; the certificates are in `certificates/` and the searches in the dataset. Run them from this
directory (they load `../certificates/psi.sage` for exact Ψ-minimization).

| script | paper | what it does |
|---|---|---|
| `geodesic_phi5.sage`, `plot_geodesic.py`, `geodesic_phi5.txt` | Figure 1, §6 | Φ₅ at 66 points of the isodual geodesic from D₄ through M\* to A₂⊗A₂: κ by branch and bound over the torus of lifting forms (`kappa_bb.py`, bracketing the global maximum to relative accuracy 10⁻⁸), η by exact Ψ-minimization of a rational Gram matrix within 10⁻¹³ of the floating-point one. `geodesic_phi5.txt` is the output used for the figure |
| `scan.sage`, `nearest_copies.py`, `nearest_local.py`, `copies_pivE6_D5s.txt` | §7 (Φ₆ landscape) | Φ₆ along the geodesic from π_vE₆ to the nearest copy of D₅\* |
| `mu_check.sage` | §9 | for m = 6, 7, 8, hill climbs over perturbations of A_m, D_m, E_m and their duals, looking for a Ψ-minimal basis with r₁ > γ_m Δ^{1/m} |
| `kappa_bb.py` | — | the floating-point branch and bound for κ used by the scripts above |

Usage, for example:
```
sage geodesic_phi5.sage > geodesic_phi5.txt && python3 plot_geodesic.py
sage scan.sage copies_pivE6_D5s.txt 0 41
sage mu_check.sage 6 200 20 1
```
