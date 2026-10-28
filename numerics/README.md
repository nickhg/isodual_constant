# Numerical computations (not used in any proof)

These scripts produce the floating-point results that the paper labels as numerical. None of them is part of the
proof of any theorem; the certificates are in `certificates/` and the searches in the dataset. Run them from this
directory (they load `../certificates/psi.sage` for exact Ψ-minimization).

| script | paper | what it does |
|---|---|---|
| `geodesic_phi5.sage`, `plot_geodesic.py`, `geodesic_phi5.txt` | Figure 1, §6 | Φ₅ at 66 points of the isodual geodesic from D₄ through M\* to A₂⊗A₂: κ by branch and bound over the torus of lifting forms (`kappa_bb.py`, bracketing the global maximum to relative accuracy 10⁻⁸), η by exact Ψ-minimization of a rational Gram matrix within 10⁻¹³ of the floating-point one. `geodesic_phi5.txt` is the output used for the figure |
| `scan.sage`, `nearest_copies.py`, `nearest_local.py`, `copies_pivE6_D5s.txt` | §7 (Φ₆ landscape) | Φ₆ along the geodesic from π_vE₆ to the nearest copy of D₅\* |
| `e6_radius.py`, `e6_radius.txt` | remark "Eutaxy radii of E₆ and E₇" | 300 polar-vertex ascents for the first-order radius r(E₆) over the C_f of the 19,440 Ψ-minimal flags (needs `../certificates/e6_A.sobj`, written by `e6_flags.sage`); `e6_radius.txt` is the output |
| `e7_dump.sage`, `e7_radius.py`, `e7_radius.txt` | remark "Eutaxy radii of E₆ and E₇" | the same for E₇: `e7_dump.sage` writes the C_f of all 544,320 Ψ-minimal flags to `e7_flags.npz` (118 MB, not in the repository; about 30 s), and `e7_radius.py 300 10 7` runs 300 ascents (about 4 hours); `e7_radius.txt` is the output |
| `mu_check.sage` | §9 | for m = 6, 7, 8, hill climbs over perturbations of A_m, D_m, E_m and their duals, looking for a Ψ-minimal basis with r₁ > γ_m Δ^{1/m} |
| `kappa_bb.py` | — | the floating-point branch and bound for κ used by the scripts above |

Usage, for example:
```
sage geodesic_phi5.sage > geodesic_phi5.txt && python3 plot_geodesic.py
sage scan.sage copies_pivE6_D5s.txt 0 41
sage mu_check.sage 6 200 20 1
(cd ../certificates && sage e6_flags.sage) && sage -python e6_radius.py 300 > e6_radius.txt
sage e7_dump.sage && python3 e7_radius.py 300 10 7 > e7_radius.txt
```
