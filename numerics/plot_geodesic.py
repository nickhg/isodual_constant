"""Figure: Phi_5 along the isodual geodesic from D4 through M_* to A2 (x) A2 (geodesic_phi5.txt), against
the distance s = t ln 2 from D4 in the metric d of Theorem D5."""
import numpy as np, matplotlib
matplotlib.use("Agg"); import matplotlib.pyplot as plt
d = np.loadtxt("geodesic_phi5.txt"); t, lo = d[:, 0], d[:, 1]; s = t * np.log(2)
fig, ax = plt.subplots(figsize=(6.0, 3.2))
ax.plot(s, lo, color="#1f4e8c", lw=1.6)
ax.axhline(6, color="0.55", lw=0.8, ls="--"); ax.text(s[-1], 6.02, r"$6$", ha="right", va="bottom", fontsize=9, color="0.4")
pts = [(0.0, 3 * (36 * np.sqrt(2) / 35) ** 2, r"$D_4$", "left"), (np.log(1.2), 2 ** 16 / (3 ** 4 * 5 ** 3), r"$M_*$", "center"),
       (np.log(2), 16 / 3, r"$A_2\otimes A_2$", "left")]
for x, y, lab, ha in pts:
    ax.plot([x], [y], "o", color="#b2182b", ms=4.5, zorder=3)
    ax.annotate(lab, (x, y), textcoords="offset points", xytext=(6 if ha == "left" else 0, 7), ha=ha, fontsize=10)
ax.set_xlabel(r"distance from $D_4$ along the geodesic")
ax.set_ylabel(r"$\Phi_5$")
ax.set_xlim(s[0], s[-1]); ax.set_ylim(4.0, 6.75)
for sp in ("top", "right"): ax.spines[sp].set_visible(False)
fig.tight_layout(); fig.savefig("fig_phi5_geodesic.pdf")
print("max on path: %.6f at s = %.4f (ln(6/5) = %.4f)" % (lo.max(), s[lo.argmax()], np.log(1.2)))
