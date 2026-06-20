#!/usr/bin/env python
"""Compare fastQpick read-level bootstrap SEs against the analytic multinomial
SE on the real yeast RNA-seq sample, and render the Application figure."""
import os
import glob

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

HERE = os.path.dirname(os.path.abspath(__file__))
FIG_DIR = os.path.join(os.path.dirname(HERE), "figures")
os.makedirs(FIG_DIR, exist_ok=True)

targets = np.loadtxt(os.path.join(HERE, "targets.txt"), dtype=str)
baseline = np.loadtxt(os.path.join(HERE, "counts", "baseline.txt"))

rep_files = sorted(glob.glob(os.path.join(HERE, "counts", "rep_*.txt")))
R = np.vstack([np.loadtxt(f) for f in rep_files])   # (B, T) est_counts
B = R.shape[0]
print(f"Loaded {B} bootstrap replicates over {R.shape[1]} transcripts")

# Per-replicate proportions among pseudoaligned reads.
P = R / R.sum(axis=1, keepdims=True)               # (B, T)
N = baseline.sum()                                  # ~ n_pseudoaligned
p_hat = baseline / N                                # point-estimate proportions

boot_se = P.std(axis=0, ddof=1)
analytic_se = np.sqrt(p_hat * (1.0 - p_hat) / N)

# Restrict the comparison to expressed transcripts (>= 100 est_counts) so that
# the SE comparison is not dominated by near-zero-count Monte Carlo noise.
expressed = baseline >= 100
print(f"Expressed transcripts (>=100 counts): {expressed.sum()}")

rel = np.abs(boot_se[expressed] - analytic_se[expressed]) / analytic_se[expressed]
ratio = boot_se[expressed] / analytic_se[expressed]
print(f"median |boot-analytic|/analytic : {np.median(rel):.1%}")
print(f"mean   boot/analytic ratio      : {ratio.mean():.4f}")
print(f"Pearson r (log SE)              : "
      f"{np.corrcoef(np.log(analytic_se[expressed]), np.log(boot_se[expressed]))[0,1]:.5f}")

# ---- Detection robustness (Panel C) -------------------------------------
# A discrete presence/absence statistic for which the Gaussian standard error
# of Panels A/B is not meaningful: a transcript is either called or not. The
# bootstrap assigns each transcript a detection probability pi_t, the fraction
# of replicates whose estimated count reaches a detection threshold theta. This
# identifies which low-abundance calls persist and which are driven by a handful
# of reads and go away under resampling -- a robustness question with no
# closed-form answer.
THETA = 5.0
detect_prob = (R >= THETA).mean(axis=0)            # (T,) bootstrap detection prob.
persist = detect_prob >= 0.95
fragile = (detect_prob > 0.05) & (detect_prob < 0.95)
n_fragile = int(fragile.sum())
n_called_fragile = int((fragile & (baseline >= THETA)).sum())
n_absent_appear = int((fragile & (baseline < THETA)).sum())
print(f"detection threshold theta         : {THETA:g} est. counts")
print(f"fragile calls (0.05<pi<0.95)      : {n_fragile}")
print(f"  called at baseline, may vanish  : {n_called_fragile}")
print(f"  absent at baseline, may appear  : {n_absent_appear}")

# ---- Figure -------------------------------------------------------------
fig, (ax1, ax2, ax3) = plt.subplots(1, 3, figsize=(16.5, 4.3))

# Panel A: standardized bootstrap deviations pooled across all expressed
# transcripts. For each transcript the deviation of each replicate's estimated
# proportion from the point estimate is divided by that transcript's analytic
# SE; if the analytic SE captures the bootstrap spread, the pooled deviations
# follow the standard normal regardless of transcript abundance.
z = ((P[:, expressed] - p_hat[expressed]) / analytic_se[expressed]).ravel()
z_std = z.std(ddof=1)
ax1.hist(z, bins=80, density=True, color="#4C72B0", alpha=0.6,
         label=f"pooled $z$  (std={z_std:.2f})")
xs = np.linspace(-4, 4, 400)
ax1.plot(xs, np.exp(-0.5 * xs ** 2) / np.sqrt(2 * np.pi),
         color="#C44E52", lw=2, label=r"$\mathcal{N}(0,1)$")
ax1.set_xlim(-4, 4)
ax1.set_xlabel(r"standardized deviation  $(\hat{p}^{(b)}-\hat{p})\,/\,$SE$_{\mathrm{analytic}}$")
ax1.set_ylabel("density")
ax1.set_title("A. Standardized bootstrap deviations vs. $\\mathcal{N}(0,1)$")
ax1.legend(frameon=False, fontsize=9)
print(f"pooled standardized-deviation std : {z_std:.4f}")

# Panel B: bootstrap SE vs analytic SE across all expressed transcripts.
ase = analytic_se[expressed]
bse = boot_se[expressed]
lo = min(ase.min(), bse.min()) * 0.8
hi = max(ase.max(), bse.max()) * 1.2
ax2.plot([lo, hi], [lo, hi], color="#888888", ls="--", lw=1, label="y = x")
sc = ax2.scatter(ase, bse, c=np.log10(baseline[expressed]), cmap="viridis",
                 s=8, alpha=0.6)
ax2.set_xscale("log"); ax2.set_yscale("log")
ax2.set_xlim(lo, hi); ax2.set_ylim(lo, hi)
ax2.set_xlabel(r"analytic SE  $\sqrt{\hat{p}(1-\hat{p})/N}$")
ax2.set_ylabel("bootstrap SE")
ax2.set_title("B. Bootstrap SE recovers analytic SE")
cb = fig.colorbar(sc, ax=ax2)
cb.set_label(r"$\log_{10}$ est. counts", fontsize=8)
ax2.legend(frameon=False, fontsize=9, loc="upper left")

# Panel C: bootstrap detection probability vs. point estimate for low-abundance
# transcripts. Each transcript is colored by whether its call persists across
# replicates, is fragile, or is essentially never made. The point estimate alone
# (a single count, x-axis) places every transcript deterministically on one side
# of the threshold, but the bootstrap reveals a band of calls whose presence
# depends on the particular reads that were sequenced.
# Persist and fragile transcripts are plotted in full so the legend counts match
# the totals reported in the text; the rarely-called class is restricted to the
# low-count regime to avoid burying the panel under the large mass of
# never-expressed transcripts. Baseline counts of zero cannot be placed on the log
# axis, so points are floored to the left edge, which keeps visible the transcripts
# that are absent at baseline yet appear under resampling.
XFLOOR = 0.3
cat = np.where(persist, 0, np.where(fragile, 1, 2))
lowab = baseline > XFLOOR                           # rarely-called bulk exclusion
colors = np.array(["#55A868", "#DD8452", "#C44E52"])   # persist / fragile / lost
labels = ["persists ($\\pi\\geq0.95$)", "fragile ($0.05<\\pi<0.95$)",
          "rarely called ($\\pi\\leq0.05$)"]
for k in (0, 1, 2):                                 # persists, fragile, rarely called
    m = (cat == k) if k != 2 else ((cat == k) & lowab)
    xk = np.clip(baseline[m], XFLOOR, None)
    ax3.scatter(xk, detect_prob[m], s=10, alpha=0.6, color=colors[k],
                label=f"{labels[k]}  (n={int(m.sum())})")
ax3.set_xlim(left=XFLOOR * 0.9)
ax3.axvline(THETA, color="#888888", ls="--", lw=1)
ax3.text(THETA * 1.1, 0.06, f"threshold $\\theta={THETA:g}$", rotation=90,
         va="bottom", ha="left", fontsize=8, color="#888888")
ax3.set_xscale("log")
ax3.set_ylim(-0.03, 1.03)
ax3.set_xlabel("point-estimate count $\\hat{c}_t$ (baseline)")
ax3.set_ylabel(r"bootstrap detection probability $\pi_t$")
ax3.set_title("C. Detection robustness of low-abundance calls")
ax3.legend(frameon=False, fontsize=8, loc="center right")

fig.tight_layout()
out = os.path.join(FIG_DIR, "bootstrap_realdata.png")
fig.savefig(out, dpi=200)
print("Saved", out)
