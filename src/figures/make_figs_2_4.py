# -*- coding: utf-8 -*-
"""Figures 2-4 for the manuscript (Nature Methods, double column 183 mm).
Output: fig2_scaling.svg, fig3_planted.svg, fig4_baselines.svg (-> EMF via soffice).
All text Arial 7pt, grayscale-safe.
Data sources: tables/*.csv (E-C01/E-D01/E-D02/E-F05/E-B09/E-A01-GPU runs).
"""
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

plt.rcParams.update({
    "font.family": "Arial",
    "font.size": 7,
    "svg.fonttype": "none",
    "axes.linewidth": 0.6,
    "xtick.major.width": 0.6,
    "ytick.major.width": 0.6,
    "xtick.major.size": 2.5,
    "ytick.major.size": 2.5,
    "legend.frameon": False,
})
MM = 1 / 25.4
C = {"cpu": "#000000", "gpu": "#666666", "aux": "#999999"}
LS = {"cooc": "-", "cond": "--", "mutex": ":"}

# ---------------- Fig 2: scaling ----------------
N = [1992, 4000, 8000]
cpu = [28.94, 46.78, 90.55]
gpu = [0.677, 0.693, 0.696]
cpu_rate = [51.5, 31.9, 16.5]
gpu_rate = [1.49e9 / 0.677 / 1e6, 1.49e9 / 0.693 / 1e6, 1.49e9 / 0.696 / 1e6]

fig, (a1, a2) = plt.subplots(1, 2, figsize=(183 * MM, 62 * MM))
a1.plot(N, cpu, "o-", color=C["cpu"], lw=1.0, ms=3, label="CPU (AVX2, 1 core)")
a1.plot(N, gpu, "s--", color=C["gpu"], lw=1.0, ms=3, label="GPU (BMMA wave kernel)")
a1.set_yscale("log")
a1.set_xlabel("Genomes N (M = 54,590 families)")
a1.set_ylabel("Wall time (s), log scale")
a1.set_xticks(N)
a1.legend(loc="center right", fontsize=6)
for x, y, s in zip(N, cpu, gpu):
    a1.annotate(f"{y:.2f} s", (x, y), textcoords="offset points", xytext=(-2, -11), fontsize=5.5, color=C["gpu"], ha="right")
for x, y in zip(N, cpu):
    a1.annotate(f"{y:.1f} s", (x, y), textcoords="offset points", xytext=(-2, 5), fontsize=5.5, color=C["cpu"], ha="right")
a1.text(0.02, 0.98, "a", transform=a1.transAxes, fontsize=9, fontweight="bold", va="top")

a2.plot(N, cpu_rate, "o-", color=C["cpu"], lw=1.0, ms=3, label="CPU throughput")
a2.plot(N, gpu_rate, "s--", color=C["gpu"], lw=1.0, ms=3, label="GPU throughput")
a2.set_yscale("log")
a2.set_xlabel("Genomes N")
a2.set_ylabel("Pairs / s (Mpair s$^{-1}$), log scale")
a2.set_xticks(N)
a2.set_xlim(1500, 8600)
a2.set_ylim(1, 1e4)
a2.legend(loc="upper right", fontsize=6)
a2.text(0.02, 0.98, "b", transform=a2.transAxes, fontsize=9, fontweight="bold", va="top")
fig.tight_layout(w_pad=2.0)
fig.savefig("fig2_scaling.svg")
plt.close(fig)

# ---------------- Fig 3: planted recall ----------------
bins = ["<=5", "6-10", "11-20", "21-40", "41-80", "81-160", ">160"]
rec = {
    "cooc":  [0.479, 0.500, 0.441, 0.475, 0.255, 0.559, 0.611],
    "cond":  [0.0, 0.0, 0.156, 0.610, 0.353, 0.609, 0.615],
    "mutex": [0.0, 0.0, 0.0, 0.0, 0.0, 0.157, 0.462],
}
fig, (a1, a2) = plt.subplots(1, 2, figsize=(183 * MM, 60 * MM))
x = np.arange(len(bins))
w = 0.26
for k, (typ, vals) in enumerate(rec.items()):
    a1.bar(x + (k - 1) * w, vals, width=w, color=["#000000", "#666666", "#aaaaaa"][k],
           edgecolor="none", label=typ)
a1.set_xticks(x)
a1.set_xticklabels(bins, fontsize=6)
a1.set_xlabel("Planted support bin (genomes)")
a1.set_ylabel("Recall at BH-FDR 0.05")
a1.legend(loc="center right", fontsize=6)
a1.text(0.02, 0.98, "a", transform=a1.transAxes, fontsize=9, fontweight="bold", va="top")

fdr = ["0.1", "0.05", "0.01", "0.001"]
overall = [0.288, 0.287, 0.279, 0.272]
a2.plot(fdr, overall, "o-", color=C["cpu"], lw=1.0, ms=3, label="FDR level (base settings)")
a2.plot([0, 1, 2], [0.287, 0.209, 0.184], "s--", color=C["gpu"], lw=1.0, ms=3,
        label="Support floors (FDR 0.05)")
a2.set_xticks([0, 1, 2])
a2.set_xticklabels(["base", "c11>=10", "minsup>=20"], fontsize=6)
a2.set_ylabel("Overall recall (n = 1,000)")
a2.set_ylim(0, 0.35)
a2.legend(loc="center right", fontsize=6)
a2.text(0.02, 0.98, "b", transform=a2.transAxes, fontsize=9, fontweight="bold", va="top")
fig.tight_layout(w_pad=2.0)
fig.savefig("fig3_planted.svg")
plt.close(fig)

# ---------------- Fig 4: approximate baselines ----------------
fig, (a1, a2) = plt.subplots(1, 2, figsize=(183 * MM, 62 * MM))
cats = ["co-occ.\npairs", "exclusion\npairs", "all pairs"]
exact = [100, 100, 100]
lsh = [15.1, 0.69, 12.75]
x = np.arange(3)
a1.bar(x - 0.18, exact, width=0.36, color="#000000", label="Exact pipeline")
a1.bar(x + 0.18, lsh, width=0.36, color="#aaaaaa", label="MinHash-LSH")
a1.set_xticks(x)
a1.set_xticklabels(cats, fontsize=6)
a1.set_ylabel("Recall of z-screened pairs (%)")
a1.set_ylim(0, 122)
a1.legend(loc="upper right", fontsize=6)
for xi, v in zip(x, lsh):
    a1.annotate(f"{v:g}%", (xi + 0.18, v), textcoords="offset points", xytext=(0, 2),
                fontsize=5.5, ha="center")
a1.text(0.02, 0.98, "a", transform=a1.transAxes, fontsize=9, fontweight="bold", va="top")

frac = [16.7, 9.8, 5.7, 2.9, 1.9]
surv = [38.5, 25.0, 17.8, 11.6, 8.1]
a2.plot(frac, surv, "o-", color=C["cpu"], lw=1.0, ms=3)
for xi, yi in zip(frac, surv):
    a2.annotate(f"{yi:g}%", (xi, yi), textcoords="offset points", xytext=(0, 4),
                fontsize=5.5, ha="center")
a2.set_xlabel("Pairs tested (% of all pairs)")
a2.set_ylabel("z-screened signal retained (%)")
a2.set_xlim(0, 19)
a2.set_ylim(0, 45)
a2.text(0.02, 0.98, "b", transform=a2.transAxes, fontsize=9, fontweight="bold", va="top")
fig.tight_layout(w_pad=2.0)
fig.subplots_adjust(bottom=0.22, top=0.94)
fig.savefig("fig4_baselines.svg")
plt.close(fig)
print("SVGs written")
