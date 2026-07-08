"""12. 能動循環の指紋 — 滞在 RNA の解離が「競合 RNA 濃度」で速くなり、やがて飽和する
       （Wagner, RNA Biology 2013 / Fender et al. 2010）。

出典: Wagner EGH (2013) "Cycling of RNAs on Hfq." RNA Biol 10(4):619-626.
（原著の生化学: Fender A et al. (2010) Genes Dev 24:2621; Hopkins et al. 2011 など）

11 は「能動交換で sRNA が速く乗れる」を扱った。12 はその**定量的な指紋**を見る:

  passive なら、滞在 RNA の解離速度は固有値（競合濃度に無関係で一定・とても遅い）。
  active cycling なら、**滞在 RNA の見かけ解離速度は競合 RNA 濃度とともに上がり、
  やがて飽和する**（サブユニット上の段階的置換の律速に達するため）。Fig.2B。

これを、競合が一時的に同居する中間体を経る置換で再現する:
    RH -> R + H            passive な自発解離（k_i、非常に遅い）
    RH + C <=> RHC         競合 C が Hfq-RH 複合体に取り付く（k_ass / k_dC）
    RHC -> CH + R          段階的置換で滞在 R が追い出される（k_rear、一次の律速）

「chase」実験: 滞在 RNA R を Hfq に付けておき、t=0 で競合 C を各濃度で加え、
Hfq 上に残る R（=RH+RHC）が半分になるまでの時間（半減期）を測る。
"""

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

from ecell4 import run_simulation
from ecell4_base.core import NetworkModel, ReactionRule, Species

RH, RHC, CH, Rf, C, H = (Species(x) for x in ["RH", "RHC", "CH", "R", "C", "H"])
K_I, K_ASS, K_DC, K_REAR = 0.005, 0.02, 1.0, 0.5   # /min
SP = ["RH", "RHC", "CH", "R", "C", "H"]


def Rr(a, b, k):
    return ReactionRule(a, b, k)


def build():
    mdl = NetworkModel()
    for x in [Rr([RH], [Rf, H], K_I),
              Rr([RH, C], [RHC], K_ASS), Rr([RHC], [RH, C], K_DC),
              Rr([RHC], [CH, Rf], K_REAR)]:
        mdl.add_reaction_rule(x)
    return mdl


def chase(c0, t_end=300.0, ndiv=6000):
    ret = run_simulation(t_end, y0={"RH": 1.0, "C": float(c0)}, model=build(),
                         solver="ode", ndiv=ndiv, species_list=SP)
    a = ret.as_array()
    return a[:, 0], a[:, 1] + a[:, 2]     # time, resident-on-Hfq (RH+RHC)


def halflife(c0):
    t, res = chase(c0)
    idx = np.where(res <= 0.5)[0]
    if len(idx) == 0:
        return np.nan
    i = idx[0]
    if i == 0:
        return t[0]
    return t[i - 1] + (0.5 - res[i - 1]) * (t[i] - t[i - 1]) / (res[i] - res[i - 1])


def main():
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(12, 4.6))

    # (左) chase: 競合が多いほど滞在 R は速く追い出される
    for c0, c in [(0, "C0"), (10, "C2"), (100, "C1"), (1000, "C3")]:
        t, res = chase(c0, t_end=30.0, ndiv=600)
        ax1.plot(t, res, "-", color=c, lw=2, label=f"competitor = {c0}")
    ax1.axhline(0.5, ls=":", c="gray", lw=1)
    ax1.set_xlabel("time after adding competitor (min)")
    ax1.set_ylabel("resident RNA still on Hfq")
    ax1.set_title("chase: more competitor → faster release\n(passive = flat)")
    ax1.legend(fontsize=8.5)

    # (右) 見かけ解離速度 vs 競合濃度: 上がって飽和（能動循環の指紋, Fig.2B）
    c_grid = np.array([0, 1, 3, 10, 30, 100, 300, 1000, 3000], dtype=float)
    kapp = np.array([np.log(2) / halflife(c) for c in c_grid])
    ax2.plot(np.log10(c_grid + 1), kapp, "o-", color="C4", lw=2, label="active cycling")
    ax2.axhline(K_I, ls="--", c="gray", lw=1.2, label=f"passive (intrinsic = {K_I}/min)")
    ax2.axhline(K_REAR, ls=":", c="C3", lw=1, label=f"rearrangement limit ({K_REAR}/min)")
    ax2.set_xlabel("competitor concentration  log10(C+1)")
    ax2.set_ylabel("apparent dissociation rate  (/min)")
    ax2.set_title("hallmark: apparent k_diss rises with\ncompetitor, then saturates")
    ax2.legend(fontsize=8.5)

    fig.tight_layout()
    out = "outputs/12_active_cycling.png"
    fig.savefig(out, dpi=120)
    print("saved:", out)
    for c0 in [0, 1, 10, 100, 1000]:
        print(f"competitor={c0:5g}  resident half-life={halflife(c0):7.2f} min")


if __name__ == "__main__":
    main()
