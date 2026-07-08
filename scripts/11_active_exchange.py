"""11. 強固に結合するのに分単位で応答できる謎 — Hfq 上の能動的 RNA 交換
       （Park, Prevost, Heideman et al., eLife 2021）。

出典: Park S, Prevost K, Heideman EM, et al. (2021) "Dynamic interactions between
the RNA chaperone Hfq, small regulatory RNAs, and mRNAs in live bacterial cells."
eLife 10:e64207.

単一分子ライブセル観察で分かったこと:
  - 通常増殖時、**Hfq はほとんど mRNA に占有**されている（distal 面）。Hfq は量が律速。
  - in vitro では sRNA-Hfq / mRNA-Hfq の寿命は >100 分（外れない）。なのにストレス時、
    sRNA は数分で標的を制御しはじめる。**どうやって占有済みの Hfq に素早く乗るのか？**
  - 答え: 空くのを待つ(passive)のではなく、**入ってきた sRNA が結合済み mRNA を直接追い出す
    能動交換(active exchange)**。これで分単位の応答が可能になる。

モデル（純粋な結合速度論、生成/分解なし）:
  H(遊離Hfq), m(遊離mRNA), s(遊離sRNA), mH(mRNA·Hfq), sH(sRNA·Hfq)
    m + H <=> mH,  s + H <=> sH     （会合 k_on / 解離 k_off。k_off=0.01/min = ~70分寿命）
    active のみ:  mH + s <=> sH + m  （能動交換 k_ex: sRNA が mRNA を置換）

事前状態: Hfq は全て mRNA に占有(mH=10)、余剰 mRNA=90。t=0 で sRNA=20 を誘導し、
機能型の sRNA·Hfq(sH) が立ち上がる速さを passive と active で比べる。
"""

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

from ecell4 import run_simulation
from ecell4_base.core import NetworkModel, ReactionRule, Species

H, M, S, MH, SH = (Species(x) for x in ["H", "m", "s", "mH", "sH"])
K_ON, K_OFF = 0.1, 0.01
SP = ["H", "m", "s", "mH", "sH"]


def R(a, b, k):
    return ReactionRule(a, b, k)


def build(k_ex):
    rules = [R([M, H], [MH], K_ON), R([MH], [M, H], K_OFF),
             R([S, H], [SH], K_ON), R([SH], [S, H], K_OFF)]
    if k_ex > 0:
        rules += [R([MH, S], [SH, M], k_ex), R([SH, M], [MH, S], k_ex)]
    mdl = NetworkModel()
    for r in rules:
        mdl.add_reaction_rule(r)
    return mdl


def traj(k_ex, t_end=60.0, ndiv=120):
    ret = run_simulation(t_end, y0={"mH": 10.0, "m": 90.0, "s": 20.0},
                         model=build(k_ex), solver="ode", ndiv=ndiv, species_list=SP)
    return ret.as_array()


def main():
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(12, 4.6))

    # (左) 機能型 sRNA·Hfq の立ち上がり: passive(遅い) vs active(速い)
    for k_ex, c, lab in [(0.0, "C0", "passive (wait for mRNA to fall off)"),
                         (0.05, "C3", "active exchange (sRNA displaces mRNA)")]:
        a = traj(k_ex)
        ax1.plot(a[:, 0], a[:, 5], "-", color=c, lw=2, label=lab)   # sH = col 5
    ax1.axvspan(0, 5, color="gray", alpha=.08)
    ax1.text(2.5, 0.15, "stress\nresponse\nwindow", ha="center", fontsize=7.5, color="gray")
    ax1.set_xlabel("time after sRNA induction (min)")
    ax1.set_ylabel("functional sRNA-Hfq  (sH)")
    ax1.set_title("active exchange loads sRNA in minutes\ndespite >100 min binding")
    ax1.legend(fontsize=8.5)

    # (右) 定常での Hfq 占有: mRNA が多数派、でも sRNA も一定分をすばやく確保
    for i, (k_ex, lab) in enumerate([(0.0, "passive\n(t=60)"), (0.05, "active\n(t=60)")]):
        v = dict(zip(SP, traj(k_ex)[-1][1:]))
        bottom = 0
        for key, col, klab in [("H", "#b8c4be", "free Hfq"),
                               ("mH", "#e0a24b", "mRNA-bound"),
                               ("sH", "#0e9e6e", "sRNA-bound")]:
            ax2.bar(i, v[key], bottom=bottom, color=col, width=.55,
                    label=klab if i == 0 else None)
            bottom += v[key]
    ax2.set_xticks([0, 1]); ax2.set_xticklabels(["passive", "active"])
    ax2.set_ylabel("Hfq (of 10 total)")
    ax2.set_title("Hfq stays mostly mRNA-occupied;\nactive lets sRNA grab a share fast")
    ax2.legend(fontsize=8.5, loc="upper left")

    fig.tight_layout()
    out = "outputs/11_active_exchange.png"
    fig.savefig(out, dpi=120)
    print("saved:", out)
    for k_ex, lab in [(0.0, "passive"), (0.05, "active")]:
        a = traj(k_ex)
        print(f"{lab:8s}: sH(5min)={a[10][5]:5.2f}  sH(15min)={a[30][5]:5.2f}  sH(60min)={a[-1][5]:5.2f}")


if __name__ == "__main__":
    main()
