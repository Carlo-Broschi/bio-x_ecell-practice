"""06. Hfq を「限られた共有資源」として陽に入れる（Adamson & Lim, PLoS Comput Biol 2011）。

出典: Adamson DN, Lim HN (2011) "Essential Requirements for Robust Signaling in
Hfq Dependent Small RNA Networks." PLoS Comput Biol 7(8):e1002138.

02〜05 では sRNA と mRNA が直接出会って壊れ合うとしていた。だが実際は、多くの sRNA は
RNA シャペロン Hfq に一度乗ってから標的 mRNA と対合する。Hfq は数が限られた共有資源。
この論文は「単一の sRNA-mRNA ペアが Hfq 経由で効率よく・頑健に duplex を作る条件」を問う。

反応網（単一ペア、独立結合＋RNA解離。Fig.2A/3D）:
    ∅ -> s, ∅ -> m                  (sRNA/mRNA 転写: a_s, a_m)
    s + H <=> sH,  m + H <=> mH      (Hfq への結合: 会合 ka / 解離 kd)
    sH + m <=> T,  mH + s <=> T      (三者複合体 T への2経路)
    T -> D + H                       (アニールして duplex 形成, Hfq 放出: k5)
    s,m,sH,mH,T は分解(β)し、複合体は分解時に Hfq を放出、D も分解
    (Hfq は合成も分解もされず、総量 H_tot = H+sH+mH+T が保存される)

見どころ:
  - **Hfq 濃度に対して duplex 形成は釣鐘状**。少なすぎ=不足、多すぎ=sRNA と mRNA が
    別々の singly-bound 複合体(sH, mH)に隔離されて三者複合体ができない。中間が最適。
  - **頻繁な RNA 解離(kd 大)は、この窓を高Hfq側へ広げる**（＝Hfqロバスト性を上げる）。
    論文の主要メッセージのひとつ。
"""

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

from ecell4 import run_simulation
from ecell4_base.core import NetworkModel, ReactionRule, Species

H, S, M, SH, MH, T, D = (Species(x) for x in ["H", "s", "m", "sH", "mH", "T", "D"])
A_S = A_M = 10.0   # 転写速度
BETA = 1.0         # 分解・希釈(全種共通)
KA = 1.0           # 会合速度(全結合共通 = 独立結合)
K5 = 10.0          # 三者複合体 -> duplex
SPECIES = ["H", "s", "m", "sH", "mH", "T", "D"]


def build(h_tot, kd):
    mdl = NetworkModel()
    rules = [
        ReactionRule([], [S], A_S), ReactionRule([], [M], A_M),
        ReactionRule([S], [], BETA), ReactionRule([M], [], BETA),
        ReactionRule([S, H], [SH], KA), ReactionRule([SH], [S, H], kd),
        ReactionRule([M, H], [MH], KA), ReactionRule([MH], [M, H], kd),
        ReactionRule([SH, M], [T], KA), ReactionRule([T], [SH, M], kd),
        ReactionRule([MH, S], [T], KA), ReactionRule([T], [MH, S], kd),
        ReactionRule([T], [D, H], K5),
        ReactionRule([SH], [H], BETA), ReactionRule([MH], [H], BETA),
        ReactionRule([T], [H], BETA), ReactionRule([D], [], BETA),
    ]
    for r in rules:
        mdl.add_reaction_rule(r)
    return mdl


def steady(h_tot, kd):
    ret = run_simulation(400.0, y0={"H": h_tot}, model=build(h_tot, kd),
                         solver="ode", ndiv=1, species_list=SPECIES)
    return dict(zip(SPECIES, ret.as_array()[-1][1:]))


def pct_duplex(h_tot, kd):
    v = steady(h_tot, kd)
    tot_m = v["m"] + v["mH"] + v["T"] + v["D"]   # 標的 mRNA の総量
    return 100 * v["D"] / tot_m if tot_m > 0 else 0.0


def main():
    rel = np.logspace(-3, 4, 43)          # relative Hfq = H_tot / (a_m/β)
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(12, 4.6))

    # (左) 釣鐘 + 解離が窓を広げる
    for kd, c, lab in [(0.3, "C0", "rare dissociation  (kd=0.3)"),
                       (10.0, "C3", "frequent dissociation  (kd=10)")]:
        y = np.array([pct_duplex(r * 10, kd) for r in rel])
        ax1.plot(np.log10(rel), y, "-", color=c, lw=2, label=lab)
    ax1.axvspan(-3, -1.3, color="gray", alpha=.07)
    ax1.axvspan(1.3, 4, color="gray", alpha=.07)
    ax1.text(-2.2, 4, "insufficient\nHfq", ha="center", fontsize=8, color="gray")
    ax1.text(2.6, 4, "Hfq\nsequestration", ha="center", fontsize=8, color="gray")
    ax1.set_xlabel("relative Hfq   log10( H_tot / target mRNA )")
    ax1.set_ylabel("% target mRNA in duplex")
    ax1.set_title("duplex formation is bell-shaped in Hfq")
    ax1.legend(loc="upper right", fontsize=8)

    # (右) なぜ高Hfqで落ちるか: 標的 mRNA の内訳
    labels = ["low\n(relHfq=0.1)", "optimal\n(relHfq=1)", "high\n(relHfq=100)"]
    parts = ["m", "mH", "T", "D"]
    part_lab = ["free m", "m·Hfq (singly)", "ternary", "duplex"]
    colors = ["#b8c4be", "#e0a24b", "#8c78e6", "#0e9e6e"]
    data = []
    for r in [0.1, 1.0, 100.0]:
        v = steady(r * 10, 0.3)
        tot = v["m"] + v["mH"] + v["T"] + v["D"]
        data.append([v[p] / tot for p in parts])
    data = np.array(data)
    bottom = np.zeros(3)
    x = np.arange(3)
    for j, (pl, col) in enumerate(zip(part_lab, colors)):
        ax2.bar(x, data[:, j], bottom=bottom, color=col, label=pl, width=.6)
        bottom += data[:, j]
    ax2.set_xticks(x); ax2.set_xticklabels(labels, fontsize=9)
    ax2.set_ylabel("fraction of target mRNA")
    ax2.set_title("why high Hfq fails: mRNA locked in singly-bound m·Hfq")
    ax2.legend(fontsize=8, loc="center left", bbox_to_anchor=(1.0, 0.5))

    fig.tight_layout()
    out = "outputs/06_hfq_robustness.png"
    fig.savefig(out, dpi=120, bbox_inches="tight")
    print("saved:", out)
    for r in [0.01, 0.1, 1, 3, 10, 100, 1000]:
        print(f"relHfq={r:7g}  %duplex(kd0.3)={pct_duplex(r*10,0.3):5.1f}  "
              f"%duplex(kd10)={pct_duplex(r*10,10.0):5.1f}")


if __name__ == "__main__":
    main()
