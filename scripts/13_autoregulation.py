"""13. Hfq は自分の翻訳を抑える — 負の自己制御による恒常性（Morita & Aiba, RNA 2019）。

出典: Morita T, Aiba H (2019) "Mechanism and physiological significance of
autoregulation of the Escherichia coli hfq gene." RNA 25(2):264-276.

Hfq タンパク質は自分の mRNA の 5'-UTR（distal 面）に結合して、自分の翻訳を止める（負の自己制御）。
論文の実験: hfq mRNA を ~3倍過剰発現しても、Hfq タンパク質は ~1.4倍しか増えず頭打ち。
生理的意義: 過剰 Hfq は rim を介して無関係な RNA を巻き込み毒性 → 一定範囲に保つ必要がある。

負フィードバックの最小モデル（Hfq が自分の mRNA を隔離＝翻訳沈黙化）:
    ∅ -> m            (hfq 転写 α。プラスミド/IPTG で可変)
    m -> ∅            (mRNA 分解 β_m)
    m -> m + p        (翻訳: 自由な mRNA からのみ)
    p -> ∅            (Hfq 希釈 β_p)
    m + p <=> mp      (Hfq が自分の mRNA UTR に結合 → 翻訳沈黙。auto のみ)
    mp -> p           (沈黙 mRNA も分解し、Hfq を放出)

負フィードバックの2大効果を見る:
  (1) 恒常性: 転写 α を上げても Hfq 量は緩衝される（no-auto は線形、auto は頭打ち）。
  (2) 応答が速い: 定常に到達する時間が短くなる（Rosenfeld/Alon の古典）。
"""

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

from ecell4 import run_simulation
from ecell4_base.core import NetworkModel, ReactionRule, Species

M, P, MP = (Species(x) for x in ["m", "p", "mp"])
B_M, K_TL, B_P, K_B, K_U = 0.3, 2.0, 0.05, 0.4, 0.1
SP = ["m", "p", "mp"]


def R(a, b, k):
    return ReactionRule(a, b, k)


def build(alpha, auto):
    rules = [R([], [M], alpha), R([M], [], B_M), R([M], [M, P], K_TL), R([P], [], B_P)]
    if auto:
        rules += [R([M, P], [MP], K_B), R([MP], [M, P], K_U), R([MP], [P], B_M)]
    mdl = NetworkModel()
    for r in rules:
        mdl.add_reaction_rule(r)
    return mdl


def steady(alpha, auto):
    v = dict(zip(SP, run_simulation(600.0, y0={}, model=build(alpha, auto),
                                    solver="ode", ndiv=1, species_list=SP).as_array()[-1][1:]))
    return v["m"] + v["mp"], v["p"] + v["mp"]     # total mRNA, total Hfq protein


def timecourse(alpha, auto, t_end=50.0, ndiv=2000):
    a = run_simulation(t_end, y0={}, model=build(alpha, auto), solver="ode",
                       ndiv=ndiv, species_list=SP).as_array()
    return a[:, 0], a[:, 2] + a[:, 3]             # time, total Hfq protein


def main():
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(12, 4.6))

    # (1) 恒常性: fold-change vs 転写 α
    a_grid = np.array([1, 2, 3, 4, 6, 8], dtype=float)
    m0a, p0a = steady(1, True)
    m0n, p0n = steady(1, False)
    m_fold = [steady(a, True)[0] / m0a for a in a_grid]
    p_auto = [steady(a, True)[1] / p0a for a in a_grid]
    p_noauto = [steady(a, False)[1] / p0n for a in a_grid]
    ax1.plot(a_grid, m_fold, "s--", color="gray", label="hfq mRNA")
    ax1.plot(a_grid, p_noauto, "o-", color="C0", label="Hfq protein, no autoreg")
    ax1.plot(a_grid, p_auto, "o-", color="C3", label="Hfq protein, autoregulated")
    ax1.set_xlabel("hfq transcription rate  (fold, ~IPTG)")
    ax1.set_ylabel("fold change vs baseline")
    ax1.set_title("homeostasis: autoregulation buffers Hfq level")
    ax1.legend(fontsize=8.5)
    ax1.text(0.03, 0.97, f"at 8x transcription:\n  mRNA {m_fold[-1]:.1f}x\n"
             f"  protein (no auto) {p_noauto[-1]:.1f}x\n  protein (auto)  {p_auto[-1]:.1f}x",
             transform=ax1.transAxes, fontsize=7.5, va="top", family="monospace",
             bbox=dict(boxstyle="round", fc="white", ec="gray", alpha=.8))

    # (2) 応答速度: 正規化した立ち上がり（自己制御ありは速い）
    def norm_tc(auto):
        t, p = timecourse(1.0, auto)
        return t, p / p[-1]
    ht = {}
    for auto, c, lab in [(False, "C0", "no autoreg"), (True, "C3", "autoregulated")]:
        t, y = norm_tc(auto)
        ax2.plot(t, y, "-", color=c, lw=2, label=lab)
        ht[auto] = t[np.where(y >= 0.5)[0][0]]
    ax2.axhline(0.5, ls=":", c="gray", lw=1)
    ax2.set_xlabel("time (min)")
    ax2.set_ylabel("Hfq protein / steady state")
    ax2.set_title("negative feedback speeds the response")
    ax2.legend(fontsize=9)
    ax2.text(0.55, 0.25, f"time to half:\n  no auto:  {ht[False]:.1f} min\n"
             f"  auto:     {ht[True]:.1f} min",
             transform=ax2.transAxes, fontsize=8, family="monospace",
             bbox=dict(boxstyle="round", fc="white", ec="gray", alpha=.8))

    fig.tight_layout()
    out = "outputs/13_autoregulation.png"
    fig.savefig(out, dpi=120)
    print("saved:", out)
    for a in [1, 3, 8]:
        au = steady(a, True); na = steady(a, False)
        print(f"alpha={a}: mRNA {au[0]/m0a:.2f}x  protein auto {au[1]/p0a:.2f}x  noauto {na[1]/p0n:.2f}x")


if __name__ == "__main__":
    main()
