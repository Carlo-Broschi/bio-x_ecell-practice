"""08. sRNA どうしが Hfq を奪い合う — 競合とクロストーク
       （Moon & Gottesman, Mol Microbiol 2011）。

出典: Moon K, Gottesman S (2011) "Competition among Hfq-binding small RNAs in
Escherichia coli." Mol Microbiol 82(6):1545-1562.

この論文は実験（Northern blot）が主で式は無いが、確立した結論は明快:
  - Hfq は量が限られた共有資源。ある sRNA(例 OxyS)を過剰発現すると Hfq を奪い、
    他の sRNA(例 DsrA)の蓄積・活性が下がる（競合＝クロストーク）。
  - Hfq を増やすと、この競合は緩和・逆転する（rescue）。

そこで 06/07 の Hfq 陽モデルを **2つの sRNA-mRNA 系（s1/m1 と s2/m2）が1つの Hfq を共有**する形に
拡張し、上の2予測を再現する。s2/m2 を「レポーター」、s1 を過剰発現する「競合者」とする。
"""

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

from ecell4 import run_simulation
from ecell4_base.core import NetworkModel, ReactionRule, Species

H = Species("H")
BETA, KA, KD, K5 = 1.0, 1.0, 1.0, 10.0


def R(a, b, k):
    return ReactionRule(a, b, k)


def pair_rules(i, a_s, a_m):
    s, m, sH, mH, T, D = (Species(f"{x}{i}") for x in ["s", "m", "sH", "mH", "T", "D"])
    return [
        R([], [s], a_s), R([], [m], a_m), R([s], [], BETA), R([m], [], BETA),
        R([s, H], [sH], KA), R([sH], [s, H], KD),
        R([m, H], [mH], KA), R([mH], [m, H], KD),
        R([sH, m], [T], KA), R([T], [sH, m], KD),
        R([mH, s], [T], KA), R([T], [mH, s], KD),
        R([T], [D, H], K5), R([sH], [H], BETA), R([mH], [H], BETA),
        R([T], [H], BETA), R([D], [], BETA),
    ]


NAMES = ["H"] + [f"{x}{i}" for i in (1, 2) for x in ["s", "m", "sH", "mH", "T", "D"]]


def build(h_tot, a_s1):
    mdl = NetworkModel()
    for r in pair_rules(1, a_s1, 10.0) + pair_rules(2, 10.0, 10.0):
        mdl.add_reaction_rule(r)
    return mdl


def activities(h_tot, a_s1):
    ret = run_simulation(600.0, y0={"H": h_tot}, model=build(h_tot, a_s1),
                         solver="ode", ndiv=1, species_list=NAMES)
    v = dict(zip(NAMES, ret.as_array()[-1][1:]))
    def pct(i):
        tot = v[f"m{i}"] + v[f"mH{i}"] + v[f"T{i}"] + v[f"D{i}"]
        return 100 * v[f"D{i}"] / tot if tot > 0 else 0.0
    return pct(1), pct(2)


def main():
    a_grid = np.logspace(0, 3, 31)   # 競合 sRNA s1 の転写量（過剰発現）
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(12, 4.6))

    # (左) レポーター s2 の活性 vs 競合 s1、Hfq 律速 vs 潤沢
    for h_tot, c, lab in [(12.0, "C3", "limiting Hfq (=12)"),
                          (40.0, "C0", "abundant Hfq (=40)")]:
        y = np.array([activities(h_tot, a)[1] for a in a_grid])
        ax1.plot(np.log10(a_grid), y, "-", color=c, lw=2, label=lab)
    ax1.set_xlabel("competitor sRNA s1 transcription  (log10)")
    ax1.set_ylabel("reporter s2 activity  (% m2 in duplex)")
    ax1.set_title("overexpressing s1 steals Hfq → suppresses s2\n(extra Hfq rescues)")
    ax1.legend(fontsize=9)

    # (右) Hfq 律速下での取り合い: s1 が勝ち s2 が負ける（ヒエラルキー）
    a1 = np.array([activities(12.0, a)[0] for a in a_grid])
    a2 = np.array([activities(12.0, a)[1] for a in a_grid])
    ax2.plot(np.log10(a_grid), a1, "-", color="C1", lw=2, label="competitor s1 activity")
    ax2.plot(np.log10(a_grid), a2, "-", color="C4", lw=2, label="reporter s2 activity")
    ax2.set_xlabel("competitor sRNA s1 transcription  (log10)")
    ax2.set_ylabel("activity  (% own mRNA in duplex)")
    ax2.set_title("limiting Hfq: s1 wins the tug-of-war, s2 loses")
    ax2.legend(fontsize=9)

    fig.tight_layout()
    out = "outputs/08_srna_competition.png"
    fig.savefig(out, dpi=120)
    print("saved:", out)
    for a in [1, 10, 100, 1000]:
        lo = activities(12.0, a); hi = activities(40.0, a)
        print(f"aS1={a:5g}  Hfq12: s1={lo[0]:5.1f} s2={lo[1]:5.1f} | Hfq40: s2={hi[1]:5.1f}")


if __name__ == "__main__":
    main()
