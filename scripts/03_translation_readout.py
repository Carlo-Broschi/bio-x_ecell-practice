"""03. sRNA サイレンシング回路に翻訳段を足し、タンパク質出力のノイズを見る。

02 の 5 反応に、翻訳とタンパク質分解を追加:

    m -> m + P    (k_p : 翻訳。mRNA は触媒として残る)
    P -> ∅        (b_p : タンパク質分解)

問い: mRNA レベルのゆらぎは、下流のタンパク質でどうなる?
  - タンパク質が mRNA より長寿命 (b_p < b_m) なら、時間平均 (低域通過フィルタ) が効き
    CV_P < CV_m。ただし閾値近傍のノイズ増幅は減衰しつつ P へ伝播する。
"""

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

from ecell4 import run_simulation
from ecell4_base.core import NetworkModel, ReactionRule, Species

M, S, PROT = Species("m"), Species("s"), Species("P")
A_M, B_M, B_S, K = 10.0, 1.0, 1.0, 100.0
K_P, B_P = 5.0, 0.5  # 翻訳速度 / タンパク質分解 (b_p < b_m: P は長寿命)


def build_model(a_s):
    mdl = NetworkModel()
    for rr in [
        ReactionRule([], [M], A_M),
        ReactionRule([], [S], a_s),
        ReactionRule([M], [], B_M),
        ReactionRule([S], [], B_S),
        ReactionRule([M, S], [], K),
        ReactionRule([M], [M, PROT], K_P),  # m -> m + P
        ReactionRule([PROT], [], B_P),       # P -> ∅
    ]:
        mdl.add_reaction_rule(rr)
    return mdl


def ode_steady(a_s):
    ret = run_simulation(150.0, y0={"m": 0, "s": 0, "P": 0}, model=build_model(a_s),
                         solver="ode", ndiv=1500, species_list=["m", "s", "P"])
    m, _, p = ret.as_array()[-1][1:4]
    return m, p


def gillespie_cv(a_s, t_end=600.0, ndiv=6000, burn=0.25, seed=0):
    ret = run_simulation(t_end, y0={"m": 0, "s": 0, "P": 0}, model=build_model(a_s),
                         solver="gillespie", ndiv=ndiv,
                         species_list=["m", "s", "P"], rndseed=seed)
    tail = ret.as_array()[int(ndiv * burn):]
    m, p = tail[:, 1], tail[:, 3]
    cv = lambda x: x.std() / x.mean() if x.mean() > 0.3 else np.nan
    return cv(m), cv(p)


def main():
    a_grid = np.array([0, 2, 4, 6, 8, 9, 10, 11, 12, 14, 16], dtype=float)
    m_ode, p_ode = np.array([ode_steady(a) for a in a_grid]).T
    cv_m, cv_p = np.array([gillespie_cv(a) for a in a_grid]).T

    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(12, 4.5))

    ax1.plot(a_grid, m_ode, "o-", color="C0", label="mRNA  <m>")
    ax1.plot(a_grid, p_ode / (K_P / B_P), "s--", color="C4",
             label="Protein  <P> / (k_p/b_p)")
    ax1.axvline(A_M, ls=":", c="gray")
    ax1.set_xlabel("sRNA transcription rate  a_s")
    ax1.set_ylabel("steady state (protein rescaled to mRNA)")
    ax1.set_title("protein inherits the threshold response")
    ax1.legend()

    ax2.plot(a_grid, cv_m, "o-", color="C1", label="CV of mRNA")
    ax2.plot(a_grid, cv_p, "s-", color="C4", label="CV of Protein")
    ax2.axvline(A_M, ls=":", c="gray")
    ax2.set_xlabel("sRNA transcription rate  a_s")
    ax2.set_ylabel("noise  CV = SD/mean")
    ax2.set_title("protein filters noise (CV_P < CV_m) but threshold propagates")
    ax2.legend()

    fig.tight_layout()
    out = "outputs/03_translation_readout.png"
    fig.savefig(out, dpi=120)
    print("saved:", out)
    for a, mo, po, cm, cp in zip(a_grid, m_ode, p_ode, cv_m, cv_p):
        print(f"a_s={a:4.0f}  m={mo:6.2f}  P={po:7.2f}  CV_m={cm:5.2f}  CV_P={cp:5.2f}")


if __name__ == "__main__":
    main()
