"""04. sRNA の作用機構が「化学量論的」か「触媒的」かで応答形状が変わる。

02/03 の共分解は化学量論的 (stoichiometric): sRNA も mRNA と一緒に消える。

    m + s -> ∅        sRNA は 1:1 で消費される  -> 鋭い閾値 (threshold-linear)

これを触媒的 (catalytic) に変える: sRNA は mRNA を壊すが自分は再利用される。

    m + s -> s        sRNA は残る (酵素的)      -> なだらかな滴定 (Michaelis-Menten 的)

触媒の定常 sRNA は s = a_s / b_s で mRNA と独立。よって
    <m> = a_m / (b_m + k_cat * s) = a_m / (b_m + k_cat * a_s / b_s)
という双曲線 (Hill 係数 1) になり、閾値スイッチは消える。

生物学的含意: 同じ「sRNA が mRNA を抑える」でも、消費されるか触媒かで
「デジタルなスイッチ」か「アナログな可変抵抗」かが決まる。
"""

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

from ecell4 import run_simulation
from ecell4_base.core import NetworkModel, ReactionRule, Species

M, S = Species("m"), Species("s")
A_M, B_M, B_S = 10.0, 1.0, 1.0
K_STO = 100.0   # 化学量論的ペアリング (強く速い)
K_CAT = 0.1     # 触媒的分解 (半減点を a_s=10 に合わせた)


def build(a_s, mode):
    mdl = NetworkModel()
    rules = [
        ReactionRule([], [M], A_M),
        ReactionRule([], [S], a_s),
        ReactionRule([M], [], B_M),
        ReactionRule([S], [], B_S),
    ]
    if mode == "stoichiometric":
        rules.append(ReactionRule([M, S], [], K_STO))    # m + s -> ∅
    else:
        rules.append(ReactionRule([M, S], [S], K_CAT))   # m + s -> s
    for rr in rules:
        mdl.add_reaction_rule(rr)
    return mdl


def ode_m(a_s, mode):
    ret = run_simulation(150.0, y0={"m": 0, "s": 0}, model=build(a_s, mode),
                         solver="ode", ndiv=1500, species_list=["m", "s"])
    return ret.as_array()[-1][1]


def main():
    a_grid = np.linspace(0, 20, 41)
    sto = np.array([ode_m(a, "stoichiometric") for a in a_grid])
    cat = np.array([ode_m(a, "catalytic") for a in a_grid])

    fig, ax = plt.subplots(figsize=(7.5, 5))
    ax.plot(a_grid, sto / A_M, "-", color="C0", lw=2,
            label="stoichiometric  m+s→∅  (threshold switch)")
    ax.plot(a_grid, cat / A_M, "-", color="C3", lw=2,
            label="catalytic  m+s→s  (graded / MM-like)")
    ax.axvline(A_M, ls=":", c="gray")
    ax.axhline(0.5, ls=":", c="gray", lw=1)
    ax.set_xlabel("sRNA transcription rate  a_s")
    ax.set_ylabel("relative mRNA  <m> / <m>(a_s=0)")
    ax.set_title("silencing mode sets response shape")
    ax.legend()

    fig.tight_layout()
    out = "outputs/04_catalytic_vs_stoichiometric.png"
    fig.savefig(out, dpi=120)
    print("saved:", out)
    for a in [0, 5, 10, 15, 20]:
        print(f"a_s={a:4.0f}  stoich m/m0={ode_m(float(a),'stoichiometric')/A_M:5.2f}"
              f"   catalytic m/m0={ode_m(float(a),'catalytic')/A_M:5.2f}")


if __name__ == "__main__":
    main()
