"""可逆反応 A <=> B を E-Cell4 で解く最初の練習。

同じモデルを 2 つのソルバで回して比較する:
  - ode      : 決定論的（常微分方程式）。濃度は滑らかな曲線。
  - gillespie: 確率論的（Gillespie 法）。分子数が少ないと揺らぐ。

平衡では正逆反応の速度が釣り合うので
    kf * [A] = kr * [B]  =>  [A]/[B] = kr / kf = 1.0 / 0.25 = 4
となり、初期 60 分子は A:B = 48:12 付近に落ち着く。
"""

import matplotlib

matplotlib.use("Agg")  # ヘッドレス（GUI なし）で PNG に保存するため
import matplotlib.pyplot as plt

from ecell4 import get_model, reaction_rules, run_simulation

KF, KR = 0.25, 1.0
T_END, N_DIV = 10.0, 200
Y0 = {"A": 60, "B": 0}


def main():
    with reaction_rules():
        A > B | KF  # noqa: F821
        B > A | KR  # noqa: F821
    model = get_model()

    fig, axes = plt.subplots(1, 2, figsize=(11, 4), sharey=True)
    for ax, solver in zip(axes, ["ode", "gillespie"]):
        ret = run_simulation(
            T_END, y0=Y0, model=model, solver=solver,
            ndiv=N_DIV, species_list=["A", "B"], rndseed=0,
        )
        arr = ret.as_array()          # 各行 = [t, A, B]
        t = [row[0] for row in arr]
        ax.plot(t, [row[1] for row in arr], label="A")
        ax.plot(t, [row[2] for row in arr], label="B")
        ax.axhline(48, ls=":", c="gray", lw=1)
        ax.axhline(12, ls=":", c="gray", lw=1)
        ax.set_title(f"A ⇌ B  ({solver})")
        ax.set_xlabel("time")
        ax.legend()
        print(f"[{solver:9s}] t=end -> {arr[-1]}")
    axes[0].set_ylabel("# molecules")
    fig.tight_layout()
    out = "outputs/01_reversible_reaction.png"
    fig.savefig(out, dpi=120)
    print("saved:", out)


if __name__ == "__main__":
    main()
