# 公式トラック — E-Cell4 公式例題の移植

このディレクトリは **E-Cell4 公式ドキュメント**（https://ecell4.e-cell.org/ ）の\
**チュートリアル／例題を、そのまま移植**したもの。`../notebooks/`（bio-d オリジナルの自作ノート）とは区別する。

各ノートは冒頭に**出典 URL** を明記し、公式コードを可能な限りそのまま使う。インストール版（ecell4 1.2.2）で走るための\
最小限の調整（重いソルバの差し替え・明らかな doc 誤植の訂正）は本文に明示する。すべて `nbconvert --execute` で実行検証済み。

| ノート | 公式ページ | 学べる公式機能 |
|---|---|---|
| `o1_simple_equilibrium.ipynb` | Examples / Simple Equilibrium | 低レベル `Factory→world→simulator→step`、`create_binding/unbinding_reaction_rule`、ソルバ差し替え（egfrd 等） |
| `o2_lotka_volterra.ipynb` | Examples / Lotka-Volterra 2D | 速度式(rate-law) ODE、素反応 Gillespie、空間 meso への拡張 |
| `o3_dual_phosphorylation.ipynb` | Examples / Dual Phosphorylation Cycle | `@species_attributes`/`@reaction_rules` デコレータ、酵素反応の連結記法、MAPK 型カスケード |

## 公式が用意している素材の全体像（参考）

公式は「解くべき練習問題」ではなく、**チュートリアル（how-to）と完成した例題モデル**を提供している。

- **チュートリアル(10)**: Brief Tour / How to Build a Model / Initial Condition / Run a Simulation /
  Log & Visualize / Solve ODEs with Rate Law / Rule-based Modeling / More about Brief Tour /
  Spatial Gillespie / Spatiocyte（単分子）。
- **例題(Examples)**: Attractors, Drosophila 概日時計, Dual Phosphorylation Cycle, Simple EGFR,
  Glycolysis & MCA, Action Potentials, **Lotka-Volterra 2D**, MinDE(meso/Spatiocyte),
  **Simple Equilibrium**, Tyson1991(細胞周期), Unit System, GPCR, sGFRD。
- **テストモデル**: Birth-Death, Homodimerization & Annihilation, Reversible(＋拡散律速), MSD。

本トラックはこのうち代表 3 例を移植した。残り（空間 Spatiocyte MinDE、Rule-based、Glycolysis/MCA、Tyson 細胞周期 など）も\
同じ要領で追加可能。
