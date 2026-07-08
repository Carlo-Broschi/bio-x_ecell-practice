# bio-d: E-Cell4 practice

[E-Cell4](https://ecell4.e-cell.org/) を使った細胞シミュレーションの練習リポジトリ。
決定論（ODE）と確率論（Gillespie）を中心に、反応ネットワークのモデル化を段階的に手を動かして学ぶ。

## 環境

- Python 3.14 / macOS arm64（`ecell4_base` の公式 wheel が対応）
- パッケージ管理は [`uv`](https://docs.astral.sh/uv/)

```bash
git clone https://github.com/Carlo-Broschi/bio-d_ecell-practice.git
cd bio-d_ecell-practice
uv sync                       # .venv を作成し依存を同期
uv run jupyter lab            # Notebook で練習する場合
```

スクリプトを直接実行する場合:

```bash
uv run python scripts/01_reversible_reaction.py   # outputs/ に PNG を保存
```

## 構成

| パス | 内容 |
|---|---|
| `notebooks/01_getting_started.ipynb` | 可逆反応 A ⇌ B を ODE と Gillespie で解く最初の練習 |
| `scripts/01_reversible_reaction.py`  | 上と同内容のスクリプト版（ヘッドレスで図を保存・動作検証用） |
| `outputs/` | 生成物（`.gitignore` 済み。再現可能なので追跡しない） |

## E-Cell4 の最小 API メモ

```python
from ecell4 import get_model, reaction_rules, run_simulation

with reaction_rules():
    A > B | 0.25      # A -> B, 速度定数 0.25
    B > A | 1.0
model = get_model()

ret = run_simulation(10.0, y0={'A': 60, 'B': 0}, model=model,
                     solver='ode',        # 'ode' | 'gillespie' | 'meso' | ...
                     species_list=['A', 'B'])
ret.as_array()        # ndarray（各行 = [t, A, B]）
ret.as_dataframe()    # pandas DataFrame（要 pandas）
ret.plot()            # 組み込みプロット
```

- 可逆反応の平衡比は `[A]/[B] = kr/kf`。上の例なら `1.0/0.25 = 4`（60 分子 → A:B = 48:12）。
- `solver` を差し替えるだけで決定論／確率論を切り替えられるのが E-Cell4 の勘所。

## 参考

- E-Cell4 ドキュメント: https://ecell4.e-cell.org/
- プロジェクトページ: https://www.e-cell.org/projects/ecell4.html
