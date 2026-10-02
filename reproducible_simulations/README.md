# Competition, Switching, and Control in Coupled Adoption-Opinion Dynamics

MATLAB code to reproduce the simulations of the paper by
M. Alutto, F. Dabbene, A. Fontan, K. H. Johansson, C. Ravazzi.

## Usage

```matlab
run_all                          % all experiments
run_all({'fig3_uniqueness'})     % selected folders only
```

Figures and logs are written to `results/output/<folder>/`.
To run a single script, call `setup_paths` first.

Requires MATLAB R2020a or later (no toolboxes).

## Contents

| Folder | Script | Paper |
|---|---|---|
| `experiments/fig1_market_share` | `realdata_ex.m` | Fig. 1 |
| `experiments/fig3_uniqueness` | `adop_multi_unique.m`, `adop_multi_nonunique.m` | Fig. 3(a), 3(b) |
| `experiments/fig4_entry` | `kingmaker.m` | Fig. 4 |
| `experiments/fig5-7_optimal_control` | `run_control_three_problems_noD.m` | Figs. 5–7 |
| `experiments/extra_kingmaker_analysis` | `kingmaker_fig1_homogeneous.m`, `kingmaker_phasemap_region.m` | not in paper |
| `experiments/extra_problem4_kingmaker_control` | `run_control_problem4_*.m` | not in paper |
| `experiments/exploratory` | `adop_multi.m` | not in paper |

- `src/`: simulation engines for model (1) and the LHS of condition (9).
- `results/reference_figures/`: the figures as originally produced (`paper/` and `extra/`).
- Random seeds are fixed in each script with `rng`.

## Expected output (Fig. 3)

| | LHS of (9) | max. difference in final mean adoption between the two runs |
|---|---|---|
| Fig. 3(a) | 0.6283 | 8.49e-06 |
| Fig. 3(b) | 2159.02 | 1.42e-02 |
