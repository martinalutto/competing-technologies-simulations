# Competing technologies — reproducible simulations

MATLAB code for the simulations in the TAC extension on **adoption–opinion
dynamics of competing technologies on networks**. There are *n* agents and
*m* technologies. Each agent *i* and technology *k* has these states:
undecided `S`, adopters `A`, dissatisfied `D` and opinion `X`. The controls
`u` act on the opinion dynamics.

## Quick start

```matlab
>> cd reproducible_simulations
>> run_all                          % run everything
>> run_all({'01_uniqueness'})       % or only some groups
```

`run_all` adds `src/` to the path and runs every script in each experiment
group, each inside its own output folder. It writes:

- `results/output/<group>/<script>_figN.png` (and `.fig` in MATLAB): every figure the script opens
- `results/output/<group>/<script>.log`: the console output (numbers to compare against below)
- `results/output/environment.txt`: the software version used

To run a single script by hand, call `setup_paths` once and then run the script.

## Requirements

- **MATLAB R2020a or later.** Scripts use `tiledlayout`, `exportgraphics`
  and local functions at the end of script files. No toolboxes are needed.
- **GNU Octave 11** runs groups `01`, `02` and `04`. In `04`, only the final
  `exportgraphics` call fails, and the figure is still saved. Octave cannot
  run group `03` because those scripts use script-local functions.

## Repository layout

```
src/                      shared simulation engines (functions)
  simulate_adoption_multi.m   S/A/D/X dynamics, m technologies, vectorized over agents
  simulate_kingmaker.m        same dynamics with per-technology entry times t_entry(k)
  build_params_assumption.m   W, tilde_W row-stochastic + beta, gamma satisfying Assumption (i),(iv)
  hyp_unique_value.m          left-hand side of eq:hyp-unique
experiments/
  01_uniqueness/          equilibrium uniqueness (eq:hyp-unique satisfied vs violated)
  02_kingmaker/           kingmaker effect: closed form, phase map, kingmaker region
  03_optimal_control/     optimal control (forward-backward sweep / PMP), Problems 1-4
  04_real_data/           smartphone market-share motivating example
  exploratory/            original exploratory scripts (not run by run_all)
results/reference_figures/  figures as produced for the paper (to compare against)
docs/Simulazioni.pdf      specification of the simulations / figures
```

## Experiments → figures

| Script | Seed | Output | Reference figure | Octave runtime |
|---|---|---|---|---|
| `01_uniqueness/adop_multi_unique.m` | `rng(1)` | Case A: eq:hyp-unique holds, so two initial conditions give the same equilibrium | `sim_unique.png/.eps` | < 1 s |
| `01_uniqueness/adop_multi_nonunique.m` | `rng(3)` | Case B: eq:hyp-unique is violated, so the equilibrium depends on A(0) | `sim_nounique.png/.eps` | ~13 s |
| `02_kingmaker/kingmaker_fig1_homogeneous.m` | none (deterministic) | Fig. 1: closed-form ratio a²*/a¹* vs z₃ (homogeneous case) | `kingmaker_fig1_homogeneous.png` | ~3 s |
| `02_kingmaker/kingmaker_phasemap_region.m` | `rng(3)` | Fig. 2: heatmap of ΔA = A²−A¹ over (γ³, ξ³); Fig. 3: kingmaker region | `kingmaker_fig2_heatmap.png`, `kingmaker_fig3_region.png` | ~4 min (30×30 sweep) |
| `03_optimal_control/run_control_three_problems_noD.m` | `rng(2)` | Problems 1–3: cumulative adoption, targeted terminal adoption with budget | `prob1.png`, `prob2.png`, `prob3_new.png`, `avg_controls.png` | MATLAB only |
| `03_optimal_control/run_control_problem4_kingmaker_linear_budget.m` | `rng(4)` | Problem 4: kingmaker control, linear budget (bang-bang PMP) | `prob4.png` | MATLAB only |
| `03_optimal_control/run_control_problem4_kingmaker_quadratic_control.m` | `rng(2)` | Problem 4: kingmaker control, quadratic cost (projected gradient) | `prob4.png` | MATLAB only |
| `04_real_data/realdata_ex.m` | none (data hard-coded) | Smartphone market share 2006–2009 (BlackBerry, Windows Mobile, iOS) | `smartphone_market_share.pdf/.eps` | < 1 s |

> ⚠️ To confirm: which `03_optimal_control` script and which panel produced
> `prob3.png` / `prob3_new.png` / `avg_controls.png` / `prob4.png`. These
> were saved by hand from the figure window, so the scripts don't record it.

## Expected numerical output (checks)

Values below come from GNU Octave 11.3 on macOS (Apple Silicon), recorded
2026-10-02. MATLAB with the same seeds should give the same random streams
(Mersenne Twister). Small floating-point differences are possible.

**Case A** (`adop_multi_unique.m`)
- eq:hyp-unique LHS = 0.6283 (< 1)
- max |ΔA equilibrium| between the two runs = 8.49e-06; max |ΔX| = 4.74e-07

**Case B** (`adop_multi_nonunique.m`)
- eq:hyp-unique LHS = 2159.02 (≥ 1)
- max |ΔA equilibrium| = 1.42e-02; max |ΔX| = 1.31e-02

**Fig. 1** (`kingmaker_fig1_homogeneous.m`)
- R(0) = δ¹/δ² = 0.6667, R(∞) = 1.6667, crossing at z₃* = 0.25

**Figs. 2–3** (`kingmaker_phasemap_region.m`)
- Baseline without tech 3: A¹ = 6.436, A² = 5.440, so tech 1 wins
- Outcomes over the grid: no reversal 27.0% | kingmaker 27.2% | tech 3 dominant 45.9%

## Modelling notes

- **Assumptions.** `build_params_assumption.m` makes `W` and `tilde_W`
  row-stochastic with strictly positive entries, so `W` is irreducible
  (Assumption i). It rescales `beta` and `gamma` so that Σₖβᵢᵏ ≤ 1 and
  Σₖγᵢᵏ ≤ 1 (Assumption iv). Opinions need X(0) > 0 (Assumption iii).
- **Uniqueness test.** eq:hyp-unique is checked on the parameters as drawn,
  with `assert(val<1)` in Case A and `assert(val>=1)` in Case B. No parameter
  is tuned to force the condition. The factor in δ is ≥ 2 and grows like
  1/δ², so satisfying the condition with non-negligible ξ needs δ close to 1.
- **What uniqueness covers.** eq:hyp-unique (λ, ξ, δ) controls whether
  **opinions X** become independent of the initial condition. Whether the
  **adoption shares A** do also depends on the A↔D recirculation through
  γ and δ.
- **Entry times.** In `simulate_kingmaker.m`, a technology that has not yet
  entered (t < t_entry(k)) is frozen and adds nothing to any coupling term.
  With `t_entry = [1 1 1]` the engine matches the per-agent loop of
  `exploratory/kingmaker.m` to 3.5e-17.
- **Closed form.** The Fig. 1 formula is exact in the homogeneous case
  (n = 1). It matches a single-node simulation to ~1e-16.
- **Vectorization check.** `simulate_adoption_multi.m` matches the original
  per-agent loop in `exploratory/adop_multi.m` to ~1e-16.
