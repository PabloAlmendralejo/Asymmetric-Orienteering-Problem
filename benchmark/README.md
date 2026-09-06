# Benchmark Suite

Parameterized asymmetric OP instances with fatigue for evaluating solver performance.
To the best of our knowledge, this is the first benchmark set for the asymmetric OP with terrain-based costs and fatigue.

## Generate instances

```bash
cd benchmark
python generate_instances.py --output-dir instances
```

This creates 66 instances (3 seeds per configuration) varying:
- Node count: 20, 30, 40, 50, 75, 100, 125, 150, 200
- Asymmetry: 0%, 10%, 25%, 50%
- Budget tightness: the budget B is set to 70% (loose), 50% (medium), or 30% (tight) of the nearest-neighbour tour cost
- Plus 6 extreme/boundary cases and 2 extra size×asymmetry combinations outside the main sweeps

All instances are solved at the same fixed `λ = 4.1e-5`, `ρ = 0` (see the top-level README's "Fatigue Model" section for the calibration sources) rather than swept per instance — an earlier version of this suite also varied a nominal fatigue rate per instance, but that dimension no longer reflects a real experimental contrast once λ/ρ are fixed, so it isn't a reported results axis.

## Solvers

Five solver variants are provided:

| Solver | File | Formulation |
|---|---|---|
| MTZ B&C | `benchmark_solver_mtz.cpp` | MTZ time propagation + McCormick linearisation |
| Flow B&C | `benchmark_solver_flow.cpp` | Single-commodity flow + tightened coupling + all cuts (B1-B4) |
| Ablation | `benchmark_solver_ablation.cpp` | Flow B&C with command-line flags to toggle each cut family |
| Classical baseline | `solver_op_classical.cpp` | Undirected, symmetric, fatigue-free OP (Fischetti et al. 1998) |
| Heuristics | `compare_heuristics.cpp` | Greedy / GA / ACO / SA comparison |

All B&C solvers share the same SA warm-start (deterministic, seed 42) and use best-first search for benchmark instances.

### Compile and run

```bash
# MTZ B&C
cl.exe /O2 /std:c++20 /EHsc benchmark_solver_mtz.cpp /I<highs_include> /link highs.lib /out:benchmark_solver_mtz.exe
benchmark_solver_mtz.exe instances --fixed-lambda=4.1e-5 --fixed-rho=0

# Flow B&C
cl.exe /O2 /std:c++20 /EHsc benchmark_solver_flow.cpp /I<highs_include> /link highs.lib /out:benchmark_solver_flow.exe
benchmark_solver_flow.exe instances --fixed-lambda=4.1e-5 --fixed-rho=0

# Ablation (with flags)
cl.exe /O2 /std:c++20 /EHsc benchmark_solver_ablation.cpp /I<highs_include> /link highs.lib /out:benchmark_solver_ablation.exe
benchmark_solver_ablation.exe instances --fixed-lambda=4.1e-5 --fixed-rho=0 --config=all --coupling=on --fatigue-covers=on --routing=on --cycle=on --path=on --csv=results.csv
```

### Ablation study

The ablation solver supports the following flags (all default to `on`):

| Flag | Controls |
|---|---|
| `--coupling=on/off` | Tightened flow-arc coupling bounds |
| `--fatigue-covers=on/off` | Fatigue-aware cover cut weights |
| `--routing=on/off` | Routing infeasibility cuts (B2) |
| `--cycle=on/off` | Directed cycle cover cuts (B3) |
| `--path=on/off` | Directed path inequality cuts (B4) |
| `--config=<name>` | Label for CSV output |
| `--csv=<file>` | Append results to CSV file |

### Key results (66 instances, 900s time limit, λ = 4.1e-5, ρ = 0)

**Formulation comparison (MTZ vs Flow):**

| Metric | MTZ | Flow |
|---|---|---|
| Mean gap | 26.2% | 14.3% |
| Instances optimal | 9/66 | 24/66 |
| Max gap | 64.1% | 57.6% |

**Ablation study (flow formulation, cumulative):**

| Config | Cuts enabled | Mean gap | Opt. |
|---|---|---|---|
| 1 | Base (SECs + conn. + tight. + fat.) | 13.29% | 24/66 |
| 2 | + B1 (lifted covers) | 14.64% | 23/66 |
| 3 | + B2 (routing infeas.) | 14.29% | 25/66 |
| 4 | + B3 (cycle covers) | 12.80% | 22/66 |
| 5 | + B4 (path ineq.) | 14.78% | 23/66 |

No cut family improves both mean gap and proven-optimal count simultaneously on this suite — see the paper's ablation section for the instance-by-instance discussion of why.

**Warm-start ablation (config 5, with vs. without the SA incumbent):**

| Warm start | Proven optimal | Total failures (0 pts found) | Total pts (66 instances) |
|---|---|---|---|
| Off | 23/66 | 42/66 | 16,550 |
| On | 24/66 | 0/66 | 70,540 |

**Heuristic comparison (equal 2s budgets, 10 seeds per stochastic method):**

| Method | Mean gap below SA | Std | Min | Max |
|---|---|---|---|---|
| Greedy | 19.0% | 8.8% | 3.4% | 40.0% |
| GA | 54.1% | 14.8% | 15.7% | 81.1% |
| ACO | 8.6% | 5.1% | 0.0% | 22.9% |
| SA | (reference) | | | |

All differences significant: exact Wilcoxon signed-rank test on the 66 paired instance means, Holm-Bonferroni corrected across the 3 comparisons ($p_{\text{Holm}} \approx 4.9\times10^{-12}$ for each), with a matched-pairs rank-biserial effect size of $r = 1.00$ in every case (SA's mean score exceeds every competitor's on all 66 instances, no ties). Run `python analyze_heuristics.py` after `compare_heuristics.exe` to reproduce these statistics from its CSV output.

## Output

Each instance produces:
- `op_output_*.json` with SA and B&C results (route, points, timing, optimality status, gap)
- Summary table printed to stdout

## Instance format

Each `op_input_*.json` contains:
- `cm`: asymmetric cost matrix (n+1 × n+1, node 0 = depot)
- `gain`/`loss`: cumulative uphill/downhill elevation per arc (metres), for `φ+`/`φ-`
- `dist`: path distance per arc (metres), for the `μ·δ` term
- `pts`: score for each node (node 0 = 0)
- `bud_raw`: base time budget
- `bud_eff`: effective budget (with fatigue)
- `fatigue_rate`: nominal λ from generation (overridden by `--fixed-lambda` at solve time; see "Generate instances" above)
- `rho_default`/`mu_default`: nominal ρ and derived μ from generation (ρ overridden by `--fixed-rho`)
