# Benchmark Suite

Parameterized asymmetric OP instances with fatigue for evaluating solver performance.
To the best of our knowledge, this is the first benchmark set for the asymmetric OP with terrain-based costs and fatigue.

## Generate instances

```bash
cd benchmark
python generate_instances.py --output-dir instances
```

This creates 21 instances varying:
- Node count: 20, 30, 40, 50, 75, 100
- Asymmetry: 0%, 10%, 25%, 50%
- Fatigue rate: 0.0, 0.1, 0.2, 0.3
- Budget tightness: the budget B is set to 70% (loose), 50% (medium), or 30% (tight) of the nearest-neighbour tour cost
- Plus 4 extreme cases (symmetric easy, asymmetric hard, high asymmetry no fatigue, symmetric high fatigue)

## Solvers

Four solver variants are provided:

| Solver | File | Formulation |
|---|---|---|
| MTZ B&C | `benchmark_solver.cpp` | MTZ time propagation + McCormick linearisation |
| Flow B&C | `benchmark_solver_flow.cpp` | Single-commodity flow + tightened coupling + all cuts (B1-B4) |
| Ablation | `benchmark_solver_ablation.cpp` | Flow B&C with command-line flags to toggle each cut family |
| Heuristics | `compare_heuristics.cpp` | Greedy / GA / ACO / SA comparison |

All B&C solvers share the same SA warm-start (deterministic, seed 42) and use best-first search for benchmark instances.

### Compile and run

```bash
# MTZ B&C
cl.exe /O2 /std:c++20 /EHsc benchmark_solver.cpp /I<highs_include> /link highs.lib /out:benchmark_solver.exe
benchmark_solver.exe instances

# Flow B&C
cl.exe /O2 /std:c++20 /EHsc benchmark_solver_flow.cpp /I<highs_include> /link highs.lib /out:benchmark_solver_flow.exe
benchmark_solver_flow.exe instances

# Ablation (with flags)
cl.exe /O2 /std:c++20 /EHsc benchmark_solver_ablation.cpp /I<highs_include> /link highs.lib /out:benchmark_solver_ablation.exe
benchmark_solver_ablation.exe instances --config=all --coupling=on --fatigue-covers=on --routing=on --cycle=on --path=on --csv=results.csv
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

### Key results (21 instances, 900s time limit)

**Formulation comparison (MTZ vs Flow):**

| Metric | MTZ | Flow |
|---|---|---|
| Mean gap | 16.0% | 5.4% |
| Instances optimal | 4/21 | 7/21 |
| Max gap | 42.0% | 10.5% |

**Ablation study (flow formulation, cumulative):**

| Config | Cuts enabled | Mean gap | Opt. |
|---|---|---|---|
| 1 | Base (SECs + conn. + tight. + fat.) | 5.54% | 6/21 |
| 2 | + B1 (lifted covers) | 5.85% | 7/21 |
| 3 | + B2 (routing infeas.) | 5.41% | 7/21 |
| 4 | + B3 (cycle covers) | 5.35% | 7/21 |
| 5 | + B4 (path ineq.) | 5.42% | 7/21 |

**Heuristic comparison (equal 2s budgets, 10 seeds):**

| Method | Mean gap below SA |
|---|---|
| Greedy | 23.0% |
| GA | 49.2% |
| ACO | 8.1% |
| SA | (reference) |

All differences significant (Wilcoxon signed-rank, p < 0.001).

## Output

Each instance produces:
- `op_output_*.json` with SA and B&C results (route, points, timing, optimality status, gap)
- Summary table printed to stdout

## Instance format

Each `op_input_*.json` contains:
- `cm`: asymmetric cost matrix (n+1 × n+1, node 0 = depot)
- `pts`: score for each node (node 0 = 0)
- `bud_raw`: base time budget
- `bud_eff`: effective budget (with fatigue)
- `fatigue_rate`: λ parameter
