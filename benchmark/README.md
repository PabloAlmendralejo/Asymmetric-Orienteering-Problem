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
- Budget tightness: loose (~70%), medium (~50%), tight (~30%) of nearest-neighbour tour cost
- Plus 4 extreme cases (symmetric easy, asymmetric hard, high asymmetry no fatigue, symmetric high fatigue)

## Solvers

Two solver variants are provided:

| Solver | File | Formulation |
|---|---|---|
| Original B&C | `benchmark_solver.cpp` | MTZ time propagation + McCormick linearisation of fatigue |
| Flow B&C | `benchmark_solver_flow.cpp` | Single-commodity flow variables, fatigue budget is linear (no McCormick) |

All solvers share the same SA warm-start (seed 42, deterministic) for fair comparison.

### Compile and run

```bash
# Original B&C
cl.exe /O2 /std:c++20 /EHsc benchmark_solver.cpp /I<highs_include> /link highs.lib /out:benchmark_solver.exe
benchmark_solver.exe instances

# Flow B&C
cl.exe /O2 /std:c++20 /EHsc benchmark_solver_flow.cpp /I<highs_include> /link highs.lib /out:benchmark_solver_flow.exe
benchmark_solver_flow.exe instances
```

## Heuristic comparison

A separate heuristic comparison evaluates four methods under equal 2-second time budgets
across all 21 instances (10 seeds per stochastic method):

```bash
cl.exe /O2 /std:c++20 /EHsc compare_heuristics.cpp /I<highs_include> /link highs.lib /out:compare_heuristics.exe
compare_heuristics.exe instances
```

| Method | Mean gap below SA | Std | Min gap | Max gap |
|---|---|---|---|---|
| Greedy | 23.0% | 9.4% | 3.7% | 40.0% |
| GA | 49.2% | 11.0% | 24.6% | 71.9% |
| ACO | 8.1% | 4.9% | 1.3% | 24.1% |
| SA | — | — | — | — |

SA significantly outperforms all alternatives (Wilcoxon signed-rank, all p < 0.001).

## Output

Each instance produces:
- `op_output_*.json` with SA and B&C results (route, points, timing, optimality status, gap)
- Summary table printed to stdout with points, nodes, time, optimality status, and LP gap

## Instance format

Each `op_input_*.json` contains:
- `cm`: asymmetric cost matrix (n+1 × n+1, node 0 = depot)
- `pts`: score for each node (node 0 = 0)
- `bud_raw`: base time budget
- `bud_eff`: effective budget (with fatigue)
- `fatigue_rate`: λ parameter
