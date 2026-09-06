# Optimal Route Planning for Orienteering on Real Terrain

End-to-end framework for solving the Asymmetric Orienteering Problem with Fatigue (AOPF) on real terrain, combining a terrain-aware asymmetric cost model with exact Branch-and-Cut optimisation.

## Overview

The pipeline has two stages:

1. **Python** — Preprocesses orienteering map files (OMAP), digital elevation models (DEM),
   and map imagery into an asymmetric cost matrix using:
   - IOF symbol-based cost raster (ACR)
   - Hypsometry Cost Raster (HCR) with IQR scaling
   - Minetti metabolic slope model (directional)
   - A state-based, order-dependent fatigue model (see below)

2. **C++** — Solves the Orienteering Problem exactly via Branch-and-Cut with HiGHS,
   warm-started by Simulated Annealing.

## Fatigue Model

The original linear model, `f(t) = 1 + λt/B`, is provably order-invariant (the
same set of arcs gives the same total fatigue cost regardless of visit order),
which collapses the problem to a much simpler one — this was the central
objection in EJOR-D-26-01706's review. It's replaced by a per-node fatigue
*state* `G_i`, propagated arc-by-arc:

```
ψ_ij = φ+_ij − ρ·φ-_ij + μ·δ_ij
```

`φ+`/`φ-` are cumulative uphill/downhill elevation (metres) and `δ` is path
distance (metres) along the cost-optimal path for arc `(i,j)`. Rather than a
sensitivity sweep, `λ` and `ρ` are fixed to literature-cited values:
`λ = 4.1×10⁻⁵`, derived from grade-adjusted-pace decline over an ultra-trail
race (Jaén-Carrillo et al. 2025), and `ρ = 0`, since post-ultramarathon
downhill running economy has been shown to *worsen*, not recover
(Vernillo et al. 2015). `μ` is derived (not swept) from the Minetti curve
already used for the base cost. Full derivation and the calibration
discussion are in `paper/orienteering_paper.tex`, Sec. 3.6 and Sec. 7.

The true (clipped) fatigue state is `F_j = max(0, F_i + ψ_ij)`; the LP
relaxation uses the unclipped `G_i` (a valid upper bound on `F_i`) for the
McCormick (MTZ) and flow-coupling constraints — see Sec. 3 of the paper for
the full formulation, including why the flow formulation needs a *two-sided*
coupling bound (`Ĝ_i` and `Ǧ_i`) where the original model only needed one.

### Known gaps

- **Terrain diversity**: validated on two Spanish terrains (Torremocha,
  La Muela). Torremocha's IQR-based cost-fusion scaling is degenerate
  (`c = 0.00`, see Sec. 8 of the paper) because its own ACR cost distribution
  has near-zero spread — HCR morphometry effectively carries the whole
  combined-cost signal there, unlike La Muela (`c = 1.09`).
- **No modern exact-OP baseline**: compared only against a classical
  (undirected, symmetric, fatigue-free) Fischetti et al. (1998) baseline
  (`benchmark/solver_op_classical.cpp`). A modern exact OP solver such as
  Kobeaga et al. (2021) would need IBM CPLEX, which has no Windows ARM64
  build (this development machine's architecture) — out of scope for now.

## Project Structure

```
├── python/
│   ├── run_pipeline.py              # Main entry point
│   ├── requirements.txt
│   ├── core/
│   │   ├── preprocessing.py         # Raw files → terrain cache (.npz)
│   │   ├── cost_functions.py        # Minetti, fatigue (incl. derive_mu), ISOM weights
│   │   ├── terrain_analysis.py      # HCR (TRI, PC, SLO), IQR scaling
│   │   ├── omap_parser.py           # OMAP XML parsing, symbol classification
│   │   ├── rasterization.py         # Line/polygon rasterization
│   │   ├── coordinate_transforms.py # Pixel ↔ UTM ↔ grid transforms
│   │   ├── pathfinding.py           # Anisotropic Dijkstra, cost/gain/loss/dist matrices
│   │   ├── route_optimizer.py       # Greedy + SA route optimization
│   │   ├── control_placement.py     # Control point generation (7 strategies)
│   │   └── visualize_results.py     # Map, cost surfaces, routes, hillshade
│   └── config/
│       ├── torremocha.py
│       └── la_muela.py
├── cpp/
│   ├── solver_flow_torremocha.cpp   # Flow B&C solver, Torremocha
│   ├── solver_flow_la_muela.cpp     # Flow B&C solver, La Muela
│   ├── solver_mtz_torremocha.cpp    # MTZ B&C solver, Torremocha
│   └── solver_mtz_la_muela.cpp      # MTZ B&C solver, La Muela
├── benchmark/
│   ├── generate_instances.py        # Parameterized instance generator
│   ├── benchmark_solver_mtz.cpp     # MTZ B&C solver for benchmarks
│   ├── benchmark_solver_flow.cpp    # Flow B&C solver for benchmarks
│   ├── benchmark_solver_ablation.cpp # Ablation study solver (cut toggles)
│   ├── solver_op_classical.cpp      # Classical (undirected, fatigue-free) OP baseline
│   ├── compare_heuristics.cpp       # Greedy / GA / ACO / SA comparison
│   ├── analyze_heuristics.py        # Exact Wilcoxon + Holm-Bonferroni + effect size
│   ├── instances/                   # 66 synthetic benchmark instances
│   └── README.md
├── paper/
│   └── orienteering_paper.tex       # LaTeX paper
├── data/                            # Input terrain data (Git LFS)
│   ├── torremocha/
│   └── la_muela/
└── .gitignore
```

## Setup

### Python
```bash
cd python
pip install -r requirements.txt
```

### C++ (Branch-and-Cut solver)
- MSVC or GCC with C++20 support
- [HiGHS](https://highs.dev/) — needs headers + a compiled `.lib`/`.a`, not just
  the runtime `.dll`/`.so`. If not already built:
  ```bash
  git clone https://github.com/ERGO-Code/HiGHS.git
  cmake -B HiGHS/build -DCMAKE_BUILD_TYPE=Release -DBUILD_SHARED_LIBS=OFF -DFAST_BUILD=ON HiGHS
  cmake --build HiGHS/build --config Release
  ```

## Usage

### Full pipeline
```bash
cd python

# First run (preprocessing from raw files):
python run_pipeline.py torremocha --preprocess

# Subsequent runs (uses cached .npz):
python run_pipeline.py torremocha

# La Muela:
python run_pipeline.py la_muela --preprocess
```

### C++ solver
```bash
# Compile (Windows/MSVC), e.g. for the Torremocha flow solver:
cl.exe /O2 /MD /std:c++20 /EHsc solver_flow_torremocha.cpp ^
  /I<highs_include> /I<highs_build_dir> ^
  /link /LIBPATH:<highs_lib_dir> highs.lib /out:solver_flow_torremocha.exe

# Run (from directory with op_input_*.json files):
solver_flow_torremocha.exe
solver_flow_torremocha.exe --fixed-lambda=4.1e-5 --fixed-rho=0  # override lambda/rho
```

### Visualization
```bash
cd python/core
python visualize_results.py torremocha --output-dir figures/
python visualize_results.py la_muela --output-dir figures/
```

Generates: map overview, cost surface panels (ACR/HCR/combined/elevation/slope),
hillshade, directional asymmetry heatmap, and route overlay with A* traced paths.

## Study Areas

| | Torremocha | La Muela |
|---|---|---|
| Location | Cáceres, Spain | Salamanca, Spain |
| Area | 2.1 × 2.9 km | 2.8 × 3.2 km |
| Resolution | 2.0 m | 2.0 m |
| Relief | 55 m | 464 m |
| Cost asymmetry | 17.5% | 36.6% |
| IQR scaling (c) | 0.00 | 1.09 |

Across both terrains (7 control-point distributions each, 14 instances
total), the MTZ formulation proves 10/14 instances optimal within the
15-minute time limit versus 9/14 for the flow formulation, though flow's
mean LP-relaxation gap on the unproven instances is lower (4.81% vs.
5.57%). See the paper's Results section for the full per-distribution
tables and discussion.

## Benchmark Suite

Parameterized instances (66 total, 3 seeds per configuration) for the
asymmetric OP with fatigue — the first benchmark set for this problem
variant. See `benchmark/README.md` for details.

Instances vary across:
- Node count: 20, 30, 40, 50, 75, 100, 125, 150, 200
- Asymmetry: 0% to 50%
- Budget tightness: loose (70%), medium (50%), tight (30%) of NN tour cost

All instances are solved at the same fixed `λ = 4.1×10⁻⁵`, `ρ = 0` (see
"Fatigue Model" above) — an earlier version of this suite also varied a
nominal fatigue rate per instance and swept `λ`/`ρ` per run, but that
dimension no longer reflects a real experimental contrast under the fixed
calibration, so it isn't reported as a separate results axis.

A classical (undirected, symmetric, fatigue-free) OP baseline
(`solver_op_classical.cpp`, per Fischetti et al. 1998) and a heuristic
comparison (Greedy / GA / ACO vs. SA, `compare_heuristics.cpp`) are also
included — the latter with exact Wilcoxon signed-rank tests,
Holm-Bonferroni correction across the 3 pairwise comparisons, and a
rank-biserial effect size (`analyze_heuristics.py`).

## B&C Solver Features

Two LP formulations:
- **MTZ formulation**: state propagation (`G_i`) with tightened big-M constants, McCormick linearisation of the bilinear fatigue term
- **Flow formulation**: single-commodity flow variables (`g_ij`) with two-sided coupling bounds, natively linear fatigue budget (no McCormick)

Valid inequalities (both formulations):
- Directed subtour elimination (Kosaraju SCC)
- Connectivity cuts (reverse BFS + Dinic max-flow)
- Lifted cover cuts with fatigue-aware arc weights (B1)
- Routing infeasibility cuts based on assignment relaxation (B2)
- Directed cycle cover cuts (B3)
- Directed path inequality cuts (B4)

Other features:
- Fatigue-aware arc elimination
- SA warm-start with auto-calibrated iterations/restarts
- DFS node selection for terrain instances, BFS for benchmarks

## Data

Input terrain data (TIF maps, DEMs, OMAP files) are included under
`data/torremocha/` and `data/la_muela/`, tracked via Git LFS.

## Citation

See `paper/orienteering_paper.tex` for the full methodology and results.

## Author

Pablo Borrego Ramos — paborrego@alumnos.unex.es
