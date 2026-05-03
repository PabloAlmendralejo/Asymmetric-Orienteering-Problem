# Optimal Route Planning for Orienteering on Real Terrain

End-to-end framework for solving the Asymmetric Orienteering Problem with Fatigue (AOPF) on real terrain, combining a terrain-aware asymmetric cost model with exact Branch-and-Cut optimisation.

## Overview

The pipeline has two stages:

1. **Python** — Preprocesses orienteering map files (OMAP), digital elevation models (DEM),
   and map imagery into an asymmetric cost matrix using:
   - IOF symbol-based cost raster (ACR)
   - Hypsometry Cost Raster (HCR) with IQR scaling
   - Minetti metabolic slope model (directional)
   - Cumulative fatigue model

2. **C++** — Solves the Orienteering Problem exactly via Branch-and-Cut with HiGHS,
   warm-started by Simulated Annealing.

## Project Structure

```
├── python/
│   ├── run_pipeline.py              # Main entry point
│   ├── requirements.txt
│   ├── core/
│   │   ├── preprocessing.py         # Raw files → terrain cache (.npz)
│   │   ├── cost_functions.py        # Minetti, fatigue, ISOM weights
│   │   ├── terrain_analysis.py      # HCR (TRI, PC, SLO), IQR scaling
│   │   ├── omap_parser.py           # OMAP XML parsing, symbol classification
│   │   ├── rasterization.py         # Line/polygon rasterization
│   │   ├── coordinate_transforms.py # Pixel ↔ UTM ↔ grid transforms
│   │   ├── pathfinding.py           # Anisotropic Dijkstra, cost matrix
│   │   ├── route_optimizer.py       # Greedy + SA route optimization
│   │   └── control_placement.py     # Control point generation (7 strategies)
│   ├── visualization/
│   │   └── visualize_results.py     # Map, cost surfaces, routes, hillshade
│   └── config/
│       ├── torremocha.py
│       └── la_muela.py
├── cpp/
│   ├── op_bnc_with_highs.cpp             # B&C solver (generic / La Muela)
│   └── op_bnc_with_highs_torremocha.cpp  # B&C solver (Torremocha)
├── benchmark/
│   ├── generate_instances.py        # Parameterized instance generator
│   ├── benchmark_solver.cpp         # MTZ B&C solver for benchmarks
│   ├── benchmark_solver_flow.cpp    # Flow B&C solver for benchmarks
│   ├── benchmark_solver_ablation.cpp # Ablation study solver
│   ├── compare_heuristics.cpp       # Greedy / GA / ACO / SA comparison
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
- [HiGHS](https://highs.dev/) LP solver library

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
# Compile (Windows/MSVC):
cl.exe /O2 /std:c++20 /EHsc op_bnc_with_highs_torremocha.cpp ^
  /I<highs_include> /link /LIBPATH:<highs_lib> highs.lib ^
  /out:solver.exe

# Run (from directory with op_input_*.json files):
solver.exe
```

### Visualization
```bash
cd python/visualization
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
| Cost asymmetry | 16.1% | 50.0% |
| IQR scaling (c) | 0.87 | 2.18 |
| Proven optimal (flow) | 4/7 instances | 5/7 instances |

## Key Results

### Real terrain instances
Across 14 instances (7 per terrain), the flow-based B&C proves optimality on
9 of 14 instances (vs. 8 for MTZ) within a 15-minute time limit.

### Formulation comparison (21 benchmark instances)
The flow formulation reduces the mean LP gap from 16.0% to 5.4% compared to
the MTZ formulation, and proves 7 instances optimal versus 4 for MTZ.

### Heuristic comparison
Under equal 2-second time budgets, SA significantly outperforms Greedy (23.0% gap),
GA (49.2% gap), and ACO (8.1% gap) across 21 benchmark instances
(Wilcoxon signed-rank, all p < 0.001).

## Benchmark Suite

Parameterized instances for the asymmetric OP with fatigue — the first
benchmark set for this problem variant. See `benchmark/README.md` for details.

Instances vary across:
- Node count: 20, 30, 40, 50, 75, 100
- Asymmetry: 0% to 50%
- Fatigue rate: 0.0 to 0.3
- Budget tightness: loose (70%), medium (50%), tight (30%) of NN tour cost

## B&C Solver Features

Two LP formulations:
- **MTZ formulation**: MTZ time propagation with tightened big-M constants, McCormick linearisation of bilinear fatigue term
- **Flow formulation**: single-commodity flow variables, natively linear fatigue budget (no McCormick), tightened flow-arc coupling bounds

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
