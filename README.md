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
│   ├── benchmark_solver.cpp         # Standalone B&C solver for benchmarks
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
| CRS | EPSG:25829 | EPSG:25830 |
| Grid | 1056 × 1463 | 1376 × 1622 |
| Resolution | 2.0 m | 2.0 m |
| Elevation | 447–501 m | 1164–1628 m |
| Relief | 55 m | 464 m |
| Cost asymmetry | 16.1% | 50.0% |
| HCR scaling (c) | 0.87 | 2.18 |
| Proven optimal | 5/7 instances | 5/7 instances |

## Key Results

### Real terrain instances
Across 14 instances (7 per terrain), the B&C incumbent matches the SA warm-start
on all instances, with optimality proven on 10 of 14. LP relaxation gaps on
unproven instances range from 8.6% to 19.5%.

### Heuristic comparison
Under equal 2-second time budgets, SA significantly outperforms Greedy (23% gap),
GA (49% gap), and ACO (8.1% gap) across 21 benchmark instances
(Wilcoxon signed-rank, all p < 0.001).

### Benchmark scaling
The B&C proves optimality for n ≤ 30 within 15 minutes. At n ≥ 40, gaps range
from 17–25% depending on asymmetry and fatigue rate. The McCormick envelope is
the primary source of LP relaxation looseness: gaps increase from 6.6% (λ = 0)
to 24.9% (λ = 0.3).

## Benchmark Suite

Parameterized instances for the asymmetric OP with fatigue — the first
benchmark set for this problem variant.

```bash
cd benchmark
python generate_instances.py --output-dir instances

# Compile and run B&C solver:
benchmark_solver.exe instances

# Heuristic comparison (no HiGHS needed):
compare_heuristics.exe instances
```

Instances vary across:
- Node count: 20, 30, 40, 50, 75, 100
- Asymmetry: 0% to 50%
- Fatigue rate: 0.0 to 0.3
- Budget tightness: loose, medium, tight

## B&C Solver Features

- MTZ time propagation with tightened big-M constants
- McCormick linearisation of bilinear fatigue term
- Directed subtour elimination (Kosaraju SCC)
- Depot-unreachable detection (reverse BFS)
- Lifted cover cuts on fatigue budget knapsack
- Fatigue-aware arc elimination
- LP gap tracking and explored/unexplored node counts
- Deterministic SA with 10x restarts for robust warm starts

## Data

Input terrain data (TIF maps, DEMs, OMAP files) are included under
`data/torremocha/` and `data/la_muela/`, tracked via Git LFS.

## Citation

See `paper/orienteering_paper.tex` for the full methodology and results.

## Author

Pablo Borrego Ramos — pabloalmendralejo@gmail.com
