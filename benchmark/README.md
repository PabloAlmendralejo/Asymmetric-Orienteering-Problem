# Benchmark Suite

Parameterized asymmetric OP instances for evaluating solver performance.

## Generate instances

```bash
cd benchmark
python generate_instances.py --output-dir instances
```

This creates ~21 instances varying:
- Node count: 20, 30, 40, 50, 75, 100
- Asymmetry: 0%, 10%, 25%, 50%
- Fatigue rate: 0.0, 0.1, 0.2, 0.3
- Budget tightness: loose, medium, tight
- Plus 4 extreme cases

## Solve

```bash
# Compile
cl.exe /O2 /std:c++20 /EHsc benchmark_solver.cpp /I<highs_include> /link highs.lib /out:benchmark_solver.exe

# Run (reads all op_input_*.json from instances/)
benchmark_solver.exe instances
```

## Output

Each instance produces:
- `op_output_*.json` with SA and B&C results
- Summary table with points, nodes, time, optimality status, and LP gap
