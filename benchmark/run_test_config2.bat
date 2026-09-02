@echo off
setlocal
set EXE=benchmark_solver_ablation.exe
set CSV=test_config2_results.csv

echo Testing config 2 (no lower bound) on bench_008 and bench_020
if exist %CSV% del %CSV%

echo === bench_008 ===
%EXE% instances\op_input_bench_008_n40_a10_f20_med.json --config=2_plus_B1 --coupling=on --fatigue-covers=on --covers=on --routing=off --cycle=off --path=off --csv=%CSV%

echo === bench_020 ===
%EXE% instances\op_input_bench_020_high_asym_no_fatigue.json --config=2_plus_B1 --coupling=on --fatigue-covers=on --covers=on --routing=off --cycle=off --path=off --csv=%CSV%

echo Done. Results in %CSV%
