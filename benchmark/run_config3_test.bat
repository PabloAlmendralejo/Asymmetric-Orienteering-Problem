@echo off
setlocal
set EXE=benchmark_solver_ablation.exe
set DIR=instances

echo === Config 3/5: Base + B1 + B2 (routing infeasibility) ===
echo Testing on bench_008 and bench_020 only
echo.

%EXE% %DIR% --config=3_plus_B1_B2 --coupling=on --fatigue-covers=on --covers=on --routing=on --cycle=off --path=off
