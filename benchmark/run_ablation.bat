@echo off
setlocal
set EXE=benchmark_solver_ablation.exe
set DIR=instances
set CSV=ablation_results.csv

echo Ablation study: 6 configurations x 21 instances
echo Results will be written to %CSV%
echo.

if exist %CSV% del %CSV%

echo === Config 1/6: Base flow B^&C (no new cuts) ===
%EXE% %DIR% --config=1_base --coupling=off --fatigue-covers=off --routing=off --cycle=off --path=off --csv=%CSV%

echo === Config 2/6: + tightened coupling ===
%EXE% %DIR% --config=2_coupling --coupling=on --fatigue-covers=off --routing=off --cycle=off --path=off --csv=%CSV%

echo === Config 3/6: + fatigue-aware covers ===
%EXE% %DIR% --config=3_fat_covers --coupling=on --fatigue-covers=on --routing=off --cycle=off --path=off --csv=%CSV%

echo === Config 4/6: + routing infeasibility ===
%EXE% %DIR% --config=4_routing --coupling=on --fatigue-covers=on --routing=on --cycle=off --path=off --csv=%CSV%

echo === Config 5/6: + cycle covers + path inequalities ===
%EXE% %DIR% --config=5_cycle_path --coupling=on --fatigue-covers=on --routing=on --cycle=on --path=on --csv=%CSV%

echo === Config 6/6: All cuts (same as 5, for verification) ===
%EXE% %DIR% --config=6_all --coupling=on --fatigue-covers=on --routing=on --cycle=on --path=on --csv=%CSV%

echo.
echo ============================================
echo Ablation study complete. Results in %CSV%
echo ============================================
