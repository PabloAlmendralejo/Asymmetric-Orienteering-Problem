@echo off
setlocal
set EXE=benchmark_solver_ablation.exe
set DIR=instances
set CSV=ablation_results.csv

echo Ablation study: 6 configurations x 21 instances
echo Results will be written to %CSV%
echo.

if exist %CSV% del %CSV%

echo === Config 1/6: Base (SECs + connectivity + tight coupling + fatigue-aware covers, no B1-B4) ===
%EXE% %DIR% --config=1_base --coupling=on --fatigue-covers=on --covers=off --routing=off --cycle=off --path=off --csv=%CSV%

echo === Config 2/6: Base + B1 only (lifted covers) ===
%EXE% %DIR% --config=2_B1_only --coupling=on --fatigue-covers=on --covers=on --routing=off --cycle=off --path=off --csv=%CSV%

echo === Config 3/6: Base + B2 only (routing infeasibility) ===
%EXE% %DIR% --config=3_B2_only --coupling=on --fatigue-covers=on --covers=off --routing=on --cycle=off --path=off --csv=%CSV%

echo === Config 4/6: Base + B3 only (cycle covers) ===
%EXE% %DIR% --config=4_B3_only --coupling=on --fatigue-covers=on --covers=off --routing=off --cycle=on --path=off --csv=%CSV%

echo === Config 5/6: Base + B4 only (path inequalities) ===
%EXE% %DIR% --config=5_B4_only --coupling=on --fatigue-covers=on --covers=off --routing=off --cycle=off --path=on --csv=%CSV%

echo === Config 6/6: Base + all (B1 + B2 + B3 + B4) ===
%EXE% %DIR% --config=6_all --coupling=on --fatigue-covers=on --covers=on --routing=on --cycle=on --path=on --csv=%CSV%

echo.
echo ============================================
echo Ablation study complete. Results in %CSV%
echo ============================================
