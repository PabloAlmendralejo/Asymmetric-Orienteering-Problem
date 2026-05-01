@echo off
setlocal

set HIGHS_INC=/I"C:\Users\borrepa\Downloadsnz_dem_jsons\HiGHS\highs" /I"C:\Users\borrepa\Downloadsnz_dem_jsons\HiGHS\build"
set HIGHS_LIB=/LIBPATH:"C:\Users\borrepa\Downloadsnz_dem_jsons\HiGHS\build\Release\bin" highs.lib
set EXE=benchmark_solver_ablation.exe
set DIR=instances
set CSV=ablation_results.csv

echo Compiling ablation solver...
cl /O2 /std:c++20 /EHsc benchmark_solver_ablation.cpp %HIGHS_INC% /link %HIGHS_LIB% /out:%EXE%
if errorlevel 1 (echo COMPILE FAILED & goto :end)

echo Ablation study (cumulative): 5 configurations x 21 instances
echo Results will be written to %CSV%
echo.

if exist %CSV% del %CSV%

echo === Config 1/5: Base (SECs + connectivity + tight coupling + fatigue-aware covers, no B1-B4) ===
%EXE% %DIR% --config=1_base --coupling=on --fatigue-covers=on --covers=off --routing=off --cycle=off --path=off --csv=%CSV%

echo === Config 2/5: Base + B1 (lifted covers) ===
%EXE% %DIR% --config=2_plus_B1 --coupling=on --fatigue-covers=on --covers=on --routing=off --cycle=off --path=off --csv=%CSV%

echo === Config 3/5: Base + B1 + B2 (routing infeasibility) ===
%EXE% %DIR% --config=3_plus_B1_B2 --coupling=on --fatigue-covers=on --covers=on --routing=on --cycle=off --path=off --csv=%CSV%

echo === Config 4/5: Base + B1 + B2 + B3 (cycle covers) ===
%EXE% %DIR% --config=4_plus_B1_B2_B3 --coupling=on --fatigue-covers=on --covers=on --routing=on --cycle=on --path=off --csv=%CSV%

echo === Config 5/5: Base + B1 + B2 + B3 + B4 (path inequalities) ===
%EXE% %DIR% --config=5_all --coupling=on --fatigue-covers=on --covers=on --routing=on --cycle=on --path=on --csv=%CSV%

echo.
echo ============================================
echo Ablation study complete. Results in %CSV%
echo ============================================

:end
