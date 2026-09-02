@echo off
setlocal

set HIGHS_INC=/I"C:\Users\borrepa\Downloadsnz_dem_jsons\HiGHS\highs" /I"C:\Users\borrepa\Downloadsnz_dem_jsons\HiGHS\build"
set HIGHS_LIB=/LIBPATH:"C:\Users\borrepa\Downloadsnz_dem_jsons\HiGHS\build\Release\bin" highs.lib
set TMPDIR=test_instances_b2

if exist %TMPDIR% rmdir /s /q %TMPDIR%
mkdir %TMPDIR%
copy instances\op_input_bench_008_n40_a10_f20_med.json %TMPDIR%\
copy instances\op_input_bench_020_high_asym_no_fatigue.json %TMPDIR%\

echo === Compiling original ablation v0 (from v2 + minimal changes) ===
cl /O2 /std:c++20 /EHsc benchmark_solver_ablation_v0.cpp %HIGHS_INC% /link %HIGHS_LIB% /out:ablation_v0.exe
if errorlevel 1 (echo COMPILE FAILED & goto :end)

echo === Running config 3 equivalent: coupling=on fatigue-covers=on routing=on ===
ablation_v0.exe %TMPDIR% --config=v0_config3 --coupling=on --fatigue-covers=on --routing=on > test_v0_log.txt 2>&1
echo Done. Results in test_v0_log.txt

:end
rmdir /s /q %TMPDIR%
