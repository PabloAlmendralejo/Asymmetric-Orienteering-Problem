@echo off
setlocal

set HIGHS_INC=/I"C:\Users\borrepa\Downloadsnz_dem_jsons\HiGHS\highs" /I"C:\Users\borrepa\Downloadsnz_dem_jsons\HiGHS\build"
set HIGHS_LIB=/LIBPATH:"C:\Users\borrepa\Downloadsnz_dem_jsons\HiGHS\build\Release\bin" highs.lib
set TMPDIR=test_instances_b2

if exist %TMPDIR% rmdir /s /q %TMPDIR%
mkdir %TMPDIR%
copy instances\op_input_bench_008_n40_a10_f20_med.json %TMPDIR%\
copy instances\op_input_bench_020_high_asym_no_fatigue.json %TMPDIR%\

echo === Test: flow solver, no lower bound, no B2 ===
cl /O2 /std:c++20 /EHsc test_flow_no_lb.cpp %HIGHS_INC% /link %HIGHS_LIB% /out:test_flow_no_lb.exe
if errorlevel 1 (echo COMPILE FAILED & goto :end)
test_flow_no_lb.exe %TMPDIR% > test_flow_no_lb_log.txt 2>&1
echo Done. Results in test_flow_no_lb_log.txt

:end
rmdir /s /q %TMPDIR%
