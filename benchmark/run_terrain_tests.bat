@echo off
setlocal

set HIGHS_INC=/I"C:\Users\borrepa\Downloadsnz_dem_jsons\HiGHS\highs" /I"C:\Users\borrepa\Downloadsnz_dem_jsons\HiGHS\build"
set HIGHS_LIB=/LIBPATH:"C:\Users\borrepa\Downloadsnz_dem_jsons\HiGHS\build\Release\bin" highs.lib
set EXE=benchmark_solver_ablation.exe
set TERRAIN=..\..\Orienteering C++

echo Compiling ablation solver...
cl /O2 /std:c++20 /EHsc benchmark_solver_ablation.cpp %HIGHS_INC% /link %HIGHS_LIB% /out:%EXE%
if errorlevel 1 (echo COMPILE FAILED & goto :end)

REM === La Muela (non-torremocha files) ===
set TMPDIR=terrain_la_muela
if exist %TMPDIR% rmdir /s /q %TMPDIR%
mkdir %TMPDIR%
copy "%TERRAIN%\op_input_clustered.json" %TMPDIR%\
copy "%TERRAIN%\op_input_elev_biased.json" %TMPDIR%\
copy "%TERRAIN%\op_input_mixed_density.json" %TMPDIR%\
copy "%TERRAIN%\op_input_path_biased.json" %TMPDIR%\
copy "%TERRAIN%\op_input_ring.json" %TMPDIR%\
copy "%TERRAIN%\op_input_sparse_far.json" %TMPDIR%\
copy "%TERRAIN%\op_input_standard.json" %TMPDIR%\

echo.
echo === Flow BnC on La Muela (7 instances) ===
%EXE% %TMPDIR% --config=flow_la_muela --coupling=on --fatigue-covers=on --covers=on --routing=on --cycle=off --path=off > terrain_flow_la_muela_log.txt 2>&1
echo Done. Results in terrain_flow_la_muela_log.txt
rmdir /s /q %TMPDIR%

REM === Torremocha ===
set TMPDIR=terrain_torremocha
if exist %TMPDIR% rmdir /s /q %TMPDIR%
mkdir %TMPDIR%
copy "%TERRAIN%\op_input_torremocha_clustered.json" %TMPDIR%\
copy "%TERRAIN%\op_input_torremocha_elev_biased.json" %TMPDIR%\
copy "%TERRAIN%\op_input_torremocha_mixed_density.json" %TMPDIR%\
copy "%TERRAIN%\op_input_torremocha_path_biased.json" %TMPDIR%\
copy "%TERRAIN%\op_input_torremocha_ring.json" %TMPDIR%\
copy "%TERRAIN%\op_input_torremocha_sparse_far.json" %TMPDIR%\
copy "%TERRAIN%\op_input_torremocha_standard.json" %TMPDIR%\

echo.
echo === Flow BnC on Torremocha (7 instances) ===
%EXE% %TMPDIR% --config=flow_torremocha --coupling=on --fatigue-covers=on --covers=on --routing=on --cycle=off --path=off > terrain_flow_torremocha_log.txt 2>&1
echo Done. Results in terrain_flow_torremocha_log.txt
rmdir /s /q %TMPDIR%

echo.
echo All terrain tests done.

:end
