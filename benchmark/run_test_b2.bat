@echo off
setlocal

set EXE=benchmark_solver_ablation.exe
set TMPDIR=test_instances_b2

if exist %TMPDIR% rmdir /s /q %TMPDIR%
mkdir %TMPDIR%
copy instances\op_input_bench_008_n40_a10_f20_med.json %TMPDIR%\
copy instances\op_input_bench_020_high_asym_no_fatigue.json %TMPDIR%\

echo === Test: coupling=on, lower-bound=OFF, fatigue-covers=on, covers=on ===
%EXE% %TMPDIR% --config=no_lb --coupling=on --lower-bound=off --fatigue-covers=on --covers=on --routing=off --cycle=off --path=off > test_b2_log.txt 2>&1
echo Done. Results in test_b2_log.txt

rmdir /s /q %TMPDIR%
