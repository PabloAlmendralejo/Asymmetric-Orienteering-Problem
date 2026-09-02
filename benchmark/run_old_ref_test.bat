@echo off
setlocal

set TMPDIR=test_instances_b2

if exist %TMPDIR% rmdir /s /q %TMPDIR%
mkdir %TMPDIR%
copy instances\op_input_bench_008_n40_a10_f20_med.json %TMPDIR%\
copy instances\op_input_bench_020_high_asym_no_fatigue.json %TMPDIR%\

echo Compiling old_ablation_ref.cpp...
cl /O2 /std:c++20 /EHsc old_ablation_ref.cpp /I. /link highs.lib /out:old_ablation_ref.exe
if errorlevel 1 (
    echo COMPILE FAILED
    goto :end
)

echo === Running old ref binary: coupling=on, fatigue-covers=on, routing=off, cycle=off, path=off ===
old_ablation_ref.exe %TMPDIR% --config=old_config3 --coupling=on --fatigue-covers=on --routing=off --cycle=off --path=off > old_ref_log.txt 2>&1
echo Done. Results in old_ref_log.txt

:end
rmdir /s /q %TMPDIR%
