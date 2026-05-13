$files = @(
    "D:\FEMM Source\femmTestFiles\probe_perp_lenz_results.txt",
    "D:\FEMM Source\femmTestFiles\probe_perp_lenz_data.csv",
    "D:\FEMM Source\femmTestFiles\probe_perp_lenz_idx.txt",
    "D:\FEMM Source\femmTestFiles\probe_pl_tmp.fem",
    "D:\FEMM Source\femmTestFiles\probe_pl_tmp.ans"
)
foreach ($f in $files) { if (Test-Path $f) { Remove-Item $f } }

$femm = "D:\FEMM Source\TestBin\femm.exe"
$script = "D:\FEMM Source\femmTestFiles\probe_perp_lenz.lua"

Write-Host "Run 1 (1 kHz)..."
Start-Process -FilePath $femm -ArgumentList "-lua-script=`"$script`"" -Wait -NoNewWindow
Write-Host "Run 2 (10 kHz)..."
Start-Process -FilePath $femm -ArgumentList "-lua-script=`"$script`"" -Wait -NoNewWindow
Write-Host "Run 3 (100 kHz)..."
Start-Process -FilePath $femm -ArgumentList "-lua-script=`"$script`"" -Wait -NoNewWindow

Write-Host ""
Write-Host "=== RESULTS ==="
Get-Content "D:\FEMM Source\femmTestFiles\probe_perp_lenz_results.txt"
