# watch_results.ps1 - Monitor test_aniso_results.txt for changes
$file = "D:\FEMM Source\femmTestFiles\test_aniso_results.txt"
$logfile = "D:\FEMM Source\TestBin\femm_debug.log"

Write-Host "Watching for results..." -ForegroundColor Cyan
Write-Host "Press Ctrl+C to stop."
Write-Host ""

$lastResultMod = [datetime]::MinValue
$lastLogMod    = [datetime]::MinValue

while ($true) {
    # Check test results file
    if (Test-Path $file) {
        $mod = (Get-Item $file).LastWriteTime
        if ($mod -gt $lastResultMod) {
            $lastResultMod = $mod
            Write-Host ("=" * 60) -ForegroundColor Green
            Write-Host "TEST RESULTS  [$mod]" -ForegroundColor Green
            Write-Host ("=" * 60) -ForegroundColor Green
            Get-Content $file | Write-Host
            Write-Host ""
        }
    }

    # Check debug log
    if (Test-Path $logfile) {
        $mod2 = (Get-Item $logfile).LastWriteTime
        if ($mod2 -gt $lastLogMod) {
            $lastLogMod = $mod2
            Write-Host ("=" * 60) -ForegroundColor Yellow
            Write-Host "FEMM DEBUG LOG  [$mod2]" -ForegroundColor Yellow
            Write-Host ("=" * 60) -ForegroundColor Yellow
            Get-Content $logfile | Select-Object -Last 30 | Write-Host
            Write-Host ""
        }
    }

    Start-Sleep -Milliseconds 500
}
