$ErrorActionPreference = 'Stop'
$taskDest = (Resolve-Path -LiteralPath (Join-Path $PSScriptRoot '..')).Path
$taskNames = @('M3_xbeta_f-100kHz_comparison', 'M3_xbeta_f-50kHz_comparison', 'M3_xbeta_f+50kHz_comparison', 'M3_xbeta_f+100kHz_comparison', 'M3_xbeta_four_frequencies', 'bluntplate_M3_xbeta_four_frequencies', 'rough_M3_xbeta_four_frequencies')
Push-Location -LiteralPath $taskDest
try {
    foreach ($taskName in $taskNames) {
        & 'E:\Program Files\Tecplot\Tecplot 360 EX 2024 R1\bin\preplot.exe' ($taskName + '.dat') ($taskName + '.plt') 2>&1 | Out-File -LiteralPath ($taskName + '_preplot.log') -Encoding utf8
        if ($LASTEXITCODE -ne 0) { throw ('Preplot failed: ' + $taskName) }
    }
    $taskProc = Start-Process -FilePath 'E:\Program Files\Tecplot\Tecplot 360 EX 2024 R1\bin\tec360.exe' -ArgumentList '-b', '-p', 'tecplot_roundtrip.mcr' -WorkingDirectory $taskDest -WindowStyle Hidden -PassThru
    $taskProc.WaitForExit()
    if ($taskProc.ExitCode -ne 0) { throw ('Tecplot exit: ' + $taskProc.ExitCode) }
    if (Select-String -LiteralPath 'batch.log' -Pattern 'ERROR' -Quiet) { throw 'See batch.log' }
} finally { Pop-Location }
