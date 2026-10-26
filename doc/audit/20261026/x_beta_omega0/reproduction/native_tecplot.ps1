$ErrorActionPreference = 'Stop'
$targetDir = Join-Path (Get-Location) 'outputs\M3_omega0_xbeta_20261026'
$names = @('M3_omega0_xbeta_comparison','bluntplate_M3_omega0_xbeta','rough_M3_omega0_xbeta')
$preplot = 'E:\Program Files\Tecplot\Tecplot 360 EX 2024 R1\bin\preplot.exe'
$tec360 = 'E:\Program Files\Tecplot\Tecplot 360 EX 2024 R1\bin\tec360.exe'
Push-Location -LiteralPath $targetDir
try {
    foreach ($name in $names) {
        $result = & $preplot ($name+'.dat') ($name+'.plt') 2>&1
        $result | Set-Content -LiteralPath ($name+'_preplot.log') -Encoding utf8
        if ($LASTEXITCODE -ne 0) { throw "Preplot failed: $name" }
    }
    $process = Start-Process -FilePath $tec360 -ArgumentList @('-b','-p','tecplot_roundtrip.mcr') -WorkingDirectory $targetDir -WindowStyle Hidden -PassThru
    if (-not $process.WaitForExit(30000)) { throw "Tecplot still running; PID $($process.Id)" }
    if ($process.ExitCode -ne 0) { throw "Tecplot failed: $($process.ExitCode)" }
    Get-Content -LiteralPath 'batch.log' -Tail 12
} finally { Pop-Location }
