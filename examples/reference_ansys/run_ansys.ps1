param([string]$AnsysExe = $env:ANSYS_EXE, [int]$Nproc = 4)
$ErrorActionPreference = 'Stop'
Set-Location -LiteralPath $PSScriptRoot
if (-not $AnsysExe) {
  $base = Join-Path ${env:ProgramFiles} 'ANSYS Inc'
  if (-not (Test-Path $base)) { throw "ANSYS directory not found: $base. Set ANSYS_EXE." }
  $candidates = Get-ChildItem -LiteralPath $base -Directory | Where-Object { $_.Name -match '^v\d+' } | Sort-Object Name -Descending
  foreach ($v in $candidates) {
    $exes = Get-ChildItem -Path (Join-Path $v.FullName 'ANSYS\bin\winx64') -Filter 'ansys*.exe' -ErrorAction SilentlyContinue
    if ($exes) { $AnsysExe = ($exes | Sort-Object Name -Descending | Select-Object -First 1).FullName; break }
  }
}
if (-not $AnsysExe -or -not (Test-Path -LiteralPath $AnsysExe)) { throw "ANSYS executable unavailable: $AnsysExe" }
Write-Host "Running $AnsysExe"
$env:ANS_CONSEC = 'YES'
& $AnsysExe -b -np $Nproc -j PIZZA -i run.inp -o ansys.out
$exitCode = $LASTEXITCODE
if ($exitCode -ne 0) { throw "ANSYS process exited $exitCode. Inspect ansys.out" }
if (-not (Test-Path -LiteralPath history.csv)) { throw 'Missing history.csv - ANSYS did not finish post-processing.' }
if (-not (Test-Path -LiteralPath wall_samples.csv)) { throw 'Missing wall_samples.csv - ANSYS did not finish post-processing.' }
Write-Host 'ANSYS_RESULT_READY'
