# Golden-run driver: copy a legacy package to build/, overlay patches, run one entry
# script under matlab -batch with a diary, and capture everything the run wrote into
# tests/golden/<slug>/<entry>_<date>/.  See README.md in this folder for the protocol.
param(
    [Parameter(Mandatory = $true)][string]$Slug,
    [Parameter(Mandatory = $true)][string]$Entry,
    [int]$TimeoutMinutes = 720,
    # Build copies live OUTSIDE the repo: the repo sits in Dropbox, whose sync client
    # locks freshly written files (breaks the fresh-copy step) and would pointlessly
    # sync large MCMC scratch output.
    [string]$BuildRoot = (Join-Path $env:LOCALAPPDATA 'bvar-toolkit\golden_runs')
)
$ErrorActionPreference = 'Stop'
$repo    = (Resolve-Path (Join-Path $PSScriptRoot '..\..')).Path
$legacy  = Join-Path $repo "replications\$Slug\legacy"
$buildTop = Join-Path $BuildRoot $Slug
$patches = Join-Path $repo "tests\golden_runs\patches\$Slug"
if (-not (Test-Path $legacy)) { throw "No legacy folder for slug '$Slug'" }

# 1. Fresh build copy
if (Test-Path $buildTop) { Remove-Item -Recurse -Force $buildTop }
New-Item -ItemType Directory -Force $buildTop | Out-Null
Copy-Item -Recurse -Force "$legacy\*" $buildTop

# 2. Overlay patches (patched files carry a header comment stating what changed and why)
if (Test-Path $patches) {
    Copy-Item -Recurse -Force "$patches\*" $buildTop
    Write-Host "Applied patches from $patches"
}

# 3. Locate the entry script (packages may nest one level, e.g. sp_code\sp_code\)
$entryFile = Get-ChildItem -Recurse -File -Path $buildTop -Filter $Entry | Select-Object -First 1
if (-not $entryFile) { throw "Entry script '$Entry' not found under $buildTop" }
$runDir = $entryFile.DirectoryName
$entryName = [System.IO.Path]::GetFileNameWithoutExtension($Entry)

# 4. Snapshot, then run MATLAB headless, capturing stdout/stderr.
# NOTE: Start-Process does NOT quote ArgumentList elements (PS 5.1) - quote the -batch
# command manually. The command must therefore contain no double quotes itself.
$before = @{}
Get-ChildItem -Recurse -File $buildTop | ForEach-Object { $before[$_.FullName] = $_.LastWriteTimeUtc }
$stamp = Get-Date -Format 'yyyyMMdd_HHmm'
$mcmd = "disp(['MATLAB ' version]); disp(char(datetime)); tic; try; $entryName; catch err; disp('=== RUN ERRORED ==='); disp(getReport(err)); end; toc"
if ($mcmd -match '"') { throw 'Internal: -batch command must not contain double quotes' }
$outLog = Join-Path $runDir "golden_log_$entryName.txt"
$errLog = Join-Path $runDir "golden_stderr_$entryName.txt"
Write-Host "Running $Entry in $runDir (timeout $TimeoutMinutes min)..."
$p = Start-Process -FilePath 'matlab' -ArgumentList @('-batch', ('"' + $mcmd + '"')) `
    -WorkingDirectory $runDir -NoNewWindow -PassThru `
    -RedirectStandardOutput $outLog -RedirectStandardError $errLog
if (-not $p.WaitForExit($TimeoutMinutes * 60 * 1000)) {
    $p.Kill(); Write-Warning "Timed out after $TimeoutMinutes minutes; partial log kept."
}
Write-Host "MATLAB exit code: $($p.ExitCode)"

# 5. Capture new/modified files as goldens
$dest = Join-Path $repo "tests\golden\$Slug\${entryName}_$stamp"
New-Item -ItemType Directory -Force $dest | Out-Null
Get-ChildItem -Recurse -File $buildTop | Where-Object {
    (-not $before.ContainsKey($_.FullName)) -or ($_.LastWriteTimeUtc -gt $before[$_.FullName])
} | ForEach-Object {
    $rel = $_.FullName.Substring($buildTop.Length + 1) -replace '\\', '__'
    Copy-Item $_.FullName (Join-Path $dest $rel)
}
Write-Host "Golden capture: $dest"
Get-ChildItem $dest | Select-Object Name, Length | Format-Table -AutoSize
