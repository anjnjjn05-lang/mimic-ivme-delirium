$ErrorActionPreference = 'Stop'

$packageRoot = (Resolve-Path -LiteralPath $PSScriptRoot).Path
$codeRoot = Join-Path $packageRoot 'R'

if (-not $env:PGPASSWORD) { throw 'Set PGPASSWORD in the current PowerShell session.' }
if (-not $env:MIMIC_DB_HOST) { $env:MIMIC_DB_HOST = 'localhost' }
if (-not $env:MIMIC_DB_PORT) { $env:MIMIC_DB_PORT = '5432' }
if (-not $env:MIMIC_DB_NAME) { $env:MIMIC_DB_NAME = 'mimiciv' }
if (-not $env:MIMIC_DB_USER) { $env:MIMIC_DB_USER = 'postgres' }

$env:ANALYSIS_ROOT = $packageRoot
$env:ANALYSIS_RAW = Join-Path $packageRoot 'outputs\restricted_patient_level'
$env:ANALYSIS_MAIN = Join-Path $packageRoot 'outputs\aggregate_review\main'
$env:ANALYSIS_SUPP = Join-Path $packageRoot 'outputs\aggregate_review\supplementary'
$env:ANALYSIS_ENV = Join-Path $packageRoot 'outputs\aggregate_review\environment'
$env:ANALYSIS_VALIDATION = Join-Path $packageRoot 'outputs\aggregate_review\validation'
$env:ANALYSIS_INTERNAL_FIG = Join-Path $packageRoot 'outputs\aggregate_review\internal'
$env:ANALYSIS_REPORT = Join-Path $packageRoot 'outputs\aggregate_review\reports'
$env:EXPECTED_PROCEDURE_ICU_PAIRS = '9662'
$env:EXPECTED_FINAL_N = '1939'
$env:EXPECTED_EVENTS = '322'

New-Item -ItemType Directory -Force -Path @(
    $env:ANALYSIS_RAW, $env:ANALYSIS_MAIN, $env:ANALYSIS_SUPP,
    $env:ANALYSIS_ENV, $env:ANALYSIS_VALIDATION,
    $env:ANALYSIS_INTERNAL_FIG, $env:ANALYSIS_REPORT
) | Out-Null

$rscript = $null
if ($env:R_SCRIPT) {
    if (-not (Test-Path -LiteralPath $env:R_SCRIPT -PathType Leaf)) {
        throw "R_SCRIPT does not point to a file: $env:R_SCRIPT"
    }
    $rscript = (Resolve-Path -LiteralPath $env:R_SCRIPT).Path
}
if (-not $rscript) {
    $rscriptCommand = Get-Command Rscript -ErrorAction SilentlyContinue
    if ($rscriptCommand) { $rscript = $rscriptCommand.Source }
}
if (-not $rscript) { throw 'Rscript was not found. Add it to PATH or set R_SCRIPT.' }

$steps = @(
    '00_check_environment.R', '01_build_analysis_dataset.R',
    '02_extract_covariates.R', '03_overlap_weighting_bootstrap.R',
    '04_spline_and_aipw.R', '05_binary_sensitivity.R',
    '06_threshold_analysis.R', '07_figures_tables.R',
    '08_validate_dose_extraction.R'
)

Push-Location $codeRoot
try {
    foreach ($step in $steps) {
        Write-Host "Running $step"
        & $rscript $step
        if ($LASTEXITCODE -ne 0) { throw "$step failed with exit code $LASTEXITCODE" }
    }
}
finally { Pop-Location }

Write-Host 'Pipeline complete. All generated files remain below the ignored outputs directory.'
