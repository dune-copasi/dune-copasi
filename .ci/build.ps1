[CmdletBinding()]
param(
  [string]$RepoRoot,
  [string]$WorkDir,
  [string]$InstallPrefix,
  [switch]$CleanWorkDir
)

Set-StrictMode -Version Latest
$ErrorActionPreference = 'Stop'

function Resolve-AbsolutePath {
  param([Parameter(Mandatory = $true)][string]$Path)

  return $ExecutionContext.SessionState.Path.GetUnresolvedProviderPathFromPSPath($Path)
}

function Invoke-ExternalCommand {
  param(
    [Parameter(Mandatory = $true)][string]$FilePath,
    [string[]]$ArgumentList = @()
  )

  Write-Host "> $FilePath $($ArgumentList -join ' ')"
  & $FilePath @ArgumentList
  if ($LASTEXITCODE -ne 0) {
    throw "Command failed with exit code ${LASTEXITCODE}: $FilePath $($ArgumentList -join ' ')"
  }
}

function Invoke-ExternalCommandWithRetry {
  param(
    [Parameter(Mandatory = $true)][string]$FilePath,
    [string[]]$ArgumentList = @(),
    [int]$MaxAttempts = 4,
    [int]$InitialDelaySeconds = 5,
    [string]$CleanupPath
  )

  for ($attempt = 1; $attempt -le $MaxAttempts; $attempt++) {
    try {
      Invoke-ExternalCommand -FilePath $FilePath -ArgumentList $ArgumentList
      return
    } catch {
      if ($attempt -ge $MaxAttempts) {
        throw
      }

      if ($CleanupPath -and (Test-Path -LiteralPath $CleanupPath)) {
        Remove-Item -LiteralPath $CleanupPath -Recurse -Force
      }

      $delaySeconds = $InitialDelaySeconds * $attempt
      Write-Warning "Attempt $attempt of $MaxAttempts failed for '$FilePath'. Retrying in $delaySeconds seconds..."
      Start-Sleep -Seconds $delaySeconds
    }
  }
}

function Repair-GitSymlinks {
  param([Parameter(Mandatory = $true)][string]$ModuleDir)

  $entries = & git -C $ModuleDir ls-files -s
  if ($LASTEXITCODE -ne 0) {
    throw "Failed to inspect git symlinks in $ModuleDir"
  }

  foreach ($entry in $entries) {
    if ($entry -notmatch '^120000 ') {
      continue
    }

    $parts = $entry -split "`t", 2
    if ($parts.Count -ne 2) {
      continue
    }

    $linkPath = Join-Path $ModuleDir $parts[1]
    if (-not (Test-Path $linkPath)) {
      continue
    }

    $item = Get-Item -LiteralPath $linkPath -Force
    if ($item.Attributes -band [System.IO.FileAttributes]::ReparsePoint) {
      continue
    }

    $relativeTarget = (Get-Content -LiteralPath $linkPath -Raw).Trim()
    if (-not $relativeTarget) {
      continue
    }

    $sourcePath = [System.IO.Path]::GetFullPath((Join-Path $item.DirectoryName $relativeTarget))
    Write-Host "Repairing git symlink placeholder: $linkPath -> $sourcePath"
    Copy-Item -LiteralPath $sourcePath -Destination $linkPath -Force
  }
}

function Normalize-CMakeValue {
  param([Parameter(Mandatory = $true)][string]$Value)

  if ($Value.Length -ge 2) {
    $first = $Value[0]
    $last = $Value[$Value.Length - 1]
    if (($first -eq "'" -and $last -eq "'") -or ($first -eq '"' -and $last -eq '"')) {
      return $Value.Substring(1, $Value.Length - 2)
    }
  }

  return $Value
}

function Get-CMakeEnvironmentArguments {
  $args = @()
  $excludedNames = [System.Collections.Generic.HashSet[string]]::new([System.StringComparer]::OrdinalIgnoreCase)
  foreach ($name in @('CMAKE_GENERATOR', 'CMAKE_INSTALL_PREFIX', 'CMAKE_PREFIX_PATH')) {
    $null = $excludedNames.Add($name)
  }

  $environment = [System.Environment]::GetEnvironmentVariables()
  foreach ($entry in $environment.GetEnumerator() | Sort-Object Key) {
    $name = [string]$entry.Key
    if ($excludedNames.Contains($name)) {
      continue
    }

    $isCMakeOrDune = $name.StartsWith('CMAKE_') -or $name.StartsWith('DUNE_')
    if (-not $isCMakeOrDune -and $name -notin @('BUILD_SHARED_LIBS', 'BUILD_TESTING')) {
      continue
    }

    $value = Normalize-CMakeValue ([string]$entry.Value)
    if ($value) {
      $args += "-D$name=$value"
    }
  }

  return $args
}

if ($IsWindows) {
  if (-not $env:CMAKE_MSVC_RUNTIME_LIBRARY) {
    $env:CMAKE_MSVC_RUNTIME_LIBRARY = 'MultiThreaded$<$<CONFIG:Debug>:Debug>'
  }
  if ($env:CMAKE_MSVC_RUNTIME_LIBRARY -and -not $env:CMAKE_POLICY_DEFAULT_CMP0091) {
    $env:CMAKE_POLICY_DEFAULT_CMP0091 = 'NEW'
  }
}

function Build-CMakeProject {
  param(
    [Parameter(Mandatory = $true)][string]$SourceDir,
    [Parameter(Mandatory = $true)][string]$BuildDir,
    [Parameter(Mandatory = $true)][string]$InstallPrefix,
    [Parameter(Mandatory = $true)][string]$PrefixPath
  )

  $configureArgs = @(
    '-S'
    $SourceDir
    '-B'
    $BuildDir
    "-DCMAKE_INSTALL_PREFIX=$InstallPrefix"
    "-DCMAKE_PREFIX_PATH=$PrefixPath"
  )
  if ($env:CMAKE_GENERATOR) {
    $configureArgs = @('-G', $env:CMAKE_GENERATOR) + $configureArgs
  }
  $configureArgs += Get-CMakeEnvironmentArguments
  Invoke-ExternalCommand cmake $configureArgs

  Invoke-ExternalCommand cmake @('--build', $BuildDir)
  Invoke-ExternalCommand cmake @('--install', $BuildDir)
}

$scriptDir = Split-Path -Parent $PSCommandPath

if (-not $RepoRoot) {
  $RepoRoot = Join-Path $scriptDir '..'
}
$RepoRoot = Resolve-AbsolutePath $RepoRoot

if (-not $InstallPrefix) {
  if ($env:CMAKE_INSTALL_PREFIX) {
    $InstallPrefix = $env:CMAKE_INSTALL_PREFIX
  } elseif ($env:INSTALL_PREFIX) {
    $InstallPrefix = $env:INSTALL_PREFIX
  } else {
    $InstallPrefix = Join-Path $RepoRoot '.local/dune'
  }
}
$InstallPrefix = Resolve-AbsolutePath $InstallPrefix

if (-not $WorkDir) {
  $workRoot = if ($env:RUNNER_TEMP) { $env:RUNNER_TEMP } elseif ($env:TEMP) { $env:TEMP } else { Join-Path $RepoRoot '.tmp' }
  $WorkDir = Join-Path $workRoot 'dune-dependencies'
}
$WorkDir = Resolve-AbsolutePath $WorkDir

if (-not (Test-Path $RepoRoot)) {
  throw "Repository root does not exist: $RepoRoot"
}

if ($CleanWorkDir -and (Test-Path $WorkDir)) {
  Write-Host "Removing existing dependency workspace: $WorkDir"
  Remove-Item -Recurse -Force $WorkDir
}

New-Item -ItemType Directory -Force -Path $InstallPrefix | Out-Null
New-Item -ItemType Directory -Force -Path $WorkDir | Out-Null
New-Item -ItemType Directory -Force -Path (Join-Path $WorkDir 'build') | Out-Null

$env:CMAKE_INSTALL_PREFIX = $InstallPrefix
if (-not $env:CMAKE_PREFIX_PATH) {
  $env:CMAKE_PREFIX_PATH = $InstallPrefix
} elseif ($env:CMAKE_PREFIX_PATH -notlike "*$InstallPrefix*") {
  $env:CMAKE_PREFIX_PATH = "$InstallPrefix;$($env:CMAKE_PREFIX_PATH)"
}

$modules = @(
  [pscustomobject]@{ Name = 'dune-common'; Repo = 'https://gitlab.dune-project.org/core/dune-common.git'; Branch = 'master' }
  [pscustomobject]@{ Name = 'dune-geometry'; Repo = 'https://gitlab.dune-project.org/core/dune-geometry.git'; Branch = 'master' }
  [pscustomobject]@{ Name = 'dune-uggrid'; Repo = 'https://gitlab.dune-project.org/staging/dune-uggrid.git'; Branch = 'master' }
  [pscustomobject]@{ Name = 'dune-grid'; Repo = 'https://gitlab.dune-project.org/core/dune-grid.git'; Branch = 'master' }
  [pscustomobject]@{ Name = 'dune-typetree'; Repo = 'https://gitlab.dune-project.org/staging/dune-typetree.git'; Branch = 'master' }
  [pscustomobject]@{ Name = 'dune-localfunctions'; Repo = 'https://gitlab.dune-project.org/core/dune-localfunctions.git'; Branch = 'master' }
  [pscustomobject]@{ Name = 'dune-istl'; Repo = 'https://gitlab.dune-project.org/core/dune-istl.git'; Branch = 'master' }
  [pscustomobject]@{ Name = 'dune-functions'; Repo = 'https://gitlab.dune-project.org/staging/dune-functions.git'; Branch = 'master' }
  [pscustomobject]@{ Name = 'dune-multidomaingrid'; Repo = 'https://gitlab.dune-project.org/liam.keegan/dune-multidomaingrid.git'; Branch = 'msvc' }
  [pscustomobject]@{ Name = 'dune-pdelab'; Repo = 'https://gitlab.dune-project.org/pdelab/dune-pdelab.git'; Branch = 'msvc' }
)

Write-Host "Repository root:  $RepoRoot"
Write-Host "Install prefix:   $InstallPrefix"
Write-Host "Dependency root:  $WorkDir"
Write-Host "MSVC runtime:     $env:CMAKE_MSVC_RUNTIME_LIBRARY"
Write-Host "Modules:"
$modules | ForEach-Object { Write-Host "  - $($_.Name) [$($_.Branch)]" }

foreach ($module in $modules) {
  $moduleDir = Join-Path $WorkDir $module.Name
  $buildDir = Join-Path (Join-Path $WorkDir 'build') $module.Name

  if (-not (Test-Path $moduleDir)) {
    Invoke-ExternalCommandWithRetry -FilePath git -CleanupPath $moduleDir -ArgumentList @(
      'clone'
      '--branch'
      $module.Branch
      '--depth'
      '1'
      '--single-branch'
      $module.Repo
      $moduleDir
    )
  } else {
    Write-Host "Reusing existing checkout: $moduleDir"
  }

  if ($IsWindows) {
    Repair-GitSymlinks $moduleDir
  }

  Build-CMakeProject -SourceDir $moduleDir -BuildDir $buildDir -InstallPrefix $InstallPrefix -PrefixPath $env:CMAKE_PREFIX_PATH
}

$duneCopasiBuildDir = Join-Path (Join-Path $WorkDir 'build') 'dune-copasi'
Write-Host "Building dune-copasi from $RepoRoot"
Build-CMakeProject -SourceDir $RepoRoot -BuildDir $duneCopasiBuildDir -InstallPrefix $InstallPrefix -PrefixPath $env:CMAKE_PREFIX_PATH

$buildTesting = $env:BUILD_TESTING
if (-not $buildTesting -or $buildTesting -notmatch '^(0|OFF|FALSE|NO)$') {
  Write-Host "Building dune-copasi test targets from $duneCopasiBuildDir"
  $buildTestArgs = @('--build', $duneCopasiBuildDir, '--target', 'build_unit_tests', 'build_system_tests', 'build_docs_tests')
  Invoke-ExternalCommand cmake $buildTestArgs

  Write-Host "Running dune-copasi tests from $duneCopasiBuildDir"
  $ctestArgs = @('--test-dir', $duneCopasiBuildDir, '--output-on-failure')
  Invoke-ExternalCommand ctest $ctestArgs
}
