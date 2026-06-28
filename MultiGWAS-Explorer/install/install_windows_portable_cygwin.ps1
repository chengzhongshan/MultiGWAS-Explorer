[CmdletBinding()]
param(
    [string]$PortableRoot = 'H:\TMP4SAS\CygwinPortablePipeline',
    [string]$RepoRoot = (Resolve-Path (Join-Path $PSScriptRoot '..')).Path,
    [string]$PortableReleaseUrl = 'https://github.com/MachinaCore/CygwinPortable/releases/download/1.4.0.0/CygwinPortable_1.4.0.0_WithDefaultCygwin.7z',
    [string]$CygwinSetupUrl = 'https://cygwin.com/setup-x86_64.exe',
    [string]$PortableArchivePath = '',
    [string]$Phase2Script = 'install/install_cygwin.sh',
    [string]$CygwinPackages = '',
    [string]$PackageCacheDir = '',
    [switch]$SkipPhase2,
    [switch]$SkipPackageRefresh,
    [switch]$ForceDownload,
    [switch]$ForceExtract
)

$ErrorActionPreference = 'Stop'
$ProgressPreference = 'SilentlyContinue'

function Write-InstallLog {
    param([string]$Message)
    Write-Host "[install] $Message"
}

function Fail {
    param([string]$Message)
    throw "[install] ERROR: $Message"
}

function Require-Command {
    param([string]$Name)
    if (-not (Get-Command $Name -ErrorAction SilentlyContinue)) {
        Fail "Required command not found: $Name"
    }
}

function Convert-ToPortableCygwinPath {
    param([Parameter(Mandatory = $true)][string]$WindowsPath)
    $resolved = [System.IO.Path]::GetFullPath($WindowsPath)
    if ($resolved -match '^([A-Za-z]):\\(.*)$') {
        $drive = $matches[1].ToLowerInvariant()
        $rest = ($matches[2] -replace '\\', '/')
        if ([string]::IsNullOrEmpty($rest)) {
            return "/mnt/$drive"
        }
        return "/mnt/$drive/$rest"
    }
    Fail "Cannot convert path to portable Cygwin form: $WindowsPath"
}

function Resolve-PortableArchivePath {
    if ($PortableArchivePath) {
        return $PortableArchivePath
    }
    $cacheDir = Join-Path $PortableRoot 'cache'
    New-Item -ItemType Directory -Force -Path $cacheDir | Out-Null
    return (Join-Path $cacheDir ([System.IO.Path]::GetFileName($PortableReleaseUrl)))
}

function Resolve-PortablePackageCacheDir {
    if ($PackageCacheDir) {
        return [System.IO.Path]::GetFullPath($PackageCacheDir)
    }
    return (Join-Path $PortableRoot 'cache\cygwin-pkg')
}

function Download-PortableArchive {
    param([string]$SourceUrl, [string]$DestinationPath)
    if ((-not $ForceDownload) -and (Test-Path $DestinationPath) -and ((Get-Item $DestinationPath).Length -gt 0)) {
        Write-InstallLog "Reusing existing archive $DestinationPath"
        return
    }
    Write-InstallLog "Downloading portable Cygwin from $SourceUrl"
    New-Item -ItemType Directory -Force -Path (Split-Path -Parent $DestinationPath) | Out-Null
    Invoke-WebRequest -Uri $SourceUrl -OutFile $DestinationPath
}

function Expand-PortableArchive {
    param([string]$ArchivePath, [string]$DestinationRoot)
    $bashPathProbe = Join-Path $DestinationRoot 'App\Runtime\Cygwin\bin\bash.exe'
    if ((-not $ForceExtract) -and (Test-Path $bashPathProbe)) {
        Write-InstallLog "Portable Cygwin already extracted under $DestinationRoot"
        return
    }
    Require-Command '7z'
    Write-InstallLog "Extracting portable Cygwin into $DestinationRoot"
    New-Item -ItemType Directory -Force -Path $DestinationRoot | Out-Null
    & 7z x $ArchivePath "-o$DestinationRoot" -y | Out-Host
}

function Resolve-PortableBash {
    param([string]$DestinationRoot)
    $preferred = @(
        (Join-Path $DestinationRoot 'App\Runtime\Cygwin\bin\bash.exe'),
        (Join-Path $DestinationRoot 'App\Runtime\cygwin\bin\bash.exe'),
        (Join-Path $DestinationRoot 'App\cygwin\bin\bash.exe'),
        (Join-Path $DestinationRoot 'App\Cygwin\bin\bash.exe')
    )
    foreach ($candidate in $preferred) {
        if (Test-Path $candidate) {
            return (Resolve-Path $candidate).Path
        }
    }
    $found = Get-ChildItem -Path $DestinationRoot -Recurse -Filter 'bash.exe' -File -ErrorAction SilentlyContinue |
        Where-Object { $_.FullName -match 'App\\Runtime\\Cygwin\\bin\\bash\.exe$|App\\Runtime\\cygwin\\bin\\bash\.exe$|App\\cygwin\\bin\\bash\.exe$|App\\Cygwin\\bin\\bash\.exe$' } |
        Select-Object -First 1
    if ($found) {
        return $found.FullName
    }
    Fail "Could not find portable Cygwin bash.exe under $DestinationRoot"
}

function Resolve-PortableCygwinRootPath {
    param([string]$DestinationRoot)
    $preferred = @(
        (Join-Path $DestinationRoot 'App\Runtime\Cygwin'),
        (Join-Path $DestinationRoot 'App\Runtime\cygwin'),
        (Join-Path $DestinationRoot 'App\Cygwin'),
        (Join-Path $DestinationRoot 'App\cygwin')
    )
    foreach ($candidate in $preferred) {
        if (Test-Path $candidate) {
            return [System.IO.Path]::GetFullPath($candidate)
        }
    }
    return [System.IO.Path]::GetFullPath($preferred[0])
}

function Reset-PortableCygwinRoot {
    param([string]$PortableCygwinRoot)
    if (Test-Path $PortableCygwinRoot) {
        Get-ChildItem -Force -Path $PortableCygwinRoot | Remove-Item -Recurse -Force
    } else {
        New-Item -ItemType Directory -Force -Path $PortableCygwinRoot | Out-Null
    }
}

function Copy-DefaultPortableData {
    param([string]$DestinationRoot)
    $defaultDataRoot = Join-Path $DestinationRoot 'App\DefaultData\cygwin'
    $portableDataRoot = Join-Path $DestinationRoot 'Data'
    $portableHome = Join-Path $portableDataRoot 'home'
    New-Item -ItemType Directory -Force -Path $portableDataRoot | Out-Null
    New-Item -ItemType Directory -Force -Path $portableHome | Out-Null

    if (Test-Path (Join-Path $defaultDataRoot 'home')) {
        Copy-Item -Path (Join-Path $defaultDataRoot 'home\*') -Destination $portableHome -Recurse -Force -ErrorAction SilentlyContinue
    }
}

function Write-PortableFstab {
    param([string]$PortableCygwinRoot)
    $portableDataHome = Join-Path $PortableRoot 'Data\home'
    New-Item -ItemType Directory -Force -Path $portableDataHome | Out-Null
    $fstabPath = Join-Path $PortableCygwinRoot 'etc\fstab'
    $fstabLines = @(
        '# /etc/fstab'
        '# Auto-generated for isolated portable pipeline bootstrap'
        "$PortableCygwinRoot\bin /usr/bin ntfs binary,auto,noacl 0 0"
        "$PortableCygwinRoot\lib /usr/lib ntfs binary,auto,noacl 0 0"
        "$PortableCygwinRoot / ntfs override,binary,auto,noacl 0 0"
        "$portableDataHome /home ntfs override,binary,auto,noacl 0 0"
        'none /mnt cygdrive binary,noacl,posix=0,user 0 0'
    )
    Set-Content -Path $fstabPath -Value $fstabLines
}

function Invoke-PortablePackageRefresh {
    param([string]$PortableCygwinRoot)
    $setupExeWindows = Join-Path $PortableRoot 'cache\setup-x86_64.exe'
    $pkgCache = Resolve-PortablePackageCacheDir
    Download-PortableArchive -SourceUrl $CygwinSetupUrl -DestinationPath $setupExeWindows
    New-Item -ItemType Directory -Force -Path $pkgCache | Out-Null
    Reset-PortableCygwinRoot -PortableCygwinRoot $PortableCygwinRoot

    if ([string]::IsNullOrWhiteSpace($CygwinPackages)) {
        $CygwinPackages = 'bash,curl,cygwin,gcc-core,gcc-g++,gnuplot-base,ImageMagick,libgd-devel,make,perl,perl-File-Which,perl-GD,perl-JSON,perl-JSON-MaybeXS,perl-Mojolicious,pkg-config,python3,python312,python312-devel,python312-imaging,python312-pip,python312-setuptools,python312-wheel,unzip,wget,which,zip'
    }

    Write-InstallLog "Refreshing portable Cygwin packages under $PortableCygwinRoot"
    $arguments = @(
        '-q',
        '-B',
        '-g',
        '-n',
        '-N',
        '--no-write-registry',
        '-R', $PortableCygwinRoot,
        '-l', $pkgCache,
        '-s', 'https://mirrors.kernel.org/sourceware/cygwin/',
        '-P', $CygwinPackages
    )
    $setupProcess = Start-Process -FilePath $setupExeWindows -ArgumentList $arguments -WorkingDirectory $pkgCache -PassThru -Wait
    if ($setupProcess.ExitCode -ne 0) {
        Fail "Portable Cygwin package refresh exited with code $($setupProcess.ExitCode)"
    }
}

function Warn-AboutExistingCygwinProcesses {
    param([string]$PortableCygwinRoot)
    $normalizedRoot = $PortableCygwinRoot.ToLowerInvariant()
    $foreign = Get-CimInstance Win32_Process -ErrorAction SilentlyContinue |
        Where-Object {
            $_.ExecutablePath -and
            $_.ExecutablePath.ToLowerInvariant() -match 'cygwin|bash\.exe|mintty\.exe|perl\.exe|python3(\.\d+)?\.exe' -and
            $_.ExecutablePath.ToLowerInvariant() -notlike "$normalizedRoot*"
        }
    if ($foreign) {
        Write-Warning "Existing Cygwin-related processes are active outside the portable root. Fully isolated mount validation may be contaminated until those processes exit."
    }
}

function Invoke-Phase2Installer {
    param([string]$PortableBashPath, [string]$PortableCygwinRoot)
    $repoRootUnix = Convert-ToPortableCygwinPath -WindowsPath $RepoRoot
    $phase2Unix = Convert-ToPortableCygwinPath -WindowsPath (Join-Path $RepoRoot $Phase2Script)
    $setupExeWindows = Join-Path $PortableRoot 'cache\setup-x86_64.exe'
    $setupExeUnix = Convert-ToPortableCygwinPath -WindowsPath $setupExeWindows
    $cygwinRootWindows = [System.IO.Path]::GetFullPath($PortableCygwinRoot)

    $commonExports = @(
        "export CYGWIN_SETUP_EXE='$setupExeUnix'",
        "export CYGWIN_ROOT_WINDOWS='$cygwinRootWindows'"
    )
    if ($CygwinPackages) {
        $escapedPackages = $CygwinPackages.Replace("'", "'\\''")
        $commonExports += "export CYGWIN_PACKAGES='$escapedPackages'"
    }

    $phase2CmdParts = @("cd '$repoRootUnix'") + $commonExports + @(
        "export CYGWIN_SKIP_PACKAGE_UPDATE=1",
        "/usr/bin/bash '$phase2Unix'"
    )
    $phase2Command = ($phase2CmdParts -join '; ')

    Write-InstallLog "Running repo-local pipeline bootstrap phase inside portable Cygwin"
    & $PortableBashPath -lc $phase2Command
    if ($LASTEXITCODE -ne 0) {
        Fail "Portable Cygwin repo-local bootstrap phase exited with code $LASTEXITCODE"
    }
}

$PortableRoot = [System.IO.Path]::GetFullPath($PortableRoot)
$RepoRoot = [System.IO.Path]::GetFullPath($RepoRoot)

Write-InstallLog "Portable root: $PortableRoot"
Write-InstallLog "Repository root: $RepoRoot"

if (-not (Test-Path $RepoRoot)) {
    Fail "Repository root does not exist: $RepoRoot"
}

$archivePath = Resolve-PortableArchivePath
Download-PortableArchive -SourceUrl $PortableReleaseUrl -DestinationPath $archivePath
Expand-PortableArchive -ArchivePath $archivePath -DestinationRoot $PortableRoot
$portableCygwinRoot = Resolve-PortableCygwinRootPath -DestinationRoot $PortableRoot

Copy-DefaultPortableData -DestinationRoot $PortableRoot
Warn-AboutExistingCygwinProcesses -PortableCygwinRoot $portableCygwinRoot

if (-not $SkipPackageRefresh) {
    Invoke-PortablePackageRefresh -PortableCygwinRoot $portableCygwinRoot
} else {
    Write-InstallLog "Skipping portable Cygwin package refresh because -SkipPackageRefresh was requested"
}

Write-PortableFstab -PortableCygwinRoot $portableCygwinRoot
$portableBash = Resolve-PortableBash -DestinationRoot $PortableRoot
Write-InstallLog "Portable Cygwin bash: $portableBash"

if (-not $SkipPhase2) {
    Invoke-Phase2Installer -PortableBashPath $portableBash -PortableCygwinRoot $portableCygwinRoot
} else {
    Write-InstallLog "Skipping phase-2 pipeline installation because -SkipPhase2 was requested"
}

Write-InstallLog "Windows portable Cygwin bootstrap completed"
