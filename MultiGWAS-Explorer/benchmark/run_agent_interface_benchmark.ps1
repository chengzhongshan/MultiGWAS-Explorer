param(
    [int]$Runs = 10,
    [int]$RecoveryRuns = 3,
    [int]$Port = 18081,
    [switch]$UseExistingServer,
    [string]$OutputDir = ''
)

$ErrorActionPreference = 'Stop'
$repo = Split-Path -Parent $PSScriptRoot
$outDir = if ([string]::IsNullOrWhiteSpace($OutputDir)) {
    Join-Path $PSScriptRoot 'agent_interface'
} else {
    [IO.Path]::GetFullPath($OutputDir)
}
New-Item -ItemType Directory -Force -Path $outDir | Out-Null
$repoFull = [IO.Path]::GetFullPath($repo).TrimEnd([IO.Path]::DirectorySeparatorChar, [IO.Path]::AltDirectorySeparatorChar)
$outDirFull = [IO.Path]::GetFullPath($outDir)
if (-not $outDirFull.StartsWith($repoFull + [IO.Path]::DirectorySeparatorChar, [StringComparison]::OrdinalIgnoreCase)) {
    throw "OutputDir must be inside the repository so MCP background logs can be resolved: $outDirFull"
}
$outDirRel = $outDirFull.Substring($repoFull.Length + 1).Replace('\', '/')
$perl = (Get-Command perl -ErrorAction Stop).Source
$specRel = 'configs/spec_pgc_scz_sex_common_automation.json'
$spec = Join-Path $repo $specRel
$artifactRel = @(
    'configs/auto_PGC_SCZ_female_vs_male_diff_effects_merge.json',
    'configs/auto_PGC_SCZ_female_vs_male_diff_effects_diff.json',
    'configs/auto_PGC_SCZ_female_vs_male_diff_effects_preset.json',
    'configs/auto_PGC_SCZ_female_vs_male_diff_effects_runner.json'
)

function Get-StringSha256([string]$Text) {
    $sha = [System.Security.Cryptography.SHA256]::Create()
    try {
        $bytes = [System.Text.Encoding]::UTF8.GetBytes($Text)
        return ([BitConverter]::ToString($sha.ComputeHash($bytes))).Replace('-', '').ToLowerInvariant()
    } finally {
        $sha.Dispose()
    }
}

function Get-ArtifactSet([string]$PathName, [int]$RunNumber) {
    $records = @()
    foreach ($rel in $artifactRel | Sort-Object) {
        $path = Join-Path $repo $rel
        if (-not (Test-Path -LiteralPath $path)) { throw "Missing expected artifact: $path" }
        $h = (Get-FileHash -Algorithm SHA256 -LiteralPath $path).Hash.ToLowerInvariant()
        $records += [pscustomobject]@{
            path = $PathName
            run = $RunNumber
            artifact = $rel.Replace('\', '/')
            sha256 = $h
            bytes = (Get-Item -LiteralPath $path).Length
        }
    }
    $canonical = ($records | ForEach-Object { "$($_.artifact)=$($_.sha256)`n" }) -join ''
    return [pscustomobject]@{ SetHash = Get-StringSha256 $canonical; Records = $records }
}

function Invoke-Mcp([hashtable]$Payload) {
    $json = $Payload | ConvertTo-Json -Depth 12 -Compress
    return Invoke-RestMethod -Uri "http://127.0.0.1:$Port/mcp" -Method Post -ContentType 'application/json' -Body $json
}

function Get-ToolText($Response) {
    if ($null -ne $Response.error) { return "MCP ERROR: $($Response.error.message)" }
    $encoded = $Response.result.content[0].text
    try {
        $inner = $encoded | ConvertFrom-Json
        if ($null -ne $inner.content) { return [string]$inner.content[0].text }
    } catch {}
    return [string]$encoded
}

function Invoke-AgentConfig([string]$SpecArgument, [string]$LogPath) {
    $start = [Diagnostics.Stopwatch]::StartNew()
    $call = Invoke-Mcp @{
        jsonrpc = '2.0'; id = 100; method = 'tools/call'
        params = @{
            name = 'auto_prepare_and_run_diff_gwas'
            arguments = @{
                spec_file = $SpecArgument
                mode = 'configs'
                skip_plots = 'true'
                output_file = $LogPath
            }
        }
    }
    $text = Get-ToolText $call
    if ($text -notmatch 'PID:\s*(\d+)') {
        $start.Stop()
        return [pscustomobject]@{ Success = $false; Seconds = $start.Elapsed.TotalSeconds; Text = $text }
    }
    $pidValue = [int]$Matches[1]
    $finalText = ''
    $localLogPath = Join-Path $repo ($LogPath.Replace('/', '\'))
    $localErrPath = $localLogPath + '.stderr.log'
    for ($poll = 0; $poll -lt 100; $poll++) {
        Start-Sleep -Milliseconds 200
        $logText = if (Test-Path -LiteralPath $localLogPath) {
            Get-Content -LiteralPath $localLogPath -Raw -ErrorAction SilentlyContinue
        } else { '' }
        $errText = if (Test-Path -LiteralPath $localErrPath) {
            Get-Content -LiteralPath $localErrPath -Raw -ErrorAction SilentlyContinue
        } else { '' }
        $finalText = [string]$logText + "`n" + [string]$errText
        if ($finalText -match 'Automation complete\.') {
            break
        }
        if ($finalText -match '(?i)(does not exist|no such file|cannot open|fatal|error:)') { break }
    }
    $start.Stop()
    $ok = $finalText -match 'Automation complete\.'
    return [pscustomobject]@{ Success = [bool]$ok; Seconds = $start.Elapsed.TotalSeconds; Text = $finalText }
}

$serverStdout = Join-Path $outDir 'server.stdout.log'
$serverStderr = Join-Path $outDir 'server.stderr.log'
$server = $null
if (-not $UseExistingServer) {
    $server = Start-Process -FilePath $perl `
        -ArgumentList @('server.pl', 'daemon', '-m', 'production', '-l', "http://127.0.0.1:$Port") `
        -WorkingDirectory $repo -WindowStyle Hidden -PassThru `
        -RedirectStandardOutput $serverStdout -RedirectStandardError $serverStderr
}

$results = @()
$manifest = @()
$recoveries = @()
Write-Host "Agent benchmark parameters: Runs=$Runs RecoveryRuns=$RecoveryRuns Port=$Port UseExistingServer=$UseExistingServer"
try {
    $ready = $false
    for ($attempt = 0; $attempt -lt 100; $attempt++) {
        Start-Sleep -Milliseconds 200
        try {
            $init = Invoke-Mcp @{ jsonrpc = '2.0'; id = 1; method = 'initialize'; params = @{} }
            if ($init.result.serverInfo.name -eq 'PerlBioServer') { $ready = $true; break }
        } catch {}
    }
    if (-not $ready) { throw "MCP server did not become ready on port $Port" }

    for ($run = 1; $run -le $Runs; $run++) {
        $cliLog = Join-Path $outDir ("cli_run_{0:D2}.log" -f $run)
        $sw = [Diagnostics.Stopwatch]::StartNew()
        Push-Location $repo
        try {
            $savedErrorPreference = $ErrorActionPreference
            $ErrorActionPreference = 'Continue'
            & $perl 'auto_prepare_and_run_diff_gwas.pl' '--spec' $specRel '--mode' 'configs' '--skip-plots' 2>&1 |
                Out-File -LiteralPath $cliLog -Encoding utf8
            $cliExit = $LASTEXITCODE
            $ErrorActionPreference = $savedErrorPreference
        } finally {
            $ErrorActionPreference = 'Stop'
            Pop-Location
        }
        $sw.Stop()
        $cliSet = Get-ArtifactSet 'CLI' $run
        $manifest += $cliSet.Records
        $results += [pscustomobject]@{
            path = 'CLI'; run = $run; success = [int]($cliExit -eq 0)
            elapsed_seconds = [math]::Round($sw.Elapsed.TotalSeconds, 6)
            artifact_set_sha256 = $cliSet.SetHash; parity_with_cli = 1
        }

        Start-Sleep -Milliseconds 1100
        $agentLogRel = "$outDirRel/agent_run_{0:D2}.log" -f $run
        $agent = Invoke-AgentConfig $specRel $agentLogRel
        $agentSet = Get-ArtifactSet 'MCP_AGENT' $run
        $manifest += $agentSet.Records
        $parity = [int]($agentSet.SetHash -eq $cliSet.SetHash)
        $results += [pscustomobject]@{
            path = 'MCP_AGENT'; run = $run; success = [int]$agent.Success
            elapsed_seconds = [math]::Round($agent.Seconds, 6)
            artifact_set_sha256 = $agentSet.SetHash; parity_with_cli = $parity
        }
        Write-Host "Completed paired CLI/MCP run $run of $Runs"
    }

    for ($run = 1; $run -le $RecoveryRuns; $run++) {
        Start-Sleep -Milliseconds 1100
        $badLog = "$outDirRel/recovery_bad_{0:D2}.log" -f $run
        $bad = Invoke-AgentConfig 'configs/INTENTIONALLY_MISSING_BENCHMARK.json' $badLog
        Start-Sleep -Milliseconds 1100
        $goodLog = "$outDirRel/recovery_good_{0:D2}.log" -f $run
        $good = Invoke-AgentConfig $specRel $goodLog
        $goodSet = Get-ArtifactSet 'MCP_AGENT_RECOVERY' $run
        $recoveries += [pscustomobject]@{
            run = $run
            injected_failure_detected = [int](-not $bad.Success)
            corrected_retry_success = [int]$good.Success
            recovery_seconds = [math]::Round($good.Seconds, 6)
            recovered_artifact_set_sha256 = $goodSet.SetHash
        }
    }
} finally {
    if (-not $UseExistingServer -and $null -ne $server -and -not $server.HasExited) {
        Stop-Process -Id $server.Id -Force
        $server.WaitForExit()
    }
}

$results | Export-Csv -Delimiter "`t" -NoTypeInformation -LiteralPath (Join-Path $outDir 'agent_interface_runs.tsv')
$manifest | Export-Csv -Delimiter "`t" -NoTypeInformation -LiteralPath (Join-Path $outDir 'artifact_manifest.tsv')
$recoveries | Export-Csv -Delimiter "`t" -NoTypeInformation -LiteralPath (Join-Path $outDir 'failure_recovery.tsv')

$cliTimes = @($results | Where-Object path -eq 'CLI' | ForEach-Object elapsed_seconds)
$agentTimes = @($results | Where-Object path -eq 'MCP_AGENT' | ForEach-Object elapsed_seconds)
$summary = @(
    '# Agent-interface benchmark',
    '',
    "Predeclared task: generate the four deterministic configuration artifacts from ``$specRel`` in ``configs`` mode with plotting disabled.",
    '',
    "- CLI successes: $(($results | Where-Object { $_.path -eq 'CLI' -and $_.success -eq 1 }).Count)/$Runs",
    "- MCP-agent successes: $(($results | Where-Object { $_.path -eq 'MCP_AGENT' -and $_.success -eq 1 }).Count)/$Runs",
    "- Exact four-artifact checksum parity: $(($results | Where-Object { $_.path -eq 'MCP_AGENT' -and $_.parity_with_cli -eq 1 }).Count)/$Runs",
    "- Median CLI elapsed seconds: $((($cliTimes | Sort-Object)[[math]::Floor(($cliTimes.Count - 1) / 2)]))",
    "- Median MCP-agent elapsed seconds: $((($agentTimes | Sort-Object)[[math]::Floor(($agentTimes.Count - 1) / 2)]))",
    "- Injected failures detected: $(($recoveries | Where-Object injected_failure_detected -eq 1).Count)/$RecoveryRuns",
    "- Corrected retries recovered: $(($recoveries | Where-Object corrected_retry_success -eq 1).Count)/$RecoveryRuns",
    '',
    'The MCP timing includes HTTP dispatch, background-process launch, and status polling; the CLI timing does not.'
)
$summary | Out-File -LiteralPath (Join-Path $outDir 'SUMMARY.md') -Encoding utf8
Write-Host "Wrote agent-interface benchmark to $outDir"
