param(
    [Parameter(Mandatory = $true)]
    [string]$CsvPath,

    [Parameter(Mandatory = $true)]
    [string]$OutputXlsx,

    [string]$SheetName = "Table_S1_Common_Loci"
)

$ErrorActionPreference = "Stop"
Add-Type -AssemblyName System.IO.Compression
Add-Type -AssemblyName System.IO.Compression.FileSystem

function Escape-XmlText {
    param([string]$Text)
    if ($null -eq $Text) { return "" }
    return [System.Security.SecurityElement]::Escape([string]$Text)
}

function ConvertTo-ExcelColumn {
    param([int]$Index)
    $name = ""
    while ($Index -gt 0) {
        $remainder = ($Index - 1) % 26
        $name = [char](65 + $remainder) + $name
        $Index = [math]::Floor(($Index - 1) / 26)
    }
    return $name
}

function New-CellXml {
    param(
        [int]$Row,
        [int]$Col,
        [object]$Value,
        [int]$Style,
        [switch]$AsText
    )

    $cellRef = "{0}{1}" -f (ConvertTo-ExcelColumn $Col), $Row
    if ($AsText) {
        $escaped = Escape-XmlText ([string]$Value)
        return "<c r=`"$cellRef`" t=`"inlineStr`" s=`"$Style`"><is><t>$escaped</t></is></c>"
    }

    if ($null -eq $Value -or [string]::IsNullOrWhiteSpace([string]$Value)) {
        return "<c r=`"$cellRef`" s=`"$Style`"/>"
    }

    $num = 0.0
    if ([double]::TryParse([string]$Value, [ref]$num)) {
        return "<c r=`"$cellRef`" s=`"$Style`"><v>$([System.Globalization.CultureInfo]::InvariantCulture.TextInfo.ToLower($num.ToString('R', [System.Globalization.CultureInfo]::InvariantCulture)))</v></c>"
    }

    $escaped = Escape-XmlText ([string]$Value)
    return "<c r=`"$cellRef`" t=`"inlineStr`" s=`"$Style`"><is><t>$escaped</t></is></c>"
}

function Add-ZipEntryFromText {
    param(
        [System.IO.Compression.ZipArchive]$Archive,
        [string]$EntryName,
        [string]$Content
    )
    $entry = $Archive.CreateEntry($EntryName)
    $writer = New-Object System.IO.StreamWriter($entry.Open(), [System.Text.UTF8Encoding]::new($false))
    try {
        $writer.Write($Content)
    }
    finally {
        $writer.Dispose()
    }
}

if (-not (Test-Path -LiteralPath $CsvPath)) {
    throw "CSV file not found: $CsvPath"
}

$rows = Import-Csv -LiteralPath $CsvPath
if (-not $rows -or $rows.Count -eq 0) {
    throw "CSV has no data rows: $CsvPath"
}

$columns = @(
    "CHR","BP","SNP","gene","Smallest_ASSOC_P",
    "ALL_STD_DIFF_P","EUR_STD_DIFF_P","ASN_STD_DIFF_P",
    "ALL_FEMALE_BETA","ALL_FEMALE_SE","ALL_FEMALE_P",
    "ALL_MALE_BETA","ALL_MALE_SE","ALL_MALE_P",
    "ALL_DIFF_BETA","ALL_DIFF_SE","ALL_DIFF_P",
    "EUR_FEMALE_BETA","EUR_FEMALE_SE","EUR_FEMALE_P",
    "EUR_MALE_BETA","EUR_MALE_SE","EUR_MALE_P",
    "EUR_DIFF_BETA","EUR_DIFF_SE","EUR_DIFF_P",
    "ASN_FEMALE_BETA","ASN_FEMALE_SE","ASN_FEMALE_P",
    "ASN_MALE_BETA","ASN_MALE_SE","ASN_MALE_P",
    "ASN_DIFF_BETA","ASN_DIFF_SE","ASN_DIFF_P"
)

$headerRow2 = @(
    "CHR","BP","SNP","Gene","Smallest_ASSOC_P",
    "ALL_STD_DIFF_P","EUR_STD_DIFF_P","ASN_STD_DIFF_P",
    "BETA","SE","P",
    "BETA","SE","P",
    "BETA","SE","P",
    "BETA","SE","P",
    "BETA","SE","P",
    "BETA","SE","P",
    "BETA","SE","P",
    "BETA","SE","P",
    "BETA","SE","P"
)

$groupHeaders = @(
    @{ Start = 1; End = 5; Label = "Locus summary" },
    @{ Start = 6; End = 8; Label = "Standardized differential association" },
    @{ Start = 9; End = 11; Label = "ALL_FEMALE" },
    @{ Start = 12; End = 14; Label = "ALL_MALE" },
    @{ Start = 15; End = 17; Label = "ALL_DIFF" },
    @{ Start = 18; End = 20; Label = "EUR_FEMALE" },
    @{ Start = 21; End = 23; Label = "EUR_MALE" },
    @{ Start = 24; End = 26; Label = "EUR_DIFF" },
    @{ Start = 27; End = 29; Label = "ASN_FEMALE" },
    @{ Start = 30; End = 32; Label = "ASN_MALE" },
    @{ Start = 33; End = 35; Label = "ASN_DIFF" }
)

$columnWidths = @(
    8,12,16,18,16,
    16,16,16,
    13,13,13,
    13,13,13,
    13,13,13,
    13,13,13,
    13,13,13,
    13,13,13,
    13,13,13,
    13,13,13
)

$sheetData = New-Object System.Collections.Generic.List[string]

$row1Cells = New-Object System.Collections.Generic.List[string]
for ($c = 1; $c -le $columns.Count; $c++) {
    $label = ""
    foreach ($group in $groupHeaders) {
        if ($group.Start -eq $c) {
            $label = $group.Label
            break
        }
    }
    $row1Cells.Add((New-CellXml -Row 1 -Col $c -Value $label -Style 1 -AsText))
}
$sheetData.Add("<row r=`"1`" ht=`"22`" customHeight=`"1`">$($row1Cells -join '')</row>")

$row2Cells = New-Object System.Collections.Generic.List[string]
for ($c = 1; $c -le $headerRow2.Count; $c++) {
    $row2Cells.Add((New-CellXml -Row 2 -Col $c -Value $headerRow2[$c - 1] -Style 2 -AsText))
}
$sheetData.Add("<row r=`"2`" ht=`"22`" customHeight=`"1`">$($row2Cells -join '')</row>")

$textColumns = @("SNP", "gene")
$integerColumns = @("CHR", "BP")

$excelRow = 3
foreach ($row in $rows) {
    $cells = New-Object System.Collections.Generic.List[string]
    for ($c = 0; $c -lt $columns.Count; $c++) {
        $name = $columns[$c]
        $value = $row.$name
        if ($textColumns -contains $name) {
            $cells.Add((New-CellXml -Row $excelRow -Col ($c + 1) -Value $value -Style 3 -AsText))
        }
        elseif ($integerColumns -contains $name) {
            $cells.Add((New-CellXml -Row $excelRow -Col ($c + 1) -Value $value -Style 4))
        }
        else {
            $cells.Add((New-CellXml -Row $excelRow -Col ($c + 1) -Value $value -Style 5))
        }
    }
    $sheetData.Add("<row r=`"$excelRow`">$($cells -join '')</row>")
    $excelRow++
}

$lastRow = $excelRow - 1
$lastCol = ConvertTo-ExcelColumn $columns.Count
$mergeRefs = @(
    "A1:E1","F1:H1","I1:K1","L1:N1","O1:Q1","R1:T1",
    "U1:W1","X1:Z1","AA1:AC1","AD1:AF1","AG1:AI1"
)
$mergeXml = ($mergeRefs | ForEach-Object { "<mergeCell ref=`"$_`"/>" }) -join ""

$colsXml = New-Object System.Collections.Generic.List[string]
for ($i = 0; $i -lt $columnWidths.Count; $i++) {
    $width = [System.Globalization.CultureInfo]::InvariantCulture.TextInfo.ToLower(([double]$columnWidths[$i]).ToString("0.##", [System.Globalization.CultureInfo]::InvariantCulture))
    $colNum = $i + 1
    $colsXml.Add("<col min=`"$colNum`" max=`"$colNum`" width=`"$width`" customWidth=`"1`"/>")
}

$safeSheetName = if ($SheetName.Length -gt 31) { $SheetName.Substring(0, 31) } else { $SheetName }
$safeSheetName = $safeSheetName -replace '[\\\/\?\*\[\]:]', '_'

$worksheetXml = @"
<?xml version="1.0" encoding="UTF-8" standalone="yes"?>
<worksheet xmlns="http://schemas.openxmlformats.org/spreadsheetml/2006/main" xmlns:r="http://schemas.openxmlformats.org/officeDocument/2006/relationships">
  <dimension ref="A1:${lastCol}${lastRow}"/>
  <sheetViews>
    <sheetView workbookViewId="0">
      <pane ySplit="2" topLeftCell="A3" activePane="bottomLeft" state="frozen"/>
      <selection pane="bottomLeft" activeCell="A3" sqref="A3"/>
    </sheetView>
  </sheetViews>
  <sheetFormatPr defaultRowHeight="15"/>
  <cols>
    $($colsXml -join "`n    ")
  </cols>
  <sheetData>
    $($sheetData -join "`n    ")
  </sheetData>
  <autoFilter ref="A2:${lastCol}${lastRow}"/>
  <mergeCells count="$($mergeRefs.Count)">
    $mergeXml
  </mergeCells>
</worksheet>
"@

$stylesXml = @"
<?xml version="1.0" encoding="UTF-8" standalone="yes"?>
<styleSheet xmlns="http://schemas.openxmlformats.org/spreadsheetml/2006/main">
  <fonts count="2">
    <font>
      <sz val="11"/>
      <color theme="1"/>
      <name val="Calibri"/>
      <family val="2"/>
    </font>
    <font>
      <b/>
      <sz val="11"/>
      <color theme="1"/>
      <name val="Calibri"/>
      <family val="2"/>
    </font>
  </fonts>
  <fills count="4">
    <fill><patternFill patternType="none"/></fill>
    <fill><patternFill patternType="gray125"/></fill>
    <fill><patternFill patternType="solid"><fgColor rgb="FFE0E0E0"/><bgColor indexed="64"/></patternFill></fill>
    <fill><patternFill patternType="solid"><fgColor rgb="FFF2F2F2"/><bgColor indexed="64"/></patternFill></fill>
  </fills>
  <borders count="2">
    <border>
      <left/><right/><top/><bottom/><diagonal/>
    </border>
    <border>
      <left style="thin"><color auto="1"/></left>
      <right style="thin"><color auto="1"/></right>
      <top style="thin"><color auto="1"/></top>
      <bottom style="thin"><color auto="1"/></bottom>
      <diagonal/>
    </border>
  </borders>
  <cellStyleXfs count="1">
    <xf numFmtId="0" fontId="0" fillId="0" borderId="0"/>
  </cellStyleXfs>
  <cellXfs count="6">
    <xf numFmtId="0" fontId="0" fillId="0" borderId="0" xfId="0"/>
    <xf numFmtId="0" fontId="1" fillId="2" borderId="1" xfId="0" applyFont="1" applyFill="1" applyBorder="1" applyAlignment="1">
      <alignment horizontal="center" vertical="center"/>
    </xf>
    <xf numFmtId="0" fontId="1" fillId="3" borderId="1" xfId="0" applyFont="1" applyFill="1" applyBorder="1" applyAlignment="1">
      <alignment horizontal="center" vertical="center" wrapText="1"/>
    </xf>
    <xf numFmtId="0" fontId="0" fillId="0" borderId="1" xfId="0" applyBorder="1"/>
    <xf numFmtId="1" fontId="0" fillId="0" borderId="1" xfId="0" applyNumberFormat="1" applyBorder="1"/>
    <xf numFmtId="11" fontId="0" fillId="0" borderId="1" xfId="0" applyNumberFormat="1" applyBorder="1"/>
  </cellXfs>
  <cellStyles count="1">
    <cellStyle name="Normal" xfId="0" builtinId="0"/>
  </cellStyles>
</styleSheet>
"@

$workbookXml = @"
<?xml version="1.0" encoding="UTF-8" standalone="yes"?>
<workbook xmlns="http://schemas.openxmlformats.org/spreadsheetml/2006/main" xmlns:r="http://schemas.openxmlformats.org/officeDocument/2006/relationships">
  <sheets>
    <sheet name="$([System.Security.SecurityElement]::Escape($safeSheetName))" sheetId="1" r:id="rId1"/>
  </sheets>
</workbook>
"@

$workbookRelsXml = @"
<?xml version="1.0" encoding="UTF-8" standalone="yes"?>
<Relationships xmlns="http://schemas.openxmlformats.org/package/2006/relationships">
  <Relationship Id="rId1" Type="http://schemas.openxmlformats.org/officeDocument/2006/relationships/worksheet" Target="worksheets/sheet1.xml"/>
  <Relationship Id="rId2" Type="http://schemas.openxmlformats.org/officeDocument/2006/relationships/styles" Target="styles.xml"/>
</Relationships>
"@

$rootRelsXml = @"
<?xml version="1.0" encoding="UTF-8" standalone="yes"?>
<Relationships xmlns="http://schemas.openxmlformats.org/package/2006/relationships">
  <Relationship Id="rId1" Type="http://schemas.openxmlformats.org/officeDocument/2006/relationships/officeDocument" Target="xl/workbook.xml"/>
</Relationships>
"@

$contentTypesXml = @"
<?xml version="1.0" encoding="UTF-8" standalone="yes"?>
<Types xmlns="http://schemas.openxmlformats.org/package/2006/content-types">
  <Default Extension="rels" ContentType="application/vnd.openxmlformats-package.relationships+xml"/>
  <Default Extension="xml" ContentType="application/xml"/>
  <Override PartName="/xl/workbook.xml" ContentType="application/vnd.openxmlformats-officedocument.spreadsheetml.sheet.main+xml"/>
  <Override PartName="/xl/worksheets/sheet1.xml" ContentType="application/vnd.openxmlformats-officedocument.spreadsheetml.worksheet+xml"/>
  <Override PartName="/xl/styles.xml" ContentType="application/vnd.openxmlformats-officedocument.spreadsheetml.styles+xml"/>
</Types>
"@

$outputDir = Split-Path -Parent $OutputXlsx
if ([string]::IsNullOrWhiteSpace($outputDir)) {
    $outputDir = (Get-Location).Path
}
if (-not (Test-Path -LiteralPath $outputDir)) {
    New-Item -ItemType Directory -Path $outputDir | Out-Null
}
$outputFull = [System.IO.Path]::GetFullPath((Join-Path $outputDir (Split-Path -Leaf $OutputXlsx)))

if (Test-Path -LiteralPath $outputFull) {
    Remove-Item -LiteralPath $outputFull -Force
}

$fs = [System.IO.File]::Open($outputFull, [System.IO.FileMode]::Create)
try {
    $archive = New-Object System.IO.Compression.ZipArchive($fs, [System.IO.Compression.ZipArchiveMode]::Create, $false)
    try {
        Add-ZipEntryFromText -Archive $archive -EntryName "[Content_Types].xml" -Content $contentTypesXml
        Add-ZipEntryFromText -Archive $archive -EntryName "_rels/.rels" -Content $rootRelsXml
        Add-ZipEntryFromText -Archive $archive -EntryName "xl/workbook.xml" -Content $workbookXml
        Add-ZipEntryFromText -Archive $archive -EntryName "xl/_rels/workbook.xml.rels" -Content $workbookRelsXml
        Add-ZipEntryFromText -Archive $archive -EntryName "xl/worksheets/sheet1.xml" -Content $worksheetXml
        Add-ZipEntryFromText -Archive $archive -EntryName "xl/styles.xml" -Content $stylesXml
    }
    finally {
        $archive.Dispose()
    }
}
finally {
    $fs.Dispose()
}

Write-Output "Wrote Excel workbook: $outputFull"
