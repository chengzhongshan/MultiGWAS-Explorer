# Tooling Strategy

## Best Default: Hybrid

Use a Codex skill for scientific workflow control and use MCP tools for atomic execution.

This is better than putting the whole GWAS analysis into `server.pl` because future GWAS tasks vary by columns, traits, cohorts, alleles, genome builds, and plot requirements. A skill lets Codex inspect the current dataset and choose or patch scripts accordingly.

## Put In The Skill

- Workflow order and scientific checks.
- Column-mapping decisions.
- Interpretation cautions, such as not treating A1/A2 as REF/ALT without documentation.
- File-size and SAS ODA quota strategy.
- Which local scripts to reuse or adapt.
- Validation requirements and expected outputs.

## Put In Perl MCP Tools

Good MCP tools are small, stable, and reusable:

- `run_perl_or_bash_cmd`
- `run_sas_codes_or_script_in_ODA`
- `list_oda_files`
- `upload_oda_file`
- `download_oda_file`
- `delete_oda_file`
- `bgzip_tabix_index`
- `tabix_query`
- `zcat_head`
- `file_size`

These tools should accept explicit paths and return machine-readable status where possible.

## Keep As Project Scripts

Use standalone Perl/Bash/SAS scripts for domain-specific transformations:

- Merging a particular GWAS release.
- Pairing female and male GWAS tags.
- Computing differential effects.
- Preparing a SAS macro-specific wide table.
- Generating a specific plot style.
- Deciding whether local GTF gene tracks should include non-protein-coding
  genes or be restricted to protein-coding genes only.

Project scripts are easier to version, test on real files, and reuse outside Codex.

## Decision Rule

- If the step needs scientific judgment or dataset-specific column mapping: skill + project script.
- If the step is a repeated mechanical operation with stable inputs/outputs: MCP tool.
- If the step is a one-off exploratory command: use `run_perl_or_bash_cmd`.
- If failure would waste SAS ODA quota or produce huge files: make it a script with explicit checks.

## Recommended MCP Server Additions

Add these only after the current workflow stabilizes:

1. `oda_file_exists(path)` returning `{exists, size}`.
2. `oda_download_verified(remote, local_dir)` returning the local path and size.
3. `oda_upload_verified(local)` returning remote basename and size.
4. `tabix_index(input, chr_col, start_col, end_col, output)` using bgzip/tabix.
5. `tabix_region(input_bgz, region)` returning a small preview/count.

Avoid embedding high-level GWAS merge or Manhattan plot logic directly in the MCP server.
