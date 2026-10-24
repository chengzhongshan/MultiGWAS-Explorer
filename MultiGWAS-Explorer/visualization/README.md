# Clustergrams and mutation tracks for MultiGWAS-Explorer

This optional SAS module adds the author's expression clustergram, branch-selection and gene/protein mutation pipelines. It carries the tested display improvements from *AI-Powered SAS*: readable full gene labels, exact inherited subtree heights, aligned mutation guides, vertical variant labels, categorical mutation colors, larger markers, compact lower clearance and thicker CNV intervals. Existing GWAS entrypoints and their defaults are unchanged.

## Run an example

From this `visualization` directory, with Perl, JSON::PP, LWP::UserAgent, Digest::SHA and a configured MCP4SAS loopback service:

```bash
perl scripts/run_visualization.pl list
perl scripts/run_visualization.pl run cluster_branch_exact
perl scripts/run_visualization.pl run tp53_protein
perl scripts/run_visualization.pl run tp53_gene
perl scripts/run_visualization.pl run tp53_isoforms
perl scripts/run_visualization.pl run cnv_intervals
```

The client uses `BOOK_MCP_URL` (default `http://127.0.0.1:18789/mcp`) for compatibility with the source reader client. It uploads the catalog's local input files and macro definitions, captures the complete submitted SAS/log/tables/images in a fresh ignored `.runs/` directory, and verifies recorded artifact hashes with `perl scripts/run_visualization.pl verify .runs/RUN_DIRECTORY`. Hash verification establishes capture integrity; numerical correctness requires the separate checks in `tests/`. No SAS account details or private server runtime are bundled. This module does not install MCP4SAS or fit GWAS models.

## Use your own SAS tables

In a separate SAS session, set the absolute local or SAS-server path to this module, then load a family:

```sas
%let MGX_VIZ_ROOT=/your/server/path/visualization;
%include "&MGX_VIZ_ROOT/load_visualization.sas";
%MGXLoadVisualization(family=clustergram);
```

Choose `family=mutation` for the real TP53 pipelines, or `family=cnv` for SNP/interval overlays. The loader includes pinned definitions and named display helpers, sets `BOOK_OUT` to WORK by default, and executes no trailing source demos. These pinned dependency definitions can replace same-named macros in memory; use an isolated session rather than loading them over a running differential-GWAS job. Supply `outdir=` to select an existing writable SAS-server directory.

For a clustergram, replace the input DATA step of `companion/research/clustergram_branches/programs/cluster_branch.sas` with your own matrix. One row represents a gene; the identifier is `gene`, and the ten named numeric columns are tissues. The example first computes log2(TPM+1), then within-gene sample Z scores, and calls `%clustergram4sas` with `stdize=0`, `perform_log2_for_numeric_vars=0`, `switch_rownames_as_colnames=1` and `cluster_type=3`. Its underlying PROC CLUSTER uses complete linkage. The data are 40 fixed GTEx v8 median profiles, not donor measurements or an inferred biomarker signature. The branch example calls `%review_tree_branches` and `%subbranches`, retains the original tree and restores its exact numeric node heights before rendering CL10; it does not recluster the eight selected genes. Read the source before changing a branch identifier. The native rounding discrepancy is preserved and checked separately.

For mutations, use `companion/research/ch17_layout/programs/ch17_gdc_tp53_protein.sas` or `ch17_gdc_tp53_gene.sas` as the complete caller. The named `%BookGDCMutationProteinTrack` and `%BookGDCMutationGeneTrack` derivatives retain the original `%Long_format_muts2ProteinTrack`/gene-track pipeline with explicit numeric query intervals and scoped warning fixes. The `Book` prefix identifies the source derivative; it is not a new upstream original. Parameters such as `dotsize`, `barthickness`, `design_height`, `text_rotate_angle`, label headroom and spacing are separate controls. `%BookChapter17Layout` is a final display helper over an already generated `BedGraph`, not a raw-data import macro. It enlarges marker strokes and tightens the lower bound while PROC COMPARE checks the original drawing columns. The protein and gene examples use 14px points and a 16px annotation body; the output height is 300 or 350 design pixels. Avoid reducing height until labels touch panel boundaries.

CNV thickness is independent of SNP size. The CNV example changes only the HIGHLOWPLOT interval stroke from 4 to 8px after the original pipeline. It retains ST/END endpoints and masks CNV midpoint circles only in the scatter layer. The underlying macro already accepts `highlow_line_cmd=`; direct callers can choose `%str(thickness=8 color=darkorange pattern=solid)` without editing the macro. Synthetic CNV data remain labeled synthetic. The genome-wide chromosome palette and intentionally omitted tick labels are untouched.

## Interpret the TP53 examples correctly

The six real GDC substitutions have ten variant-project presences across TCGA-BRCA and TCGA-LUAD; these are not patients, mutation frequencies or a complete catalogue. Genomic coordinates use GRCh38, and protein coordinates use UniProt P04637. Dashed guides and circles share the original residue/genomic coordinate, rather than the reshaped interval endpoint. Nonsense and missense consequences retain explicit colors. The nine-isoform display uses independently checked UniProt sequences, with blue shared sequence and orange alternative C-terminal sequence. Its canonical-site projections do not establish isoform-specific tumor mutations or expression. Coincident full-length intervals exposed a generic annotation-transformer edge case; the final isoform display reads verified sequence endpoints directly instead of inheriting the collapsed placeholders.

## Provenance and limits

`SOURCE_MANIFEST.json` pins each included SAS/data file by SHA-256 and names the source book path. Original source notices remain attached. The original author-library archive matches COVID19_GWAS_Analyzer commit `4fa6fa710b8506b4aa410c0dfed304d447e41c8a` for the files identified in its manifests; four local tree helpers and named book derivatives have separate hashes, not invented upstream commit identities. MultiGWAS-Explorer's GPL-3 license applies to this module's distribution; preserve attribution and the existing data-source conditions. GTEx v8 medians, GDC/Ensembl GRCh38 and UniProt source manifests accompany the small inputs.

Only the bounded matrix, branch, selected variants and numeric query paths are claimed tested. Large matrices, arbitrary branch combinations, legacy internal-standardization label failures and other SAS platforms remain outside that evidence. A successful plotting run is not a new GWAS, variant-calling or isoform-expression analysis. Review the machine-readable validation and complete logs before generalizing a result.

## Reusable oncoprint macro

Run `perl scripts/run_visualization.pl run oncoprint` for an invented 12-sample/eight-gene assay fixture. The generic `%MGXOncoPrint` is in `companion/research/multigwas_port/MGXOncoPrint.sas`; its template preserves the book's multi-event glyphs, burden bars, clinical strip and positive/assessable counts. It is independent of the GWAS macro library. Include that macro, then supply four input tables and the provided `oncoprint_template.sas` path:

```sas
%MGXOncoPrint(samples=work.my_samples,genes=work.my_genes,
 callability=work.my_callability,variants=work.my_events,
 prefix=my_onco,outlib=work,outdir=/your/output/path,image=my_oncoprint,
 template_file=/your/module/path/companion/research/multigwas_port/oncoprint_template.sas,
 width=7.6,height=5.4,dpi=450);
```

Input schema: samples have character `sample`/`cohort` and unique numeric `sample_index` 1..N; genes have character `gene` and unique numeric `gene_index` 1..G. Callability has `sample`, `gene_index` and numeric `callable` (1 assessed, 0 unavailable); omitted pairs display as not assessable. Events have character `sample`, `gene`, `variant_id`, `mutation_class`, with classes Missense, Truncating, Splice or Inframe. Map source consequence names explicitly before calling. Gene indices run from bottom to top. Up to eight cohort palette colors are recycled; for many categories, provide a revised explicit cohort map rather than assuming all colors are unique.

The macro does not infer cohort membership, patient counts or clinical significance. It deduplicates exact event records, rejects conflicting classifications and events without assessable pairs, and retains each distinct event as a separate colored glyph. Frequencies count positive gene-sample pairs once and divide by assessed pairs; zero-denominator percentages are missing. Outputs are PREFIX_tiles/frequency/glyphs/burden/render. `summary_thickness=12` gives the vertical burden bars and horizontal gene bars the same nominal pixel stroke. Top bars use the cohort strip's colors; `gene_summary_color=CX6B5B95` supplies a distinct gene-summary color. These affect appearance, not counts. `MGX_ONCOPRINT_RC=0` means the macro reached rendering; inspect the SAS log and outputs as well. Input datasets are not replaced. WORK._MGXO_* scratch tables and numeric MGX formats are reserved; use an isolated session for simultaneous jobs. The `oncoprint_invalid` lesson intentionally reports FAILED while recording return-code checks.
