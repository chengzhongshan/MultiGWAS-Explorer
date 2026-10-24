# Optional SAS visualization module validation

The standalone module was executed through MCP4SAS runtime af05aed11a74a2146823023f3bf4b6132aaa4a1c against SAS OnDemand. Five current captures are retained: cluster_branch_exact, tp53_protein, oncoprint, oncoprint_edges and the intentionally FAILED oncoprint_invalid. Exact source/input/dependency hashes, submitted SAS, complete redacted logs, numerical outputs, generated templates and PNG/SVG files are in EXECUTION.json and its referenced directories.

Run from the visualization directory:

    perl tests/validate_visualization.pl

Actual local result: **2426 passed; 0 failed**. This checks the recorded executions; it does not submit fresh SAS jobs. Run the example client separately to reproduce a live job in your environment. New runs go under ignored .runs/ and are not automatically accepted as replacements for the archived evidence.

Independent tests check 400 official GTEx medians, log and sample-Z transformations, all 400 heatmap cells, every greedy complete-link merge in both trees and their mean-normalized heights. The CL10 branch retains eight leaves and exact original heights; actual final PROC COMPARE is zero, while the native rounding discrepancy remains recorded as 4096. Every displayed full and selected gene name is checked in SVG.

The protein example retains six GDC sites and ten project presences, with aligned dashed guides, complete vertical labels and an actual zero PROC COMPARE status. Scientific CSVs match the protected previous outputs. Neither a project-presence dot nor a blank combination estimates mutation prevalence.

The reusable oncoprint retains 32 distinct events in 31 positive pairs, 94 assessed opportunities out of 96, and gene-specific denominators. The 2x2 case checks duplicate removal, two categories in one cell, a 100% result and a missing percentage with visible 0/0. Missing datasets, unsupported classes and events in unassessable pairs yield three named errors and return code one; their job remains FAILED. Top sample-summary strokes use the cohort strip's colors, gene summaries use a separate purple, and emitted vertical/horizontal strokes are both 12px. Numerical and graphical checks are separate.

The first oncoprint upload trial failed because an unquoted template filename was not rewritten by the client. A quoted argument plus macro dequoting and home-path normalization fixes the actual cause. Earlier styling trials and expected-error captures remain diagnostics. Native SAS group styling initially made the second group's top bars dashed; explicit solid line attributes correct it. A zero-denominator label initially inherited a missing x coordinate; a defined label position now displays 0/0 without inventing a percentage. No old validated source or output was silently rewritten.

The module adds callable clustergram/branch definitions, complete original-pipeline mutation examples, scoped display helpers and MGXOncoPrint. The existing differential-GWAS entrypoints, chromosome palette and label omissions are unchanged. Source snapshots are intentionally isolated from a running GWAS session because same-named dependency macros can replace definitions in memory.

## Scope and remaining limits

The matrix, selected branch, numeric mutation queries and assay fixtures are bounded tests. Arbitrary branches, huge matrices, internal-standardization label failures, other SAS platforms and every class mapping remain unverified. The generic full-length protein interval transformer edge case is handled explicitly in the isoform example, not declared fixed for every annotation input. The module does not fit new GWAS models, infer isoform expression or rerun donor-level/variant-calling pipelines. The original warning lines remain in the logs. Source licensing and public data provenance are preserved; no credentials or user account configuration are distributed.

### cluster_branch_exact warning lines

- `WARNING: A very large output size of (3958, 4375) is in effect. This could make Java VM run out of memory and result in some Java `
- `WARNING: A very large output size of (2560, 2840) is in effect. This could make Java VM run out of memory and result in some Java `
- `WARNING: This CREATE TABLE statement recursively references the target table. A consequence of this is a possible data integrity `
- `WARNING: A very large output size of (4000, 1920) is in effect. This could make Java VM run out of memory and result in some Java `
- `WARNING: The variable 'new_y_:'n in the DROP, KEEP, or RENAME list has never been referenced.`

### oncoprint warning lines

- `WARNING: A very large output size of (3420, 2430) is in effect. This could make Java VM run out of memory and result in some Java `

### oncoprint_edges warning lines

No SAS WARNING lines.

### oncoprint_invalid warning lines

No SAS WARNING lines.

### tp53_protein warning lines

No SAS WARNING lines.

