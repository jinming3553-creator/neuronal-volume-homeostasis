# LRRC8A study — analysis code

This repository contains the custom R scripts used for the human genetic, epigenetic,
transcriptomic, and visualization analyses supporting the study "Disrupted neuronal volume homeostasis transforms electrical dysfunction into structural pathology."

## Scope

The scripts reproduce the custom analysis steps used to integrate schizophrenia-associated
DMRs/DMPs, Brain-mMeta mQTL data, PGC3 GWAS/SMR results, BrainxQTL mQTL–eQTL effects,
BrainSeq expression data, and selected genomic visualizations.

The study is primarily experimental; these scripts document the custom computational analyses
used in the human discovery component.

## Scripts

- `00_merge_PGC3_SMR_results.R` — merges PGC3 SMR result sets.
- `01_DMR_GWAS_overlap.R` — annotates schizophrenia-associated DMRs by overlap with PGC3 loci.
- `02_select_DMR_probes.R` — selects Brain-mMeta probes within ±2.5 kb of DMRs.
- `03_merge_mQTL_SMR_DMP.R` — merges mQTL–SMR results with DMP results.
- `04_filter_DMPs_by_significant_DMR.R` — retains DMPs near significant DMRs.
- `04b_DMR_DMP_overlap.R` — constructs DMR–DMP overlap tables.
- `05_mQTL_eQTL_correlation.R` — correlates cg13630871 mQTL and nearby-gene eQTL effect sizes.
- `06_CpG_gene_expression_correlation.R` — evaluates CpG–gene expression correlations.
- `07_BrainSeq_differential_expression.R` — summarizes LRRC8-family expression differences in BrainSeq data.
- `08_SMR_plots.R` — generates SMR locus/effect plots using the SMR plotting functions.
- `09_RCircos_plot.R` — generates circular genomic visualization of SMR results.

## Data

Large or third-party input datasets are not redistributed in this code archive. They should be
obtained from the original resources described in the manuscript and its Materials and Methods,
including Brain-mMeta/GSE74193, PGC3, BrainxQTL, BrainSeq, PsychENCODE, ENIGMA2, and
1000 Genomes where applicable.

The scripts retain the original analysis logic and expected input structures. Local absolute paths
from the working analysis environment have been removed.

## Path configuration

By default, scripts expect data under:

- `../data/LRRC8A`
- `../data/LRRC8A_update`
- `../external_code`

Alternatively set:

```bash
export LRRC8A_DATA_ROOT=/path/to/LRRC8A
export LRRC8A_UPDATE_ROOT=/path/to/LRRC8A_update
export LRRC8A_EXTERNAL_CODE_ROOT=/path/to/external_code
```

## Software

The analyses use R. `09_RCircos_plot.R` additionally uses the `RCircos` package.
`08_SMR_plots.R` calls the SMR plotting functions in `plot_SMR.r`; obtain the compatible
SMR plotting script/software from the original SMR distribution/source used for the analysis.

## Reproducibility note

These are the custom analysis scripts used in the study, cleaned for public release by replacing
machine-specific file paths and renaming files descriptively. No scientific thresholds or analysis
logic were intentionally changed during this cleanup.
