# Multi-cohort landscape of dysregulated lactylation-related genes in heart failure

Companion code for: *Multi-cohort landscape of dysregulated lactylation-related genes in heart failure: chromatin regulation, immune remodeling, and Mendelian randomization appraisal of causality.*

## Study design

- **Discovery**: GSE57338 (177 HF / 136 non-failing donor myocardium, Affymetrix HuGene 1.1 ST)
- **Replication**: GSE5406 (194 HF / 16 NF, HG-U133A) and GSE116250 (50 HF / 14 NF, RNA-seq)
- **Analyses**: limma differential expression of 332 curated lactylation-related genes (LRGs) → 13-gene LASSO discriminative signature with four validation scenarios → WGCNA co-expression → GO/GSEA/KEGG enrichment → CIBERSORT/MCP-counter/EPIC immune deconvolution with lineage-level concordance → two-sample Mendelian randomization (eQTLGen blood cis-eQTL × HERMES HF GWAS, 47,309 cases / 930,014 controls) with GTEx left-ventricle tissue-source sensitivity, nested cross-validation, bootstrap selection stability, and Meff-adjusted multiple testing; GTEx v8 left-ventricle expression boundary of the LRG set; donor-inflammation-stratified immune sensitivity analysis.

## Repository layout

```
R_scripts/    Full analysis chain, one script per stage (run in order):
              M0   environment setup (package dependencies)
              M1   stages 1-3: GEO parsing + preprocessing + limma (GSE57338)
              M2   stage 4: LASSO + four-scenario ROC validation
              M3   stage 5: WGCNA + module enrichment
              M4   stage 7: functional enrichment + immune infiltration
              M5   stage 8: immune three-method sensitivity (MCP-counter/EPIC/subset)
              M6   stage 6: two-sample MR (31 candidate genes)
              M7   GTEx left-ventricle cardiac-eQTL MR sensitivity (local reproduction)
              M8   Meff-adjusted multiple-testing correction (26 genes)
              M9   nested CV + bootstrap stability + DeLong CI summary
              M10  optional: local KEGG rerun (offline fallback for M4)
              M11  GTEx v8 left-ventricle expression boundary (Suppl. Table S8)
              M12  donor-inflammation-stratified immune sensitivity (Suppl. Table S9 / Fig. S1)
              Appendix1_replot_figures.R   optional: base versions of Figures 4B/5C/7
              CIBERSORT_v1.04.R           official CIBERSORT source (Newman et al. 2015)
              opengwas_config.R           OpenGWAS API configuration (MR modules)
data/         lrg332.txt (332 curated LRG symbols) + LM22.txt (CIBERSORT reference
              matrix); the GEO series-matrix files are downloaded from GEO per script M1
results/      key result tables reproduced by the pipeline (Table 3 and MR
              supplementary tables), CSV format
```

## Requirements

- R 4.6.1 (all numerical results in the paper were produced with this version)
- Key packages: `limma`, `glmnet`, `WGCNA`, `pROC`, `clusterProfiler`, `TwoSampleMR` (see `M0_environment_setup.R` for the full installation list)
- Public summary statistics: eQTLGen cis-eQTL; HERMES heart-failure GWAS (GWAS Catalog GCST009541); eQTL Catalogue (imported GTEx v8, heart left ventricle); 1000 Genomes phase 3 EUR (LD reference)
- An OpenGWAS token (free registration at https://api.opengwas.io) for the MR modules

## Usage

Each script is self-contained and idempotent. Set the `DIR` variable at the top of each script to your working directory, then run in order (M0 → M12). Every script prints the **expected values** (numbers reported in the paper) at the end so that a successful reproduction can be verified at a glance. The MR script caches per-gene results and resumes automatically after an interruption.

## Citation

If you use this code, please cite the paper and the underlying public resources (GEO, eQTLGen, HERMES/GWAS Catalog, eQTL Catalogue, GTEx).

## License

MIT
