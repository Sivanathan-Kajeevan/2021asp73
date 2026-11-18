 
ANALYSIS OF GENE EXPRESSION CHANGES DURING PINEAPPLE FRUIT RIPENING USING R
======================================================================

1. PROJECT OVERVIEW
-------------------
This project compares expression signals between mature green and mature
yellow pineapple fruit using the public microarray dataset GSE38521.
The analysis was performed in R using RStudio.

Research question:
Which measured pineapple clones show expression differences between mature
green and mature yellow fruit?

A clone is a microarray measurement target. Several clones can represent
the same gene, so clone counts must not be described as unique gene counts.

2. OBJECTIVES
-------------
- Download and document a real public biological dataset.
- Inspect and clean the original microarray intensity measurements.
- Apply background correction and within-array normalisation.
- Account for the dye-swap arrangement.
- Identify clone-level expression differences using limma.
- Visualise results and explain their biological and technical limitations.

3. DATA SOURCE
--------------
Repository: NCBI Gene Expression Omnibus (GEO)
Accession: GSE38521
Organism: Ananas comosus (pineapple)
Experiment: Two-colour microarray expression profiling
Comparison: Mature yellow fruit (pine9) versus mature green fruit (pine3)

Dataset URL:
https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE38521

Arrays and dye assignments:
GSM944528: Yellow = Cy3; green = Cy5
GSM944529: Yellow = Cy5; green = Cy3 (dye swap)
GSM944530: Yellow = Cy3; green = Cy5

Each array compares both stages. The three columns are three array
comparisons, not separate green and yellow sample groups.

The analysis uses the original supplementary intensity files, excluding
files labelled "2xSD+new". The initially downloaded processed GEO matrix
was inspected but was not the input to the final raw-intensity workflow.

Data provenance and reuse:
Cite the GEO accession and original publication. Check the source records
for applicable reuse terms; this README does not assign a new licence to
the source data. This is a secondary analysis, with no new field or
laboratory sampling performed for this project.

4. SOFTWARE AND PACKAGES
------------------------
Software: R and RStudio
Packages: BiocManager, GEOquery, Biobase, limma, ggplot2, pheatmap

Example package installation in the R Console:

if (!requireNamespace("BiocManager", quietly = TRUE)) {
  install.packages("BiocManager")
}
BiocManager::install(
  c("GEOquery", "Biobase", "limma"),
  ask = FALSE,
  update = FALSE
)
install.packages(c("ggplot2", "pheatmap"))

An internet connection is required for package installation and the initial
downloads. Exact R and package versions should be recorded using
sessionInfo(); they were not supplied with the reported outputs.

5. ANALYSIS WORKFLOW
--------------------
1. Download GEO metadata, the processed matrix and supplementary archive.
2. Inspect sample definitions, dye assignments and raw-file columns.
3. Import the three original intensity tables.
4. Retain NORMAL spots with valid Clone IDs and finite, non-negative
   foreground/background measurements. Exclude zero foreground values and
   values at or above the assumed 16-bit scanner ceiling of 65,535.
5. Apply limma normexp background correction with offset = 50.
6. Apply within-array loess normalisation separately to each array.
7. Orient all log-ratios as log2(yellow / green). Because limma initially
   calculates log2(Cy5 / Cy3), reverse the sign for GSM944528 and GSM944530.
8. Combine repeated spots for each Clone ID using median normalised
   log-ratios and median average log-intensities.
9. Match Clone IDs across arrays. Do not match by spot number or row order.
10. Retain clones with average log2 intensity above log2(100) in at least
    two arrays. This is an exploratory threshold, including the offset.
11. Combine supplied annotations and inspect possible control descriptions.
12. Fit a single-coefficient limma model to the three oriented log-ratios.
    Set fit$Amean to the actual mean average log-intensity and apply
    eBayes(trend = TRUE).
13. Adjust p-values across tested clones using Benjamini-Hochberg (BH).
14. Select clones with adjusted p-value < 0.05 and estimated absolute
    log2 fold change >= log2(1.5).
15. Export result tables, quality diagnostics, volcano plot and heatmap.

The statistical test assesses a difference from zero. The additional
1.5-fold screen is an effect-size selection rule; these p-values do not
test whether the true change exceeds 1.5-fold.

6. MAIN RESULTS
---------------
Original distinct Clone IDs shared across all files: 9,303
Matched clones after spot-quality processing:        9,100
Clones removed by the intensity filter:              2,645
Clones tested:                                      6,455
Clones with BH-adjusted p-value < 0.05:                 301
Clones also meeting the 1.5-fold criterion:             300
Selected clones higher in green fruit:                 233
Selected clones higher in yellow fruit:                 67

Positive log2(yellow / green): Higher expression in yellow fruit.
Negative log2(yellow / green): Higher expression in green fruit.

Selected clones higher in green account for 77.7% of the selected set;
those higher in yellow account for 22.3%.

Correlation after intensity filtering:
Arrays 1 and 2: 0.804
Arrays 1 and 3: 0.900
Arrays 2 and 3: 0.780

Log-ratio standard deviation after filtering:
Array 1: 0.966
Array 2: 1.055
Array 3: 0.722

The arrays broadly agree in pattern, but Array 3 shows smaller changes.
Correlation does not demonstrate identical scale or quality.

Example top result:
Clone ID: JBW017C07
Description: unsequenced
Mean log2(yellow / green): approximately -5.188
BH-adjusted p-value: approximately 0.000745
Yellow/green ratio: approximately 0.0274
Interpretation: The estimated signal is about 36.5 times higher in green
fruit. Its function cannot be assigned from the supplied annotation.

These numbers reflect the supplied R outputs. They are exploratory
clone-level results, not independently validated gene discoveries.

7. FIGURES AND THEIR INTERPRETATION
----------------------------------
GSM944528_MA_plots.png
GSM944529_MA_plots.png
GSM944530_MA_plots.png
  Show signal differences before and after within-array loess.
  "Before" means after background correction but before loess.
  These panels use scanner log2(Cy5 / Cy3), before stage-direction changes.

normalisation_quality_checks.png
  Shows boxplots and pairwise array agreement for the 9,100 matched clones
  before intensity filtering. Array 3 has a narrower distribution.

variance_trend.png
  Shows residual-spread behaviour across intensity values and the fitted
  limma prior trend. The vertical axis is sqrt(sigma), not fold change.

volcano_plot.png
  Each point represents one clone.
  X-axis: Mean log2(yellow / green).
  Y-axis: -log10(BH-adjusted p-value).
  Green: Selected clones higher in green fruit.
  Yellow: Selected clones higher in yellow fruit.
  Grey: Clones not meeting both selection rules.
  Dashed lines mark the adjusted p-value and fold-change cutoffs.

top_30_clone_heatmap.png
  Shows up to 30 selected clones ranked by adjusted p-value, with raw
  p-value used to break ties. Rows are not standardised.
  Green values indicate higher expression in green fruit; yellow values
  indicate higher expression in yellow fruit.
  All displayed top-ranked clones are green-associated. This selected
  subset does not represent the full set of 300 selected clones.
  The row tree groups similar expression profiles, not evolutionary history.

8. EXPECTED PROJECT STRUCTURE
-----------------------------
README.txt
code/
  01_setup.R
  02_download_data.R
  03_inspect_data.R
  04_inspect_raw.R
  05_import_raw.R
  06_normalise.R
  07_quality_checks.R
  08_filter_expression.R
  09_annotations.R
  10_differential_expression.R
  11_results_plots.R
data/
  pineapple_GSE38521.rds
  sample_metadata.csv
  GSE38521/                  Downloaded archive and extracted files
  raw_intensity_tables.rds
  normalised_pineapple.rds
  filtered_pineapple.rds
  clone_annotations.rds
  limma_fit.rds
figures/
  GSM944528_MA_plots.png
  GSM944529_MA_plots.png
  GSM944530_MA_plots.png
  normalisation_quality_checks.png
  variance_trend.png
  volcano_plot.png
  top_30_clone_heatmap.png
results/
  array_correlations.csv
  intensity_filter_decisions.csv
  clone_annotations.csv
  all_clone_results.csv
  selected_clone_results.csv
  selected_described_clones.csv
reports/
  Pineapple_Ripening_Project_Report.pdf
  Pineapple_Ripening_Results_Explained.pdf

This structure is a guide based on the step-by-step workflow. Confirm that
the files exist locally and move the two PDFs into reports/ if desired.
The setup script may currently be at the project root; place a copy in
code/ if using the structure above.

9. HOW TO RUN
-------------
1. Download or clone the complete project repository.
2. Open its RStudio project, if included, or create an RStudio project in
   the repository root. Relative paths must resolve from that root.
3. Install the required packages using the commands in Section 4.
4. Run the numbered scripts in order, from 01_setup.R through
   11_results_plots.R. Save the corrected compressed-file preview code in
   04_inspect_raw.R so it recognises .txt.gz files.
5. The first download requires internet access. GEOquery may reuse cached
   downloads on subsequent runs.
6. Review Console messages, result tables and quality plots. Do not assume
   that successfully saving a plot confirms biological validity.
7. Record the software environment from the project root:

writeLines(capture.output(sessionInfo()), "results/sessionInfo.txt")

To run one script from the R Console, for example:
source("code/10_differential_expression.R")

The complete workflow should be rerun from a fresh R session before a
claim of start-to-finish reproducibility is made. Exact package versions
and an independently checked full rerun are not provided in this README.

10. LIMITATIONS
---------------
- Only three array comparisons were available. Confirm independent
  biological sampling in the original experimental records before making
  broad population-level claims.
- Arrays 1 and 2 supplied median intensities; Array 3 supplied mean
  intensities. Standardising names does not remove this difference.
- Array 3 has a smaller log-ratio spread.
- Many top clones are unsequenced or poorly annotated. Several clones may
  correspond to the same gene, and cross-hybridisation is possible.
- Large branches remain in the MA plots. Their biological or technical
  origin has not been fully resolved.
- The intensity threshold is exploratory. Thresholds of 75, 100 and 150
  retained 8,095, 6,455 and 3,094 clones, respectively. Differential-expression
  stability across these alternatives has not yet been tested.
- No independent laboratory validation, functional enrichment analysis or
  causal experiment was completed in the reported workflow.

11. CONCLUSION
--------------
This analysis identified 300 clones meeting both the statistical and
estimated fold-change criteria when comparing mature green and mature
yellow pineapple fruit. Of these, 233 were higher in green and 67 were
higher in yellow. These expression signals are candidates for further
annotation and validation. The findings do not establish 300 unique genes
or prove which genes cause pineapple ripening.

12. FUTURE WORK
---------------
- Reconcile clone descriptions with updated supplementary annotations.
- Map reliable sequences to unique pineapple genes.
- Inspect extreme signals and potential probe artefacts.
- Rerun models under alternative intensity filters.
- Assess sensitivity to individual arrays.
- Validate promising candidates using independent biological material.

13. REFERENCES
--------------
1. NCBI Gene Expression Omnibus. GSE38521.
   https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE38521

2. Koia JH, Moyle RL, Botella JR. (2012).
   Microarray analysis of gene expression profiles in ripening pineapple
   fruits. BMC Plant Biology, 12, 240.
   https://doi.org/10.1186/1471-2229-12-240

3. Bioconductor. limma User's Guide.
   https://bioconductor.org/packages/release/bioc/vignettes/limma/inst/doc/usersguide.pdf

4. Ritchie ME et al. (2015).
   limma powers differential expression analyses for RNA-sequencing and
   microarray studies. Nucleic Acids Research, 43(7), e47.
   https://doi.org/10.1093/nar/gkv007