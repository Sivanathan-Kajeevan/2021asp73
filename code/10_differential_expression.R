library(limma)

# Load expression data and annotations
filtered <- readRDS("data/filtered_pineapple.rds")
annotation <- readRDS("data/clone_annotations.rds")

M <- filtered$M
A <- filtered$A

stopifnot(
  identical(rownames(M), annotation$Clone_ID),
  identical(dimnames(M), dimnames(A)),
  all(is.finite(M)),
  all(is.finite(A))
)

# Each column already contains log2(yellow / green).
# A single coefficient estimates the average stage difference.
design <- matrix(
  1,
  nrow = ncol(M),
  ncol = 1,
  dimnames = list(colnames(M), "Yellow_vs_Green")
)

# Fit clone-level models
fit <- lmFit(M, design)

# Supply actual average intensities for intensity-dependent moderation
fit$Amean <- rowMeans(A)

fit <- eBayes(fit, trend = TRUE)

# Retrieve all results with Benjamini-Hochberg adjustment
results <- topTable(
  fit,
  coef = "Yellow_vs_Green",
  number = Inf,
  adjust.method = "BH",
  sort.by = "P"
)

# Attach annotations in the correct order
results$Clone_ID <- rownames(results)

matched <- match(results$Clone_ID, annotation$Clone_ID)

results$Contig_ID <- annotation$Contig_ID[matched]
results$Gene_name <- annotation$Gene_name[matched]

# Yellow-to-green expression ratio
results$Yellow_Green_ratio <- 2^results$logFC

# Statistical significance and an additional effect-size screen
results$FDR_significant <- results$adj.P.Val < 0.05

results$Selected <- results$FDR_significant &
  abs(results$logFC) >= log2(1.5)

results$Direction <- ifelse(
  results$logFC > 0,
  "Higher in yellow",
  ifelse(results$logFC < 0, "Higher in green", "No change")
)

# Save results
dir.create("results", showWarnings = FALSE)

write.csv(
  results,
  "results/all_clone_results.csv",
  row.names = FALSE
)

write.csv(
  results[results$Selected, ],
  "results/selected_clone_results.csv",
  row.names = FALSE
)

saveRDS(fit, "data/limma_fit.rds")

# Save the variance-trend diagnostic
png(
  "figures/variance_trend.png",
  width = 1000, height = 800, res = 120
)
plotSA(fit, main = "Variance trend across clones")
dev.off()

# Print a summary
cat("\nClones tested:", nrow(results), "\n")

cat(
  "FDR < 0.05:",
  sum(results$FDR_significant), "\n"
)

cat(
  "FDR < 0.05 and at least 1.5-fold change:",
  sum(results$Selected), "\n"
)

cat("\nSelected clones by direction:\n")
print(table(results$Direction[results$Selected]))

cat("\nTop 10 results:\n")
print(
  head(results[, c(
    "Clone_ID", "Contig_ID", "Gene_name",
    "logFC", "P.Value", "adj.P.Val",
    "Yellow_Green_ratio"
  )], 10),
  row.names = FALSE
)