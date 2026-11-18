# Load normalised log-ratios and average intensities
normalised_data <- readRDS(
  "data/normalised_pineapple.rds"
)

M <- normalised_data$M
A <- normalised_data$A

stopifnot(
  identical(dimnames(M), dimnames(A)),
  all(is.finite(M)),
  all(is.finite(A))
)

# Exploratory intensity threshold.
# Intensities include the offset of 50 used in background correction.
intensity_cutoff <- log2(100)

# Require sufficient intensity in at least two arrays
keep <- rowSums(A > intensity_cutoff) >= 2

M_filtered <- M[keep, , drop = FALSE]
A_filtered <- A[keep, , drop = FALSE]

# Save all filtering decisions
filter_details <- data.frame(
  Clone_ID = rownames(M),
  Arrays_above_threshold = rowSums(A > intensity_cutoff),
  Retained = keep
)

dir.create("results", showWarnings = FALSE)

write.csv(
  filter_details,
  "results/intensity_filter_decisions.csv",
  row.names = FALSE
)

saveRDS(
  list(
    M = M_filtered,
    A = A_filtered,
    intensity_cutoff = intensity_cutoff,
    minimum_arrays = 2
  ),
  "data/filtered_pineapple.rds"
)

cat("Clones before filtering:", nrow(M), "\n")
cat("Clones retained:", nrow(M_filtered), "\n")
cat("Clones removed:", sum(!keep), "\n")

cat("\nCorrelations after filtering:\n")
print(round(cor(M_filtered), 3))

cat("\nLog-ratio spread per array (standard deviation):\n")
print(round(apply(M_filtered, 2, sd), 3))

# Show how alternative thresholds affect the retained count
cat("\nThreshold sensitivity — number of retained clones:\n")

for (threshold in c(75, 100, 150)) {
  retained <- rowSums(A > log2(threshold)) >= 2
  cat("Intensity >", threshold, ":", sum(retained), "\n")
}