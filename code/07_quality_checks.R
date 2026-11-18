# Load normalised data
pineapple_normalised <- readRDS(
  "data/normalised_pineapple.rds"
)

M <- pineapple_normalised$M
A <- pineapple_normalised$A

stopifnot(
  all(is.finite(M)),
  all(is.finite(A))
)

dir.create("figures", showWarnings = FALSE)
dir.create("results", showWarnings = FALSE)

# Save correlations
write.csv(
  cor(M),
  "results/array_correlations.csv"
)

# Create one figure containing four quality checks
png(
  "figures/normalisation_quality_checks.png",
  width = 1600,
  height = 1200,
  res = 150
)

par(mfrow = c(2, 2), mar = c(5, 5, 3, 1))

# 1. Distribution of log-ratios
boxplot(
  M,
  names = c("Array 1", "Array 2", "Array 3"),
  col = c("#36A6A0", "#E9A23B", "#8172B3"),
  outline = FALSE,
  ylab = "Normalised log2(yellow / green)",
  main = "Log-ratio distributions"
)
abline(h = 0, col = "red", lty = 2)

# 2–4. Compare every pair of arrays
pairs_to_plot <- list(c(1, 2), c(1, 3), c(2, 3))

axis_limits <- range(M)

for (pair in pairs_to_plot) {
  i <- pair[1]
  j <- pair[2]
  
  smoothScatter(
    M[, i], M[, j],
    xlim = axis_limits,
    ylim = axis_limits,
    xlab = paste("Array", i, "log2(yellow / green)"),
    ylab = paste("Array", j, "log2(yellow / green)"),
    main = paste(
      "Arrays", i, "and", j,
      "| r =", round(cor(M[, i], M[, j]), 3)
    )
  )
  
  abline(a = 0, b = 1, col = "red", lty = 2)
}

dev.off()

cat(
  "\nSaved: figures/normalisation_quality_checks.png\n"
)