library(limma)

raw_tables <- readRDS("data/raw_intensity_tables.rds")

dir.create("figures", showWarnings = FALSE)

# Direction for log2(yellow / green).
# limma calculates log2(R / G), where:
# R = W685 (Cy5), G = W595 (Cy3).
direction <- c(
  GSM944528 = -1,
  GSM944529 =  1,
  GSM944530 = -1
)

normalised <- list()

for (sample in names(raw_tables)) {
  
  x <- raw_tables[[sample]]
  x[["Clone ID"]] <- trimws(x[["Clone ID"]])
  
  intensity_columns <- c(
    "Spot Intensity (W595)",
    "Spot Intensity (W685)",
    "Background Intensity (W595)",
    "Background Intensity (W685)"
  )
  
  # Keep normal spots with valid IDs and numerical measurements
  valid <- trimws(x[["Flagged As"]]) == "NORMAL" &
    !is.na(x[["Clone ID"]]) &
    nzchar(x[["Clone ID"]])
  
  for (column in intensity_columns) {
    valid <- valid &
      is.finite(x[[column]]) &
      x[[column]] >= 0
  }
  
  # Exclude zero foreground and values at the 16-bit scanner ceiling
  valid <- valid &
    x[["Spot Intensity (W595)"]] > 0 &
    x[["Spot Intensity (W685)"]] > 0 &
    x[["Spot Intensity (W595)"]] < 65535 &
    x[["Spot Intensity (W685)"]] < 65535
  
  valid[is.na(valid)] <- FALSE
  x <- x[valid, ]
  
  # Build a two-colour intensity object for this slide
  RG <- new("RGList", list(
    R = matrix(x[["Spot Intensity (W685)"]], ncol = 1),
    G = matrix(x[["Spot Intensity (W595)"]], ncol = 1),
    Rb = matrix(x[["Background Intensity (W685)"]], ncol = 1),
    Gb = matrix(x[["Background Intensity (W595)"]], ncol = 1)
  ))
  
  # Background correction produces positive adjusted intensities
  RG_corrected <- backgroundCorrect(
    RG,
    method = "normexp",
    offset = 50
  )
  
  # Ratios before and after within-slide normalisation
  M_before <- log2(RG_corrected$R / RG_corrected$G)
  A_before <- (
    log2(RG_corrected$R) + log2(RG_corrected$G)
  ) / 2
  
  MA <- normalizeWithinArrays(
    RG_corrected,
    method = "loess"
  )
  
  # Save diagnostic plots in the scanner's R/G direction
  png(
    file.path("figures", paste0(sample, "_MA_plots.png")),
    width = 1400, height = 650, res = 120
  )
  
  par(mfrow = c(1, 2))
  
  plot(
    A_before[, 1], M_before[, 1],
    pch = 16, cex = 0.25,
    xlab = "Average log2 intensity",
    ylab = "log2(Cy5 / Cy3)",
    main = paste(sample, "before loess")
  )
  abline(h = 0, col = "red")
  
  plot(
    MA$A[, 1], MA$M[, 1],
    pch = 16, cex = 0.25,
    xlab = "Average log2 intensity",
    ylab = "Normalised log2(Cy5 / Cy3)",
    main = paste(sample, "after loess")
  )
  abline(h = 0, col = "red")
  
  dev.off()
  
  # Orient every slide as yellow / green
  M_oriented <- MA$M[, 1] * direction[[sample]]
  
  # Combine repeated spots using their median normalised log-ratio
  clone_M <- tapply(
    M_oriented, x[["Clone ID"]], median
  )
  
  clone_A <- tapply(
    MA$A[, 1], x[["Clone ID"]], median
  )
  
  normalised[[sample]] <- list(
    M = clone_M,
    A = clone_A
  )
  
  cat(
    "\n", sample,
    ": retained", nrow(x), "spots;",
    length(clone_M), "clones\n"
  )
}

# Match clones by ID across the three slides
common_clones <- Reduce(
  intersect,
  lapply(normalised, function(x) names(x$M))
)

M_matrix <- do.call(cbind, lapply(
  normalised, function(x) x$M[common_clones]
))

A_matrix <- do.call(cbind, lapply(
  normalised, function(x) x$A[common_clones]
))

rownames(M_matrix) <- common_clones
rownames(A_matrix) <- common_clones

saveRDS(
  list(
    M = M_matrix,
    A = A_matrix,
    slides = normalised
  ),
  "data/normalised_pineapple.rds"
)

cat("\nFinal matrix dimensions:\n")
print(dim(M_matrix))

cat("\nNormalised log-ratio summary:\n")
print(summary(as.vector(M_matrix)))

cat("\nCorrelation between slides:\n")
print(round(cor(M_matrix), 3))