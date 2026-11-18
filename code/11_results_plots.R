library(ggplot2)
library(pheatmap)

results <- read.csv(
  "results/all_clone_results.csv",
  stringsAsFactors = FALSE
)

filtered <- readRDS("data/filtered_pineapple.rds")

dir.create("figures", showWarnings = FALSE)

# Recalculate selection explicitly
results$Selected <- (
  results$adj.P.Val < 0.05 &
    abs(results$logFC) >= log2(1.5)
)

results$Plot_group <- "Other"

results$Plot_group[
  results$Selected & results$logFC > 0
] <- "Higher in yellow"

results$Plot_group[
  results$Selected & results$logFC < 0
] <- "Higher in green"

# Avoid infinite plotting values if a p-value is zero
results$Minus_log10_FDR <- -log10(
  pmax(results$adj.P.Val, .Machine$double.xmin)
)

# Volcano plot: adjusted p-values on the vertical axis
volcano <- ggplot(
  results,
  aes(
    x = logFC,
    y = Minus_log10_FDR,
    colour = Plot_group
  )
) +
  geom_point(alpha = 0.65, size = 1.5) +
  geom_vline(
    xintercept = c(-log2(1.5), log2(1.5)),
    linetype = "dashed"
  ) +
  geom_hline(
    yintercept = -log10(0.05),
    linetype = "dashed"
  ) +
  scale_colour_manual(values = c(
    "Higher in green" = "#16856B",
    "Higher in yellow" = "#D59400",
    "Other" = "grey75"
  )) +
  labs(
    title = "Pineapple fruit ripening: clone-level changes",
    x = "Mean log2(yellow / green)",
    y = "-log10(BH-adjusted p-value)",
    colour = NULL
  ) +
  theme_minimal(base_size = 12) +
  theme(legend.position = "bottom")

ggsave(
  "figures/volcano_plot.png",
  plot = volcano,
  width = 8,
  height = 6,
  dpi = 300
)

# Select up to 30 clones with the smallest adjusted p-values
selected <- results[results$Selected, ]
selected <- selected[
  order(selected$adj.P.Val, selected$P.Value),
]

if (nrow(selected) >= 2) {
  top_clones <- head(selected, 30)
  
  heatmap_matrix <- filtered$M[
    match(top_clones$Clone_ID, rownames(filtered$M)),
    ,
    drop = FALSE
  ]
  
  # Keep actual log-ratios rather than standardising each row
  colour_limit <- max(abs(heatmap_matrix))
  
  pheatmap(
    heatmap_matrix,
    scale = "none",
    cluster_rows = TRUE,
    cluster_cols = FALSE,
    labels_col = c("Array 1", "Array 2", "Array 3"),
    color = colorRampPalette(
      c("#16856B", "white", "#D59400")
    )(100),
    breaks = seq(
      -colour_limit, colour_limit,
      length.out = 101
    ),
    main = "Top selected clones: log2(yellow / green)",
    fontsize_row = 7,
    border_color = NA,
    filename = "figures/top_30_clone_heatmap.png",
    width = 8,
    height = 10
  )
}

# Export selected clones with potentially informative descriptions.
# This screens annotation labels; it does not validate gene functions.
uninformative <- (
  is.na(selected$Gene_name) |
    !nzchar(trimws(selected$Gene_name)) |
    grepl(
      "unsequenced|undiscovered|unknown|poor seq|<150",
      selected$Gene_name,
      ignore.case = TRUE
    )
)

described <- selected[!uninformative, ]

write.csv(
  described,
  "results/selected_described_clones.csv",
  row.names = FALSE
)

cat("\nSelected clones with informative descriptions:",
    nrow(described), "\n")

cat("\nTop described results:\n")
print(
  head(described[, c(
    "Clone_ID", "Contig_ID", "Gene_name",
    "logFC", "adj.P.Val"
  )], 15),
  row.names = FALSE
)

cat("\nPlots saved in the figures folder.\n")
