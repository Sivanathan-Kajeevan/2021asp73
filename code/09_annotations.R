# Load filtered data and original annotations
filtered <- readRDS("data/filtered_pineapple.rds")
raw_tables <- readRDS("data/raw_intensity_tables.rds")

# Collect annotations from all three slides
annotation_rows <- do.call(
  rbind,
  lapply(raw_tables, function(x) {
    data.frame(
      Clone_ID = trimws(x[["Clone ID"]]),
      Contig_ID = trimws(x[["UniGene cluster ID"]]),
      Gene_name = trimws(x[["Gene name"]]),
      stringsAsFactors = FALSE
    )
  })
)

# Combine distinct annotations for each clone
combine_labels <- function(values) {
  values <- unique(values[
    !is.na(values) & nzchar(values)
  ])
  
  paste(sort(values), collapse = "; ")
}

annotation <- aggregate(
  annotation_rows[, c("Contig_ID", "Gene_name")],
  by = list(Clone_ID = annotation_rows$Clone_ID),
  FUN = combine_labels
)

# Match annotation order to the filtered expression matrix
annotation <- annotation[
  match(rownames(filtered$M), annotation$Clone_ID),
]

stopifnot(
  identical(annotation$Clone_ID, rownames(filtered$M))
)

# Flag possible controls for inspection.
# This is a screening rule, not automatic removal.
annotation$Possible_control <- grepl(
  "control|blank|empty|spike| Arabidopsis|luciferase|GFP",
  paste(
    annotation$Clone_ID,
    annotation$Contig_ID,
    annotation$Gene_name
  ),
  ignore.case = TRUE
)

write.csv(
  annotation,
  "results/clone_annotations.csv",
  row.names = FALSE
)

saveRDS(
  annotation,
  "data/clone_annotations.rds"
)

cat("\nFirst 15 annotations:\n")
print(head(annotation, 15), row.names = FALSE)

cat("\nPossible controls:\n")
print(
  annotation[annotation$Possible_control, ],
  row.names = FALSE
)

cat("\nMost frequent gene descriptions:\n")
print(head(
  sort(table(annotation$Gene_name), decreasing = TRUE),
  15
))