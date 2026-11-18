# Locate original files, excluding the filtered 2xSD versions
raw_folder <- file.path("data", "GSE38521", "extracted")

files <- sort(list.files(
  raw_folder,
  pattern = "_raw\\.txt\\.gz$",
  recursive = TRUE,
  full.names = TRUE
))

files <- files[
  !grepl("2xSD", basename(files), fixed = TRUE)
]

stopifnot(length(files) == 3)

# Read the compressed tables
read_raw <- function(file) {
  connection <- gzfile(file, open = "rt")
  on.exit(close(connection))
  
  read.delim(
    connection,
    check.names = FALSE,
    stringsAsFactors = FALSE,
    quote = "",
    comment.char = ""
  )
}

raw_tables <- lapply(files, read_raw)

names(raw_tables) <- sub(
  "_.*", "", basename(files)
)

# Standardise the intensity column names
for (i in seq_along(raw_tables)) {
  x <- raw_tables[[i]]
  
  names(x) <- sub(
    "Median Intensity", "Intensity",
    names(x), fixed = TRUE
  )
  
  names(x) <- sub(
    "Mean Intensity", "Intensity",
    names(x), fixed = TRUE
  )
  
  raw_tables[[i]] <- x
}

# Check each slide
for (sample in names(raw_tables)) {
  x <- raw_tables[[sample]]
  
  cat("\n---", sample, "---\n")
  
  cat("Number of spots:", nrow(x), "\n")
  cat("Unique Clone IDs:",
      length(unique(x[["Clone ID"]])), "\n")
  
  cat("\nSpot flags:\n")
  print(table(x[["Flagged As"]], useNA = "ifany"))
  
  cat("\nIntensity summary:\n")
  print(summary(x[, c(
    "Spot Intensity (W595)",
    "Spot Intensity (W685)",
    "Background Intensity (W595)",
    "Background Intensity (W685)"
  )]))
}

# Check the number of Clone IDs shared across all slides
clone_sets <- lapply(raw_tables, function(x) {
  ids <- trimws(x[["Clone ID"]])
  unique(ids[!is.na(ids) & nzchar(ids)])
})

shared_clones <- Reduce(intersect, clone_sets)

cat("\nClone IDs shared across all three slides:",
    length(shared_clones), "\n")

# Save imported data without modifying the original files
saveRDS(raw_tables, "data/raw_intensity_tables.rds")