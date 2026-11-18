library(GEOquery)

# Download supplementary files
getGEOSuppFiles(
  "GSE38521",
  baseDir = "data",
  makeDirectory = TRUE
)

# Find the downloaded archive
raw_folder <- file.path("data", "GSE38521")

archives <- list.files(
  raw_folder,
  pattern = "\\.tar$",
  full.names = TRUE
)

if (length(archives) != 1) {
  stop("Expected one TAR archive. Check the downloaded files.")
}

# Extract the files
extract_folder <- file.path(raw_folder, "extracted")
dir.create(extract_folder, showWarnings = FALSE)

untar(archives[1], exdir = extract_folder)

# Display the filenames
raw_files <- list.files(
  extract_folder,
  recursive = TRUE,
  full.names = TRUE
)

cat("\n--- Extracted files ---\n")
print(basename(raw_files))

# Inspect the first lines of each text file
text_files <- raw_files[
  grepl("\\.txt$", raw_files, ignore.case = TRUE)
]

for (file in text_files) {
  cat("\n---", basename(file), "---\n")
  cat(readLines(file, n = 8, warn = FALSE), sep = "\n")
  cat("\n")
}

# Find compressed text files
text_files <- list.files(
  file.path("data", "GSE38521", "extracted"),
  pattern = "\\.txt\\.gz$",
  recursive = TRUE,
  full.names = TRUE
)

# Read the first 12 lines of each compressed file
for (file in text_files) {
  cat("\n---", basename(file), "---\n")
  
  connection <- gzfile(file, open = "rt")
  preview <- readLines(connection, n = 12, warn = FALSE)
  close(connection)
  l
  cat(preview, sep = "\n")
  cat("\n")
}