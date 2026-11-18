# Create project folders
folders <- c("code", "data", "figures", "results")

for (folder in folders) {
  dir.create(folder, showWarnings = FALSE)
}

# Install the Bioconductor package manager
if (!requireNamespace("BiocManager", quietly = TRUE)) {
  install.packages("BiocManager")
}

# Install packages for downloading and analysing expression data
BiocManager::install(
  c("GEOquery", "limma"),
  ask = FALSE,
  update = FALSE
)

# Install plotting packages
install.packages(c("ggplot2", "pheatmap"))

# Check that packages load successfully
library(GEOquery)
library(limma)
library(ggplot2)
library(pheatmap)

# Display your R version
R.version.string

message("Step 1 completed: project and packages are ready.")
