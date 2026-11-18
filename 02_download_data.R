library(GEOquery)
library(Biobase)

dir.create("data", showWarnings = FALSE)

# Download processed expression data and sample information
geo_data <- getGEO(
  "GSE38521",
  GSEMatrix = TRUE,
  getGPL = FALSE,
  destdir = "data"
)

# Check how many expression datasets were returned
length(geo_data)

# Select the first dataset
pineapple <- geo_data[[1]]

# Save a local copy for future steps
saveRDS(pineapple, "data/pineapple_GSE38521.rds")

# Extract expression values and sample information
expression_data <- exprs(pineapple)
sample_info <- pData(pineapple)

# Inspect dimensions and sample names
dim(expression_data)
colnames(expression_data)

# Inspect the first few expression values
head(expression_data)

# Check the numerical scale and missing values
summary(as.vector(expression_data))
sum(is.na(expression_data))

# Display metadata fields relevant to the experiment
metadata_fields <- grep(
  "title|source_name|characteristics|label|data_processing",
  names(sample_info),
  value = TRUE,
  ignore.case = TRUE
)

print(sample_info[, metadata_fields, drop = FALSE])

# Export sample information for reference
write.csv(
  sample_info,
  "data/sample_metadata.csv",
  row.names = TRUE
)

