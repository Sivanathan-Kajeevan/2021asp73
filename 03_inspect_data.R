library(GEOquery)
library(Biobase)

# Load the saved dataset
pineapple <- readRDS("data/pineapple_GSE38521.rds")
sample_info <- pData(pineapple)

# Print relevant metadata clearly, one field at a time
fields <- grep(
  "title|source_name|characteristics|label|data_processing|description",
  names(sample_info),
  value = TRUE,
  ignore.case = TRUE
)

for (field in fields) {
  cat("\n---", field, "---\n")
  print(sample_info[[field]])
}

# Download one individual GEO sample
sample_one <- getGEO(
  "GSM944528",
  destdir = "data"
)

# Show how GEO defines the sample table columns
print(Columns(sample_one))

# Check missing values in each array
expression_data <- exprs(pineapple)
colSums(is.na(expression_data))

# Count probes with 0, 1, 2, or 3 available measurements
table(
  available_measurements = rowSums(!is.na(expression_data))
)