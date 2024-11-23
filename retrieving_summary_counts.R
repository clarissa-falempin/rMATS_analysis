# In this script we are going to retrieve the summary counts data for each gene. the original output was rendered using the "incorrectly" merged GTF MANE clinical annotations, and the second run is with the correct GTF file

base_dir <- '~/cafalempin/rmats'

# List of genes
gene_list <- c('snapin', 'ptpn11', 'palb2', 'ifg1r', 'sdha', 'eftud2', 'phf8', 'msh2', 'mlh1')

# Define subdirectories for 'before' and 'after' analyses
output_directories_before <- c(
  'first_run/snapin/output/',
  'first_run/ptpn11/output/',
  'first_run/palb2/output/',
  'first_run/ifg1r/output/',
  'first_run/sdha/output/',
  'first_run/eftud2/output/',
  'first_run/phf8/output/',
  'first_run/msh2/output/',
  'first_run/mlh1/output/'
)

output_directories_after <- c(
  'second_run/snapin/output/',
  'second_run/ptpn11/output/',
  'second_run/palb2/output/',
  'second_run/ifg1r/output/',
  'second_run/sdha/output/',
  'second_run/eftud2/output/',
  'second_run/phf8/output/',
  'second_run/msh2/output/',
  'second_run/mlh1/output/'
)

# Construct the full paths for 'before' and 'after' directories
output_dirs_before <- file.path(base_dir, output_directories_before)
output_dirs_after <- file.path(base_dir, output_directories_after)

# Function to preview the first few lines of a file
preview_file <- function(file_path, num_lines = 6) {
  if (file.exists(file_path)) {
    lines <- readLines(file_path, n = num_lines)
    cat(lines, sep = "\n")
  } else {
    warning(paste("File does not exist:", file_path))
  }
}

# Preview summary files from 'before' directories
for (i in seq_along(output_dirs_before)) {
  dir_path <- output_dirs_before[i]
  gene <- gene_list[i]
  summary_file <- file.path(dir_path, "summary.txt")
  cat("Previewing 'before' file for gene:", gene, "\n")
  cat("File path:", summary_file, "\n")
  preview_file(summary_file)
  cat("\n")
}

# Preview summary files from 'after' directories
for (i in seq_along(output_dirs_after)) {
  dir_path <- output_dirs_after[i]
  gene <- gene_list[i]
  summary_file <- file.path(dir_path, "summary.txt")
  cat("Previewing 'after' file for gene:", gene, "\n")
  cat("File path:", summary_file, "\n")
  preview_file(summary_file)
  cat("\n")
}

# Function to read event data from a directory
read_event_data <- function(dir_path) {
  file_path <- file.path(dir_path, "summary.txt")
  if (file.exists(file_path)) {
    data <- read.table(file_path, header = TRUE, sep = "\t", stringsAsFactors = FALSE)
    return(data)
  } else {
    warning(paste("File does not exist:", file_path))
    return(NULL)
  }
}

# Function to compare before and after data
compare_event_data <- function(event_data_before, event_data_after, columns_to_compare) {
  if (is.null(event_data_before) || is.null(event_data_after)) {
    return(NULL)
  }
  
  # Ensure both datasets contain the columns to compare
  missing_columns_before <- setdiff(columns_to_compare, colnames(event_data_before))
  missing_columns_after <- setdiff(columns_to_compare, colnames(event_data_after))
  
  if (length(missing_columns_before) > 0) {
    warning(paste("Missing columns in before data:", paste(missing_columns_before, collapse=", ")))
    return(NULL)
  }
  if (length(missing_columns_after) > 0) {
    warning(paste("Missing columns in after data:", paste(missing_columns_after, collapse=", ")))
    return(NULL)
  }
  
  # Ensure both datasets have the same number of rows
  if (nrow(event_data_before) != nrow(event_data_after)) {
    warning("Number of rows in before and after data do not match")
    return(NULL)
  }
  
  comparison <- data.frame(RowIndex = 1:nrow(event_data_before))
  comparison$EventType <- event_data_before$EventType
  comparison$EventTypeDescription <- event_data_before$EventTypeDescription
  
  for (column in columns_to_compare) {
    comparison[[paste0(column, "_before")]] <- event_data_before[[column]]
    comparison[[paste0(column, "_after")]] <- event_data_after[[column]]
    comparison[[paste0(column, "_diff")]] <- event_data_after[[column]] - event_data_before[[column]]
  }
  
  return(comparison)
}

# Columns to compare
columns_to_compare <- c(
  "TotalEventsJC", "TotalEventsJCEC", "SignificantEventsJC", 
  "SigEventsJCSample1HigherInclusion", "SigEventsJCSample2HigherInclusion", 
  "SignificantEventsJCEC", "SigEventsJCECSample1HigherInclusion", "SigEventsJCECSample2HigherInclusion"
)

# Initialize a list to store comparison results
summary_comparison <- list()

# Loop over each gene and compare data
for (i in seq_along(output_dirs_before)) {
  dir_before <- output_dirs_before[i]
  dir_after <- output_dirs_after[i]
  
  gene_name <- gene_list[i]
  
  event_data_before <- read_event_data(dir_before)
  event_data_after <- read_event_data(dir_after)
  
  if (is.null(event_data_before) || is.null(event_data_after)) {
    warning(paste("Event data missing for gene:", gene_name))
    next
  }
  
  comparison_result <- compare_event_data(event_data_before, event_data_after, columns_to_compare)
  
  if (!is.null(comparison_result)) {
    summary_comparison[[gene_name]] <- comparison_result
  }
}

# Display all data tables in summary_comparison
for (gene_name in names(summary_comparison)) {
  cat("\nGene:", gene_name, "\n")
  print(summary_comparison[[gene_name]])
}

# Optionally, view the comparison results for each gene
for (gene_name in gene_list) {
  if (!is.null(summary_comparison[[gene_name]])) {
    View(summary_comparison[[gene_name]], title = paste("Comparison for", gene_name))
  }
}