#This R script contains code that aims to determine a common cutoff on the output columns rendered once the --individual-counts parameter is enabled in rMATS
  #each event type has its own specific columns... the latest information about this parameter can be found here: https://github.com/Xinglab/rmats-turbo 

#this code is not yet fully developed so please feel free to use any interesting aspects of it to further analyze your rMATS data 

library(dplyr)
library(ggplot2)

# Set the base directory
base_dir <- '~/cafalempin/rmats/'

# List of genes
gene_list <- c('snapin', 'ptpn11', 'palb2', 'ifg1r',
               'sdha', 'eftud2', 'phf8', 'msh2', 'mlh1')

# Define subdirectories for analyses
output_subdirs <- c(
  'snapin/output/',
  'ptpn11/output/',
  'palb2/output/',
  'ifg1r/output/',
  'sdha/output/',
  'eftud2/output/',
  'phf8/output/',
  'msh2/output/',
  'mlh1/output/'
)

# Construct the full paths
output_directories <- file.path(base_dir, output_subdirs)
names(output_directories) <- gene_list  # Assign gene names to the paths

# All event types and their corresponding columns
event_types <- list(
  SE = c("upstream_to_target_count", "upstream_to_downstream_count"),
  MXE = c("upstream_to_first_count", "upstream_to_second_count", "first_count",
          "second_count", "first_to_downstream_count", "second_to_downstream_count"),
  A3SS = c("long_to_flanking_count", "short_to_flanking_count"),
  A5SS = c("long_to_flanking_count", "short_to_flanking_count"),
  RI = c("upstream_to_intron_count", "upstream_to_downstream_count", "intron_count")
)


### CREATING THE HELPER FUNCTIONS ###

# Function to read through the event data
read_event_data <- function(directory, event_type) {
  file_path <- file.path(directory, paste0(event_type, ".MATS.JC.txt"))
  if (file.exists(file_path)) {
    read.table(file_path, header = TRUE, sep = "\t", stringsAsFactors = FALSE)
  } else {
    warning(paste("File does not exist:", file_path))
    return(NULL)
  }
}

# Function to filter event data according to significance of p-value and FDR values
filter_event_data <- function(event_data, p_threshold = 0.05, fdr_threshold = 0.05) {
  event_data <- event_data %>% 
    mutate(new_FDR = p.adjust(PValue, method = "BH"))
  
  filtered_data <- event_data %>%
    filter(PValue < p_threshold, new_FDR < fdr_threshold)
  return(filtered_data)
}

# Function to extract the counts in each of the columns
extract_column_replicate_counts <- function(data, columns_to_be_split) {
  data <- data %>%
    rowwise() %>%
    mutate(
      coverage_meets_cutoff = sum(sapply(columns_to_be_split, function(col) {
        if (!is.na(get(col))) {
          # Split the comma-separated values into lists
          replicate_list <- strsplit(as.character(get(col)), ",")[[1]]
          replicate_values <- as.numeric(replicate_list)
          # Counting the number of replicates that meet the cutoff of at least 5
          sum(replicate_values >= 5, na.rm = TRUE)
        } else {
          # If the column has NA, return 0
          0
        }
      }))
    ) %>%
    ungroup() %>%
    filter(coverage_meets_cutoff >= 2) %>%  # Adjust the cutoff value as needed
    select(-coverage_meets_cutoff)
  return(data)
}

# Function to create histogram plots
create_histogram <- function(data, column, gene_name, event_type) {
  replicate_values <- unlist(strsplit(as.character(data[[column]]), ","))
  replicate_values <- as.numeric(replicate_values)
  
  p <- ggplot() + 
    geom_histogram(aes(x = replicate_values), binwidth = 1) +
    labs(title = paste("Histogram of", column, "for gene", gene_name, "event type", event_type),
         x = "Replicate values", y = "Frequency") +
    theme_minimal() +
    xlim(0, 300) +
    ylim(0, 500)
  
  return(p)
}

### PROCESSING THE DATA ###

# Initialize lists to store results for each gene and event type
all_data <- list()
all_plots <- list()

for (gene in gene_list) {
  dir <- output_directories[[gene]]
  if (!dir.exists(dir)) {
    warning(paste("Directory does not exist:", dir))
    next
  }
  
  all_data[[gene]] <- list()
  all_plots[[gene]] <- list()
  
  for (event_type in names(event_types)) {
    event_data <- read_event_data(dir, event_type)
    if (is.null(event_data)) {
      next
    }
    filtered_data <- filter_event_data(event_data)
    
    # Prepare data for each event type
    columns_to_be_split <- event_types[[event_type]]
    processed_data <- extract_column_replicate_counts(filtered_data, columns_to_be_split)
    
    # Store processed data
    all_data[[gene]][[event_type]] <- processed_data
    
    # Store plots
    all_plots[[gene]][[event_type]] <- list()
    for (column in columns_to_be_split) {
      p <- create_histogram(processed_data, column, gene, event_type)
      all_plots[[gene]][[event_type]][[column]] <- p
    }
  }
}

### EXAMPLE USAGE ###

# View the processed data for a specific gene and event type
gene_of_interest <- 'snapin'
event_type_of_interest <- 'SE'

if (!is.null(all_data[[gene_of_interest]][[event_type_of_interest]])) {
  View(all_data[[gene_of_interest]][[event_type_of_interest]])
}

# Print the plots for a specific gene, event type, and column
columns_to_plot <- event_types[[event_type_of_interest]]

for (column in columns_to_plot) {
  if (!is.null(all_plots[[gene_of_interest]][[event_type_of_interest]][[column]])) {
    print(all_plots[[gene_of_interest]][[event_type_of_interest]][[column]])
  }
}

