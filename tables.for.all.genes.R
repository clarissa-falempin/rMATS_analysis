#This R script is possibly the main and more important one of this directory
#it extracts, wrangles and parses the rMATS output data files and merges them into one big table onto which a series of cutoffs are applied
#first a cutoff on the count coverage from the RNA-seq data is applied (where at least one of the read count columns needs to be equal to 20 or above, otherwise the event is filtered out)
#new FDR values are determine from the p-values using the Benjamini Hochberg method and a cutoff on significance is applied (PValue < 0.05 and FDR < 0.05)
library(dplyr)
library(tidyr)
library(ggplot2)

# List of the gene directories to the rMATS data (including output ran with the MANE clinical selects and  Ensembl GTF files 
#The output from individual counts data is also present (this was out of curiosity to see what exactly this parameter does from the most recent rMATS turbo release March 2024)

# Set the base directory
base_dir <- '~/cafalempin/rmats/'

# With either Ensembl or MANE clinical annotations
output_directories <- c(
  'snapin/ensembl/output/',
  'ptpn11/ensembl/output/',
  'palb2/ensembl/output/',
  'ifg1r/ensembl/output/',
  'sdha/ensembl/output/',
  'eftud2/ensembl/output/',
  'phf8/ensembl/output/',
  'msh2/ensembl/output/',
  'mlh1/ensembl/output/',
  
  'snapin/MANE/output/',
  'ptpn11/MANE/output/',
  'palb2/MANE/output/',
  'ifg1r/MANE/output/',
  'sdha/MANE/output/',
  'eftud2/MANE/output/',
  'phf8/MANE/output/',
  'msh2/MANE/output/',
  'mlh1/MANE/output/'
)

# All event types
event_types <- c("SE", "MXE", "A3SS", "A5SS", "RI")

###DEFINING THE HELPER FUNCTIONS
# Function to read through the event data for JC files
read_event_data <- function(directory, event_type) {
  file_path <- paste0(directory, event_type, ".MATS.JC.txt")
  read.table(file_path, header = TRUE, sep = "\t")
}


# Function to combine all event data into an EventType column
combine_event_data <- function(event_data_list) {
  bind_rows(lapply(names(event_data_list), function(event_type) {
    mutate(event_data_list[[event_type]], EventType = event_type)
  }))
}

# Function to filter event data according to significance of p-value and FDR values (new_FDRs determined from the p-values according to Benjamini Hochberg calculations)
filter_event_data <- function(event_data, p_threshold = 0.05) {
  event_data %>% filter(PValue < p_threshold, new_FDR < 0.05)
}

#Function to replace each FDR p-value with value of 0 to <2.2e-16 (see https://groups.google.com/g/rmats-user-group/c/TW534af62fg)
replace_zero_pvalues <- function(event_data, replacement_value = 2.2e-16) {
    event_data %>%
      mutate(
        PValue = if_else(PValue == 0, replacement_value, PValue),
        new_FDR = if_else(new_FDR == 0, replacement_value, new_FDR)
      )
  }

#Function to separate annotated and unannotated exons
#By looping through each directory, we are retrieving all unannotated events present in the fromGTF.SpliceSite[AS_event].txt files into their own data table
unannotated_list <- function(directories, event_types) {
  unannotated_data <- lapply(event_types, function(event_type) {
    for (directory in directories) {
      file_path <-
        paste0(directory,
               "fromGTF.novelSpliceSite.",
               event_type,
               ".txt")
    }
  })
  names(unannotated_data) <- event_types
  unannotated_data
}

add_ID_novel_column <-
  function(coverage_values_tables,
           unannotated_list) {
    if (!is.data.frame(coverage_values_tables)) {
      stop("coverage_values_tables must be a data frame")
    }
    # Loop through each event type and add ID_novel column
    for (event_type in names(unannotated_list)) {
      unannotated_data <- unannotated_list[[event_type]]
      coverage_values_tables <- coverage_values_tables %>%
        mutate(ID_novel = ifelse(ID %in% unannotated_data$ID, 'yes', 'no'))
    }
    coverage_values_tables
  }

# Function to extract coverage values
extract_coverage_values <- function(data) {
  data %>%
    rowwise() %>%
    mutate(coverage1 = IJC_SAMPLE_1 + SJC_SAMPLE_1) %>%
    ungroup() %>%
    mutate(
      first_IJC = sapply(IJC_SAMPLE_2, function(x)
        as.numeric(strsplit(x, ",")[[1]][1])),
      first_SJC = sapply(SJC_SAMPLE_2, function(x)
        as.numeric(strsplit(x, ",")[[1]][1])),
      coverage2a = first_IJC + first_SJC,
      second_IJC = sapply(IJC_SAMPLE_2, function(x)
        as.numeric(strsplit(x, ",")[[1]][2])),
      second_SJC = sapply(SJC_SAMPLE_2, function(x)
        as.numeric(strsplit(x, ",")[[1]][2])),
      coverage2b = second_IJC + second_SJC,
      third_IJC = sapply(IJC_SAMPLE_2, function(x)
        as.numeric(strsplit(x, ",")[[1]][3])),
      third_SJC = sapply(SJC_SAMPLE_2, function(x)
        as.numeric(strsplit(x, ",")[[1]][3])),
      coverage2c = third_IJC + third_SJC,
      fourth_IJC = sapply(IJC_SAMPLE_2, function(x)
        as.numeric(strsplit(x, ",")[[1]][4])),
      fourth_SJC = sapply(SJC_SAMPLE_2, function(x)
        as.numeric(strsplit(x, ",")[[1]][4])),
      coverage2d = fourth_IJC + fourth_SJC
    ) %>%
    select(
      -first_IJC,
      -first_SJC,
      -second_IJC,
      -second_SJC,
      -third_IJC,
      -third_SJC,
      -fourth_IJC,
      -fourth_SJC
    )
}

# Function to plot the distribution of the inclusion level difference data
plot_distribution_data <- function(filtered_cutoff) {
  coverage_long <- filtered_cutoff %>%
    select(IncLevelDifference) %>%
    pivot_longer(cols = IncLevelDifference,
                 names_to = "IncLevelDifference",
                 values_to = "Value") %>%
    mutate(Row = row_number())
  ggplot(coverage_long, aes(x = Value)) +
    geom_histogram(binwidth = 0.02) +  # Adjust binwidth as necessary
    labs(title = "Inclusion level Histogram", x = "Inclusion Level Differences", y = "Count") +
    theme_minimal() +
    geom_vline(
      xintercept = c(-0.01, 0.01),
      linetype = "dashed",
      color = "red"
    )
}

# Function to filter data based on coverage cutoff
apply_coverage_cutoff <- function(coverage_data) {
  coverage_data %>%
    filter(coverage1 >= 20 |
             coverage2a >= 20 |
             coverage2b >= 20 |
             coverage2c >= 20 |
             coverage2d >= 20)
}

# Function to filter inclusion level differences over 0.15 and under -0.15
filter_inclusion_cutoffs <- function(filtered_cutoff) {
  filtered_cutoff %>%
    filter(IncLevelDifference < -0.15 | IncLevelDifference > 0.15)
}

#new graphs which should include the data from the gathered tables but also including a small stripe with the control event
gene_name <-
  c("SNAPIN",
    "PTPN11",
    "PALB2",
    "EFTUD2",
    "IGF1R",
    "SDHA",
    "MSH2",
    "MLH1",
    "PHF8")
events_per_gene <-
  function(filtered_cutoff, gene_name) {
    #had to modify the background data to be the filtered for significance so that i would see all counts displayed it the graph
    all_inclusion <- filtered_cutoff %>%
      select(IncLevelDifference, geneSymbol) %>%
      pivot_longer(cols = IncLevelDifference,
                   names_to = "IncLevelDifference",
                   values_to = "Value")
    inclusion_long <- all_inclusion %>%
      filter(geneSymbol == gene_name)
    p <- ggplot() +
      geom_histogram(
        data = all_inclusion,
        aes(x = Value),
        binwidth = 0.02,
        alpha = 0.1 ,
        color = "#e9ecef" ,
        position = 'identity'
      ) +
      geom_vline(
        data = inclusion_long,
        aes(xintercept = Value),
        color = "#404080",
        alpha = 1.0
      ) +
      labs(
        title = paste(
          "Inclusion Level Histogram for filtered data including events for",
          gene_name
        ),
        x = "Inclusion Level Differences",
        y = "Count"
      ) +
      theme_minimal() +
      geom_vline(
        xintercept = c(-0.15 , 0.15) ,
        linetype = "dashed",
        color = "red"
      )
    return(p)
  }

# Initialize lists to store results
unfiltered_tables <- list()
coverage_values_tables <- list()
unannotated_data_tables <- list()
coverage_val_annotation_tables <- list()
coverage_cutoff_tables <- list()
filtered_cutoff_tables <- list()
changed_PValues_under_0 <- list()
filtered_inclusion_tables <- list()
coverage_plots <- list()
gene_plots <- list()

# Loop through each directory
for (dir in output_directories) {
  cat("\nProcessing directory:", dir, "\n")
  
  # Read event data
  event_data_list <-
    lapply(event_types, read_event_data, directory = dir)
  names(event_data_list) <- event_types
  
  # Combine event data without filtering
  unparsed_table <- combine_event_data(event_data_list)
  unparsed_table <-
    unparsed_table %>% mutate(new_FDR = p.adjust(PValue, method = "BH"))
  
  # Select relevant columns if they exist
  relevant_columns <-
    c(
      "ID",
      "geneSymbol",
      "EventType",
      "chr",
      "strand",
      "IJC_SAMPLE_1",
      "SJC_SAMPLE_1",
      "IJC_SAMPLE_2",
      "SJC_SAMPLE_2",
      "PValue",
      "new_FDR",
      "IncLevelDifference",
      "exonStart_0base",
      "exonEnd",
      "upstreamES" ,
      "upstreamEE" ,
      "riExonStart_0base",
      "riExonEnd",
      "downstreamES",
      "downstreamEE",
      "X1stExonStart_0base",
      "X1stExonEnd",
      "X2ndExonStart_0base",
      "X2ndExonEnd",
      "longExonStart_0base",
      "longExonEnd",
      "shortES",
      "shortEE",
      "flankingES",
      "flankingEE",
      "1stExonStart_0base",
      "1stExonEnd",
      "2ndExonStart_0base",
      "2ndExonEnd"
    )
  
  
  unfiltered_table <- unparsed_table %>% select(any_of(relevant_columns))
  
  # Filter event data
  filtered_table <- filter_event_data(unfiltered_table, 0.05)
  
  # Extract coverage values
  coverage_values <- extract_coverage_values(unparsed_table)
  
  #adding the column for annotated and unannotated events
  unannotated_data <-
    unannotated_list(output_directories, event_types)
  coverage_vals_annotations <-
    add_ID_novel_column(coverage_values, unannotated_data)
  
  # Apply coverage cutoff
  coverage_cutoff <-
    apply_coverage_cutoff(coverage_vals_annotations)
  
  #filter the coverage for significance
  filtered_cutoff <- filter_event_data(coverage_cutoff, 0.05)
  corrected_event_data <- replace_zero_pvalues(filtered_cutoff)
  filtered_inclusion <-
    filter_inclusion_cutoffs(corrected_event_data)
  
  # Ploting the coverage data
  coverage_plot <- plot_distribution_data(filtered_inclusion)
  
  #saving the histograms as well
  plot_per_gene <- events_per_gene(filtered_inclusion, gene_name)
  
  # Store results
  unfiltered_tables[[dir]] <- unfiltered_table
  coverage_values_tables[[dir]] <- coverage_values
  unannotated_data_tables[[dir]] <- unannotated_data
  coverage_val_annotation_tables[[dir]] <-
    coverage_vals_annotations
  coverage_cutoff_tables[[dir]] <- coverage_cutoff
  filtered_cutoff_tables[[dir]] <- filtered_cutoff
  changed_PValues_under_0[[dir]] <- corrected_event_data
  filtered_inclusion_tables[[dir]] <- filtered_inclusion
  coverage_plots[[dir]] <- coverage_plot
  gene_plots[[dir]] <- plot_per_gene
}

# To access the results for each directory, change the value in the brackets which also depends on the order of your initial output directories list (above)
unfiltered_tables[[output_directories[1]]] #snapin unfiltered table
coverage_cutoff_tables[[output_directories[2]]] #ptpn11 gene
filtered_inclusion_tables[[output_directories[2]]]
filtered_cutoff_tables[[output_directories[2]]]
coverage_plots[[output_directories[2]]]
coverage_val_annotation_tables[[output_directories[2]]]

View(filtered_cutoff_tables[[output_directories[1]]])
# events_per_gene(snapin_filtered, "SNAPIN") #expect 4 counts for snapin
# events_per_gene(mlh1_filtered, "MLH1") #expect 10 counts for MLH1
# events_per_gene(msh2_filtered, "MSH2") #expect 3 counts (values for inclevediff are the following; -0.138, 0.307 and 0.501)
# events_per_gene(sdha_filtered, "SDHA") # expect 3 counts (values  -0.017, 0.605 and 0.526)
# events_per_gene(eftud2_filtered, "EFTUD2") # 5 counts (values -0.311, -0.54, -0.569, -0.903 and 0.036)
# events_per_gene(phf8_filtered, "PHF8") # 8 counts
# events_per_gene(ifg1r_filtered, "IGF1R") # 1 count (value 0.325)
# events_per_gene(palb2_filtered, "PALB2") # 2 counts
# events_per_gene(ptpn11_filtered, "PTPN11") # 11 counts


#preparing the Venn diagram between JC locations and JCEC genomic locations
#already have the JC event data and simply need to retrieve the JCEC data for one gene (or multiple)
library(VennDiagram)
library(grid)


directory_1 <-
  '~/cafalempin/SNAPIN/ensembl/output/' #you can replace with another directory here to retrieve a table of the JCEC output for a particular gene

read_events_JCEC <- function(directory_1, event_type) {
  file_path <- paste0(directory_1, event_type, ".MATS.JCEC.txt")
  read.table(file_path, header = TRUE, sep = "\t")
}

event_data_list_JCEC <-
  lapply(event_types, read_events_JCEC, directory = directory_1)
names(event_data_list_JCEC) <- event_types

#the MXE event has 25 columns whereas all the other events has 23 columns, we need to account for this in the code
all_columns <-
  unique(unlist(lapply(event_data_list_JCEC, colnames)))
# Add missing columns to each data frame, filling them with NA
event_data_list_JCEC <- lapply(event_data_list_JCEC, function(df) {
  missing_columns <- setdiff(all_columns, colnames(df))
  if (length(missing_columns) > 0) {
    df[missing_columns] <- NA
  }
  return(df)
})
combine_event_data <- function(event_data_list) {
  do.call(rbind, event_data_list)
}

#one of the data frames that will be used in the Venn diagram
set_JCEC_unfiltered  <- combine_event_data(event_data_list_JCEC)

#retrieving the JCEC filtered file (have to run the functions above onto the current dataset)
# Process and filter data for the JCEC files instead of the JC files above
process_data <- function(data) {
  data %>%
    mutate(new_FDR = p.adjust(PValue, method = "BH")) %>%
    filter_event_data(0.05) %>%
    extract_coverage_values() %>%
    apply_coverage_cutoff() %>%
    filter_event_data(0.05) %>%
    replace_zero_pvalues() %>%
    filter_inclusion_cutoffs()
}

#another of the data frames that will be used in the Venn diagram
set_JCEC_filtered  <- process_data(set_JCEC_unfiltered)

#set for the unfiltered data JC
set_JC_unfiltered <- unfiltered_tables[[output_directories[4]]]
set_JC_filtered <-
  filtered_inclusion_tables[[output_directories[4]]]

View(set_JCEC_unfiltered)
View(set_JCEC_filtered)
View(set_JC_filtered)
View(set_JC_unfiltered)

# Extract genomic coordinates
extract_genomic_coordinates <- function(data) {
  data %>%
    select(
      "ID",
      "geneSymbol",
      "exonStart_0base",
      "exonEnd",
      "upstreamES",
      "upstreamEE",
      "downstreamES",
      "downstreamEE",
      "X1stExonStart_0base",
      "X1stExonEnd",
      "X2ndExonStart_0base",
      "X2ndExonEnd",
      "longExonStart_0base",
      "longExonEnd",
      "shortES",
      "shortEE",
      "flankingES",
      "flankingEE",
      "riExonStart_0base",
      "riExonEnd"
    )
}

# Generate unique coordinates
generate_coordinates <- function(genomic_coordinates) {
  unique(
    c(
      na.omit(genomic_coordinates$exonStart_0base),
      na.omit(genomic_coordinates$exonEnd),
      na.omit(genomic_coordinates$upstreamES),
      na.omit(genomic_coordinates$upstreamEE),
      na.omit(genomic_coordinates$downstreamES),
      na.omit(genomic_coordinates$downstreamEE),
      na.omit(genomic_coordinates$X1stExonStart_0base),
      na.omit(genomic_coordinates$X1stExonEnd),
      na.omit(genomic_coordinates$X2ndExonStart_0base),
      na.omit(genomic_coordinates$X2ndExonEnd),
      na.omit(genomic_coordinates$longExonStart_0base),
      na.omit(genomic_coordinates$longExonEnd),
      na.omit(genomic_coordinates$shortES),
      na.omit(genomic_coordinates$shortEE),
      na.omit(genomic_coordinates$flankingES),
      na.omit(genomic_coordinates$flankingEE),
      na.omit(genomic_coordinates$riExonStart_0base),
      na.omit(genomic_coordinates$riExonEnd)
    )
  )
}

set_JCEC_unfiltered <-
  extract_genomic_coordinates(set_JCEC_unfiltered)
coordinatesJCEC <- generate_coordinates(set_JCEC_unfiltered)

set_JC_unfiltered <- extract_genomic_coordinates(set_JC_unfiltered)
coordinatesJC <- generate_coordinates(set_JC_unfiltered)

set_JCEC_filtered <- extract_genomic_coordinates(set_JCEC_filtered)
coordinatesJCEC_filtered <- generate_coordinates(set_JCEC_filtered)

set_JC_filtered <- extract_genomic_coordinates(set_JC_filtered)
coordinatesJC_filtered <- generate_coordinates(set_JC_filtered)

# Create a list with your two sets of coordinates
venn_data <-
  list(DatasetJCEC = coordinatesJCEC, DatasetJC = coordinatesJC)
venn_data_filtered <-
  list(DatasetJCEC = coordinatesJCEC_filtered, DatasetJC = coordinatesJC_filtered)

# Function to create Venn diagrams
create_venn_plot <- function(venn_data) {
  venn_plot <- venn.diagram(
    x = venn_data,
    category.names = c("Dataset JCEC", "Dataset JC"),
    filename = NULL,
    output = TRUE,
    col = c("#440154ff", '#21908dff'),
    fill = c(alpha("#440154ff", 0.3), alpha('#21908dff', 0.3)),
    height = 3000,
    width = 3000,
    resolution = 300,
    compression = "lzw",
    cex = 1.1,
    fontface = "bold",
    fontfamily = "sans",
    cat.cex = 1,
    cat.fontface = "bold",
    cat.default.pos = "outer",
    cat.fontfamily = "sans"
  )
  grid.draw(venn_plot)
}
create_venn_plot(venn_data)
create_venn_plot(venn_data_filtered)
