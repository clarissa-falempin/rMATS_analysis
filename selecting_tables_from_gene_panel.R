#need to run the tables.for.all.genes.R file before running this one first 
library(purrr)

#in this R script, we are retrieving the number of genes that can be found in each of the datasets for the potential investigation of off-target events tha tmight be detected by rMATS 

#function to get the unique gene names and their counts from each data table 
get_gene_counts <- function(data_tables, gene_symbol_column) {
  
  unique_genes_per_dataset <- list()
  all_unique_genes <- list()
  
  for (i in seq_along(data_tables)) {
    unique_genes <- data_tables[[i]] %>%
      pull(gene_symbol_column) %>%
      unique()
    unique_genes_per_dataset[[paste0("dataset_", i)]] <- unique_genes
    
    all_unique_genes <- c(all_unique_genes, unique_genes)
  }
  total_unique_genes <- length(unique(all_unique_genes))
  cat("Unique genes per dataset:\n")
  print(unique_genes_per_dataset)
  
  cat("\nTotal number of unique genes across all datasets:", total_unique_genes, "\n")
  return(list(unique_genes_per_dataset = unique_genes_per_dataset,
              total_unique_genes_count = total_unique_genes))
}


#list of all data tables to parse from 
data_tables <- list(snapin_inclusion, ptpn11_inclusion, palb2_inclusion,
                    igf1r_inclusion, sdha_inclusion, eftud2_inclusion, 
                    phf8_inclusion, msh2_inclusion, mlh1_inclusion)


gene_counts <- get_gene_counts(data_tables, "geneSymbol")
print(gene_counts[[1]])

options(max.print=100000)
sna <- gene_counts$unique_genes_per_dataset$dataset_1
ptp <- gene_counts$unique_genes_per_dataset$dataset_2
palb <- gene_counts$unique_genes_per_dataset$dataset_3
igf <- gene_counts$unique_genes_per_dataset$dataset_4
sdh <- gene_counts$unique_genes_per_dataset$dataset_5
eft <- gene_counts$unique_genes_per_dataset$dataset_6
phf <- gene_counts$unique_genes_per_dataset$dataset_7
msh <- gene_counts$unique_genes_per_dataset$dataset_8
mlh <- gene_counts$unique_genes_per_dataset$dataset_9
kdm <- gene_counts$unique_genes_per_dataset$dataset_10

#should be 68 total
genes_from_panel <- c("APC","ATM","AXIN2","BAP1","BARD1","BLM","BMPR1A","BRCA1","BRCA2","BRIP1","CDC73","CDH1","CDK4","CDKN2A","CHEK2",
                      "DICER1","FH","FLCN","HOXB13","MAX","MEN1","MITF","MLH1","MLH3","MSH2","MSH3","MSH6","MUTYH","NF1","NF2","NSD1","NTHL1",
                      "PALB2","PMS2","PMS2CL","POLD1","POLE","POT1", "PRKAR1A", "PTCH1", "PTCH2", "PTEN", "RAD51B", "RAD51C", "RAD51D", "RB1",
                      "RET", "SDHA", "SDHAF2", "SDHB", "SDHC", "SDHD", "SMAD4", "SMARCA4", "SMARCB1", "SMARCE1", "SPRED1", "STK11", "SUFU", "TERT", 
                      "TGFBR1", "TGFBR2", "TMEM127", "TP53", "TSC1", "TSC2", "VHL","WT1")


output_panel <- list()
output_profiler <- list()

parse_for_panel <- function(filtered_inclusion){
  data <- filtered_inclusion %>%
    filter(geneSymbol %in% genes_from_panel) %>%
    select("ID","GeneID", "geneSymbol", "chr","strand", "EventType" ,"PValue", "new_FDR" , "IncLevelDifference" , "exonStart_0base", "exonEnd" ,"upstreamES",  "upstreamEE" ,   "downstreamES",
           "downstreamEE",  "X1stExonStart_0base", "X1stExonEnd" , "X2ndExonStart_0base",
           "X2ndExonEnd", "longExonStart_0base", "longExonEnd", "shortES", "shortEE", "flankingES", "flankingEE","riExonStart_0base", 
           "riExonEnd" ,  "coverage1" , "coverage2a",  "coverage2b",  "coverage2c", "coverage2d")
  return(data)
}

most_significant_25 <- function(filtered_inclusion){
  table <- filtered_inclusion %>%
    arrange(new_FDR) %>% 
    slice_head(n = 25) %>%
    select("ID","GeneID", "geneSymbol", "chr","strand", "EventType" ,"PValue", "new_FDR" , "IncLevelDifference" , "exonStart_0base", "exonEnd" ,"upstreamES",  "upstreamEE" ,   "downstreamES",
           "downstreamEE",  "X1stExonStart_0base", "X1stExonEnd" , "X2ndExonStart_0base", "X2ndExonEnd", "longExonStart_0base", "longExonEnd", "shortES", "shortEE", "flankingES", "flankingEE","riExonStart_0base", 
           "riExonEnd" ,  "coverage1" , "coverage2a",  "coverage2b",  "coverage2c", "coverage2d")
}


for (i in seq_along(data_tables)) {
  cat("\nProcessing table:", i, "\n")
  data <- data_tables[[i]]
  parsed_tables <- parse_for_panel(data)
  gprofiler_lists <- most_significant_25(table)
  # Store results
  output_panel[[paste0("table", i)]] <- parsed_tables
  output_profiler[[paste0("table", i)]]<- gprofiler_lists
}

View(output_panel[[(1)]])

dim(output_panel[[(1)]]) #snapin
dim(output_panel[[(2)]]) #ptpn11
dim(output_panel[[(3)]]) #palb2
dim(output_panel[[(4)]]) #igf1r
dim(output_panel[[(5)]]) #sdha
dim(output_panel[[(6)]]) #eftud2
dim(output_panel[[(7)]]) #phf8
dim(output_panel[[(8)]]) #msh2
dim(output_panel[[(9)]]) #mlh1

n <- snapin_filtered %>%
  arrange(new_FDR) %>% 
  slice_head(n = 25) %>%
  select("ID","GeneID", "geneSymbol", "chr","strand", "EventType" ,"PValue", "new_FDR" , "IncLevelDifference" , "exonStart_0base", "exonEnd" ,"upstreamES",  "upstreamEE" ,   "downstreamES",
         "downstreamEE",  "X1stExonStart_0base", "X1stExonEnd" , "X2ndExonStart_0base", "X2ndExonEnd", "longExonStart_0base", "longExonEnd", "shortES", "shortEE", "flankingES", "flankingEE","riExonStart_0base", 
         "riExonEnd" ,  "coverage1" , "coverage2a",  "coverage2b",  "coverage2c", "coverage2d")
  
 View(n) 
  
#creating a loop to get the 20 first off targets to look into further 
get_top_significant_results <- function(data, top_n = 25) {
  data %>%
    #group_by(EventType) %>%
    filter(new_FDR < 0.05) %>%
    arrange(new_FDR) %>% 
    slice_head(n = top_n) %>%
    ungroup() %>%
}

