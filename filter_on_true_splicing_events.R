#This R script is aiming to parse and render the rMATS output for "known" or standard splicing events to use this list as a form of filtering to keep only truly diferentially spliced events from rMATS
# this goal was mentioned in the future perspectives of my discussion

#this code is not yet fully complete also feel free to also further develop onto it

#Overview: make a list of all the true splicing events found it AT LEAST 3 of the replicate control samples FOR EACH EVENT TYPE (at least 3 because some of our control event genes only contained 3 or 4 replicates)
#exception for EFTUD2 gene because it only has 2 replicate samples so there will be a separate list for this gene

librayr(dplyr)

unfiltered_tables <-
  list(
    snapin_annotations,
    ptpn11,
    palb2_dir_coverage_values,
    ifg1r_dir_coverage_values,
    sdha_dir_coverage_values,
    phf8_dir_coverage_values,
    msh2_dir_coverage_values,
    mlh1_dir_coverage_values
  )


event_types <- c("SE", "MXE", "A3SS", "A5SS", "RI")
genes <-
  c("SNAPIN",
    "PTPN11",
    "PALB2",
    "IGF1R",
    "SDHA",
    "MSH2",
    "MLH1",
    "PHF8")

# Function to extract potential true events
potential_true_events <-
  function(unfiltered_tables,
           genes,
           event_types,
           eftud2_dir_coverage_values = NULL) {
    potential_true_events <-
      function(unfiltered_tables,
               genes,
               event_types,
               eftud2_data = NULL) {
        # Define the column mappings for each event type
        column_mappings <- list(
          SE = c(
            "EventType",
            "exonStart_0base",
            "exonEnd",
            "upstreamES",
            "upstreamEE",
            "downstreamES",
            "downstreamEE"
          ),
          RI = c(
            "EventType",
            "riExonStart_0base",
            "riExonEnd",
            "upstreamES",
            "upstreamEE",
            "downstreamES",
            "downstreamEE"
          ),
          MXE = c(
            "EventType",
            "1stExonStart_0base",
            "1stExonEnd",
            "2ndExonStart_0base",
            "2ndExonEnd",
            "upstreamES",
            "upstreamEE",
            "downstreamES",
            "downstreamEE"
          ),
          A3SS = c(
            "EventType",
            "longExonStart_0base",
            "longExonEnd",
            "shortES",
            "shortEE",
            "flankingES",
            "flankingEE"
          ),
          A5SS = c(
            "EventType",
            "longExonStart_0base",
            "longExonEnd",
            "shortES",
            "shortEE",
            "flankingES",
            "flankingEE"
          )
        )
        
        true_splicing_events <- list()
        eftud2_true_splicing_events <- list()
        
        # Process each gene
        for (i in seq_along(unfiltered_tables)) {
          gene <- genes[i]
          data_table <- unfiltered_tables[[i]]
          gene_results <- list()
          
          for (event_type in event_types) {
            columns <- column_mappings[[event_type]]
            event_data <-
              data_table[data_table$EventType == event_type, columns, drop = FALSE]
            
            if (nrow(event_data) == 0) {
              next  # Skip if no data for this event type
            }
            
            # Count the occurrence of each splicing event based on the relevant columns
            event_string <- apply(event_data, 1, paste, collapse = "-")
            splicing_event_counts <- table(event_string)
            
            # Filter events that appear in at least 3 replicates
            filtered_events <-
              names(splicing_event_counts[splicing_event_counts >= 3])
            
            # Convert filtered events back to data frame format
            filtered_event_data <-
              event_data[event_string %in% filtered_events, , drop = FALSE]
            
            gene_results[[event_type]] <- filtered_event_data
          }
          true_splicing_events[[gene]] <- gene_results
        }
        
        
        
        # Process the EFTUD2 gene separately if provided
        if (!is.null(eftud2_data)) {
          gene_results <- list()
          
          for (event_type in event_types) {
            columns <- column_mappings[[event_type]]
            
            event_data <-
              eftud2_data[eftud2_data$EventType == event_type, columns, drop = FALSE]
            
            if (nrow(event_data) == 0) {
              next  # Skip if no data for this event type
            }
            
            event_string <- apply(event_data, 1, paste, collapse = "-")
            splicing_event_counts <- table(event_string)
            
            # Filter events that appear in at least 2 replicates for EFTUD2
            filtered_events <-
              names(splicing_event_counts[splicing_event_counts >= 2])
            filtered_event_data <-
              event_data[event_string %in% filtered_events, , drop = FALSE]
            
            gene_results[[event_type]] <- filtered_event_data
          }
          
          eftud2_true_splicing_events <- gene_results
        } else {
          eftud2_true_splicing_events <- NULL
        }
        
        return(
          list(
            true_splicing_events = true_splicing_events,
            eftud2_true_splicing_events = eftud2_true_splicing_events
          )
        )
      }
    
    # Call the function with your data
    result <- potential_true_events(
      unfiltered_tables = unfiltered_tables,
      genes = genes,
      event_types = event_types,
      eftud2_data = eftud2_dir_coverage_values  # Provide this if you have EFTUD2 data
    )
    
    # Access the results
    true_splicing_events <- result$true_splicing_events
    eftud2_true_splicing_events <- result$eftud2_true_splicing_events
    
    # View true splicing events for SNAPIN gene and SE event type
    snapin_se_events <- true_splicing_events[["SNAPIN"]][["SE"]]
    View(snapin_se_events)