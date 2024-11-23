# rMATS analysis

## About

Method elaborated to assess the performance fo rMATS in accurately detecting control splicing events from RNA sequencing data from patients of the Clinical Genetics group of the LUMC.

## Description

With the following R files, rMATS output files were assembled into one table and further wrangling of the data was performed. The files that are contained in this project are the following;

-   tables.for.all.genes : this R script is possibly the main and more important one of this directory since it extracts, wrangles and parses the rMATS output data files and merges them into one big table onto which a series of cutoffs are applied. The cutoffs were the following: (1) on the count coverage from the RNA-seq data (where at least one of the read count columns needs to be equal to 20 or above, otherwise the event is filtered out). (2) on the new FDR values determined from the p-values using the Benjamini Hochberg method and a cutoff on significance is applied (PValue \< 0.05 and FDR \< 0.05).

-   filter_on_true_splicing_events : from this R script the aim is to parse and render the rMATS output for "known" or standard splicing events to use this list as a form of filtering to keep only truly diferentially spliced events from rMATS. The overall aim in this file is to make a list of all the true splicing events found it at least 3 of the replicate control samples FOR EACH EVENT TYPE (at least 3 because some of our control event genes only contained 3 or 4 replicates and the EFTUD2 gene only had 2 replicates).

    *This goal was mentioned in the future perspectives of my discussion*, and this code is not yet fully elaborated so feel free to also further develop it.

-   selecting_tables_from_gene_panel : in this R script, we are retrieving the number of genes that can be found in each of the datasets for the potential investigation of off-target events that might be detected by rMATS

-   retrieving_summary_counts : this script was used to retrieve the summary counts data for each of the genes with known splicing events run with rMATS. The original output was rendered using the "incorrectly" merged GTF MANE clinical annotations, and the second run is with the correct GTF file (which is why the name of the destination files is R2_vs_all/output_151var/ and R2_vs_all/new_run/new_output_151/). The knowledge of the summary counts was a way for us to understand whether the GTF file annotations could influence the output rendered especially when it comes to the significant events reported (and if the control event we were looking for could be found in this data).

-   selecting_cutoff_individual_counts : this R script contains code that aims to determine a common cutoff on the output columns rendered once the \--individual-counts parameter is enabled in rMATS (from the latest update released in March 2024). Each event type has its own specific columns... the latest information about this parameter can be found here: <https://github.com/Xinglab/rmats-turbo>. This code is not yet fully developed so please feel free to use any interesting aspects of it to further analyze your rMATS data.

------------------------------------------------------------------------

# Commiting new changes to the repository

    cd /exports/sascstudent/cafalempin/R.scripts
    git status
    git add . (git commit -a)
    git remote add origin https://git.lumc.nl/sasc/rmats_analysis.git
    git push origin main

## Support

In the case that you would need any help with understanding the parameters, or how to set each of your parameters, there is a user group for rMATS: <https://groups.google.com/g/rmats-user-group>

## Authors and acknowledgment

I would like to thank give a special thanks to Leon and Davy and my colleagues of the SASC core for their invaluable guidance all throughout this internship project.

## License

This project was a part of a masters internship at the SASC core group of the LUMC. All scripts can be used openly and modified for further projects in the group.
