library(dplyr)
library(stringr)

#read in OTU table
otu_counts_table <- read.csv("~/Desktop/Spiromics_PiPhillin/Data/Master_OTU_Table_11_13_20.csv")

#transpose table
otu_counts_table <- t(otu_counts_table)

write.csv(otu_counts_table, file = "~/Desktop/Spiromics_PiPhillin/Data/BAL_Decontam_OTU_Table_Counts_for_Piphillin_DESeq2.csv")
