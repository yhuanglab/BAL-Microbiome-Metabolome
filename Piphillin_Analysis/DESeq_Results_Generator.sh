#!/usr/bin/bash

Rscript ~/Desktop/Spiromics_Piphillin/Code/run_deseq2.R

python ~/Desktop/Spiromics_Piphillin/Code/parse_keg.py ~/Desktop/Spiromics_Piphillin/DESeq2_Results/PFV45_DERV_ctb_DESeq2_Results.tsv > ~/Desktop/Spiromics_Piphillin/DESeq2_Results/PFV45_DERV_ctb_DESeq2_Results_Mapped_Pathways.tsv