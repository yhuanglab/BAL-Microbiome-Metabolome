#!/bin/sh

#wrapper script for downstream statistical analyses after modeling
for i in ~/Desktop/Spiromics_BAL_Microbiome_Metabolome_V2/Spiromics_Colorado_Stratified/Spiromics_No_Severe_No_Healthy/Output/DIABLO/*
do
		echo "Starting $i"

		cd $i

		Rscript ~/Desktop/Spiromics_BAL_Microbiome_Metabolome_V2/Spiromics_Colorado_Stratified/Spiromics_No_Severe_No_Healthy/Code/5_Microbiome_Metabolome_DIABLO_Metabolites_Importance_Stripchart.R
		Rscript ~/Desktop/Spiromics_BAL_Microbiome_Metabolome_V2/Spiromics_Colorado_Stratified/Spiromics_No_Severe_No_Healthy/Code/5_Microbiome_Metabolome_DIABLO_OTU_Importance_Stripchart.R

		#********#
		#Boxplots#
		#********#

		#create directories to store boxplots
		mkdir $i/Metabolite_Boxplots
		mkdir $i/OTU_Boxplots

		#generate boxplots
		Rscript ~/Desktop/Spiromics_BAL_Microbiome_Metabolome_V2/Spiromics_Colorado_Stratified/Spiromics_No_Severe_No_Healthy/Code/5_Microbiome_Metabolome_DIABLO_Metabolites_Relative_Abundance_Boxplots.R
		Rscript ~/Desktop/Spiromics_BAL_Microbiome_Metabolome_V2/Spiromics_Colorado_Stratified/Spiromics_No_Severe_No_Healthy/Code/5_Microbiome_Metabolome_DIABLO_OTU_Relative_Abundance_Boxplots.R

		echo "$i Completed"

		cd .. #exit folder
done
