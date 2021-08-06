#!/bin/bash

#clear existing output file if present
rm ~/Desktop/Spiromics_BAL_Microbiome_Metabolome_V2/Spiromics_UNC_Stratified/Spiromics_No_Severe_No_Healthy/Output/Heatmap/DIABLO_Num_Models_Heatmap/Top_Features_by_Num_Sig_Models_No_Unknowns.csv
rm ~/Desktop/Spiromics_BAL_Microbiome_Metabolome_V2/Spiromics_UNC_Stratified/Spiromics_No_Severe_No_Healthy/Output/Heatmap/DIABLO_Num_Models_Heatmap/Demographic_Features_by_Num_Sig_Models.csv


#create new output file
touch ~/Desktop/Spiromics_BAL_Microbiome_Metabolome_V2/Spiromics_UNC_Stratified/Spiromics_No_Severe_No_Healthy/Output/Heatmap/DIABLO_Num_Models_Heatmap/Top_Features_by_Num_Sig_Models_No_Unknowns.csv

for i in ~/Desktop/Spiromics_BAL_Microbiome_Metabolome_V2/Spiromics_UNC_Stratified/Spiromics_No_Severe_No_Healthy/Output/DIABLO/*
do
	echo "***********"
	echo "Starting $i"
	cd $i #Enter folder
	Rscript ~/Desktop/Spiromics_BAL_Microbiome_Metabolome_V2/Spiromics_UNC_Stratified/Spiromics_No_Severe_No_Healthy/Code/6_DIABLO_Top_Features_Extraction.R
	echo "***********"
	cd .. #Exit folder
done