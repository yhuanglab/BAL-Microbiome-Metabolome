#!/bin/bash

#clear existing output file if present
rm ~/Desktop/Spiromics_BAL_Microbiome_Metabolome_V2/Spiromics_Colorado_Stratified/Spiromics_No_Severe_No_Healthy/Output/DIABLO/DIABLO_Model_Accuracies_Groups.csv

#create new output file
touch ~/Desktop/Spiromics_BAL_Microbiome_Metabolome_V2/Spiromics_Colorado_Stratified/Spiromics_No_Severe_No_Healthy/Output/DIABLO/DIABLO_Model_Accuracies_Groups.csv

#iterate through each outcome modeled and extract mean testing out-of-sample accuracy
for i in ~/Desktop/Spiromics_BAL_Microbiome_Metabolome_V2/Spiromics_Colorado_Stratified/Spiromics_No_Severe_No_Healthy/Output/DIABLO/*
do
	echo "***********"
	echo "Starting $i"
	cd $i #Enter folder
	Rscript ~/Desktop/Spiromics_BAL_Microbiome_Metabolome_V2/Spiromics_Colorado_Stratified/Spiromics_No_Severe_No_Healthy/Code/4_DIABLO_Model_Accuracy_Extraction.R
	echo "***********"
	cd .. #Exit folder
done

#drop the repeated column headers