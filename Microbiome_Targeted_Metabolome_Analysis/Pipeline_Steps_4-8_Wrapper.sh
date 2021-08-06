#!/bin/bash

#Shell script to run all analyses downstream of DIABLO model generation (barcharts, heatmaps, networks)

#***************************************************************************************************************#
#Step 4 - Extracting DIABLO model out-of-sample accuracies for each outcome and creating barchart for comparison
echo "************************************************"
echo "4_DIABLO_Model_Accuracy_Extraction.R --- RUNNING"

#clear existing output file if present
rm ~/Desktop/Spiromics_BAL_Microbiome_Metabolome_V2/Spiromics_UNC_Stratified/Spiromics_No_Severe_No_Healthy/Output/DIABLO/DIABLO_Model_Accuracies_Groups.csv

#create new output file
touch ~/Desktop/Spiromics_BAL_Microbiome_Metabolome_V2/Spiromics_UNC_Stratified/Spiromics_No_Severe_No_Healthy/Output/DIABLO/DIABLO_Model_Accuracies_Groups.csv

#iterate through each outcome modeled and extract mean testing out-of-sample accuracy
for i in ~/Desktop/Spiromics_BAL_Microbiome_Metabolome_V2/Spiromics_UNC_Stratified/Spiromics_No_Severe_No_Healthy/Output/DIABLO/*
do
	#echo "Starting $i"
	cd $i #Enter folder
	Rscript ~/Desktop/Spiromics_BAL_Microbiome_Metabolome_V2/Spiromics_UNC_Stratified/Spiromics_No_Severe_No_Healthy/Code/4_DIABLO_Model_Accuracy_Extraction.R
	cd .. #Exit folder
done

echo "4_DIABLO_Model_Accuracy_Extraction.R --- DONE"
echo "************************************************"
echo "4_DIABLO_Model_Accuracy_Graph.R --- RUNNING"

Rscript ~/Desktop/Spiromics_BAL_Microbiome_Metabolome_V2/Spiromics_UNC_Stratified/Spiromics_No_Severe_No_Healthy/Code/4_DIABLO_Model_Accuracy_Graph.R

echo "4_DIABLO_Model_Accuracy_Graph.R --- DONE"
echo "************************************************"
#***************************************************************************************************************#


#***************************************************************************************************************#
#Step 5 - Feature importance stripcharts and boxplots for each DIABLO outcome model
echo "************************************************"
echo "5_Microbiome_Metabolome_DIABLO_OTU_Relative_Abundance_Boxplots.R --- RUNNING"
for i in ~/Desktop/Spiromics_BAL_Microbiome_Metabolome_V2/Spiromics_UNC_Stratified/Spiromics_No_Severe_No_Healthy/Output/DIABLO/*
do
		#echo "Starting $i"

		cd $i

		Rscript ~/Desktop/Spiromics_BAL_Microbiome_Metabolome_V2/Spiromics_UNC_Stratified/Spiromics_No_Severe_No_Healthy/Code/5_Microbiome_Metabolome_DIABLO_Metabolites_Importance_Stripchart.R
		Rscript ~/Desktop/Spiromics_BAL_Microbiome_Metabolome_V2/Spiromics_UNC_Stratified/Spiromics_No_Severe_No_Healthy/Code/5_Microbiome_Metabolome_DIABLO_OTU_Importance_Stripchart.R

		#********#
		#Boxplots#
		#********#

		#create directories to store boxplots
		mkdir $i/Metabolite_Boxplots
		mkdir $i/OTU_Boxplots

		#generate boxplots
		Rscript ~/Desktop/Spiromics_BAL_Microbiome_Metabolome_V2/Spiromics_UNC_Stratified/Spiromics_No_Severe_No_Healthy/Code/5_Microbiome_Metabolome_DIABLO_Metabolites_Relative_Abundance_Boxplots.R
		Rscript ~/Desktop/Spiromics_BAL_Microbiome_Metabolome_V2/Spiromics_UNC_Stratified/Spiromics_No_Severe_No_Healthy/Code/5_Microbiome_Metabolome_DIABLO_OTU_Relative_Abundance_Boxplots.R

		#echo "$i Completed"

		cd .. #exit folder
done

echo "5_Microbiome_Metabolome_DIABLO_OTU_Relative_Abundance_Boxplots.R --- DONE"
echo "************************************************"


#***************************************************************************************************************#
#Step 6 - Extract Features in 500+ models
echo "************************************************"
echo "6_DIABLO_Top_Features_Extraction.R --- RUNNING"

#clear existing output file if present
rm ~/Desktop/Spiromics_BAL_Microbiome_Metabolome_V2/Spiromics_UNC_Stratified/Spiromics_No_Severe_No_Healthy/Output/Heatmap/DIABLO_Num_Models_Heatmap/Top_Features_by_Num_Sig_Models_No_Unknowns.csv

#create new output file
touch ~/Desktop/Spiromics_BAL_Microbiome_Metabolome_V2/Spiromics_UNC_Stratified/Spiromics_No_Severe_No_Healthy/Output/Heatmap/DIABLO_Num_Models_Heatmap/Top_Features_by_Num_Sig_Models_No_Unknowns.csv

for i in ~/Desktop/Spiromics_BAL_Microbiome_Metabolome_V2/Spiromics_UNC_Stratified/Spiromics_No_Severe_No_Healthy/Output/DIABLO/*
do
	#echo "***********"
	#echo "Starting $i"
	cd $i #Enter folder
	Rscript ~/Desktop/Spiromics_BAL_Microbiome_Metabolome_V2/Spiromics_UNC_Stratified/Spiromics_No_Severe_No_Healthy/Code/6_DIABLO_Top_Features_Extraction.R
	#echo "***********"
	cd .. #Exit folder
done

echo "6_DIABLO_Top_Features_Extraction.R --- DONE"
echo "************************************************"
#***************************************************************************************************************#


#***************************************************************************************************************#
#Step 7 - Calculate direction of effect for each feature with its outcome pair and create heatmaps showing # of models
echo "************************************************"
echo "7_DIABLO_Feature_Direction_of_Effects.R --- RUNNING"

Rscript ~/Desktop/Spiromics_BAL_Microbiome_Metabolome_V2/Spiromics_UNC_Stratified/Spiromics_No_Severe_No_Healthy/Code/7_DIABLO_Feature_Direction_of_Effects.R

echo "7_DIABLO_Feature_Direction_of_Effects.R --- DONE"
echo "************************************************"
echo "7_DIABLO_Feature_Num_Models_Heatmap.R--- RUNNING"

Rscript ~/Desktop/Spiromics_BAL_Microbiome_Metabolome_V2/Spiromics_UNC_Stratified/Spiromics_No_Severe_No_Healthy/Code/7_DIABLO_Feature_Num_Models_Heatmap.R

echo "7_DIABLO_Feature_Num_Models_Heatmap.R--- DONE"
echo "************************************************"
#***************************************************************************************************************#


#***************************************************************************************************************#
#Step 8 - Pearson correlation heatmaps and network models
echo "************************************************"
echo "8_DIABLO_Feature_Pearson_Network.R --- RUNNING"

Rscript ~/Desktop/Spiromics_BAL_Microbiome_Metabolome_V2/Spiromics_UNC_Stratified/Spiromics_No_Severe_No_Healthy/Code/8_DIABLO_Feature_Pearson_Network.R

echo "8_DIABLO_Feature_Pearson_Network.R --- DONE"
echo "************************************************"
echo "8_DIABLO_Feature_Pearson_Correlation_Heatmap.R --- RUNNING"

Rscript ~/Desktop/Spiromics_BAL_Microbiome_Metabolome_V2/Spiromics_UNC_Stratified/Spiromics_No_Severe_No_Healthy/Code/8_DIABLO_Feature_Pearson_Correlation_Heatmap.R

echo "8_DIABLO_Feature_Pearson_Correlation_Heatmap.R --- DONE"
echo "************************************************"
#***************************************************************************************************************#

echo "ALL ANALYSES PERFORMED AND FIGURES GENERATED"









