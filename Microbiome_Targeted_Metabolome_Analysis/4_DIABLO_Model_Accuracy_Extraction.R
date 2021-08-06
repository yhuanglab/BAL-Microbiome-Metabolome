suppressMessages(library(dplyr))
suppressMessages(library(snow))
suppressMessages(library(tidyr))
suppressMessages(library(stringr))
suppressMessages(library(reshape2))
options(warn = -1)

target_variable <- str_remove(getwd(), "/Users/Siddharth/Desktop/Spiromics_BAL_Microbiome_Metabolome_V2/Spiromics_UNC_Stratified/Spiromics_No_Severe_No_Healthy/Output/DIABLO/") #remove path and just get target variable name

parameters <- read.csv("DIABLO_Model_Parameters.csv") #read file with model statistics

crude_accuracy <- parameters[which(parameters$"Model.Parameter" == "Testing Overall BAR"), "Mean.Value.Over.100.Resampling.Iterations"]

p_value <- parameters[which(parameters$"Model.Parameter" == "Testing Overall BAR Significance Test P-Value"), "Mean.Value.Over.100.Resampling.Iterations"]

#identify the category of the outcome modeled
outcome_names_categories_for_figures <- read.csv("~/Desktop/Spiromics_BAL_Microbiome_Metabolome_V2/Outcome_Names_Categories_for_Figures.csv")
category <- outcome_names_categories_for_figures[which(outcome_names_categories_for_figures$Model == target_variable), "Category"] %>% as.character()

#identify the name we're using for the outcome in the figure
model_name <- outcome_names_categories_for_figures[which(outcome_names_categories_for_figures$Model == target_variable), "Model_Name"] %>% as.character()

#create data frame with outcome name, figure name, category, crude error, and crude accuracy for each outcome modeled
labels <- c("Model", "Model_Name", "Category", "Accuracy", "P_Value")
values <- c(target_variable, model_name, category, crude_accuracy, p_value)

if (target_variable == "BRON_PEX_HCUDRUG_365_A.y_1_Oneormore_0_None") {
	df <- rbind(labels, values) %>% data.frame()
	colnames(df) <- NULL
	rownames(df) <- NULL

	#append to table of all outcomes
	write.table(df, file = "~/Desktop/Spiromics_BAL_Microbiome_Metabolome_V2/Spiromics_UNC_Stratified/Spiromics_No_Severe_No_Healthy/Output/DIABLO/DIABLO_Model_Accuracies_Groups.csv", 
					sep = ",", 
					append = TRUE, 
					row.names = FALSE, 
					col.names = FALSE)
} else {
	df <- rbind(labels, values) %>% data.frame()
	df <- df[-1, , drop = FALSE]
	write.table(df, file = "~/Desktop/Spiromics_BAL_Microbiome_Metabolome_V2/Spiromics_UNC_Stratified/Spiromics_No_Severe_No_Healthy/Output/DIABLO/DIABLO_Model_Accuracies_Groups.csv", 
					sep = ",", 
					append = TRUE, 
					row.names = FALSE, 
					col.names = FALSE)
}