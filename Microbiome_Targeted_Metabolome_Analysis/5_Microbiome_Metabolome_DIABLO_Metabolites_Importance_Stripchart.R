suppressMessages(library(dplyr))
suppressMessages(library(ggplot2))
suppressMessages(library(stringr))
options(warn = -1)

#*****************************************#
#Automatically iterates over all variables#
#*****************************************#

target_variable <- str_remove(getwd(), "/Users/Siddharth/Desktop/Spiromics_BAL_Microbiome_Metabolome_V2/Spiromics_UNC_Stratified/Spiromics_No_Severe_No_Healthy/Output/DIABLO/") #remove path and just get target variable name

df <- read.csv(paste("~/Desktop/Spiromics_BAL_Microbiome_Metabolome_V2/Spiromics_UNC_Stratified/Spiromics_No_Severe_No_Healthy/Output/DIABLO/", target_variable, "/DIABLO_Model_Features_ORs_P_Values.csv", sep = ""), stringsAsFactors = FALSE) #table with Metabolites, loadings, and number of models in which they were selected
df <- df %>% arrange(desc(Number_of_Final_Models_Where_Feature_Selected))
formatted_name <- target_variable #title of chart

demographic_variables <- c("(Intercept)", "GENDER2", "RACE_Nonwhite1", "AGE_DERV_ctb", "CURRENT_SMOKER_ctb1", "INHALEDSTEROIDS_ctb1", "ANTI_USE_3MTH_BEF_BRON1") 

#subset to just the metabolite features
metabolites_list <- df$Feature[which(str_sub(df$Feature, start = 1, end = 3) != "Otu"), drop = FALSE]
if (length(which(metabolites_list %in% demographic_variables)) > 0) {
	metabolites_list <- metabolites_list[-which(metabolites_list %in% demographic_variables)]
}
df <- df[which(df$Feature %in% metabolites_list), , drop = FALSE]

#*******************#
#DO NOT MODIFY BELOW#
#*******************#

otu <- df$"Feature" #extract labels for the features
importance <- df$"Number_of_Final_Models_Where_Feature_Selected" #extract importance metric (number of models out of 100 in which the feature was significant) 

#Generate Stripchart
plot <- ggplot(df, aes(x=reorder(otu, -importance), y=importance)) + 
geom_point() + 
	geom_smooth(method = "lm") + 
	#manually draw rectangle around selected otus
	theme(text = element_text(size = 10)) + 
	theme(axis.text.x = element_text(angle = 60, hjust = 1, size = 6)) + 
	labs(title = paste(formatted_name, "Metabolite Importance Stripchart", #add label of what chart shows to the name of the variable
		         sep = " "), 
		 		 x = "Metabolite", 
		 		 y="Number of Final Models in Which Feature Was Selected (Out of 100)") + 
	theme(plot.title = element_text(face="bold", hjust = 0.5, size = 10))
pdf(paste("~/Desktop/Spiromics_BAL_Microbiome_Metabolome_V2/Spiromics_UNC_Stratified/Spiromics_No_Severe_No_Healthy/Output/DIABLO/", target_variable, "/DIABLO_Model_Metabolite_Feature_Importance_Stripchart.pdf", sep = ""))
plot(plot)
dev.off()

