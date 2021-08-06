library(ggplot2)

#******************************************#
#CRUDE MODEL ACCURACY BY CLINICAL PARAMETER#
#******************************************#

#read file with model names and accuracies
df <- read.csv('~/Desktop/Spiromics_BAL_Microbiome_Metabolome_V2/Spiromics_Colorado_Stratified/Spiromics_No_Severe_No_Healthy/Output/DIABLO/DIABLO_Model_Accuracies_Groups.csv')
#manually order columns of barchart
df$Model_Name <- factor(df$Model_Name, levels = c("FEV1/FVC", 
												  "FEV1 % of predicted",
												  "FEV1 (L)",
												  "FEF25-75 (L/s)",
												  "PEFR (L/s)",
												  "FEV1 BDR (%)",
												  "FVC BDR (%)",
												  "FVC BDR (mL)",
												  "# Exacerbations in Past 12 mos. at Baseline",
												  "# Exacerbations (HCU + AB/S) in Year 1",
												  "# Exacerbations (HCU + AB/S) Post-Bronchoscopy",
												  "Mild/Moderate COPD",
												  "CAT Score",
												  "Chronic Bronchitis (SGRQ)"))
df$Sig <- ifelse(df$P_Value < 0.05, "*", "ns")

#plot
pdf("~/Desktop/Spiromics_BAL_Microbiome_Metabolome_V2/Spiromics_Colorado_Stratified/Spiromics_No_Severe_No_Healthy/Output/DIABLO/DIABLO_Classification_Accuracy_by_Clinical_Parameter_Groups.pdf", 
	width = 10, 
	height = 10)
ggplot(data = df) + 
geom_col(aes(x = Model_Name, y = Accuracy, fill = Category)) +
scale_x_discrete(labels = function(x) stringr::str_wrap(x, width = 25)) + 
geom_text(aes(x = Model_Name, 
  y = Accuracy + 0.025, 
  	label = paste(round(Accuracy, 2), " (", Sig, ")", sep = "")), size = 5) +
xlab("Clinical Measure") +
ylab("Mean Out-of-Sample Classification Accuracy Over 100 Models") +
theme(axis.text.x = element_text(angle=45, hjust=1, size = 12, color = 'black'), 
	  axis.text.y = element_text(size = 12),
	  axis.title.x = element_text(size = 15),
	  axis.title.y = element_text (size = 15), 
	  legend.position = "none") +
geom_hline(yintercept = 0.5, linetype = "dashed", color = "red")
#ggtitle("DIABLO (BAL Microbiome + Untargeted Metabolites) Classification Accuracy by Clinical Measure, All Ever-Smokers")
dev.off()