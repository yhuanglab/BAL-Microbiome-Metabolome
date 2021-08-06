library(dplyr)
library(tidyverse)

all_metadata <- read.csv("~/Desktop/Spiromics_BAL_Microbiome_Metabolome_V2/Spiromics_All_137_Subjects/Master_Metadata_11_13_20.csv")
all_metadata$GOLD_STAGE_COPD_SEVERITY_1_12_0_Not <- ifelse(all_metadata$GOLD_STAGE_COPD_SEVERITY %in% c(1, 2), 1, 0)
all_metadata$BAL_Neutrophils_Aboveequal_Below_Mean <- ifelse(all_metadata$BAL_Neutrophils >= mean(na.omit(all_metadata$BAL_Neutrophils)), 1, 0)
colnames(all_metadata)[which(colnames(all_metadata) == "BAL_Neutrophils")] <- "BAL_Neutrophils.y"
framework_1_metadata <- read.csv("~/Desktop/Spiromics_BAL_Microbiome_Metabolome_V2/Spiromics_Colorado_Stratified/Spiromics_No_Severe_No_Healthy/Data/Clinical_Parameters_04_21_20_No_Severe_No_Healthy.csv", stringsAsFactors = FALSE)
framework_2_metadata <- read.csv("~/Desktop/Spiromics_BAL_Microbiome_Metabolome_V2/Spiromics_UNC_Stratified/Spiromics_No_Severe_No_Healthy/Data/Esther_Clinical_04_21_20_No_Severe_No_Healthy.csv", stringsAsFactors = FALSE)

#test white vs. non-white
framework_1_metadata$RACE <- ifelse(framework_1_metadata$RACE == 1, 1, 0)
framework_2_metadata$RACE <- ifelse(framework_2_metadata$RACE == 1, 1, 0)

#continuous variables
continuous_variables <- c("AGE_DERV_ctb",
						  "SMOKING_PACK_YEARS_ctb",
						  "POST_FEV1FVC_DERV_ctb",
						  "PCT_POST_FEV1_ctb",
						  "PCT_POST_FVC_ctb",
						  "PFV53_DERV_ctb",
						  "PFV62_DERV_ctb",
						  "FVC_BDRESPONSE_PCT_ctb",
						  "FEV1_BDRESPONSE_HANKINSON_ctb",
						  "COPDSCORE_ctb",
						  "BAL_Neutrophils.y",
						  "V1_PERCENT_FSAD_TOTAL.y",
						  "PI10_WHOLE_TREE_LEQ20_V2")
#categorical variables
categorical_variables <- c("RACE", 
						   "GENDER", 
						   "CURRENT_SMOKER_ctb",
						   "GOLD_STAGE_COPD_SEVERITY_1_12_0_Not")

#*************#
#Summary Stats#
#*************#

summary_stats <- data.frame(colnames(framework_2_metadata[, continuous_variables]), #empty data frame to store mean and SD
							rep(0, (length(continuous_variables))),
							rep(0, (length(continuous_variables))),
							rep(0, (length(continuous_variables))),
							rep(0, (length(continuous_variables))),
							rep(0, (length(continuous_variables))),
							rep(0, (length(continuous_variables))))
rownames(summary_stats) <- summary_stats[, 1]
summary_stats <- summary_stats[, -1]
colnames(summary_stats) <- c("All_Mean", "All_SD", "Framework_1_Mean", "Framework_1_SD", "Framework_2_Mean", "Framework_2_SD")

#continuous variables
for (i in continuous_variables) {
	print(i)
	hist(all_metadata[, i])
	summary_stats[i, 1] <- median(na.omit(all_metadata[, i]))
	summary_stats[i, 2] <- IQR(na.omit(all_metadata[, i]))
	hist(framework_1_metadata[, i])
	summary_stats[i, 3] <- median(na.omit(framework_1_metadata[, i]))
	summary_stats[i, 4] <- IQR(na.omit(framework_1_metadata[, i]))
	hist(framework_2_metadata[, i])
	summary_stats[i, 5] <- median(na.omit(framework_2_metadata[, i]))
	summary_stats[i, 6] <- IQR(na.omit(framework_2_metadata[, i]))
}

summary_stats <- rownames_to_column(summary_stats, var = "Clinical_Measure") #convert rownames to column to work in tibble format
 
#merge mean and SD in the same column separated by parentheses
summary_stats <- summary_stats %>% mutate(All_Mean_SD = paste(signif(All_Mean, digits = 3), " (", signif(All_SD, digits = 3), ")", sep = ""),
										  Framework_1_Mean_SD = paste(signif(Framework_1_Mean, digits = 3), " (", signif(Framework_1_SD, digits = 3), ")", sep = ""),
										  Framework_2_Mean_SD = paste(signif(Framework_2_Mean, digits = 3), " (", signif(Framework_2_SD, digits = 3), ")", sep = ""))
summary_stats <- summary_stats %>% select(Clinical_Measure, All_Mean_SD, Framework_1_Mean_SD, Framework_2_Mean_SD)

#write summary stats for continuous variables to table to copy paste into figure
write.csv(summary_stats, file = "~/Desktop/Spiromics_BAL_Microbiome_Metabolome_V2/Spiromics_All_137_Subjects/Table_1_Continuous_Summary_Statistics.csv", row.names = FALSE)

summary_stats_2 <- data.frame(colnames(framework_2_metadata[, categorical_variables]),
							  rep(0, length(categorical_variables)),
							   rep(0, length(categorical_variables)),
							    rep(0, length(categorical_variables)))
rownames(summary_stats_2) <- summary_stats_2[, 1]
summary_stats_2 <- summary_stats_2[, -1]
colnames(summary_stats_2) <- c("All %", "Framework_1 %", "Framework_2 %")

#categorical variables - get % white, male, etc.
for (i in categorical_variables) {
	summary_stats_2[i, 1] <- paste(signif(length(which(all_metadata[, i] == "1"))/length(na.omit(all_metadata[, i])), digits = 3)*100, "%", sep = "")
	summary_stats_2[i, 2] <- paste(signif(length(which(framework_1_metadata[, i] == "1"))/length(na.omit(framework_1_metadata[, i])), digits = 3)*100, "%", sep = "")
	summary_stats_2[i, 3] <- paste(signif(length(which(framework_2_metadata[, i] == "1"))/length(na.omit(framework_2_metadata[, i])), digits = 3)*100, "%", sep = "")
}

write.csv(summary_stats_2, file = "~/Desktop/Spiromics_BAL_Microbiome_Metabolome_V2/Spiromics_All_137_Subjects/Table_1_Categorical_Summary_Statistics.csv")
