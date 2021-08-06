library(dplyr)
library(tidyverse)

#read clinical files (overall files used for modeling)
framework_1 <- read.csv("~/Desktop/Spiromics_BAL_Microbiome_Metabolome_V2/Spiromics_Colorado_Stratified/Spiromics_No_Severe_No_Healthy/Data/Clinical_Parameters_04_21_20_No_Severe_No_Healthy.csv", stringsAsFactors = FALSE)
framework_2 <- read.csv("~/Desktop/Spiromics_BAL_Microbiome_Metabolome_V2/Spiromics_UNC_Stratified/Spiromics_No_Severe_No_Healthy/Data/Esther_Clinical_04_21_20_No_Severe_No_Healthy.csv", stringsAsFactors = FALSE)

#test white vs. non-white
framework_2$RACE <- ifelse(framework_2$RACE == 1, 1, 0)
framework_1$RACE <- ifelse(framework_1$RACE == 1, 1, 0)

print(nrow(all))
print(nrow(framework_1))
print(nrow(framework_2))

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
						  "SIX_MINUTE_WALK_DISTANCE01",
						  "COPDSCORE_ctb",
						  "SGR_TOTALSCORE01",
						  "BAL_Neutrophils.y",
						  "V1_PERCENT_FSAD_TOTAL.y",
						  "PI10_WHOLE_TREE_LEQ20_V2")
#categorical variables
categorical_variables <- c("RACE", 
						   "GENDER", 
						   "CURRENT_SMOKER_ctb",
						   "GOLD_STAGE_COPD_SEVERITY_1_12_0_Not")

framework_1 <- framework_1[, which(colnames(framework_1) %in% c(continuous_variables, categorical_variables))]
framework_2 <- framework_2[, which(colnames(framework_2) %in% c(continuous_variables, categorical_variables))]

list_of_non_normal_variables <- c("SMOKING_PACK_YEARS_ctb", "FVC_BDRESPONSE_PCT_ctb", "COPDSCORE_ctb", "BAL_Neutrophils.y")

#*******************#
#Statistical Testing#
#*******************#

#data frame to store p-values from Wilcoxon, T-Test, or Fisher's Exact Test
test_results <- data.frame(colnames(framework_2[, continuous_variables]), 
						   rep(0, (length(continuous_variables))), 
						   rep(0, (length(continuous_variables))))
rownames(test_results) <- test_results[, 1]
test_results <- test_results[, -1] #remov variable names since they're rownames

colnames(test_results) <- c("Test", "P-Value")

#test shapiro-wilk and wilcoxon or t-test based on SW result
for (i in 1:length(continuous_variables)) {
	variable <- continuous_variables[i]
	if (variable %in% list_of_non_normal_variables) { #normality assumption violated
		test_results[variable, 1] <- "Wilcoxon"
		test_results[variable, 2] <- wilcox.test(framework_2[, variable], framework_1[, variable])$p.value
	}
	else {
		test_results[variable, 1] <- "T-Test"
		test_results[variable, 2] <- t.test(framework_2[, variable], framework_1[, variable])$p.value
	}
}

#create strata for Fisher Test
framework_2$stratum <- rep(1, nrow(framework_2))
framework_1$stratum <- rep(0, nrow(framework_1))
bound <- rbind(framework_2, framework_1)

#conduct test
for (i in 1:length(categorical_variables)) {
	#conduct test
	variable <- categorical_variables[i]
	test_results[variable, 1] <- "Fisher"
}
test_results["RACE", 2] <- fisher.test(bound$stratum, bound$RACE)$p.value
test_results["GENDER", 2] <- fisher.test(bound$stratum, bound$GENDER)$p.value
test_results["CURRENT_SMOKER_ctb", 2] <- fisher.test(bound$stratum, bound$CURRENT_SMOKER_ctb)$p.value
test_results["GOLD_STAGE_COPD_SEVERITY_1_12_0_Not", 2] <- fisher.test(bound$stratum, bound$GOLD_STAGE_COPD_SEVERITY_1_12_0_Not)$p.value

#adjust p-values using BH-correction for each test type
test_results$"Adjusted_P-Value" <- p.adjust(test_results$"P-Value", method = "BH")
#test_results[which(test_results$"Test" == "Wilcoxon"), "Adjusted_P-Value"] <- p.adjust(test_results[which(test_results$"Test" == "Wilcoxon"), "P-Value"], method = "BH")
#test_results[which(test_results$"Test" == "T-Test"), "Adjusted_P-Value"] <- p.adjust(test_results[which(test_results$"Test" == "T-Test"), "P-Value"], method = "BH")
#test_results[which(test_results$"Test" == "Fisher"), "Adjusted_P-Value"] <- p.adjust(test_results[which(test_results$"Test" == "Fisher"), "P-Value"], method = "BH")

#test_results$"Adjusted_P-Value" <- p.adjust(test_results[, "P-Value"], method = "BH")

test_results[, "P-Value"] <- signif(test_results[, "P-Value"], digits = 3)
test_results[, "Adjusted_P-Value"] <- signif(test_results[, "Adjusted_P-Value"], digits = 3)

highly_significant <- which(test_results[, "Adjusted_P-Value"] < 0.001)
test_results[highly_significant, "Adjusted_P-Value"] <- "<0.001"

#output to csv
write.csv(test_results, file = "~/Desktop/Spiromics_BAL_Microbiome_Metabolome_V2/Spiromics_Non_Overlapping_50_Subjects/Table_1_Statistical_Tests.csv")


