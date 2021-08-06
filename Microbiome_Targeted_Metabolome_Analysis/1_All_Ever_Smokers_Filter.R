library(dplyr)

clin <- read.csv('~/Desktop/Spiromics_BAL_Microbiome_Metabolome_V2/Spiromics_UNC_Stratified/Spiromics_No_Severe_No_Healthy/Data/Esther_Clinical_04_21_20.csv')
rownames(clin) <- clin[, 1]

nsns <- rownames(clin)[which(clin$GOLD_STAGE_COPD_SEVERITY %in% c(0, 1, 2))]

clin <- clin[nsns, ]

clin$MUCIN_CONCENTRATION_V1.y_1_Aboveequal_0_Below_Median <- ifelse(clin$MUCIN_CONCENTRATION_V1.y > median(na.omit(clin$MUCIN_CONCENTRATION_V1.y)), 1, 0) 
clin$PFV45_DERV_ctb_1_Aboveequal_0_Below_Median <- ifelse(clin$PFV45_DERV_ctb >= median(na.omit(clin$PFV45_DERV_ctb)), 1, 0) 
clin$PFV53_DERV_ctb_1_Aboveequal_0_Below_Median <- ifelse(clin$PFV53_DERV_ctb >= median(na.omit(clin$PFV53_DERV_ctb)), 1, 0) 
clin$PFV62_DERV_ctb_1_Aboveequal_0_Below_Median <- ifelse(clin$PFV62_DERV_ctb >= median(na.omit(clin$PFV62_DERV_ctb)), 1, 0) 
clin$BRON_PEX_HCUDRUG_365_A.y_1_Oneormore_0_None <- ifelse(clin$BRON_PEX_HCUDRUG_365_A.y >= 1, 1, 0)
clin$GOLD_STAGE_COPD_SEVERITY_1_12_0_Not <- ifelse(clin$GOLD_STAGE_COPD_SEVERITY %in% c(1, 2), 1, 0)
clin$PCT_POST_FEV1_ctb_1_Atabove_80_0_Below_80 <- ifelse(clin$PCT_POST_FEV1_ctb >= 80, 1, 0)
clin$POST_FEV1FVC_DERV_ctb_1_Atabove_0.7_0_Below_0.7 <- ifelse(clin$POST_FEV1FVC_DERV_ctb >= 0.7, 1, 0)
clin$FVC_BDRESPONSE_VOL_ctb_1_Atabove_200_0_Below_200 <- ifelse(clin$FVC_BDRESPONSE_VOL_ctb >= 200, 1, 0)
clin$COPDSCORE_12 <- ifelse(clin$COPDSCORE_ctb >= 10, 1, 0)

#merge inhaled steroids ctb to clin
inhaled_steroids <- read.csv("~/Desktop/Spiromics_BAL_Microbiome_Metabolome_V2/Spiromics_All_137_Subjects/INHALEDSTEROIDS_ctb.csv", stringsAsFactors = FALSE)
clin <- left_join(clin, inhaled_steroids, by = "SUBJID")

#merge antibiotic use in the 3 months before bronchoscopy to clin
antibiotic_use <- read.csv("~/Desktop/Spiromics_BAL_Microbiome_Metabolome_V2/Spiromics_All_137_Subjects/SP0959_ANTIUSE_BEF_BRON_210204.csv", stringsAsFactors = FALSE)
antibiotic_use <- antibiotic_use %>% select(SUBJID, ANTI_USE_3MTH_BEF_BRON)
clin <- left_join(clin, antibiotic_use, by = "SUBJID")

#merge one of the CT columns to to clin
ct_data_1 <- read.csv("~/Desktop/Spiromics_BAL_Microbiome_Metabolome_V2/Spiromics_All_137_Subjects/Metadata_06132019_new_Faith.csv", stringsAsFactors = FALSE)
ct_data_1 <- ct_data_1 %>% select(SUBJID, V1_PERCENT_FSAD_TOTAL)
clin <- left_join(clin, ct_data_1, by = "SUBJID")
clin$V1_PERCENT_FSAD_TOTAL_Aboveequal_Below_Median <- ifelse(clin$V1_PERCENT_FSAD_TOTAL.y >= median(na.omit(clin$V1_PERCENT_FSAD_TOTAL.y)), 1, 0)

#merge another CT column to clin
ct_data <- read.csv("~/Desktop/Spiromics_BAL_Microbiome_Metabolome_V2/Spiromics_All_137_Subjects/Core6_CT_190903.csv", stringsAsFactors = FALSE)
ct_data <- ct_data %>% select(SUBJID, PI10_WHOLE_TREE_LEQ20_V2)
clin <- left_join(clin, ct_data, by = "SUBJID")
clin$PI10_WHOLE_TREE_LEQ20_V2_Aboveequal_Below_Median <- ifelse(clin$PI10_WHOLE_TREE_LEQ20_V2 >= median(na.omit(clin$PI10_WHOLE_TREE_LEQ20_V2)), 1, 0)

#merge BAL neutrophils and eosinophils to clin
bal_neutrophils <- read.csv("~/Desktop/Spiromics_BAL_Microbiome_Metabolome_V2/Spiromics_UNC_Stratified/Spiromics_No_Severe_No_Healthy/Data/BAL_Neutrophils.csv", stringsAsFactors = FALSE)
clin <- left_join(clin, bal_neutrophils, by = "SUBJID")
clin$BAL_Neutrophils_Aboveequal_Below_Median <- ifelse(clin$BAL_Neutrophils.y >= median(na.omit(clin$BAL_Neutrophils.y)), 1, 0)
clin$BAL_Eosinophils_Aboveequal_Below_Median <- ifelse(clin$BAL_Eosinophils.y >= median(na.omit(clin$BAL_Eosinophils.y)), 1, 0)

write.csv(clin, file = '~/Desktop/Spiromics_BAL_Microbiome_Metabolome_V2/Spiromics_UNC_Stratified/Spiromics_No_Severe_No_Healthy/Data/Esther_Clinical_04_21_20_No_Severe_No_Healthy.csv')

bac <- read.csv('~/Desktop/Spiromics_BAL_Microbiome_Metabolome_V2/Spiromics_All_137_Subjects/All_137_Subjects_OTU_Relative_Abundance_Table_11_13_20.csv', stringsAsFactors = FALSE)
rownames(bac) <- bac[, 1]
bac <- bac[, -1]
bac <- bac[intersect(rownames(bac), clin$X), ]
write.csv(bac, file = '~/Desktop/Spiromics_BAL_Microbiome_Metabolome_V2/Spiromics_UNC_Stratified/Spiromics_No_Severe_No_Healthy/Data/Esther_Microbiota_No_Severe_No_Healthy.csv')

met <- read.csv('~/Desktop/Spiromics_BAL_Microbiome_Metabolome_V2/Spiromics_UNC_Stratified/Spiromics_No_Severe_No_Healthy/Data/Esther_Metabolites.csv', check.names = FALSE)
rownames(met) <- met[, 1]
met <- met[, -1]
met <- met[intersect(rownames(met), clin$X), ]
#replace missing values in metabolites with the colmin of the metabolite
for (i in 1:ncol(met)) {
  if (length(which(is.na(met[, i]))) > 0) {
    met[which(is.na(met[, i])), i] <- min(na.omit(met[, i]))
  }
}
met <- log2(met)
write.csv(met, file = '~/Desktop/Spiromics_BAL_Microbiome_Metabolome_V2/Spiromics_UNC_Stratified/Spiromics_No_Severe_No_Healthy/Data/Esther_Metabolites_No_Severe_No_Healthy.csv')
