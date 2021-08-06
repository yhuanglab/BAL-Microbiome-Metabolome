#manually create a csv file with all features that appear in at least 500 models and are on the scree plot
suppressMessages(library(tidyverse))
suppressMessages(library(stringr))
suppressMessages(library(mixOmics))
suppressMessages(library(dplyr))

#aggregate demographics to show average number of models in which they are significant
df <- read.csv("~/Desktop/Spiromics_BAL_Microbiome_Metabolome_V2/Spiromics_UNC_Stratified/Spiromics_No_Severe_No_Healthy/Output/Heatmap/DIABLO_Num_Models_Heatmap/Demographic_Features_by_Num_Sig_Models.csv", stringsAsFactors = FALSE)
df <- df %>% dplyr::select(Feature, Number_of_Models) %>% 
                      dplyr::group_by(Feature) %>%
                        dplyr::summarize(Mean_Number_of_Models = mean(Number_of_Models))
colnames(df)[2] <- "Mean_Number_of_Models_Where_Feature_Significant_Across_All_Outcomes"
write.csv(df, file = "~/Desktop/Spiromics_BAL_Microbiome_Metabolome_V2/Spiromics_UNC_Stratified/Spiromics_No_Severe_No_Healthy/Output/Heatmap/DIABLO_Num_Models_Heatmap/Demographic_Features_Mean_Num_Sig_Models.csv")

df_raw <- read.csv('~/Desktop/Spiromics_BAL_Microbiome_Metabolome_V2/Spiromics_UNC_Stratified/Spiromics_No_Severe_No_Healthy/Output/Heatmap/DIABLO_Num_Models_Heatmap/Top_Features_by_Num_Sig_Models_No_Unknowns.csv', stringsAsFactors = FALSE)
met_raw <- read.csv("~/Desktop/Spiromics_BAL_Microbiome_Metabolome_V2/Spiromics_UNC_Stratified/Spiromics_No_Severe_No_Healthy/Data/Esther_Metabolites_No_Severe_No_Healthy.csv")
rownames(met_raw) <- met_raw[, 1]
met_raw <- met_raw[, -1]

bac_raw <- read.csv("~/Desktop/Spiromics_BAL_Microbiome_Metabolome_V2/Spiromics_UNC_Stratified/Spiromics_No_Severe_No_Healthy/Data/Esther_Microbiota_No_Severe_No_Healthy.csv")
rownames(bac_raw) <- bac_raw[, 1]
bac_raw <- bac_raw[, -1]
bac_raw <- logratio.transfo(as.matrix(bac_raw), logratio = 'CLR', offset = 1e-10)
bac_raw <- data.frame(bac_raw[1:nrow(bac_raw), ])

clin <- read.csv("~/Desktop/Spiromics_BAL_Microbiome_Metabolome_V2/Spiromics_UNC_Stratified/Spiromics_No_Severe_No_Healthy/Data/Esther_Clinical_04_21_20_No_Severe_No_Healthy.csv")
rownames(clin) <- clin[, 1]
clin <- clin[, -1]
clin <- clin[intersect(rownames(clin), rownames(bac_raw)), ] #subset to only clinical data with matching microbiome data

#direction of effect based on OR (OR > 1 implies positive direction of effect)
df_raw$Direction_of_Effect <- ifelse(df_raw$Median_Log_OR_In_Final_Models_Where_Feature_Selected > 0, 1, -1)

#direction of association between OTU or metabolite and the outcome based on what you would see in a boxplot
#i.e. do people with exacerbations have higher levels of the outcome 
#df_raw$Direction_of_Effect <- NULL
#for (i in 1:nrow(df_raw)) {
#  outcome_of_interest <- df_raw[i, "Phenotype"]
#  feature_of_interest <- df_raw[i, "Feature"]
#  print(outcome_of_interest)
#  clinical <- clin[, as.character(outcome_of_interest)]

  #select outcome from met/otu table
#  if (str_sub(df_raw[i, "Feature"], start = 1, end = 3) == "Otu") {
#   bac <- bac_raw[, as.character(feature_of_interest)]
#   merged <- data.frame(bac, clinical)
#   summarized <- merged %>% arrange(clinical) %>% group_by(clinical) %>% summarize(mean = median(bac))
#    if (summarized[[1, 2]] < summarized[[2, 2]]) {
#      df_raw[i, "Direction_of_Effect"] <- 1
#    }
#    else {
#      df_raw[i, "Direction_of_Effect"] <- -1
#    }
#  }
#  else {
#   met <- met_raw[, as.character(feature_of_interest)]
#   merged <- data.frame(met, clinical)
#   summarized <- merged %>% arrange(clinical) %>% group_by(clinical) %>% summarize(mean = median(na.omit(met)))
#    if (summarized[[1, 2]] < summarized[[2, 2]]) {
#      df_raw[i, "Direction_of_Effect"] <- 1
#    }
#    else {
#      df_raw[i, "Direction_of_Effect"] <- -1
#    }
#  }
#}
#merge rows that are the same measure (i.e. FVC BDR % and mL merged)
df_raw$Phenotype_Figure_Name[which(df_raw$Phenotype_Figure_Name == "FVC BDR (%)")] <- "FVC BDR (% or mL)"
df_raw$Phenotype_Figure_Name[which(df_raw$Phenotype_Figure_Name == "FVC BDR (mL)")] <- "FVC BDR (% or mL)"
df_raw$Phenotype_Figure_Name[which(df_raw$Phenotype_Figure_Name == "FEV1 % of predicted")] <- "FEV1 (%predicted or L)"
df_raw$Phenotype_Figure_Name[which(df_raw$Phenotype_Figure_Name == "FEV1 (L)")] <- "FEV1 (%predicted or L)"

#for each merged phenotype, 
#check for features that are duplicated 
#(i.e. if a specific OTU is predictive of both FVC BDR % and mL) and pick highest magnitude one
if (nrow(df_raw[which(df_raw$Phenotype_Figure_Name == "FVC BDR (% or mL)"), ]) > 0) {
  fvc_bdr_pct_ml <- df_raw[which(df_raw$Phenotype_Figure_Name == "FVC BDR (% or mL)"), ]
  fvc_bdr_pct_ml <- fvc_bdr_pct_ml[order(fvc_bdr_pct_ml$Feature, -abs(as.numeric(fvc_bdr_pct_ml$Number_of_Models))), ]
  fvc_bdr_pct_ml <- fvc_bdr_pct_ml[!duplicated(fvc_bdr_pct_ml$Feature), ] #drop the row(s) with lower # of models
  df_raw <- df_raw[-which(df_raw$Phenotype_Figure_Name == "FVC BDR (% or mL)"), ] #drop original rows with duplicates in main data table
  df_raw <- rbind(df_raw, fvc_bdr_pct_ml) #append cleaned table to original data table
} 
if (nrow(df_raw[which(df_raw$Phenotype_Figure_Name == "FEV1 (%predicted or L)"), ]) > 0) {
  fvc_bdr_pct_ml <- df_raw[which(df_raw$Phenotype_Figure_Name == "FEV1 (%predicted or L)"), ]
  fvc_bdr_pct_ml <- fvc_bdr_pct_ml[order(fvc_bdr_pct_ml$Feature, -abs(as.numeric(fvc_bdr_pct_ml$Number_of_Models))), ]
  fvc_bdr_pct_ml <- fvc_bdr_pct_ml[!duplicated(fvc_bdr_pct_ml$Feature), ] #drop the row(s) with lower # of models
  df_raw <- df_raw[-which(df_raw$Phenotype_Figure_Name == "FEV1 (%predicted or L)"), ] #drop original rows with duplicates in main data table
  df_raw <- rbind(df_raw, fvc_bdr_pct_ml) #append cleaned table to original data table
}

#output to CSV file for heatmaps and network models
write.csv(df_raw, file = "~/Desktop/Spiromics_BAL_Microbiome_Metabolome_V2/Spiromics_UNC_Stratified/Spiromics_No_Severe_No_Healthy/Output/Heatmap/DIABLO_Num_Models_Heatmap/Top_Features_by_Num_Sig_Models_No_Unknowns_Dir_Effect_Added.csv", row.names = FALSE)

#threshold features at 500 models or more
#df_raw_500 <- df_raw %>% filter(Number_of_Models >= 500)
#write.csv(df_raw_500, file = "~/Desktop/Spiromics_BAL_Microbiome_Metabolome_V2/Spiromics_UNC_Stratified/Spiromics_No_Severe_No_Healthy/Output/Heatmap/DIABLO_Num_Models_Heatmap/DIABLO_Results_Num_Models_Heatmap_Data_Dir_Effect_Added_500_Above.csv", row.names = FALSE)

#threshold features at 250 models or more
#df_raw_250 <- df_raw %>% filter(Number_of_Models >= 250)
#write.csv(df_raw_250, file = "~/Desktop/Spiromics_BAL_Microbiome_Metabolome_V2/Spiromics_UNC_Stratified/Spiromics_No_Severe_No_Healthy/Output/Heatmap/DIABLO_Num_Models_Heatmap/DIABLO_Results_Num_Models_Heatmap_Data_Dir_Effect_Added_250_Above.csv", row.names = FALSE)


