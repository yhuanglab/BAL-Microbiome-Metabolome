
clin <- read.csv(file="~/Desktop/Spiromics_PiPhillin/Data/Master_Metadata_11_13_20.csv", row.names = 1)

clin$PFV45_DERV_ctb_1_Aboveequal_0_Below_Median <- ifelse(clin$PFV45_DERV_ctb >= median(na.omit(clin$PFV45_DERV_ctb)), 1, 0) 
clin$PFV53_DERV_ctb_1_Aboveequal_0_Below_Median <- ifelse(clin$PFV53_DERV_ctb >= median(na.omit(clin$PFV53_DERV_ctb)), 1, 0) 
clin$PFV62_DERV_ctb_1_Aboveequal_0_Below_Median <- ifelse(clin$PFV62_DERV_ctb >= median(na.omit(clin$PFV62_DERV_ctb)), 1, 0) 
clin$BRON_PEX_HCUDRUG_365_A_1_Oneormore_0_None <- ifelse(clin$BRON_PEX_HCUDRUG_365_A >= 1, 1, 0)
clin$GOLD_STAGE_COPD_SEVERITY_1_12_0_Not <- ifelse(clin$GOLD_STAGE_COPD_SEVERITY %in% c(1, 2), 1, 0)
clin$PCT_POST_FEV1_ctb_1_Atabove_80_0_Below_80 <- ifelse(clin$PCT_POST_FEV1_ctb >= 80, 1, 0)
clin$POST_FEV1FVC_DERV_ctb_1_Atabove_0.7_0_Below_0.7 <- ifelse(clin$POST_FEV1FVC_DERV_ctb >= 0.7, 1, 0)
clin$COPDSCORE_12_ctb <- ifelse(clin$COPDSCORE_ctb >= 10, 1, 0)
clin$PEX_HCUDRUG_CAT_365_1_Oneormore_0_None <- ifelse(clin$PEX_HCUDRUG_CAT_365 >= 1, 1, 0)
clin$PEX_TOT0101 <- ifelse(clin$PEX_TOT0101 >= 1, 1, 0)
clin$FEV1_BDRESPONSE_HANKINSON_ctb_1_12above <- ifelse(clin$FEV1_BDRESPONSE_HANKINSON_ctb >= 12, 1, 0)
clin$FVC_BDRESPONSE_PCT_ctb_1_12above <- ifelse(clin$FVC_BDRESPONSE_PCT_ctb >= 12, 1, 0)
clin$GOLD_STAGE_COPD_SEVERITY_1_12_0_Not <- ifelse(clin$FVC_BDRESPONSE_VOL_ctb >= 200, 1, 0)

write.csv(clin, file = "~/Desktop/Spiromics_PiPhillin/Data/Master_Metadata_11_13_20_Updated.csv")