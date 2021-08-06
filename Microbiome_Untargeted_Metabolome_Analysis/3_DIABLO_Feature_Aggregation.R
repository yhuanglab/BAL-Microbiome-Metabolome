#applying TCGA case study from mixOmics website to our OTU table and metabolite data in SPIROMICS
#creates dynamic training and testing sets and generates DIABLO models as well as output files containing model parameters
#used to model a clinical parameter of interest using metabolome and microbiome data (can add in additional datasets)
suppressMessages(library(mixOmics))
suppressMessages(library(dplyr))
suppressMessages(library(snow))
suppressMessages(library(tidyr))
suppressMessages(library(stringr))
suppressMessages(library(reshape2))
suppressMessages(library(caret))
#suppressMessages(library(car))
options(warn = -1)

response_variable <- str_remove(getwd(), "/home/madapoos/Spiromics_No_Severe_No_Healthy_Reg/Output/") #remove path and just get target variable name

coef_pred_df <- read.csv(paste("/home/madapoos/Spiromics_No_Severe_No_Healthy_Reg/Output/", response_variable, "/DIABLO_Model_Features_ORs_P_Values_Raw.csv", sep = ""))

#list of features that were significant in at least one logistic reg model 
#(Wald p < 0.05; measure of model importance)
coef_pred_df_sig <- coef_pred_df %>% 
						select(Feature) %>% 
								group_by(Feature) %>%
									summarize(Number_of_Final_Models_Where_Feature_Selected = n())

#for most important features (see above), aggregate coefficients, ORs (convert from log-odds),
#and p-values from the models in which that feature was significant
coef_pred_df_sig_aggregated <- coef_pred_df %>% 
									group_by(Feature) %>% 
										summarize(Median_Log_OR_In_Final_Models_Where_Feature_Selected = median(Log_Odds_Ratio),
									  			  Median_Univariate_Log_OR_In_Final_Models_Where_Feature_Selected = median(Univariate.Demographic.Adjusted.Log.OR),
									  			  Median_Univariate_95_CI_Lower_Bound_In_Final_Models_Where_Feature_Selected = median(Univariate.Demographic.Adjusted.Log.OR.95..CI.Lower.Bound),
									  			  Median_Univariate_95_CI_Upper_Bound_In_Final_Models_Where_Feature_Selected = median(Univariate.Demographic.Adjusted.Log.OR.95..CI.Upper.Bound),
									  			  Median_Univariate_P_Value_In_Final_Models_Where_Feature_Selected = median(Univariate.Demographic.Adjusted.Log.OR.Wald.Test.P.Value))
coef_pred_df_sig_aggregated <-  coef_pred_df_sig_aggregated %>% left_join(coef_pred_df_sig)

write.csv(coef_pred_df_sig_aggregated, paste("/home/madapoos/Spiromics_No_Severe_No_Healthy_Reg/Output/", response_variable, "/DIABLO_Model_Features_ORs_P_Values.csv", sep = ""), row.names = F)
