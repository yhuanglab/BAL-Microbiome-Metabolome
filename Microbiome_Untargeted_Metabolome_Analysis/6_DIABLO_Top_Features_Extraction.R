#extract features above the first 25% drop
#for all models that were significant (average p < 0.05)
suppressMessages(library(dplyr))
suppressMessages(library(snow))
suppressMessages(library(tidyr))
suppressMessages(library(stringr))
suppressMessages(library(reshape2))
options(warn = -1)

cut_point <- 25 #features need to be selected by at least this many models to be shown in the heatmaps

target_variable <- str_remove(getwd(), "/Users/Siddharth/Desktop/Spiromics_BAL_Microbiome_Metabolome_V2/Spiromics_Colorado_Stratified/Spiromics_No_Severe_No_Healthy/Output/DIABLO/") #remove path and just get target variable name

#check if overall model is significantly predictive out-of-sample
model_sig <- read.csv("DIABLO_Model_Parameters.csv", stringsAsFactors = FALSE)

if (model_sig[12, "Mean.Value.Over.100.Resampling.Iterations"] > 0.05) {
	print("Model is not significantly predictive, so the top features are not extracted")
} else {
	#list of OTUs and metabolites significant in at least one final logistic regression model
	top_features <- read.csv("DIABLO_Model_Features_ORs_P_Values.csv", stringsAsFactors = FALSE)

	#top_features <- read.csv('~/Desktop/Spiromics_BAL_Microbiome_Metabolome_V2/Spiromics_Colorado_Stratified/Spiromics_No_Severe_No_Healthy/Output/DIABLO/FVC_BDRESPONSE_PCT_ctb_1_12above/DIABLO_Model_Features_ORs_P_Values.csv')

	#order features by descending order with respect to number of final logistic regression models
	#in which they were significant - importance metric
	top_features <- top_features %>% arrange(desc(Number_of_Final_Models_Where_Feature_Selected))

	#create spreadsheet with demographic variables and number of models they were significant in for comparison
	demographic_variables <- c("(Intercept)", "GENDER2", "RACE_Nonwhite1", "AGE_DERV_ctb", "CURRENT_SMOKER_ctb1", "INHALEDSTEROIDS_ctb1", "ANTI_USE_3MTH_BEF_BRON1") 
	top_features_demographics <- top_features[which(top_features$Feature %in% demographic_variables), ]

	#aggregate
	top_features_demographics <- top_features_demographics %>% select(Feature, 
																		Number_of_Final_Models_Where_Feature_Selected,
																		Median_Log_OR_In_Final_Models_Where_Feature_Selected)
	
	labels <- c("Group", "Phenotype", "Phenotype_Figure_Name", "Feature", "Number_of_Models", "Median_Log_OR_In_Final_Models_Where_Feature_Significant") #column names for final data table
	#identify the category of the outcome modeled
	outcome_names_categories_for_figures <- read.csv("~/Desktop/Spiromics_BAL_Microbiome_Metabolome_V2/Outcome_Names_Categories_for_Figures.csv")
	group <- outcome_names_categories_for_figures[which(outcome_names_categories_for_figures$Model == target_variable), "Category"] %>% as.character()

	#identify the name we're using for the outcome in the figure
	phenotype_figure_name <- outcome_names_categories_for_figures[which(outcome_names_categories_for_figures$Model == target_variable), "Model_Name"] %>% as.character()

	#create data frame with outcome name, figure name, category, crude error, and crude accuracy for each outcome modeled
	group_column <- rep(group, nrow(top_features_demographics))
	phenotype_column <- rep(target_variable, nrow(top_features_demographics))
	phenotype_figure_name_column <- rep(phenotype_figure_name, nrow(top_features_demographics))
	frame <- cbind(group_column, phenotype_column, phenotype_figure_name_column, top_features_demographics)

	if (target_variable == "BRON_PEX_HCUDRUG_365_A.y_1_Oneormore_0_None") {
		df <- frame
		colnames(df) <- labels	
		rownames(df) <- NULL

		#append to table of all outcomes
		write.table(df, file = "~/Desktop/Spiromics_BAL_Microbiome_Metabolome_V2/Spiromics_Colorado_Stratified/Spiromics_No_Severe_No_Healthy/Output/Heatmap/DIABLO_Num_Models_Heatmap/Demographic_Features_by_Num_Sig_Models.csv", 
						sep = ",", 
						append = TRUE, 
						row.names = FALSE)
	} else {
		df <- rbind(labels, frame) %>% data.frame()
		df <- df[-1, , drop = FALSE]
		write.table(df, file = "~/Desktop/Spiromics_BAL_Microbiome_Metabolome_V2/Spiromics_Colorado_Stratified/Spiromics_No_Severe_No_Healthy/Output/Heatmap/DIABLO_Num_Models_Heatmap/Demographic_Features_by_Num_Sig_Models.csv",
						sep = ",", 
						append = TRUE, 
						row.names = FALSE, 
						col.names = FALSE)
	}

	#drop demographic variables
	top_features <- top_features[-which(top_features$Feature %in% demographic_variables), ]

	#drop rows with unidentified matabolites (X###.####.####)
	if (nrow(top_features[grep('^X\\d{3}.', top_features$Feature), ])) {
		top_features <- top_features[-grep('^X\\d{3}.', top_features$Feature), ]
	}

	#filter to OTU features first, create sigmoid, and select top features above the cutpoint
	otus_list <- top_features$Feature[which(str_sub(top_features$Feature, start = 1, end = 3) == "Otu"), drop = FALSE]
	top_features_otus <- top_features %>% 
							filter(Feature %in% otus_list) %>%
								filter(Number_of_Final_Models_Where_Feature_Selected >= cut_point)
	#top_features_otus <- top_features %>% 
	#						filter(Feature %in% otu_list) %>% 
	#							mutate(pct_change = (Number_of_Models_Where_Feature_Significant - lag(Number_of_Models_Where_Feature_Significant, 
	#									default = first(Number_of_Models_Where_Feature_Significant))) / lag(Number_of_Models_Where_Feature_Significant, default = first(Number_of_Models_Where_Feature_Significant)))
	#top_features_otus$pct_change <- abs(top_features_otus$pct_change * 100)
	#identify the features above the first 25%+ drop in feature importance
	#these are the most significant features and on the plateau of the sigmoid
	#drop_in_features <- which(top_features_otus$pct_change >= 33) #filter to cutpoint or more
	#if (length(drop_in_features) == 0) { #if the features never drop by 25% or more (no significant drop)
	#	if (length(which(top_features_otus$pct_change == 0)) == 0) { #check if the features never drop at all (all at 1)
	#		features_above_cut_point <- NULL #don't include features from this clinical outcome if all features are only found in 1 model each
	#	} else {
	#		drop_in_features <- which(top_features_otus$pct_change > 0) #identify the cutpoints, albeit less than 25%
	#		first_drop <- drop_in_features[1] #filter to the first cutpoint
	#		features_above_cut_point <- first_drop - 1
	#	}
	#} else { #if the features drop by 25% or more
		#if (length(drop_in_features > 1)) { # if there are multiple drops (ideal), take features above the second drop of 25%+ to maximize features
		#	second_drop <- drop_in_features[2] 
		#	features_above_cut_point <- second_drop - 1
		# #else { #if there is only one drop, look at all features above that
		#	first_drop <- drop_in_features[1]
		#	features_above_cut_point <- first_drop - 1
		#}
	#	first_drop <- drop_in_features[1]
	#	features_above_cut_point <- first_drop - 1
	#}
	labels <- c("Group", "Phenotype", "Phenotype_Figure_Name", "Feature", "Number_of_Models", "Median_Log_OR_In_Final_Models_Where_Feature_Selected") #column names for final data table

	if (is.null(top_features_otus)) {
		print(paste("No OTU features predictive of ", target_variable, "were selected in 25% of models or more so this outcome is not included"))
		if (target_variable == "BRON_PEX_HCUDRUG_365_A.y_1_Oneormore_0_None") { #write column names in case first model has no features
			df <- data.frame(t(labels))
			write.table(df, file = "~/Desktop/Spiromics_BAL_Microbiome_Metabolome_V2/Spiromics_Colorado_Stratified/Spiromics_No_Severe_No_Healthy/Output/Heatmap/DIABLO_Num_Models_Heatmap/Top_Features_by_Num_Sig_Models_No_Unknowns.csv", 
							sep = ",", 
							append = TRUE, 
							row.names = FALSE, 
							col.names = FALSE)
		}
	} else { #if there is a sigmoidal structure and we selected features, append features to the csv file
		#top_features_above_cutpoint <- top_features_otus[1:features_above_cut_point, ] #filter feature data to those above cut point from previous step
		
		#select just the importance metric column
		#top_features_above_cutpoint <- top_features_above_cutpoint %>% select(Feature, 
		#																	  Number_of_Models_Where_Feature_Significant,
		#																	  Median_Log_OR_In_Models_Where_Feature_Significant)
		top_features_above_cutpoint <- top_features_otus %>% 
											select(Feature, Number_of_Final_Models_Where_Feature_Selected, Median_Log_OR_In_Final_Models_Where_Feature_Selected)

		#identify the category of the outcome modeled
		outcome_names_categories_for_figures <- read.csv("~/Desktop/Spiromics_BAL_Microbiome_Metabolome_V2/Outcome_Names_Categories_for_Figures.csv")
		group <- outcome_names_categories_for_figures[which(outcome_names_categories_for_figures$Model == target_variable), "Category"] %>% as.character()

		#identify the name we're using for the outcome in the figure
		phenotype_figure_name <- outcome_names_categories_for_figures[which(outcome_names_categories_for_figures$Model == target_variable), "Model_Name"] %>% as.character()

		#create data frame with outcome name, figure name, category, crude error, and crude accuracy for each outcome modeled
		group_column <- rep(group, nrow(top_features_above_cutpoint))
		phenotype_column <- rep(target_variable, nrow(top_features_above_cutpoint))
		phenotype_figure_name_column <- rep(phenotype_figure_name, nrow(top_features_above_cutpoint))
		frame <- cbind(group_column, phenotype_column, phenotype_figure_name_column, top_features_above_cutpoint)

		if (target_variable == "BRON_PEX_HCUDRUG_365_A.y_1_Oneormore_0_None") {
			df <- frame
			colnames(df) <- labels	
			rownames(df) <- NULL

			#append to table of all outcomes
			write.table(df, file = "~/Desktop/Spiromics_BAL_Microbiome_Metabolome_V2/Spiromics_Colorado_Stratified/Spiromics_No_Severe_No_Healthy/Output/Heatmap/DIABLO_Num_Models_Heatmap/Top_Features_by_Num_Sig_Models_No_Unknowns.csv", 
							sep = ",", 
							append = TRUE, 
							row.names = FALSE)
		} else {
			df <- rbind(labels, frame) %>% data.frame()
			df <- df[-1, , drop = FALSE]
			write.table(df, file = "~/Desktop/Spiromics_BAL_Microbiome_Metabolome_V2/Spiromics_Colorado_Stratified/Spiromics_No_Severe_No_Healthy/Output/Heatmap/DIABLO_Num_Models_Heatmap/Top_Features_by_Num_Sig_Models_No_Unknowns.csv",
							sep = ",", 
							append = TRUE, 
							row.names = FALSE, 
							col.names = FALSE)
		}
	}

	#repeat for metabolite features separately
	metabolites_list <- top_features$Feature[which(str_sub(top_features$Feature, start = 1, end = 3) != "Otu"), drop = FALSE]
	top_features_metabolites <- top_features %>% 
							filter(Feature %in% metabolites_list) %>%
								filter(Number_of_Final_Models_Where_Feature_Selected >= cut_point)
	#top_features_metabolites <- top_features %>% 
	#							filter(Feature %in% metabolites_list) %>% 
	#								mutate(pct_change = (Number_of_Models_Where_Feature_Significant - lag(Number_of_Models_Where_Feature_Significant, 
	#										default = first(Number_of_Models_Where_Feature_Significant))) / lag(Number_of_Models_Where_Feature_Significant, default = first(Number_of_Models_Where_Feature_Significant)))
	#top_features_metabolites$pct_change <- abs(top_features_metabolites$pct_change * 100)
	#identify the features above the first 25%+ drop in feature importance
	#these are the most significant features and on the plateau of the sigmoid
	#drop_in_features <- which(top_features_metabolites$pct_change >= 33) #filter to cutpoint or more
	#if (length(drop_in_features) == 0) { #if the features never drop by 25% or more (no significant drop)
	#	if (length(which(top_features_metabolites$pct_change == 0)) == 0) { #check if the features never drop at all (all at 1)
	#		features_above_cut_point <- NULL #don't include features from this clinical outcome if all features are only found in 1 model each
	#	} else {
	#		drop_in_features <- which(top_features_metabolites$pct_change > 0) #identify the cutpoints, albeit less than 25%
	#		first_drop <- drop_in_features[1] #filter to the first cutpoint
	#		features_above_cut_point <- first_drop - 1
	#	}
	#} else { #if the features drop by 25% or more
		#if (length(drop_in_features > 1)) { # if there are multiple drops (ideal), take features above the second drop of 25%+ to maximize features
		#	second_drop <- drop_in_features[2] 
		#	features_above_cut_point <- second_drop - 1
		# #else { #if there is only one drop, look at all features above that
		#	first_drop <- drop_in_features[1]
		#	features_above_cut_point <- first_drop - 1
		#}
	#	first_drop <- drop_in_features[1]
	#	features_above_cut_point <- first_drop - 1
	#}
	labels <- c("Group", "Phenotype", "Phenotype_Figure_Name", "Feature", "Number_of_Models", "Median_Log_OR_In_Final_Models_Where_Feature_Selected") #column names for final data table

	#top_features_above_cutpoint <- top_features_metabolites[1:features_above_cut_point, ] #filter feature data to those above cut point from previous step
	
	if (is.null(top_features_metabolites)) {
		print(paste("No Metabolites features predictive of ", target_variable, "were selected in 25% of models or more so this outcome is not included"))
		if (target_variable == "BRON_PEX_HCUDRUG_365_A.y_1_Oneormore_0_None") { #write column names in case first model has no features
			df <- data.frame(t(labels))
			write.table(df, file = "~/Desktop/Spiromics_BAL_Microbiome_Metabolome_V2/Spiromics_Colorado_Stratified/Spiromics_No_Severe_No_Healthy/Output/Heatmap/DIABLO_Num_Models_Heatmap/Top_Features_by_Num_Sig_Models_No_Unknowns.csv", 
							sep = ",", 
							append = TRUE, 
							row.names = FALSE, 
							col.names = FALSE)
		}
	} else {
		#select just the importance metric column
		#top_features_above_cutpoint <- top_features_above_cutpoint %>% select(Feature, 
		#																	  Number_of_Models_Where_Feature_Significant,
		#																	  Median_Log_OR_In_Models_Where_Feature_Significant)
		top_features_above_cutpoint <- top_features_metabolites %>% 
											select(Feature, Number_of_Final_Models_Where_Feature_Selected, Median_Log_OR_In_Final_Models_Where_Feature_Selected)

		#identify the category of the outcome modeled
		outcome_names_categories_for_figures <- read.csv("~/Desktop/Spiromics_BAL_Microbiome_Metabolome_V2/Outcome_Names_Categories_for_Figures.csv")
		group <- outcome_names_categories_for_figures[which(outcome_names_categories_for_figures$Model == target_variable), "Category"] %>% as.character()

		#identify the name we're using for the outcome in the figure
		phenotype_figure_name <- outcome_names_categories_for_figures[which(outcome_names_categories_for_figures$Model == target_variable), "Model_Name"] %>% as.character()

		#create data frame with outcome name, figure name, category, crude error, and crude accuracy for each outcome modeled
		group_column <- rep(group, nrow(top_features_above_cutpoint))
		phenotype_column <- rep(target_variable, nrow(top_features_above_cutpoint))
		phenotype_figure_name_column <- rep(phenotype_figure_name, nrow(top_features_above_cutpoint))
		frame <- cbind(group_column, phenotype_column, phenotype_figure_name_column, top_features_above_cutpoint)

		if (target_variable == "BRON_PEX_HCUDRUG_365_A.y_1_Oneormore_0_None") {
			df <- frame
			colnames(df) <- labels	
			rownames(df) <- NULL

			#append to table of all outcomes
			write.table(df, file = "~/Desktop/Spiromics_BAL_Microbiome_Metabolome_V2/Spiromics_Colorado_Stratified/Spiromics_No_Severe_No_Healthy/Output/Heatmap/DIABLO_Num_Models_Heatmap/Top_Features_by_Num_Sig_Models_No_Unknowns.csv", 
							sep = ",", 
							append = TRUE, 
							row.names = FALSE)
		} else {
			df <- rbind(labels, frame) %>% data.frame()
			df <- df[-1, , drop = FALSE]
			write.table(df, file = "~/Desktop/Spiromics_BAL_Microbiome_Metabolome_V2/Spiromics_Colorado_Stratified/Spiromics_No_Severe_No_Healthy/Output/Heatmap/DIABLO_Num_Models_Heatmap/Top_Features_by_Num_Sig_Models_No_Unknowns.csv",
							sep = ",", 
							append = TRUE, 
							row.names = FALSE, 
							col.names = FALSE)
		}
	}
}

