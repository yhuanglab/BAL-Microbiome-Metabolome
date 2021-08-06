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

#*****************#
#Tuning Parameters#
#*****************#

diablo_pct <- 0.1 #% of top DIABLO features based on absolute value of loading score to use for regression analysis

#****************#
#Data Preparation#
#****************#

response_variable <- str_remove(getwd(), "/home/madapoos/Spiromics_No_Severe_No_Healthy_UNC_Reg/Output/") #remove path and just get target variable name
#response_variable <- "FEV1_BDRESPONSE_HANKINSON_ctb_1_12above"

num_iterations <- 100 #number of bootstrapped training and testing sets to run model on

#read in data sets
env <- read.csv("/home/madapoos/Spiromics_No_Severe_No_Healthy_UNC_Reg/Data/Esther_Clinical_04_21_20_No_Severe_No_Healthy.csv", header=T)
met <- read.csv("/home/madapoos/Spiromics_No_Severe_No_Healthy_UNC_Reg/Data/Esther_Metabolites_No_Severe_No_Healthy.csv", header = T, stringsAsFactors = FALSE)
bac <- read.csv("/home/madapoos/Spiromics_No_Severe_No_Healthy_UNC_Reg/Data/Esther_Microbiota_No_Severe_No_Healthy.csv", header = T, stringsAsFactors = FALSE)

#CLR transform OTUs
rownames(bac) <- bac[, 1] #set names to batch number
bac <- bac[, -1]
bac <- logratio.transfo(as.matrix(bac), logratio = 'CLR', offset = 1e-10)

#Due to low sample size, create a new RACE variable as white vs. nonwhite
env$RACE_Nonwhite <- ifelse(env$RACE == 1, 0, 1)

#impute missing values for age as median of all subjects
env[which(is.na(env$AGE_DERV_ctb)), "AGE_DERV_ctb"] <- median(na.omit(env$AGE_DERV_ctb))

env <- env[, -1]
rownames(env) <- env[, 1] #set names to batch number
env <- env[, -1] 
env <- env[intersect(rownames(env), rownames(bac)), ] #subset to only clinical data with matching microbiome data
target_variable_index <- which(colnames(env) == response_variable)
colnames(env)[target_variable_index] <- "target_variable" #rename COPDSCORE12 to allow pass by reference

bac <- bac[intersect(rownames(bac), rownames(env)), ]

rownames(met) <- met[, 1] #set names to batch number
met <- met[, -1] 

met <- met[intersect(rownames(met), rownames(bac)), ] #subset to only clinical data with matching microbiome data

#drop subjects with missing data
missing_subjects <- which(is.na(env[, "target_variable"]))
if (length(missing_subjects > 0)) {
	env <- env[-missing_subjects, , drop = FALSE]
	bac <- bac[-missing_subjects, , drop = FALSE]
	met <- met[-missing_subjects, , drop = FALSE]
}

#************#
#ML Algorithm#
#************#

#Parameters to aggregate over resampling iterations
filter_num_optimal_met_list <- list(NULL) #stores optimal number of metabolite features from t-tests
filter_num_optimal_bac_list <- list(NULL) #stores optimal number of bacterial features from t-tests
num_optimal_met_list <- list(NULL) #stores optimal number of features from metabolite list
num_optimal_bac_list <- list(NULL) #stores optimal number of bacterial features
diablo_num_optimal_met_list <- list(NULL) #stores optimal number of features from metabolite list
diablo_num_optimal_bac_list <- list(NULL) #stores optimal number of bacterial features
optimal_met_list <- list(NULL) #optimal features from metabolite list
optimal_bac_list <- list(NULL) #optimal features from bac list
diablo_training_majority_vote_list <- list(NULL)
diablo_training_weighted_vote_list <- list(NULL)
overall_bar_list <- list(NULL)
coef_pred_list <- list(NULL)
alpha_list <- list(NULL)
lambda_list <- list(NULL)

predicted <- list(NULL)
confusion_matrix_list <- NULL

#*************#
#Training Data#
#*************#

train_rows <- list(NULL) #list of samples included in each training set

#Bootstrap n training data sets of size 70% from the total data set
set.seed(100)
for (k in 1:(num_iterations + 200)) { #create twice as many training sets as needed in case model fails
	training_prop <- 0.7 #proportion of total data set to be used for training
	training_samples_count <- training_prop*(nrow(env)) #number of samples to be used in training data set

	categories <- unique(env[, "target_variable"]) #categories of response variable
	training_samples <- NULL #list of samples to include in training data
	for (i in categories) { #automatic selection of observations for training data, regardless of the number of categories
		prop_cat <- length(which(env[, "target_variable"] == i))/nrow(env) #finds proportion of patients that are from that category
		training_cat_count <- prop_cat*training_samples_count #number of observations from that category to be included in training set
		cat <- which(env[, "target_variable"] == i) #identifies observations with that category in original data set
		cat_sample <- sample(cat, training_cat_count) #samples category to make up identical distribution of each category in the training set as in the original data set
		training_samples <- c(training_samples, cat_sample) #add to list of samples to extract for training data
	}
	train_rows[[k]] <- training_samples #add to list of training data sets
}

#DEFINE PARAMETERS
j <- 1
error_counter <- 0 #count number of times modeling yields an error
broken <- FALSE #record if broken from while loop due to data sparsity
error_counter_max <- num_iterations + 200 #number of times model is allowed to fail before determining that data is too sparse to improve computational efficiency for user

while (j <= num_iterations) { #model generation over all resampling iterations
	print(paste("Model: ", j, sep = ""))
	print(j + error_counter)

	#*********************#
	#Initial Training Data#
	#*********************#

	training_data <- train_rows[[j + error_counter]] #select training data set matching iteration (or future one if errors)

	#subset all datasets to the training data observations
	env_train <- env[training_data, ]
	met_train <- met[training_data, ]
	bac_train <- bac[training_data, ]

	env_test <- env[-training_data, ] #samples not used in training data
	met_test <- met[intersect(rownames(met), rownames(env_test)), ]
	bac_test <- bac[intersect(rownames(bac), rownames(env_test)), ]

	#****************************************************************************************************************************#
	#use T-Tests to filter OTUs and metabolites to subset that are univariate associated with the outcome on the training dataset#
	#****************************************************************************************************************************#

	#read and preprocess metabolite data
	#met_raw <- read.csv("/home/madapoos/Spiromics_No_Severe_No_Healthy_UNC_Reg/Data/Esther_Metabolites_No_Severe_No_Healthy.csv", header = T)
	#rownames(met_raw) <- met_raw[, 1] #set compound names to rows
	#met_raw <- met_raw[, -1]
	#met_raw <- met_raw[intersect(rownames(met_raw), rownames(bac_train)), ] #subset to only training data

	#met_env_df <- data.frame(cbind(env_train[, "target_variable"], met_raw))

	#run individual t-tests between each metabolite and outcome
	#met_t_test <- list(NULL)
	#for (i in 2:ncol(met_env_df)) {
	#	if ((shapiro.test(met_env_df[, i])$p.value > 0.05) == TRUE) {
	#		met_t_test[i] <- t.test(met_env_df[, i], as.numeric(met_env_df[, 1]))$p.value #t-test on normally distributed metabolites
	#	} else {
	#		met_t_test[i] <- wilcox.test(met_env_df[, i] ~ met_env_df[, 1])$p.value #wilcoxon rank sum test on non-normally distributed metabolites
	#	}
	#}

	#met_t_test_adjusted <- p.adjust(unlist(met_t_test), method = "BH") #adjust all p-values using BH

	#metabolite_associations <- cbind(colnames(met_env_df)[2:ncol(met_env_df)], met_t_test_adjusted)
	#colnames(metabolite_associations) <- c("Metabolites", "P-Value")

	#associated_metabolites <- colnames(met_env_df)[which(met_t_test_adjusted < 0.05) + 1] #identify list of metabolites significantly associated with outcome of interest

	#met_train <- met_train[, associated_metabolites] #subset metabolite data to just those that are significantly associated with outcme
	#met_test <- met_test[, associated_metabolites]

	otu_raw <- read.csv("/home/madapoos/Spiromics_No_Severe_No_Healthy_UNC_Reg/Data/Esther_Microbiota_No_Severe_No_Healthy.csv", header = T, stringsAsFactors = FALSE)
	rownames(otu_raw) <- otu_raw[, 1] #set compound names to rows
	otu_raw <- otu_raw[, -1]
	otu_raw <- logratio.transfo(as.matrix(otu_raw), logratio = 'CLR', offset = 1e-10) #CLR transform OTU table
	otu_raw <- otu_raw[intersect(rownames(otu_raw), rownames(bac_train)), ] #subset to only otus data with matching microbiome data

	otu_env_df <- data.frame(cbind(env_train[, "target_variable"], otu_raw))

	#run individual t-tests between each otus and outcome
	otu_t_test <- list(NULL)
	for (i in 2:ncol(otu_env_df)) {
		if ((shapiro.test(otu_env_df[, i])$p.value > 0.05) == TRUE) {
			otu_t_test[i] <- t.test(otu_env_df[, i], as.numeric(otu_env_df[, 1]))$p.value #t-test on normally distributed otus
		} else {
			otu_t_test[i] <- wilcox.test(otu_env_df[, i] ~ otu_env_df[, 1])$p.value #wilcoxon rank sum test on non-normally distributed otus
		}
	}

	otus_associations <- cbind(colnames(otu_env_df)[2:ncol(otu_env_df)], unlist(otu_t_test))
	colnames(otus_associations) <- c("OTUs", "P-Value")

	associated_otus <- colnames(otu_env_df)[which(unlist(otu_t_test) < 0.2) + 1] #identify list of otus somewhat significantly associated with outcome of interest

	filter_num_optimal_bac_list[[j]] <- length(associated_otus)

	bac_train <- bac_train[, associated_otus, drop = FALSE] #subset otus data to just those that are significantly associated with outcme
	bac_test <- bac_test[, associated_otus, drop = FALSE]

	#**********************#
	#Final Training Dataset#
	#**********************#

	#merge datasets together for DIABLO
	overall <- list(env = env_train, met = met_train, bac = bac_train) #merge C18 met metabolites with bacteria and clinical data
	met_bac <- list(met = met_train, bac = bac_train)

	#response variable
	Y <- overall$env[, "target_variable"]

	#*******************************#
	#DIABLO-Based Feature Filtration#
	#*******************************#

	#prepare design matrix
	data <- met_bac
	design = matrix(0.1, ncol = length(data), nrow = length(data), 
	                dimnames = list(names(data), names(data))) #determines which blocks should be connected to maximize correlation or covariance between components (ex. mRNA and miRNA)
	diag(design) = 0

	#list of feature subsets to test
	ncomp <- 2
	test.keepX <- list(met = c(round(0.2*ncol(met_train)), 
							   round(0.4*ncol(met_train)),
							   round(0.6*ncol(met_train)), 
							   round(0.8*ncol(met_train)),
							   round(ncol(met_train))),
					   bac = c(round(0.2*ncol(bac_train)),
					   		   round(0.4*ncol(met_train)), 
							   round(0.6*ncol(bac_train)), 
							   round(0.8*ncol(bac_train)),
							   ncol(bac_train))) #number of features to try from each data set for optimal number of features
	tune.TCGA = try(tune.block.splsda(X = data, 
								  Y = Y, 
								  ncomp = ncomp, 
	                              test.keepX = test.keepX, 
	                              design = design,
	                              validation = 'Mfold', 
	                              folds = 5, #10-fold cross-validation
	                              nrepeat = 5,
	                              cpus = 2,
	                              dist = "centroids.dist"))
	if (class(tune.TCGA) == "try-error") {
		error_counter <- error_counter + 1 #count number of errors to break in case we have continuous errors, in which case the data is too sparse
		if (error_counter < error_counter_max) { #n/20
			next
		}
		else {
			broken <- TRUE #change to true as data is too sparse
			break #break from while loop if data too sparse to do gCCA
		}
	} else {
		list.keepX <- tune.TCGA$choice.keepX #print number of optimal features for each data set
	}
	print("b")
	#DIABLO model generation and OTU-metabolite data integration
	sgccda.res <- try(block.splsda(X = data, Y = Y, ncomp = ncomp, keepX = list.keepX, design = design))
	if (class(sgccda.res) == "try-error") {
		error_counter <- error_counter + 1 #count number of errors to break in case we have continuous errors, in which case the data is too sparse
		if (error_counter < error_counter_max) { #n/20
			next
		}
		else {
			broken <- TRUE #change to true as data is too sparse
			break #break from while loop if data too sparse to do gCCA
		}
	}
	#check DIABLO model performance in-sample and append parameters to list
	perf.diablo = try(perf(sgccda.res, validation = 'Mfold', M = 5, nrepeat = 5, dist = 'max.dist'))
	diablo_training_majority_vote_list[[j]] <- data.frame(perf.diablo$MajorityVote.error.rate)["Overall.BER", ] #training majority vote error rate
	diablo_training_weighted_vote_list[[j]] <- data.frame(perf.diablo$WeightedVote.error.rate)["Overall.BER", ] #training weighted vote error rate

	diablo_optimal_met_list <- rbind(data.frame(selectVar(sgccda.res, block = 'met', comp = 1)$met),
												data.frame(selectVar(sgccda.res, block = 'met', comp = 2)$met)) #optimal set of metabolites ranked in descending order by absolute value
	diablo_optimal_bac_list <- rbind(data.frame(selectVar(sgccda.res, block = 'bac', comp = 1)$bac),
												data.frame(selectVar(sgccda.res, block = 'bac', comp = 2)$bac)) #optimal set of OTUs ranked in descending order by absolute value

	#**************************************************#
	#Using DIABLO, select features for regression model#
	#**************************************************#

	#average loading scores of OTU and metabolite features that were selected in both components
	diablo_optimal_met_list <- aggregate(value.var ~ ., mean, data = diablo_optimal_met_list)
	rownames(diablo_optimal_met_list) <- diablo_optimal_met_list[, 1]
	diablo_optimal_bac_list <- aggregate(value.var ~ ., mean, data = diablo_optimal_bac_list)
	rownames(diablo_optimal_bac_list) <- diablo_optimal_bac_list[, 1]

	print("c")
	#order list of OTU and metabolite features identified by DIABLO in descending order by absolute value of loading
	diablo_optimal_met_list <- diablo_optimal_met_list %>% select(value.var) %>% arrange(desc(abs(value.var))) 
	diablo_optimal_bac_list <- diablo_optimal_bac_list %>% select(value.var) %>% arrange(desc(abs(value.var))) 
	diablo_num_optimal_met_list[[j]] <- nrow(diablo_optimal_met_list)
	diablo_num_optimal_bac_list[[j]] <- nrow(diablo_optimal_bac_list)

	#filter metabolite features to those in the top 10% of absolute value of loading, rounded up so we have
	#at least 1 feature from each dataset
	#diablo_optimal_met_list <- diablo_optimal_met_list[1:(ceiling(diablo_pct*nrow(diablo_optimal_met_list))), , drop = FALSE] 
	#diablo_optimal_bac_list <- diablo_optimal_bac_list[1:(ceiling(diablo_pct*nrow(diablo_optimal_bac_list))), , drop = FALSE]

	#filter training data to just those metabolites and OTUs selected by DIABLO
	met_train <- met_train[, which(colnames(met_train) %in% rownames(diablo_optimal_met_list)), drop = FALSE]
	bac_train <- bac_train[, which(colnames(bac_train) %in% rownames(diablo_optimal_bac_list)), drop = FALSE]

	#drop features with high correlations (>90%)
	#high_met_cor <- findCorrelation(cor(met_train), cutoff = 0.7)
	#if (length(high_met_cor) > 0) {
	#	met_train <- met_train[, -high_met_cor, drop = FALSE]
	#}
	#if (ncol(bac_train) > 1) {
	#	high_bac_cor <- findCorrelation(cor(bac_train), cutoff = 0.7) 
	#	if (length(high_bac_cor) > 0) {
	#		bac_train <- bac_train[, -high_bac_cor, drop = FALSE]
	#	}
	#}

	#*******************************#
	#Build Adjusted Regression Model#
	#*******************************#

#*******************************#
	#Build Adjusted Regression Model#
	#*******************************#

	print("a")	

	#data frame for demographics-adjusted model
	train_df <- cbind(env_train[, c("target_variable", "GENDER", "RACE_Nonwhite", "AGE_DERV_ctb", "CURRENT_SMOKER_ctb", "INHALEDSTEROIDS_ctb", "ANTI_USE_3MTH_BEF_BRON"), drop = FALSE], 
						bac_train, 
						met_train)

	#assume people without smoking status are in class 2 (NA)
	train_df$CURRENT_SMOKER_ctb[which(is.na(train_df$CURRENT_SMOKER_ctb))] <- 0
	train_df$INHALEDSTEROIDS_ctb[which(is.na(train_df$INHALEDSTEROIDS_ctb))] <- 0
	train_df$ANTI_USE_3MTH_BEF_BRON[which(is.na(train_df$ANTI_USE_3MTH_BEF_BRON))] <- 0

	#convert nominal demographics to separate variables and drop original nominal variables
	train_df$CURRENT_SMOKER_ctb <- factor(train_df$CURRENT_SMOKER_ctb)
	train_df$RACE_Nonwhite <- factor(train_df$RACE_Nonwhite)
	train_df$GENDER <- factor(train_df$GENDER)
	train_df$INHALEDSTEROIDS_ctb <- factor(train_df$INHALEDSTEROIDS_ctb)
	train_df$ANTI_USE_3MTH_BEF_BRON <- factor(train_df$ANTI_USE_3MTH_BEF_BRON)

	tc <- trainControl(method = "cv", number = 10, summaryFunction = twoClassSummary, classProbs = TRUE)
	train_df$target_variable <- make.names(as.factor(train_df$target_variable))

	#run penalized logistic regression on training data
	if (class(try(caret::train(target_variable ~ ., 
								data = train_df,
								method = "glmnet",
								family = "binomial", 
								trControl = tc,
								standardize = TRUE,
								#metric = "ROC", 
								tuneLength = 10))) == "try-error") { #try 100 combinations of alpha and lambda
		error_counter <- error_counter + 1
		if (error_counter < error_counter_max) { #n/20
			next
		}
		else {
			broken <- TRUE #change to true as data is too sparse
			print("Error in Penalized Logistic Regression Model")
			break #break from while loop if data too sparse
		}
	}
	mod_fit <- caret::train(target_variable ~ ., 
								data = train_df,
								method = "glmnet",
								family = "binomial", 
								trControl = tc,
								standardize = TRUE,
								#metric = "ROC", 
								tuneLength = 10)

	print("d")
	
	train_df$target_variable <- as.numeric(str_sub(train_df$target_variable, start = 2))

	final_model <- mod_fit$finalModel #extract model

	#elastic net-selected features
	coefficients <- coef(final_model, mod_fit$bestTune$lambda)
	coefficients <- as.data.frame(as.matrix(coefficients))
	coefficients <- coefficients[-which(coefficients[, 1] == 0), , drop = FALSE]
	coef_pred <- coefficients
	coef_pred$Feature <- rownames(coef_pred)
	colnames(coef_pred)[1] <- "Log_Odds_Ratio"
	coef_pred <- coef_pred %>% select(Feature, Log_Odds_Ratio)

	#extract optimal tuning parameters
	alpha <- mod_fit$bestTune$alpha
	lambda <- mod_fit$bestTune$lambda

	#subset training data to just those variables with no collinearity
	#train_df <- train_df[, which(colnames(train_df) %in% rownames(coefficients)), drop = FALSE]
	#bac_train <- bac_train[, which(colnames(bac_train) %in% colnames(train_df)), drop = FALSE]
	#met_train <- met_train[, which(colnames(met_train) %in% colnames(train_df)), drop = FALSE]


	#*************#
	#Model Testing#
	#*************#

	met_test <- met_test[, which(colnames(met_test) %in% colnames(met_train)), drop = FALSE]
	bac_test <- bac_test[, which(colnames(bac_test) %in% colnames(bac_train)), drop = FALSE]

	#data frame for model
	test_df <- cbind(env_test[, c("target_variable", "GENDER", "RACE_Nonwhite", "AGE_DERV_ctb", "CURRENT_SMOKER_ctb", "INHALEDSTEROIDS_ctb", "ANTI_USE_3MTH_BEF_BRON"), drop = FALSE], 
						bac_test, 
						met_test)

	#assume people without smoking status are in class 2 (NA)
	test_df$CURRENT_SMOKER_ctb[which(is.na(test_df$CURRENT_SMOKER_ctb))] <- 0
	test_df$INHALEDSTEROIDS_ctb[which(is.na(test_df$INHALEDSTEROIDS_ctb))] <- 0
	test_df$ANTI_USE_3MTH_BEF_BRON[which(is.na(test_df$ANTI_USE_3MTH_BEF_BRON))] <- 0

	#convert nominal demographics to separate variables and drop original nominal variables
	test_df$CURRENT_SMOKER_ctb <- factor(test_df$CURRENT_SMOKER_ctb)
	test_df$RACE_Nonwhite <- factor(test_df$RACE_Nonwhite)
	test_df$GENDER <- factor(test_df$GENDER)
	test_df$INHALEDSTEROIDS_ctb <- factor(test_df$INHALEDSTEROIDS_ctb)
	test_df$ANTI_USE_3MTH_BEF_BRON <- factor(test_df$ANTI_USE_3MTH_BEF_BRON)

	predictions <- predict(mod_fit, newdata = test_df)
	predicted[[j]] <- predictions

	#calculate balanced error rate training data
	confusion_matrix <- confusionMatrix(predictions, 
										as.factor(make.names(as.factor(test_df$target_variable))),
										positive = 'X1')
	bar <- confusion_matrix$byClass[[11]] #out-of-sample balanced accuracy

	print("f")

	#************************************************#
	#Univariate logistic regressions for each feature#
	#************************************************#

	#run univariate logistic regressions with each of the variables to get unadjusted odds ratios and p-values
	coef_pred$"Univariate Demographic-Adjusted Log OR" <- NULL
	coef_pred$"Univariate Demographic-Adjusted Log OR 95% CI Lower Bound" <- NULL
	coef_pred$"Univariate Demographic-Adjusted Log OR 95% CI Upper Bound" <- NULL
	coef_pred$"Univariate Demographic-Adjusted Log OR Wald Test P-Value" <- NULL
	coef_pred$"Univariate Demographic-Adjusted Model LRT P-Value" <- NULL
	coef_pred$"Univariate Demographic-Adjusted Model BAR" <- NULL

	for (m in rownames(coef_pred)) { #iterate through each predictor variable
		if (m == "(Intercept)") { #skip the intercept row
			next
		}
		demographics_list <- c("GENDER2", "RACE_Nonwhite1", "AGE_DERV_ctb", "CURRENT_SMOKER_ctb1", "INHALEDSTEROIDS_ctb1", "ANTI_USE_3MTH_BEF_BRON1") 
		#demographics-only model
		if (m %in% demographics_list) { #check if the variable of interest is one of the demographics (can't adjust for itself)
			df_subset <- env_train[, c("target_variable", "GENDER", "RACE_Nonwhite", "AGE_DERV_ctb", "CURRENT_SMOKER_ctb", "INHALEDSTEROIDS_ctb", "ANTI_USE_3MTH_BEF_BRON") , drop = FALSE]
			df_subset$CURRENT_SMOKER_ctb[which(is.na(df_subset$CURRENT_SMOKER_ctb))] <- 0
			df_subset$INHALEDSTEROIDS_ctb[which(is.na(df_subset$INHALEDSTEROIDS_ctb))] <- 0
			df_subset$ANTI_USE_3MTH_BEF_BRON[which(is.na(df_subset$ANTI_USE_3MTH_BEF_BRON))] <- 0
			#convert nominal demographics to separate variables and drop original nominal variables
			df_subset$CURRENT_SMOKER_ctb <- factor(df_subset$CURRENT_SMOKER_ctb)
			df_subset$RACE_Nonwhite <- factor(df_subset$RACE_Nonwhite)
			df_subset$GENDER <- factor(df_subset$GENDER)
			df_subset$INHALEDSTEROIDS_ctb <- factor(df_subset$INHALEDSTEROIDS_ctb)
			df_subset$ANTI_USE_3MTH_BEF_BRON <- factor(df_subset$ANTI_USE_3MTH_BEF_BRON)
			#model each variable individually, adjusted for demographics
			tc <- trainControl(method = "cv", number = 10)
			univariate_model <- train(as.factor(target_variable) ~ ., 
										data = df_subset,
										method = "glm",
										family = "binomial",
										trControl = tc,
										maxit = 100)

			univariate_final_model <- univariate_model$finalModel
			#likelihood ratio test for model significance
			univariate_null_deviance <- summary(univariate_final_model)$null.deviance
			univariate_residual_deviance <- summary(univariate_final_model)$deviance
			univariate_null_df <- summary(univariate_final_model)$df.null
			univariate_residual_df <- summary(univariate_final_model)$df.residual
			univariate_lrt_pval <- 1 - pchisq(univariate_null_deviance - univariate_residual_deviance, univariate_null_df - univariate_residual_df)
		
			#save coefficients to final table
			univariate_coefficient <- coef(univariate_final_model)[2]
			univariate_coef_df <- data.frame(univariate_coefficient)

			#get confidence intervals
			critval <- 1.96 #t critical value
			univariate_standard_errors <- summary(univariate_final_model)$coefficients[2, 2] #std errors for coefficients
			univariate_upper_bound <- coef(univariate_final_model)[2] + (critval * univariate_standard_errors)
			univariate_lower_bound <- coef(univariate_final_model)[2] - (critval * univariate_standard_errors)
			univariate_coef_confint <- cbind(data.frame(univariate_lower_bound), data.frame(univariate_upper_bound))

			#predictions by model in-sample
			univariate_predictions <- predict(univariate_model, newdata = df_subset)

			#calculate balanced error rate on all data
			univariate_confusion_matrix <- confusionMatrix(univariate_predictions, 
															as.factor(df_subset$target_variable), 
															positive = '1')
			univariate_bar <- univariate_confusion_matrix$byClass[[11]]

			#using Wald test from model output, determine which predictors are significant (Type 3 SS; partial)
			univariate_pred_p <- summary(univariate_final_model)$coefficients[2, 4]
			univariate_pred_p <- data.frame(univariate_pred_p)

			#append univariate ORs, confidence intervals, p-values, LRT p-values, and balanced accuracies to coefficient data frame
			coef_pred[m, "Univariate Demographic-Adjusted Log OR"] <- univariate_coefficient[1]
			coef_pred[m, "Univariate Demographic-Adjusted Log OR 95% CI Lower Bound"] <- univariate_lower_bound[1]
			coef_pred[m, "Univariate Demographic-Adjusted Log OR 95% CI Upper Bound"] <- univariate_upper_bound[1]
			coef_pred[m, "Univariate Demographic-Adjusted Log OR Wald Test P-Value"] <- univariate_pred_p[1]
			coef_pred[m, "Univariate Demographic-Adjusted Model LRT P-Value"] <- univariate_lrt_pval
			coef_pred[m, "Univariate Demographic-Adjusted Model BAR"] <- univariate_bar
			next #go to next iteration
		#for all non-intercept, non-demographic predictors
		} else {
			df_subset <- cbind(env_train[, c("target_variable", "GENDER", "RACE_Nonwhite", "AGE_DERV_ctb", "CURRENT_SMOKER_ctb", "INHALEDSTEROIDS_ctb", "ANTI_USE_3MTH_BEF_BRON") , drop = FALSE], 
									train_df[, m, drop = FALSE]) #create 2 column dataset

			df_subset$CURRENT_SMOKER_ctb[which(is.na(df_subset$CURRENT_SMOKER_ctb))] <- 0
			df_subset$INHALEDSTEROIDS_ctb[which(is.na(df_subset$INHALEDSTEROIDS_ctb))] <- 0
			df_subset$ANTI_USE_3MTH_BEF_BRON[which(is.na(df_subset$ANTI_USE_3MTH_BEF_BRON))] <- 0

			#convert nominal demographics to separate variables and drop original nominal variables
			df_subset$CURRENT_SMOKER_ctb <- factor(df_subset$CURRENT_SMOKER_ctb)
			df_subset$RACE_Nonwhite <- factor(df_subset$RACE_Nonwhite)
			df_subset$GENDER <- factor(df_subset$GENDER)
			df_subset$INHALEDSTEROIDS_ctb <- factor(df_subset$INHALEDSTEROIDS_ctb)
			df_subset$ANTI_USE_3MTH_BEF_BRON <- factor(df_subset$ANTI_USE_3MTH_BEF_BRON)

			#model each variable individually, adjusted for demographics
			tc <- trainControl(method = "cv", number = 10)
			univariate_model <- train(as.factor(target_variable) ~ ., 
										data = df_subset,
										method = "glm",
										family = "binomial",
										trControl = tc,
										maxit = 100)

			univariate_final_model <- univariate_model$finalModel
			#likelihood ratio test for model significance
			univariate_null_deviance <- summary(univariate_final_model)$null.deviance
			univariate_residual_deviance <- summary(univariate_final_model)$deviance
			univariate_null_df <- summary(univariate_final_model)$df.null
			univariate_residual_df <- summary(univariate_final_model)$df.residual
			univariate_lrt_pval <- 1 - pchisq(univariate_null_deviance - univariate_residual_deviance, univariate_null_df - univariate_residual_df)
		
			#save coefficients to final table
			univariate_coefficient <- coef(univariate_final_model)[2]
			univariate_coef_df <- data.frame(univariate_coefficient)

			#get confidence intervals
			critval <- 1.96 #t critical value
			univariate_standard_errors <- summary(univariate_final_model)$coefficients[2, 2] #std errors for coefficients
			univariate_upper_bound <- coef(univariate_final_model)[2] + (critval * univariate_standard_errors)
			univariate_lower_bound <- coef(univariate_final_model)[2] - (critval * univariate_standard_errors)
			univariate_coef_confint <- cbind(data.frame(univariate_lower_bound), data.frame(univariate_upper_bound))

			#predictions by model in-sample
			univariate_predictions <- predict(univariate_model, newdata = df_subset)

			#calculate balanced error rate on all data
			univariate_confusion_matrix <- confusionMatrix(univariate_predictions, 
															as.factor(df_subset$target_variable), 
															positive = '1')
			univariate_bar <- univariate_confusion_matrix$byClass[[11]]

			#using Wald test from model output, determine which predictors are significant (Type 3 SS; partial)
			univariate_pred_p <- summary(univariate_final_model)$coefficients[2, 4]
			univariate_pred_p <- data.frame(univariate_pred_p)

			#append univariate ORs, confidence intervals, p-values, LRT p-values, and balanced accuracies to coefficient data frame
			coef_pred[m, "Univariate Demographic-Adjusted Log OR"] <- univariate_coefficient[1]
			coef_pred[m, "Univariate Demographic-Adjusted Log OR 95% CI Lower Bound"] <- univariate_lower_bound[1]
			coef_pred[m, "Univariate Demographic-Adjusted Log OR 95% CI Upper Bound"] <- univariate_upper_bound[1]
			coef_pred[m, "Univariate Demographic-Adjusted Log OR Wald Test P-Value"] <- univariate_pred_p[1]
			coef_pred[m, "Univariate Demographic-Adjusted Model LRT P-Value"] <- univariate_lrt_pval
			coef_pred[m, "Univariate Demographic-Adjusted Model BAR"] <- univariate_bar
		}
	}

	#identify number of OTU and metabolite features that were selected by the penalized regression
	train_df <- train_df[, which(colnames(train_df) %in% rownames(coef_pred)), drop = FALSE]
	bac_train <- bac_train[, which(colnames(bac_train) %in% colnames(train_df)), drop = FALSE]
	met_train <- met_train[, which(colnames(met_train) %in% colnames(train_df)), drop = FALSE]

	print("g")

	#**************************************************************#
	#Append adjusted model parameters to list for later aggregation#
	#**************************************************************#

	num_optimal_met_list[[j]] <- ncol(met_train) #number of features from metabolite dataset
	num_optimal_bac_list[[j]] <- ncol(bac_train) #number of features from bacterial dataset
	optimal_met_list[[j]] <- colnames(met_train)
	optimal_bac_list[[j]] <- colnames(bac_train)
	overall_bar_list[[j]] <- bar #overall BAR - used to determine whether this model is the best
	confusion_matrix_list <- rbind(confusion_matrix_list, confusion_matrix) #append to data frame containing previous confusion matrices
	alpha_list[[j]] <- alpha
	lambda_list[[j]] <- lambda

	#check if this model is better than the previous best model based on BAR
	if (j == 1) { #if it's the first model it's automatically the best model
		best_model <- final_model
		best_num_filter_optimal_bac <- length(associated_otus)
		best_diablo_training_majority_vote <- data.frame(perf.diablo$MajorityVote.error.rate)["Overall.BER", 1]
		best_diablo_training_weighted_vote <- data.frame(perf.diablo$WeightedVote.error.rate)["Overall.BER", 1]	
		best_diablo_num_optimal_met <- nrow(diablo_optimal_met_list)
		best_diablo_num_optimal_bac <- nrow(diablo_optimal_bac_list)
		best_num_optimal_met <- ncol(met_train)
		best_num_optimal_bac <- ncol(bac_train)
		best_optimal_met_list <- colnames(met_train)
		best_optimal_bac_list <- colnames(bac_train)
		best_confusion_matrix <- confusion_matrix
		best_alpha <- alpha
		best_lambda <- lambda
		best_overall_bar <- bar
		saveRDS(best_model, file = paste("/home/madapoos/Spiromics_No_Severe_No_Healthy_UNC_Reg/Output/", response_variable, "/", response_variable, "_Best_Model.rds", sep = "")) #save model as R object to read in later
	} else if (bar >= best_overall_bar) {
		best_model <- final_model
		best_num_filter_optimal_bac <- length(associated_otus)
		best_diablo_training_majority_vote <- data.frame(perf.diablo$MajorityVote.error.rate)["Overall.BER", 1]
		best_diablo_testing_majority_vote <- data.frame(perf.diablo$WeightedVote.error.rate)["Overall.BER", 1]	
		best_diablo_num_optimal_met <- nrow(diablo_optimal_met_list)
		best_diablo_num_optimal_bac <- nrow(diablo_optimal_bac_list)
		best_num_optimal_met <- ncol(met_train)
		best_num_optimal_bac <- ncol(bac_train)
		best_optimal_met_list <- colnames(met_train)
		best_optimal_bac_list <- colnames(bac_train)
		best_confusion_matrix <- confusion_matrix
		best_alpha <- alpha
		best_lambda <- lambda
		best_overall_bar <- bar
		saveRDS(best_model, file = paste("/home/madapoos/Spiromics_No_Severe_No_Healthy_UNC_Reg/Output/", response_variable, "/", response_variable, "_Best_Model.rds", sep = "")) #save model as R object to read in later
	}
	coef_pred_list[[j]] <- data.frame(coef_pred) #coefficients, CIs, p-values, and univariate metrics

	j <- j + 1 #increment counter

}

#********************#
#Aggregate Parameters#
#********************#
mean_num_optimal_met <- mean(unlist(num_optimal_met_list)) #mean number of metabolite features
mean_num_optimal_bac <- mean(unlist(num_optimal_bac_list)) #mean number of OTU features
mean_diablo_training_majority_vote <- mean(unlist(diablo_training_majority_vote_list))
mean_diablo_training_weighted_vote <- mean(unlist(diablo_training_weighted_vote_list))
mean_num_diablo_optimal_met <- mean(unlist(diablo_num_optimal_met_list))
mean_num_diablo_optimal_bac <- mean(unlist(diablo_num_optimal_bac_list))
mean_num_filter_optimal_bac <- mean(unlist(filter_num_optimal_bac_list))
mean_alpha <- mean(unlist(alpha_list))
mean_lambda <- mean(unlist(lambda_list))
mean_overall_bar <- mean(unlist(overall_bar_list))
mean_overall_bar_sig_test_p_value <- t.test(unlist(overall_bar_list), alternative = "greater", mu = 0.5)$p.value #test if model accuracy significantly greater than 0.5, on average
best_overall_bar_sig_test_p_value <- NA #placeholder since you can't test if the best model is significant

model_parameters_names <- c("Number of OTUs Inputted Into DIABLO",
							"Number of DIABLO-Identified Optimal Metabolites", 
							"Number of DIABLO-Identified Optimal OTUs",
							"DIABLO Training Majority Vote BAR",
							"DIABLO Training Weighted Vote BAR",
							"Number of Metabolites in Final Penalized Regression Model",
							"Number of OTUs in Final Penalized Regression Model",
							"Mean Alpha",
							"Mean Lambda",
							"Testing Overall BAR",
							"Testing Overall BAR Significance Test P-Value")
mean_model_parameters <- c(mean_num_filter_optimal_bac,
						   mean_num_diablo_optimal_met,
						   mean_num_diablo_optimal_bac,
						   mean_diablo_training_majority_vote,
						   mean_diablo_training_weighted_vote,
						   mean_num_optimal_met,
						   mean_num_optimal_bac,
						   mean_alpha,
						   mean_lambda,
						   mean_overall_bar,
						   mean_overall_bar_sig_test_p_value)
best_model_parameters <- c(best_num_filter_optimal_bac,
						   best_diablo_num_optimal_met,
						   best_diablo_num_optimal_bac,
						   best_diablo_training_majority_vote,
						   best_diablo_training_weighted_vote,
						   best_num_optimal_met,
						   best_num_optimal_bac,
						   best_alpha,
						   best_lambda,
						   best_overall_bar,
						   best_overall_bar_sig_test_p_value)
model_parameters <- cbind(model_parameters_names, 
						  mean_model_parameters, 
						  best_model_parameters)
colnames(model_parameters) <- c("Model Parameter", 
								"Mean Value Over 100 Resampling Iterations", 
								"Value from Best Model (Model with Lowest Testing Overall BER)")
write.csv(model_parameters, file = paste("/home/madapoos/Spiromics_No_Severe_No_Healthy_UNC_Reg/Output/", response_variable, "/DIABLO_Model_Parameters.csv", sep = ""), row.names = F)


#list of features and the number of models in which each feature was selected (one measure of importance)
#optimal_met_list <- unlist(optimal_met_list) %>% data.frame()
#colnames(optimal_met_list) <- "Metabolite_Feature"
#optimal_met_list <- optimal_met_list %>% 
#						group_by(Metabolite_Feature) %>% 
#						summarize(n = n()) %>% 
#						arrange(desc(n)) #aggregate mean loading for selected variable and how many times they appear in 100 models 
#write.csv(optimal_met_list, file = paste("/home/madapoos/Spiromics_No_Severe_No_Healthy_UNC_Reg/Output/", target_variable, "/DIABLO_Model_Metabolite_Features.csv", sep = ""), row.names = F)
#optimal_bac_list <- unlist(optimal_bac_list) %>% data.frame()
#colnames(optimal_bac_list) <- "OTU_Feature"
#optimal_bac_list <- optimal_bac_list %>% 
#						group_by(OTU_Feature) %>% 
#						summarize(n = n()) #aggregate mean loading for selected variable and how many times they appear in 100 models 
#write.csv(optimal_bac_list, paste("/home/madapoos/Spiromics_No_Severe_No_Healthy_UNC_Reg/Output/", target_variable, "/DIABLO_Model_OTU_Features.csv", sep = ""), row.names = F)

#data frame of coefficients, ORs, etc. across all n models
coef_pred_df <- coef_pred_list[[1]]
for (i in 2:length(coef_pred_list)) { #create single data frame with all n coef_preds
	coef_pred_df <- dplyr::bind_rows(coef_pred_df, data.frame(coef_pred_list[[i]]))
}

write.csv(coef_pred_df, file = paste("/home/madapoos/Spiromics_No_Severe_No_Healthy_UNC_Reg/Output/", response_variable, "/DIABLO_Model_Features_ORs_P_Values_Raw.csv", sep = ""), row.names = F)


