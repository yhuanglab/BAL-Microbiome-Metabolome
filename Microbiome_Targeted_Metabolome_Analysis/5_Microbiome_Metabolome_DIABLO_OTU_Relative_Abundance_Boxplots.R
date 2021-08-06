#select number of Metabolites to generate for based on inflection point of scree plot
suppressMessages(library(dplyr))
suppressMessages(library(BiodiversityR)) #for generating rank abundance plots
suppressMessages(library(ggplot2))
suppressMessages(library(gridExtra)) #for arranging grids
suppressMessages(library(stringr))
suppressMessages(library(mixOmics))
options(warn = -1)

#***************************************#
#USER UPDATES THIS FOR EACH NEW VARIABLE#
#***************************************#

target_variable <- str_remove(getwd(), "/Users/Siddharth/Desktop/Spiromics_BAL_Microbiome_Metabolome_V2/Spiromics_UNC_Stratified/Spiromics_No_Severe_No_Healthy/Output/DIABLO/") #remove path and just get target variable name

number_of_otus_to_plot <- 10 #based on inflection point from stripchart if there is one; otherwise plot top 10 OTUs
response_variable <- target_variable #variable whose levels against which we plot OTU relative abundance
response_variable_levels <- c("0", "1+") #manually set levels of the response variable to show on the boxplot
response_variable_levels_label <- ""

#*******************#
#DO NOT MODIFY BELOW#
#*******************#
env <- read.csv("~/Desktop/Spiromics_BAL_Microbiome_Metabolome_V2/Spiromics_UNC_Stratified/Spiromics_No_Severe_No_Healthy/Data/Esther_Clinical_04_21_20_No_Severe_No_Healthy.csv")
bac <- read.csv('~/Desktop/Spiromics_BAL_Microbiome_Metabolome_V2/Spiromics_UNC_Stratified/Spiromics_No_Severe_No_Healthy/Data/Esther_Microbiota_No_Severe_No_Healthy.csv') #OTU relative abundance data
#ids <- read.csv("~/Desktop/Spiromics_BAL_Microbiome_Metabolome_V2/Spiromics_UNC_Stratified/Spiromics_No_Severe_No_Healthy/Data/Sample_IDs_C18_Pos_BALF.csv", header=T)

df <- read.csv(paste("~/Desktop/Spiromics_BAL_Microbiome_Metabolome_V2/Spiromics_UNC_Stratified/Spiromics_No_Severe_No_Healthy/Output/DIABLO/", target_variable, "/DIABLO_Model_Features_ORs_P_Values.csv", sep = ""), stringsAsFactors = FALSE) #table with Metabolites, loadings, and number of models in which they were selected
df <- df %>% arrange(desc(Number_of_Final_Models_Where_Feature_Selected))
formatted_name <- target_variable #title of chart

demographic_variables <- c("(Intercept)", "GENDER2", "RACE_Nonwhite1", "AGE_DERV_ctb", "CURRENT_SMOKER_ctb1", "INHALEDSTEROIDS_ctb1") 

#subset to just the metabolite features
otus_list <- df$Feature[which(str_sub(df$Feature, start = 1, end = 3) == "Otu"), drop = FALSE]
if (length(which(otus_list %in% demographic_variables)) > 0) {
	otus_list <- otus_list[-which(otus_list %in% demographic_variables)]
}
df <- df[which(df$Feature %in% otus_list), , drop = FALSE]
bac_list <- df[1:number_of_otus_to_plot, ]

#prepare OTU data
rownames(bac) <- bac[, 1] #set BALF number to rownames
bac <- bac[, -1]
bac <- logratio.transfo(as.matrix(bac), logratio = 'CLR', offset = 1e-10) 
bac <- data.frame(bac[1:nrow(bac), ])

#prepare clinical data
#env <- left_join(ids, env, by = 'Reisdorph..') #add in lab ids to match to bac2terial data
rownames(env) <- env[, 1] #set names to reisdorph number
env <- env[, -1] 
env <- env[intersect(rownames(env), rownames(bac)), ] #subset to only clinical data with matching microbiome data

#######################################################################################
#Analysis of each group for significance of differences in all OTU relative abundances#
#######################################################################################

bac_all <- cbind(bac, env[, response_variable]) #bind group levels to bac2terial data
bac_all <- na.omit(bac_all) #remove subjects with missing values for the response variable of interest

colnames(bac_all)[ncol(bac_all)] <- response_variable #label clinical column appropriately

group_1 <- bac_all[which(bac_all[, target_variable] == 0), ] #subset to members of group 1
group_2 <- bac_all[which(bac_all[, target_variable] == 1), ] #subset to members of group 2

p_value_list <- NULL
p_value_list_corrected <- NULL #t-tests comparing mean relative abundance of each bac2terium
for (i in 1:(ncol(bac_all)-1)) { #run on all OTUs except the last column, which is k-means cluster
	t_test <- t.test(group_1[, i], group_2[, i]) #two sample t test comparing mean relative abundance by cluster
	p_value_list <- c(p_value_list, t_test$p.value)
	p_value_list_corrected <- c(p_value_list_corrected, t_test$p.value)
}

p_value_list_corrected <- p.adjust(p_value_list_corrected, method = "BH") #adjust p-values for FDR across all OTUs
p_value_list <- cbind(p_value_list, p_value_list_corrected) #add corrected p-values to original p-values
names(p_value_list) <- colnames(group_1[1:(ncol(group_1)-1)])

group_1_means <- colMeans(group_1[1:(ncol(group_1) - 1)]) %>% data.frame()
group_2_means <- colMeans(group_2[1:(ncol(group_2) - 1)]) %>% data.frame()

df_raw <- cbind(group_1_means, group_2_means, p_value_list) #append p-values to groups
colnames(df_raw) <- c(paste("CLR-Transformed Relative Abundance in", response_variable_levels[[1]]),
				  paste("CLR-Transformed Relative Abundance in", response_variable_levels[[2]]),
				  "Uncorrected_P-Value", 
				  "BH-Adjusted_T-Test_P-Value")
write.csv(df_raw, file = paste("~/Desktop/Spiromics_BAL_Microbiome_Metabolome_V2/Spiromics_UNC_Stratified/Spiromics_No_Severe_No_Healthy/Output/DIABLO/", target_variable, "/OTU_Significance_Tests_by_Levels_of_", target_variable, ".csv", sep = ""))

#**************#
#SUBSET OF OTUS#
#**************#

bac2 <- bac[, which(colnames(bac) %in% bac_list$"Feature")] #subset to those OTUs selected by the most models

bac2 <- cbind(bac2, env[, response_variable]) #bind group levels to bac2terial data
bac2 <- na.omit(bac2) #remove subjects with missing values for the response variable of interest

colnames(bac2)[number_of_otus_to_plot + 1] <- response_variable #label clinical column appropriately

df <- df_raw[which(rownames(df_raw) %in% bac_list$"Feature"), ] #subset to those OTUs selected by the most models

bac2[, target_variable] <- as.factor(bac2[, target_variable])

myplot <- function(i) { #plots individual boxplots using ggplot
	dz <- data.frame(bac2[ target_variable], bac2[, i]) #create new data frame with just variable and otu
	colnames(dz) <- c(target_variable, "otu_relative_abundance")
	p_val_corrected <- df[i, "BH-Adjusted_T-Test_P-Value"] #gather adjusted p-value from t-test
	p_val <- df[i, "Uncorrected_P-Value"] #gather unadjusted p-value from t-test
	asterisk_selection <- function(n) {
		label <- NULL
		if(n > 0.05) {
			label <- "ns" #relative abundance difference not significant
		}
		else if (n <= 0.05 && n > 0.01) {
			label <- "*"
		}
		else if (n <= 0.01 && n > 0.001) {
			label <- "**"
		}
		else if (n <= 0.001 && n > 0.0001) {
			label <- "***"
		}
		else if (n <= 0.0001) {
			label <- "****"
		}
		return (label)
	}
	asterisk_corrected <- asterisk_selection(p_val_corrected) #determine which significance indicator to be selected
	asterisk <- asterisk_selection(p_val)
	#Generate boxplot
	pdf(paste("~/Desktop/Spiromics_BAL_Microbiome_Metabolome_V2/Spiromics_UNC_Stratified/Spiromics_No_Severe_No_Healthy/Output/DIABLO/", target_variable, "/OTU_Boxplots/", i, "_Relative_Abundance_by", target_variable, ".pdf", sep=""))
	p <- ggplot(data = dz) +
		geom_boxplot(mapping = aes(x = dz[, target_variable], y = otu_relative_abundance, fill = dz[, target_variable])) +
		annotate("text", -Inf, Inf, label = paste("BH-Corrected P-Value = ", signif(p_val_corrected, digits = 3), " (", asterisk_corrected, ")", sep = ''), hjust = 0, vjust = 3, size = 5, color = 'red') + 
		annotate("text", -Inf, Inf, label = paste("Uncorrected P-Value = ", signif(p_val, digits = 3), " (", asterisk, ")", sep = ''), hjust = 0, vjust = 1.5, size = 5, color = 'red') + 
		labs(x = "Level", y = "CLR-Transformed Relative Abundance") +
		ggtitle(str_wrap(paste(i, "CLR-Transformed Relative Abundance by Level of", response_variable, response_variable_levels_label, sep=" "), width = 40)) + #plot title
		scale_fill_manual(values = c("blue", "yellow"), guide = FALSE) +
		scale_x_discrete(labels = response_variable_levels) +
		theme(plot.title = element_text(hjust = 0.5, size = 15), #modify plot title
			  axis.title = element_text(face="bold", size = 15), #set axes to bold
			  axis.text.x = element_text(size = 12),
			  axis.text.y = element_text(size = 12),
			  legend.position = "none") #remove legend
	plot(p) #plot individual graph
	return(p)
}
p <- lapply(colnames(bac2)[1:number_of_otus_to_plot], myplot) #call function on all otus of interest
pdf(paste("~/Desktop/Spiromics_BAL_Microbiome_Metabolome_V2/Spiromics_UNC_Stratified/Spiromics_No_Severe_No_Healthy/Output/DIABLO/", target_variable, "/OTU_Boxplots/All_OTUs_Relative_Abundance_.pdf", sep=""), width = 20, height = 20)
do.call(grid.arrange, c(p, ncol = 3)) #arrange plots in a grid
dev.off()
