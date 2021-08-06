#derived from http://www.sthda.com/english/wiki/ggplot2-quick-correlation-matrix-heatmap-r-software-and-data-visualization
library(reshape2) #for correlation matrix plotting
library(ggplot2)
library(stringr)
library(tidyverse)
library(Hmisc)
library(corrplot)
library(RColorBrewer)
suppressMessages(library(circlize))
library(mixOmics)
library(gplots)
library(ComplexHeatmap)

#significance threshold for what feature correlations are shown
sig_threshold <- 0.05

#just plot features predictive of BAL neutrophils?
only_bal_neutrophils <- FALSE
only_ct <- TRUE

#relative_abundance tables
otu_relative_abundance <- read.csv('~/Desktop/Spiromics_BAL_Microbiome_Metabolome_V2/Spiromics_Colorado_Stratified/Spiromics_No_Severe_No_Healthy/Data/BALF_Decontam_Only_OTU_Relative_Abundances_Genus_04_06_19_No_Severe_No_Healthy.csv', stringsAsFactors = FALSE)
metabolite_relative_abundance <- read.csv('~/Desktop/Spiromics_BAL_Microbiome_Metabolome_V2/Spiromics_Colorado_Stratified/Spiromics_No_Severe_No_Healthy/Data/BALF_C18_Norm_Log2_Transposed_No_Severe_No_Healthy.csv', stringsAsFactors = FALSE)

#CLR transform OTU table before doing Pearson to make it follow a normal distribution
rownames(otu_relative_abundance) <- otu_relative_abundance[, 1]
otu_relative_abundance <- otu_relative_abundance[, -1]
otu_relative_abundance <- logratio.transfo(as.matrix(otu_relative_abundance), logratio = 'CLR', offset = 1e-10)
otu_relative_abundance <- data.frame(otu_relative_abundance[1:nrow(otu_relative_abundance), ])
otu_relative_abundance$X <- rownames(otu_relative_abundance)

#features for network
features <- read.csv('~/Desktop/Spiromics_BAL_Microbiome_Metabolome_V2/Spiromics_Colorado_Stratified/Spiromics_No_Severe_No_Healthy/Output/Heatmap/DIABLO_Num_Models_Heatmap/Top_Features_by_Num_Sig_Models_No_Unknowns_Dir_Effect_Added.csv', stringsAsFactors = FALSE)

#if we are just plotting features predictive of BAL neutrophils or the CT data, subset features to just those ones
if (only_bal_neutrophils == TRUE) {
    features <- features[which(features$Phenotype == "BAL_Neutrophils_Aboveequal_Below_Median"), ]
    sig_threshold <- 0.05 #increase significance threshold to pull more features out
} else if (only_ct == TRUE) {
  features <- features[which(features$Phenotype %in% c("V1_PERCENT_FSAD_TOTAL_Aboveequal_Below_Median", "PI10_WHOLE_TREE_LEQ20_V2_Aboveequal_Below_Median")), ]
  sig_threshold <- 0.05 #increase significance threshold to pull more features out
} else {
    features <- features[which(features$Phenotype != "BAL_Neutrophils_Aboveequal_Below_Median"), ]
    features <- features[which(features$Phenotype != "V1_PERCENT_FSAD_TOTAL_Aboveequal_Below_Median"), ]
    features <- features[which(features$Phenotype != "PI10_WHOLE_TREE_LEQ20_V2_Aboveequal_Below_Median"), ]
}

#features_list <- read.csv('~/Desktop/Spiromics_BAL_Microbiome_Metabolome_V2/Spiromics_Colorado_Stratified/Spiromics_No_Severe_No_Healthy/Output/Heatmap/DIABLO_Correlation_Heatmap/BALF_Heatmap_Group_Selected_Features_OTU_vs_Met_No_Unknowns_Significant.csv', stringsAsFactors = FALSE)
#features <- colnames(data.frame(features_list))
#features <- c(features, features_list$X)
#features <- features[-1]
#features <- data.frame(features)
#colnames(features)[1] <- "Feature"

#subset relative_abundance tables to features
otu_features <- features$Feature[which(str_sub(features$Feature, 1, 3) == "Otu")] #OTUs that are most predictive and in the majority of models
otu_relative_abundance_subset <- otu_relative_abundance[, which(colnames(otu_relative_abundance) %in% otu_features), drop = FALSE] #count data for predictive OTUs
if (length(which(str_sub(colnames(otu_relative_abundance_subset), -1) %in% 1:10)) > 0) {
    otu_relative_abundance_subset_nodup <- otu_relative_abundance_subset[, -which(str_sub(colnames(otu_relative_abundance_subset), -1) %in% 1:10)] #remove duplicate columns (same OTU predictive of multiple measures)
} else {
    otu_relative_abundance_subset_nodup <- otu_relative_abundance_subset
}

duplicate_keys <- c(".1", ".2", ".3", ".4", ".5", ".6", ".7", ".8", ".9", ".0") #R automatically adds this to duplicated names
metabolite_features <- features$Feature[which(str_sub(features$Feature, 1, 3) != "Otu")] #metabolites that are most predictive and in the majority of models
metabolite_relative_abundance_subset <- metabolite_relative_abundance[, which(colnames(metabolite_relative_abundance) %in% metabolite_features), drop = FALSE] #count data for predictive metabolites
if (length(which(str_sub(colnames(metabolite_relative_abundance_subset), -2) %in% duplicate_keys)) > 0) {
  metabolite_relative_abundance_subset_nodup <- metabolite_relative_abundance_subset[, -which(str_sub(colnames(metabolite_relative_abundance_subset), -2) %in% duplicate_keys)] #remove duplicate columns (same metabolite predictive of multiple measures)  
} else {
  metabolite_relative_abundance_subset_nodup <- metabolite_relative_abundance_subset
}
metabolite_relative_abundance_subset_nodup_ct <- metabolite_relative_abundance_subset_nodup

#label metabolites with correct name
met_raw <- read.csv("~/Desktop/Spiromics_BAL_Microbiome_Metabolome_V2/Spiromics_Colorado_Stratified/Spiromics_No_Severe_No_Healthy/Data/BALF_C18_Norm_Log2_Transposed_No_Severe_No_Healthy.csv")
rownames(met_raw) <- met_raw[, 1]
met_raw <- met_raw[, -1]
met_indices <- which(colnames(met_raw) %in% colnames(metabolite_relative_abundance_subset_nodup_ct))
met_raw <- met_raw[, met_indices]
met <- read.csv("~/Desktop/Spiromics_BAL_Microbiome_Metabolome_V2/Spiromics_Colorado_Stratified/Spiromics_No_Severe_No_Healthy/Data/BALF_C18_Norm_Log2_Transposed_No_Severe_No_Healthy.csv", check.names = FALSE)
rownames(met) <- met[, 1]
met <- met[, -1]
met <- met[, met_indices]
mapping <- cbind(colnames(met_raw), colnames(met))
for (i in 1:ncol(metabolite_relative_abundance_subset_nodup_ct)) {
  column_name <- colnames(metabolite_relative_abundance_subset_nodup_ct)[i]
  mapped_index <- which(colnames(met_raw) == column_name)
  mapped_name <- colnames(met)[mapped_index]
  colnames(metabolite_relative_abundance_subset_nodup_ct)[i] <- mapped_name
}

#sort metabolites by class
metabolite_classes <- read.csv("~/Desktop/Spiromics_BAL_Microbiome_Metabolome_V2/04_04_21_DIABLO_Significant_Metabolites_List_Updated_KAD_NR.csv", stringsAsFactors = FALSE)
col_order_df <- data.frame(colnames(metabolite_relative_abundance_subset_nodup_ct))
colnames(col_order_df) <- "Original_Compound_Name"
col_order_df <- left_join(col_order_df, metabolite_classes)
col_order_df <- col_order_df[, c("Original_Compound_Name", "Super_Class", "Heatmap_Color_Superclass")]
col_order_df <- col_order_df[order(col_order_df$Super_Class), ]
col_order_df[which(is.na(col_order_df$Super_Class)), "Super_Class"] <- "Other"
metabolite_relative_abundance_subset_nodup_ct <- metabolite_relative_abundance_subset_nodup_ct[, col_order_df$"Original_Compound_Name"]
metabolite_relative_abundance_subset_nodup_ct <- rbind(metabolite_relative_abundance_subset_nodup_ct, col_order_df$Super_Class, col_order_df$Heatmap_Color_Superclass)

#abbreviate long metabolites by replacing with class identifier
if (length(which(str_length(colnames(metabolite_relative_abundance_subset_nodup_ct)) > 0))) { #check if any metabolite names are over 50 characters
    long_metabolite_names <- which(str_length(colnames(metabolite_relative_abundance_subset_nodup_ct)) >= 50)
    output_abbreviation_key <- colnames(metabolite_relative_abundance_subset_nodup_ct)[long_metabolite_names] #maps original annotation to updated annotation
    for (p in 1:length(long_metabolite_names)) {
      index <- which(metabolite_classes$Original_Compound_Name == output_abbreviation_key[p])
      if (length(index > 0)) {
        if (is.na(metabolite_classes$Super_Class[index])) { #long-named metabolite was not given a class name
            colnames(metabolite_relative_abundance_subset_nodup_ct)[long_metabolite_names[p]] <- str_sub(colnames(metabolite_relative_abundance_subset_nodup_ct)[long_metabolite_names[p]], start = 1, end = 47)
            colnames(metabolite_relative_abundance_subset_nodup_ct)[long_metabolite_names[p]] <- paste0(colnames(metabolite_relative_abundance_subset_nodup_ct)[long_metabolite_names[p]], "...")
        } else {
          colnames(metabolite_relative_abundance_subset_nodup_ct)[long_metabolite_names[p]] <- metabolite_classes$Class_Name_Heatmap[index]
        }
     }
        output_abbreviation_key <- cbind(output_abbreviation_key, colnames(metabolite_relative_abundance_subset_nodup_ct)[long_metabolite_names])
    }
}

metabolite_relative_abundance_subset_nodup_ct <- as.data.frame(t(metabolite_relative_abundance_subset_nodup_ct))
metabolite_relative_abundance_subset_nodup_ct$name <- rownames(metabolite_relative_abundance_subset_nodup_ct)

#remove "Esi+" and asterisks from the endings of metabolite names
if (length(grep('Esi+', metabolite_relative_abundance_subset_nodup_ct$name))) {
    esis <- grep('Esi+', metabolite_relative_abundance_subset_nodup_ct$name)
    metabolite_relative_abundance_subset_nodup_ct$name[esis] <- gsub('Esi+', '\\1?\\2', metabolite_relative_abundance_subset_nodup_ct$name[esis]) #replace "Esi" with a "?" for easier regex
    metabolite_relative_abundance_subset_nodup_ct$name[esis] <- sub('(?<=\\?).*$', '', metabolite_relative_abundance_subset_nodup_ct$name[esis], perl=TRUE) #remove everything after question mark
    metabolite_relative_abundance_subset_nodup_ct$name[esis] <- str_sub(metabolite_relative_abundance_subset_nodup_ct$name[esis], end = -2) #remove the ? and space at the end of those metabolite names
}

#drop the asterisks from the names of some metabolites iteratively
    if(length(which(str_sub(metabolite_relative_abundance_subset_nodup_ct$name, start = -1) == "*")) > 0) {
    asterisk_names <- which(str_sub(metabolite_relative_abundance_subset_nodup_ct$name, start = -1) == "*")
    metabolite_relative_abundance_subset_nodup_ct$name[asterisk_names] <- str_sub(metabolite_relative_abundance_subset_nodup_ct$name[asterisk_names], end = -2) #drop the asterisk
}
if(length(which(str_sub(metabolite_relative_abundance_subset_nodup_ct$name, start = -1) == "*")) > 0) {
    asterisk_names <- which(str_sub(metabolite_relative_abundance_subset_nodup_ct$name, start = -1) == "*")
    metabolite_relative_abundance_subset_nodup_ct$name[asterisk_names] <- str_sub(metabolite_relative_abundance_subset_nodup_ct$name[asterisk_names], end = -2) #drop the asterisk
}

#mapper for metabolite final name to class and color
color_function <- metabolite_relative_abundance_subset_nodup_ct[, (ncol(metabolite_relative_abundance_subset_nodup_ct)-2):ncol(metabolite_relative_abundance_subset_nodup_ct)]

#drop the class and color columns from the metabolite dataset
metabolite_relative_abundance_subset_nodup_ct <- metabolite_relative_abundance_subset_nodup_ct[, -(ncol(metabolite_relative_abundance_subset_nodup_ct) - 1)]
metabolite_relative_abundance_subset_nodup_ct <- metabolite_relative_abundance_subset_nodup_ct[, -(ncol(metabolite_relative_abundance_subset_nodup_ct) - 1)]

#convert back into a numeric data frame and set column names appropriately
column_names <- metabolite_relative_abundance_subset_nodup_ct$name
metabolite_relative_abundance_subset_nodup_ct <- data.frame(t(metabolite_relative_abundance_subset_nodup_ct))
metabolite_relative_abundance_subset_nodup_ct <- metabolite_relative_abundance_subset_nodup_ct[-nrow(metabolite_relative_abundance_subset_nodup_ct), ]
colnames(metabolite_relative_abundance_subset_nodup_ct) <- column_names
metabolite_relative_abundance_subset_nodup_ct_copy <- metabolite_relative_abundance_subset_nodup_ct
metabolite_relative_abundance_subset_nodup_ct <- data.frame(lapply(metabolite_relative_abundance_subset_nodup_ct, function(x) as.numeric(as.character(x))))
colnames(metabolite_relative_abundance_subset_nodup_ct) <- colnames(metabolite_relative_abundance_subset_nodup_ct_copy)

#Average the intensity measures for metabolites that overlap in their new annotation 
#i.e. 2 different metabolite features went into DIABLO, but their new annotation is the same 
#so we combine them
if (length(unique(colnames(metabolite_relative_abundance_subset_nodup_ct))) != length(colnames(metabolite_relative_abundance_subset_nodup_ct))) {
    duplicated_metabolites <- unique(colnames(metabolite_relative_abundance_subset_nodup_ct)[which(duplicated(colnames(metabolite_relative_abundance_subset_nodup_ct)))])
    for (i in 1:length(duplicated_metabolites)) { #for each metabolite that is duplicated...
        specific_duplicated_indices <- which(colnames(metabolite_relative_abundance_subset_nodup_ct) == duplicated_metabolites[i]) #get nth metabolite that is duplicated
        first_index <- specific_duplicated_indices[1]
        subsequent_indices <- specific_duplicated_indices[2:length(specific_duplicated_indices)]
        for (j in subsequent_indices) { #iterate through subsequent duplications
          metabolite_relative_abundance_subset_nodup_ct[, first_index] <- metabolite_relative_abundance_subset_nodup_ct[, first_index] + metabolite_relative_abundance_subset_nodup_ct[, j] #add to first instance of duplicated column
        }
        metabolite_relative_abundance_subset_nodup_ct[, first_index] <- metabolite_relative_abundance_subset_nodup_ct[, first_index] / length(specific_duplicated_indices) #average all duplicated columns
        metabolite_relative_abundance_subset_nodup_ct <- metabolite_relative_abundance_subset_nodup_ct[, -subsequent_indices] #remove subsequent instances
    }
}

#generate list of omics data for spiec-easi
omics_list <- list(as.matrix(otu_relative_abundance_subset_nodup), as.matrix(metabolite_relative_abundance_subset_nodup_ct))
df <- cbind(otu_relative_abundance_subset_nodup, metabolite_relative_abundance_subset_nodup_ct)

#Generate pearson correlation matrix
cormat <- cor(df, method = "pearson")

#Get p-values for correlations and adjust to create p-value matrix
res1 <- cor.mtest(as.matrix(df), conf.level = .95) #return a list of p-values from significance test for these correlations
pAdj <- p.adjust(c(res1[[1]]), method = "BH") #correct p-values using Benjamini-Hochberg correction
resAdj <- matrix(pAdj, ncol = dim(res1[[1]])[1]) #p-value matrix

#Create lower triangle heatmap showing pairwise feature correlations
cormat_table <- data.frame(Feature_1 = rownames(cormat)[row(cormat)[upper.tri(cormat)]], 
                           Feature_2 = colnames(cormat)[col(cormat)[upper.tri(cormat)]], 
                           Pearson_Correlation = cormat[upper.tri(cormat)],
                           BH_Corrected_P_Value = resAdj[upper.tri(resAdj)]) #convert to 3-column data table

#Output heatmap correlations and p-values to table
write.csv(cormat_table, file = "~/Desktop/Spiromics_BAL_Microbiome_Metabolome_V2/Spiromics_Colorado_Stratified/Spiromics_No_Severe_No_Healthy/Output/Heatmap/DIABLO_Correlation_Heatmap/BALF_Heatmap_Group_Selected_Features_CT_Data.csv", row.names = FALSE)

#Generate heatmap
pdf("~/Desktop/Spiromics_BAL_Microbiome_Metabolome_V2/Spiromics_Colorado_Stratified/Spiromics_No_Severe_No_Healthy/Output/Heatmap/DIABLO_Correlation_Heatmap/BALF_Heatmap_Group_Selected_Features_No_Unknowns_CT_Data_.pdf", width = 15, height = 15)
corrplot(cormat, order = "hclust",
                 addrect = 3,
                 p.mat = resAdj, #p-values from correlation test, Benjamini-Hochberg corrected
                 #title = "All Pairwise Pearson Correlations for Top OTU and Metabolite Features",
                 method = "color", #only color used to signify corrleation
                 type = "lower", #lower triangle
                 tl.col = "black",
                 tl.cex = 0.5,
                 outline = FALSE, 
                 tl.srt = 45, #rotate axes labels 45 degrees
                 insig = "blank", #leave blank on no significant correlations,
                 pch.col = "red",
                 sig.level = c(.001, .01, .05),
                 pch.cex = 0.75,
                 diag = TRUE)
                 #mar=c(0,0,2,0)) #adjust margins so title is in plot
dev.off()

#Subset to just OTU-metabolite pairwise comparisons
cormat <- cormat[1:ncol(otu_relative_abundance_subset_nodup), (ncol(otu_relative_abundance_subset_nodup) + 1):ncol(cormat)] #OTU vs. Metabolite comparisons only
resAdj <- resAdj[1:ncol(otu_relative_abundance_subset_nodup), (ncol(otu_relative_abundance_subset_nodup) + 1):ncol(resAdj)]

#Generate subsetted heatmap
pdf("~/Desktop/Spiromics_BAL_Microbiome_Metabolome_V2/Spiromics_Colorado_Stratified/Spiromics_No_Severe_No_Healthy/Output/Heatmap/DIABLO_Correlation_Heatmap/BALF_Heatmap_Group_Selected_Features_OTU_vs_Met_No_Unknowns_CT_Data_.pdf", width = 15, height = 15)
corrplot(cormat, p.mat = resAdj, #p-values from correlation test, Benjamini-Hochberg corrected
                 title = "Pairwise Pearson Correlations for Top OTU and Metabolite Features",
                 method = "color", #only color 0.1used to signify corrleation
                 tl.col = "black",
                 tl.cex = 0.75,
                 tl.srt = 60, #rotate axes labels 45 degrees
                 insig = "label_sig", #leave blank on no significant correlations,
                 pch.col = "red",
                 sig.level = c(.001, .01, .05), 
                 pch.cex = 1,
                 mar=c(0,0,2,0)) #adjust margins so title is in plot
dev.off()

#Subset to just those OTUs with at least 1 significant correlation after BH correction to shrink the plot
sig_otus <- NULL
sig_otu_columns <- NULL #column indices
for (i in 1:ncol(cormat)) {
    if(min(resAdj[, i]) < sig_threshold) {
        sig_otus <- c(sig_otus, colnames(cormat)[i]) #append OTU to list of OTUs with at least 1 sig correlation with a metabolite
        sig_otu_columns <- c(sig_otu_columns, i)
    }
}
cormat <- cormat[, sig_otu_columns, drop = FALSE] #subset correlation matrix
resAdj <- resAdj[, sig_otu_columns, drop = FALSE]

#Subset to just those metabolites with at least 1 significant correlation after BH correction to shrink the plot
sig_met <- NULL
sig_met_rows <- NULL #column indices
for (i in 1:nrow(cormat)) {
    if(min(resAdj[i, ]) < sig_threshold) {
        sig_met <- c(sig_met, rownames(cormat)[i]) #append OTU to list of OTUs with at least 1 sig correlation with a metabolite
        sig_met_rows <- c(sig_met_rows, i)
    }
}

#color metabolites by class
colnames(color_function) <- c("Super_Class", "Heatmap_Color_Superclass", "name")
color_function <- color_function[unique(color_function$name), ] #drop duplicated metabolites
color_function <- color_function[which(color_function$name %in% colnames(cormat)), ] #filter to the metabolites in the correlation matrix
color_function$Heatmap_Color_Superclass[which(is.na(color_function$Heatmap_Color_Superclass))] <- "goldenrod"
cormat <- cormat[, -73]
resAdj <- resAdj[, -73]
list_of_classes <- unique(color_function$Super_Class)
list_of_colors <- col2hex(unique(color_function$Heatmap_Color_Superclass))
names(list_of_colors) <- list_of_classes
col_ha <- HeatmapAnnotation(Super_Class = color_function$Super_Class, 
                            col = list(Super_Class = list_of_colors), 
                            annotation_legend_param = list(Super_Class = list(nrow = 2, labels_gp = gpar(fontsize = 12))),
                            show_annotation_name = FALSE)

#if any OTUs have a significant correlation with metabolites, plot only significant OTU-metabolite correlations
if (length(sig_met_rows) > 0) {
    cormat <- cormat[sig_met_rows, , drop = FALSE] #subset correlation matrix
    resAdj <- resAdj[sig_met_rows, , drop = FALSE]
    write.csv(cormat, file = "~/Desktop/Spiromics_BAL_Microbiome_Metabolome_V2/Spiromics_Colorado_Stratified/Spiromics_No_Severe_No_Healthy/Output/Heatmap/DIABLO_Correlation_Heatmap/BALF_Heatmap_Group_Selected_Features_OTU_vs_Met_No_Unknowns_Significant_CT_Data.csv", row.names = TRUE)

    pdf("~/Desktop/Spiromics_BAL_Microbiome_Metabolome_V2/Spiromics_Colorado_Stratified/Spiromics_No_Severe_No_Healthy/Output/Heatmap/DIABLO_Correlation_Heatmap/BALF_Heatmap_Group_Selected_Features_OTU_vs_Met_No_Unknowns_Significant_CT_Data_.pdf", width = 30, height = 10)
    final_heatmap <- Heatmap(cormat, 
                            #add asterisks based on appropriate significance levels
                            cell_fun = function(j, i, x, y, width, height, fill) {
                                if(resAdj[i, j] < 0.001) {
                                    grid.text("***", x, y)
                                } else if (resAdj[i, j] < 0.01) {
                                    grid.text("**", x, y)
                                } else if (resAdj[i, j] < 0.05) {
                                    grid.text("*", x, y)
                                }
                            },
                            rect_gp = gpar(col = "black", lwd = 1), #color the border of squares
                            show_column_dend = FALSE,
                            show_row_dend = FALSE,
                            row_order = order(rownames(cormat)), #sort in numeric order
                            column_order = colnames(cormat), #retain order of metabolites by class as in the matrix
                            column_names_rot = 45, 
                            column_names_side = "top",
                            row_names_side = "left",
                            heatmap_legend_param = list(title = "Pearson Correlation", #legend showing correlations
                                                        at = c(-1, -0.8, -0.6, -0.4, -0.2, 0, 0.2, 0.4, 0.6, 0.8, 1),
                                                        legend_height = unit(5, "cm"),
                                                        labels_gp = gpar(fontsize = 12),
                                                        direction = "horizontal"),
                            border = TRUE,
                            row_names_gp = gpar(fontsize = 12),
                            column_names_gp = gpar(fontsize = 12),
                            top_annotation = col_ha,
                            height = (nrow(cormat)*unit(6, "mm")),
                            width = (ncol(cormat)*unit(6, "mm")))
    draw(final_heatmap, heatmap_legend_side = "bottom", annotation_legend_side = "bottom")
    dev.off()
    #corrplot(cormat, p.mat = resAdj, #p-values from correlation test, Benjamini-Hochberg corrected
    #             #title = "Pairwise Pearson Correlations for Top OTU and Metabolite Features with 1+ Significant Correlation",
    #             method = "color", #only color used to signify corrleation
    #             cl.lim=c(-1,1), 
    #             col=colorRampPalette(c("blue","white","red"))(100),
    #             tl.col = "black",
    #             tl.cex = 0.8,
    #             tl.srt = 45, #rotate axes labels 45 degrees
    #             insig = "label_sig", #leave blank on no significant correlations,
    #             pch.col = "black",
    #             sig.level = c(.001, .01, .05), 
    #             pch.cex = 1,
    #             cl.pos = "r") #put legend to the right of plot
                # cl.ratio = 0.3,
                #cl.align.text = "r",
                # cl.offset = -0.5) #put labels for legend to the right of legend
} else {
    print("No Significant OTU-Metabolite Correlations")
}

