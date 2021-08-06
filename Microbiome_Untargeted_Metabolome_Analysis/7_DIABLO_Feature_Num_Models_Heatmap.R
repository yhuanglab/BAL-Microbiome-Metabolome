suppressMessages(library(ComplexHeatmap))
suppressMessages(library(tidyverse))
suppressMessages(library(stringr))
suppressMessages(library(circlize))
suppressMessages(library(gplots))

df_raw <- read.csv('~/Desktop/Spiromics_BAL_Microbiome_Metabolome_V2/Spiromics_Colorado_Stratified/Spiromics_No_Severe_No_Healthy/Output/Heatmap/DIABLO_Num_Models_Heatmap/Top_Features_by_Num_Sig_Models_No_Unknowns_Dir_Effect_Added.csv')

#df_raw$Number_of_Models_times_Dir_Effect <- df_raw$Number_of_Models * df_raw$Direction_of_Effect #code new column so intensity of square is number of models and color is direction of effect

categories <- unique(df_raw$Group) #iterate through categories of outcomes

#order of variables for heatmap
lf_order <- c("FEV1/FVC", 
              "FEV1 (%predicted or L)", 
              "FEF25-75 (L/s)",
              "PEFR (L/s)")
bdr_order <- c("FEV1 BDR (%)",
              "FVC BDR (mL)")
ex_order <- c(#"# Exacerbations in Past 12 mos. at Baseline")
              "# Exacerbations (HCU + AB/S) in Year 1",
              "# Exacerbations (HCU + AB/S) Post-Bronchoscopy") 

for (cat in categories) {
  print(cat)
  #text sizes for the heatmap
  if (cat == "Exacerbations") {
    text_size <- 10
  } else if (cat == "Lung Function") {
    text_size <- 9
  } else if (cat == "Bronchodilator Response") {
    text_size <- 10
  } else {
    text_size <- 10
  }

  #set up data frame
  df <- df_raw[which(df_raw$Group == cat), ] #subset to category of outcomes

  #exponentiate log odds ratio column to get odds ratios
  #df$Median_Log_OR_In_Final_Models_Where_Feature_Selected <- exp(df$Median_Log_OR_In_Final_Models_Where_Feature_Selected)

  #narrow data so convert to a wide matrix with dplyr
  df_dir_effect <- df %>% 
  					select(Phenotype_Figure_Name, Feature, Median_Log_OR_In_Final_Models_Where_Feature_Selected) %>% 
  					spread(Feature, Median_Log_OR_In_Final_Models_Where_Feature_Selected)

  #convert matrix to heatmap and set row and column names
  rownames(df_dir_effect) <- df_dir_effect[, 1]
  df_dir_effect <- df_dir_effect[, -1]
  df_dir_effect <- as.matrix(t(df_dir_effect))
  #df_dir_effect <- replace(df_dir_effect, is.na(df_dir_effect), 0)

  #discretize the median log odds ratios as (-1, -0.5], (-0.5, 0], (0, 0.5], (0.5, 1]
  for (i in 1:ncol(df_dir_effect)) {
      df_dir_effect[, i] <- cut(df_dir_effect[, i], seq(-1, 1, 0.5))
  }

  #***********************************************************#
  #OTUs alone
  otus <- which(str_sub(rownames(df_dir_effect), start = 1, end = 3) == "Otu")
  df_otus <- df_dir_effect[otus, , drop = FALSE]

  #column orders for the heatmap
  if (cat == "Bronchodilator Response") {
    column_order <- bdr_order
    colnames(df_otus) <- c("FEV1 BDR (%)", "FVC BDR (mL)")
  } else if (cat == "Lung Function") {
    column_order <- lf_order
    colnames(df_otus) <- c("FEV1/FVC", "FEV1 (%predicted or L)", "FEF25-75 (L/s)", "PEFR (L/s)")
  } else if (cat == "Exacerbations") {
    column_order <- ex_order
  } else {
    column_order <- sort(colnames(df_otus))
  }

  #heatmap of otus alone
  if (length(otus > 0)) {
      pdf(paste("~/Desktop/Spiromics_BAL_Microbiome_Metabolome_V2/Spiromics_Colorado_Stratified/Spiromics_No_Severe_No_Healthy/Output/Heatmap/DIABLO_Num_Models_Heatmap/DIABLO_Results_Num_Models_Heatmap_", cat, "_OTUs.pdf", sep = ""), 
            width = 10, 
            height = 10)
      otu_heatmap <- Heatmap(df_otus, 
                              na_col = 'black',
      						            col = c("1" = "blue", "2" = "lightskyblue", "3" = "lightsalmon", "4" = "red"), #color based on log OR magnitude
      						            #name = "Number of Models (Out of 1000)",
      						            rect_gp = gpar(col = "black", lwd = 1), #color the border of squares
                          		show_column_dend = FALSE,
                        			show_row_dend = FALSE,
                              #row_order = order(as.numeric(gsub("Otu", "", rownames(df_otus)))), #order OTU number
                              row_order = order(str_sub(rownames(df_otus), start = 9)), #order by genus name in alphabetical order
                              column_order = column_order,
                              #lgd = Legend(at = seq(1, 4, by = 1), 
                              #              title = "Median Log Odds Ratio", 
                              #               labels = c("(-1, -0.5]", "(-0.5, 0]", "(0, 0.5]", "(0.5, 1]")),
                        		  heatmap_legend_param = list(title = "Median Log Odds Ratio",
                                                          at = c("1", "2", "3", "4"),
                                                          labels = c("(-1, -0.5]", 
                                                                     "(-0.5, 0]",
                                                                     "(0, 0.5]",
                                                                     "(0.5, 1]"),
                                                          direction = "vertical",
                                                          legend_height = unit(10, "cm")),
                              #heatmap_legend_param = list(legend_direction = "horizontal", legend_width = unit(5, "cm")),
                        			column_names_rot = 45, column_names_side = "top", row_names_side = "left",
                        			border = TRUE,
                        			row_names_gp = gpar(fontsize = text_size, lwd = 20),
                        			column_names_gp = gpar(fontsize = text_size),
                        			width = unit(3, "cm"))
      draw(otu_heatmap,
           heatmap_legend_side = "right")
           #show_heatmap_legend = FALSE)
      dev.off()
  #}

  #***********************************************************#
  #Metabolites alone
  if (length(otus) > 0) {
    df_metabolites <- df_dir_effect[-otus, , drop = FALSE] #remove OTUs, if any
  } else {
    df_metabolites <- df_dir_effect #otherwise, use entire table
  }

  #end if no metabolites
  #if (nrow(df_otus) == nrow(df_dir_effect)) {
  #  next
  #}
  
  #remove R formatting for metabolite names and replace with original name with correct special characters
  met_raw <- read.csv("~/Desktop/Spiromics_BAL_Microbiome_Metabolome_V2/Spiromics_Colorado_Stratified/Spiromics_No_Severe_No_Healthy/Data/BALF_C18_Norm_Log2_Transposed_No_Severe_No_Healthy.csv")
  rownames(met_raw) <- met_raw[, 1]
  met_raw <- met_raw[, -1]
  met_indices <- which(colnames(met_raw) %in% rownames(df_metabolites))
  met_raw <- met_raw[, met_indices, drop = FALSE]
  met <- read.csv("~/Desktop/Spiromics_BAL_Microbiome_Metabolome_V2/Spiromics_Colorado_Stratified/Spiromics_No_Severe_No_Healthy/Data/BALF_C18_Norm_Log2_Transposed_No_Severe_No_Healthy.csv", check.names = FALSE)
  rownames(met) <- met[, 1]
  met <- met[, -1]
  met <- met[, met_indices, drop = FALSE]
  mapping <- cbind(colnames(met_raw), colnames(met))
  for (i in 1:nrow(df_metabolites)) {
    column_name <- rownames(df_metabolites)[i]
    mapped_index <- which(colnames(met_raw) == column_name)
    mapped_name <- colnames(met)[mapped_index]
    rownames(df_metabolites)[i] <- mapped_name
  }

  metabolite_classes <- read.csv("~/Desktop/Spiromics_BAL_Microbiome_Metabolome_V2/04_04_21_DIABLO_Significant_Metabolites_List_Updated_KAD_NR.csv", stringsAsFactors = FALSE)
  #map name of each metabolite to class
  row_order_df <- data.frame(df_metabolites)
  row_order_df$Original_Compound_Name <- rownames(row_order_df)
  row_order_df <- left_join(row_order_df, metabolite_classes, by = "Original_Compound_Name")
  row_order_df <- row_order_df[, 1:(ncol(df_metabolites) + 7)] #select all columns up to Heatmap_Color column
  class_index <- which(colnames(row_order_df) == "Super_Class")
  row_order_df <- row_order_df[order(row_order_df$Super_Class), ]
  rownames(row_order_df) <- row_order_df$"Original_Compound_Name"
  row_order_df <- row_order_df[, c(1:ncol(df_metabolites), class_index, class_index + 2)]
  if (cat == "BAL Neutrophils") {
    row_order_df <- row_order_df[-nrow(row_order_df), ]
  }
  df_metabolites <- row_order_df

  #for metabolites whose names are longer than 50 characters, replace with the name of the metabolite's class
  if (length(which(str_length(rownames(df_metabolites)) > 0))) { #check if any metabolite names are over 50 characters
    long_metabolite_names <- which(str_length(rownames(df_metabolites)) >= 50)
    output_abbreviation_key <- rownames(df_metabolites)[long_metabolite_names] #maps original annotation to updated annotation
    for (p in 1:length(long_metabolite_names)) {
      index <- which(metabolite_classes$Original_Compound_Name == output_abbreviation_key[p])
      if (length(index > 0)) {
        if (is.na(metabolite_classes$Class[index])) { #long-named metabolite was not given a class name
          rownames(df_metabolites)[long_metabolite_names[p]] <- str_sub(rownames(df_metabolites)[long_metabolite_names[p]], start = 1, end = 47)
          rownames(df_metabolites)[long_metabolite_names[p]] <- paste0(rownames(df_metabolites)[long_metabolite_names[p]], "...")
        } else {
          rownames(df_metabolites)[long_metabolite_names[p]] <- metabolite_classes$Class_Name_Heatmap[index]
        }
      }
      output_abbreviation_key <- cbind(output_abbreviation_key, rownames(df_metabolites)[long_metabolite_names])
    }
  }

  #for metabolites less than 50 characters in length, replace with updated annotation if present
  #short_metabolite_names <- which(str_length(rownames(df_metabolites)) < 50)
  #short_output_abbreviation_key <- rownames(df_metabolites)[short_metabolite_names]
  #class_names <- which(str_sub(output_abbreviation_key, start = 1, end = 6) == "Class-")
  #short_output_abbreviation_key <- output_abbreviation_key[-class_names]
  #short_metabolite_names <- short_metabolite_names[-class_names]
  #for (u in 1:length(short_metabolite_names)) { #for short metabolite names, replace with updated annotation if present
  #    index <- which(metabolite_classes$Original_Compound_Name == short_output_abbreviation_key[u])
  #    if (length(index > 0)) {
  #      rownames(df_metabolites)[short_metabolite_names[u]] <- metabolite_classes$Updated_Compound_05_20_20[index]
  #    }
  }
  #print(df_metabolites)

  df_metabolites$name <- rownames(df_metabolites)

  #remove "Esi+" and asterisks from the endings of metabolite names
  if (length(grep('Esi+', df_metabolites$name))) {
    esis <- grep('Esi+', df_metabolites$name)
    df_metabolites$name[esis] <- gsub('Esi+', '\\1?\\2', df_metabolites$name[esis]) #replace "Esi" with a "?" for easier regex
    df_metabolites$name[esis] <- sub('(?<=\\?).*$', '', df_metabolites$name[esis], perl=TRUE) #remove everything after question mark
    df_metabolites$name[esis] <- str_sub(df_metabolites$name[esis], end = -2) #remove the ? and space at the end of those metabolite names
  }

  #drop the asterisks from the names of some metabolites iteratively
  if(length(which(str_sub(df_metabolites$name, start = -1) == "*")) > 0) {
    asterisk_names <- which(str_sub(df_metabolites$name, start = -1) == "*")
    df_metabolites$name[asterisk_names] <- str_sub(df_metabolites$name[asterisk_names], end = -2) #drop the asterisk
  }
  if(length(which(str_sub(df_metabolites$name, start = -1) == "*")) > 0) {
    asterisk_names <- which(str_sub(df_metabolites$name, start = -1) == "*")
    df_metabolites$name[asterisk_names] <- str_sub(df_metabolites$name[asterisk_names], end = -2) #drop the asterisk
  }

  #if new mapping overlaps across multiple clinical measures, average together in one row
  if (length(which(duplicated(df_metabolites$name) == TRUE) > 0)) {
    duplicated_indices <- which(duplicated(df_metabolites$name) == TRUE)
    duplicated_metabolites <- df_metabolites$name[duplicated_indices]
    for (a in 1:length(duplicated_metabolites)) {
      dup_ind <- which(df_metabolites$name == duplicated_metabolites[a])
      first_instance <- dup_ind[1]
      for (a in dup_ind[2:length(dup_ind)]) {
        df_metabolites[first_instance, 1:(ncol(df_metabolites) - 3)] <- (df_metabolites[first_instance, 1:(ncol(df_metabolites) - 3)] + df_metabolites[a, 1:(ncol(df_metabolites) - 3)])/2 #collapse duplicated rows
      }
      df_metabolites <- df_metabolites[-dup_ind[2:length(dup_ind)], ] #remove duplicated rows
    }
  }

  rownames(df_metabolites) <- df_metabolites$name
  df_metabolites <- df_metabolites[, 1:(ncol(df_metabolites) - 1)]

  #set column and row order for heatmap
  row_order <- df_metabolites[, c("Super_Class", "Heatmap_Color_Superclass"), drop = FALSE]
  row_order[which(is.na(row_order$Super_Class)), "Super_Class"] <- "Other" #set all metabolites with no class annotation to "Other"
  row_order[which(is.na(row_order$Heatmap_Color_Superclass)), "Heatmap_Color_Superclass"] <- "goldenrod" #set all no-class metabolites to black color
  #class annotations for the heatmap
  list_of_classes <- unique(row_order$Super_Class)
  print(list_of_classes)
  print(unique(row_order$Heatmap_Color_Superclass))
  list_of_colors <- col2hex(unique(row_order$Heatmap_Color_Superclass)) #returns hexadecimal code for each color
  names(list_of_colors) <- list_of_classes

  row_ha <- rowAnnotation(Super_Class = row_order$Super_Class, col = list(Super_Class = list_of_colors), show_annotation_name = FALSE)
  #set row order for the heatmap
  row_order <- rownames(row_order)
  write.csv(df_metabolites, file = paste("~/Desktop/Spiromics_BAL_Microbiome_Metabolome_V2/Spiromics_Colorado_Stratified/Spiromics_No_Severe_No_Healthy/Output/Heatmap/DIABLO_Num_Models_Heatmap/DIABLO_Results_Num_Models_Heatmap_", cat, "_List_of_Metabolites.csv", sep = ""))
  df_metabolites <- df_metabolites[, 1:(ncol(df_metabolites) - 2), drop = FALSE]
  #missing_class <- which(is.na(df_metabolites$Class))
  #df_metabolites$Class[missing_class] <- "None"

  if (cat == "Bronchodilator Response") {
    column_order <- bdr_order
    colnames(df_metabolites) <- c("FEV1 BDR (%)", "FVC BDR (mL)")
  } else if (cat == "Lung Function") {
    column_order <- lf_order
    colnames(df_metabolites) <- c("FEV1/FVC", "FEV1 (%predicted or L)", "FEF25-75 (L/s)", "PEFR (L/s)")
  } else if (cat == "BAL Neutrophils") {
    colnames(df_metabolites) <- c("BAL Neutrophil %")
  } else {
    column_order <- sort(colnames(df_metabolites))
  }

  pdf(paste("~/Desktop/Spiromics_BAL_Microbiome_Metabolome_V2/Spiromics_Colorado_Stratified/Spiromics_No_Severe_No_Healthy/Output/Heatmap/DIABLO_Num_Models_Heatmap/DIABLO_Results_Num_Models_Heatmap_", cat, "_Metabolites.pdf", sep = ""), width = 10, height = 10)
  metabolites_heatmap <- Heatmap(df_metabolites, 
              na_col = 'black',
              #col = col_fun, 
              #name = "Number of Models (Out of 1000)",
              rect_gp = gpar(col = "black", lwd = 1), #color the border of squares
              show_column_dend = FALSE,
              show_row_dend = FALSE,
              #row_labels = str_wrap(rownames(df_metabolites), width = 50),
              #row_order = order(rownames(df_metabolites)),
              #row_order = sort(df_metabolites$Class),
              row_order = row_order,
              column_order = column_order,
              #column_labels = str_wrap(colnames(df_metabolites), width = 30),
              col = c("1" = "blue", "2" = "lightskyblue", "3" = "lightsalmon", "4" = "red"), #color based on log OR magnitude
              heatmap_legend_param = list(title = "Median Log Odds Ratio",
                                          at = c("1", "2", "3", "4"),
                                          labels = c("(-1, -0.5]", 
                                                     "(-0.5, 0]",
                                                     "(0, 0.5]",
                                                     "(0.5, 1]"),
                                          direction = "vertical",
                                          legend_height = unit(10, "cm")),
              column_names_rot = 45, column_names_side = "top", row_names_side = "left",
              border = TRUE,
              row_names_gp = gpar(fontsize = text_size, lwd = 20),
              column_names_gp = gpar(fontsize = text_size),
              left_annotation = row_ha,
              width = unit(3, "cm"))
  #metabolites_heatmap <- metabolites_heatmap + rowAnnotation(df = color_function, col = )
  draw(metabolites_heatmap,
       heatmap_legend_side = "right")
       #show_heatmap_legend = FALSE)
  dev.off()

  if (nrow(df_otus) != 0) {
    ht_list <- otu_heatmap %v% metabolites_heatmap
    pdf(paste("~/Desktop/Spiromics_BAL_Microbiome_Metabolome_V2/Spiromics_Colorado_Stratified/Spiromics_No_Severe_No_Healthy/Output/Heatmap/DIABLO_Num_Models_Heatmap/DIABLO_Results_Num_Models_Heatmap_", cat, ".pdf", sep = ""), width = 10, height = 10)
    draw(ht_list,
          show_heatmap_legend = FALSE)
    dev.off()
  }
}
