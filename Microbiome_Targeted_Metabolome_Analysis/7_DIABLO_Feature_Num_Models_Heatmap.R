suppressMessages(library(ComplexHeatmap))
suppressMessages(library(tidyverse))
suppressMessages(library(stringr))
suppressMessages(library(circlize))

df_raw <- read.csv('~/Desktop/Spiromics_BAL_Microbiome_Metabolome_V2/Spiromics_UNC_Stratified/Spiromics_No_Severe_No_Healthy/Output/Heatmap/DIABLO_Num_Models_Heatmap/Top_Features_by_Num_Sig_Models_No_Unknowns_Dir_Effect_Added.csv')

#df_raw$Number_of_Models_times_Dir_Effect <- df_raw$Number_of_Models * df_raw$Direction_of_Effect #code new column so intensity of square is number of models and color is direction of effect

categories <- unique(df_raw$Group) #iterate through categories of outcomes

#order of variables for heatmap
lf_order <- c("FEV1/FVC", 
              "FEV1 (%predicted or L)", 
              #"FEF25-75 (L/s)", 
              "PEFR (L/s)")
bdr_order <- c("FEV1 BDR (%)", "FVC BDR (% or mL)")
ex_order <- c("# Exacerbations in Past 12 mos. at Baseline",
              "# Exacerbations (HCU + AB/S) in Year 1")
              #"# Exacerbations (HCU + AB/S) Post-Bronchoscopy")

for (cat in categories) {
  print(cat)
  #text sizes for the heatmap
  if (cat == "Exacerbations") {
    text_size <- 10
  } else if (cat == "Lung Function") {
    text_size <- 10
  } else if (cat == "Bronchodilator Response") {
    text_size <- 10
  } else {
    text_size <- 10
  }

  #set up data frame
  df <- df_raw[which(df_raw$Group == cat), ] #subset to category of outcomes

  #narrow data so convert to a wide matrix with dplyr
  df_dir_effect <- df %>% 
                select(Phenotype_Figure_Name, Feature, Median_Log_OR_In_Final_Models_Where_Feature_Selected) %>% 
                spread(Feature, Median_Log_OR_In_Final_Models_Where_Feature_Selected)

  #replace positive and negative in direction of effect with + and -
  rownames(df_dir_effect) <- df_dir_effect[, 1]
  df_dir_effect <- df_dir_effect[, -1]
  df_dir_effect <- as.matrix(t(df_dir_effect))
  #df_dir_effect <- replace(df_dir_effect, is.na(df_dir_effect), 0)

  #discretize the median log odds ratios as (-1, -0.5], (-0.5, 0], (0, 0.5], (0.5, 1]
  for (i in 1:ncol(df_dir_effect)) {
      df_dir_effect[, i] <- cut(df_dir_effect[, i], seq(-1, 1, 0.5))
  }

  #Heatmap functions
  #col_fun <- colorRamp2(c(-1, 0, 1), 
  #                      c("blue", "white", "red")) #map 0 to white and 1000 to red with range of colors
  #colors = structure(c(-1, 0, 1), names = c("blue", "white", "red")) # black, red, green, blue

  #***********************************************************#
  #OTUs alone
  otus <- which(str_sub(rownames(df_dir_effect), start = 1, end = 3) == "Otu")
  df_otus <- df_dir_effect[otus, , drop = FALSE]

  #column orders for the heatmap
  if (cat == "Bronchodilator Response") {
    column_order <- bdr_order
  } else if (cat == "Lung Function") {
    column_order <- lf_order
  } else if (cat == "Exacerbations") {
    column_order <- ex_order
  } else {
    column_order <- sort(colnames(df_otus))
  }

#legend as annotation
#ha <- HeatmapAnnotation(col = list(bar = c("-1" = "blue", "0" = "black", "1" = "red")))

#heatmap of otus alone
  if (length(otus > 0)) {
      pdf(paste("~/Desktop/Spiromics_BAL_Microbiome_Metabolome_V2/Spiromics_UNC_Stratified/Spiromics_No_Severe_No_Healthy/Output/Heatmap/DIABLO_Num_Models_Heatmap/DIABLO_Results_Num_Models_Heatmap_", cat, "_OTUs.pdf", sep = ""), 
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
  }

  #***********************************************************#
  #Metabolites alone
  if (length(otus) > 0) {
    df_metabolites <- df_dir_effect[-otus, , drop = FALSE] #remove OTUs, if any
  } else {
    df_metabolites <- df_dir_effect #otherwise, use entire table
  }

  #end if no metabolites
  if (nrow(df_otus) == nrow(df_dir_effect)) {
    next
  }
  #read and preprocess metabolite data
  met_raw <- read.csv("~/Desktop/Spiromics_BAL_Microbiome_Metabolome_V2/Spiromics_UNC_Stratified/Spiromics_No_Severe_No_Healthy/Data/Esther_Metabolites_No_Severe_No_Healthy.csv")
  rownames(met_raw) <- met_raw[, 1]
  met_raw <- met_raw[, -1]

  met_indices <- which(colnames(met_raw) %in% rownames(df_metabolites))

  met_raw <- met_raw[, met_indices, drop = FALSE]

  met <- read.csv("~/Desktop/Spiromics_BAL_Microbiome_Metabolome_V2/Spiromics_UNC_Stratified/Spiromics_No_Severe_No_Healthy/Data/Esther_Metabolites_No_Severe_No_Healthy.csv", check.names = FALSE)
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

  if (cat == "Bronchodilator Response") {
    column_order <- bdr_order
  } else if (cat == "Lung Function") {
    column_order <- lf_order
  } else {
    column_order <- sort(colnames(df_metabolites))
  }

  pdf(paste("~/Desktop/Spiromics_BAL_Microbiome_Metabolome_V2/Spiromics_UNC_Stratified/Spiromics_No_Severe_No_Healthy/Output/Heatmap/DIABLO_Num_Models_Heatmap/DIABLO_Results_Num_Models_Heatmap_", cat, "_Metabolites.pdf", sep = ""), width = 10, height = 10)
  metabolites_heatmap <- Heatmap(df_metabolites, 
              na_col = 'black',
              #col = col_fun, 
              #name = "Number of Models (Out of 1000)",
              rect_gp = gpar(col = "black", lwd = 1), #color the border of squares
              show_column_dend = FALSE,
              show_row_dend = FALSE,
              row_labels = str_wrap(rownames(df_metabolites), width = 50),
              row_order = order(rownames(df_metabolites)),
              column_order = column_order,
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
              width = unit(3, "cm"))
  draw(metabolites_heatmap,
       heatmap_legend_side = "right")
       #show_heatmap_legend = FALSE)
  dev.off()

  if (nrow(df_otus) != 0) {
    ht_list <- otu_heatmap %v% metabolites_heatmap
    pdf(paste("~/Desktop/Spiromics_BAL_Microbiome_Metabolome_V2/Spiromics_UNC_Stratified/Spiromics_No_Severe_No_Healthy/Output/Heatmap/DIABLO_Num_Models_Heatmap/DIABLO_Results_Num_Models_Heatmap_", cat, ".pdf", sep = ""), width = 10, height = 10)
    draw(ht_list,
          show_heatmap_legend = FALSE)
    dev.off()
  }
}
