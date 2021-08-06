#load libraries
library(devtools)
library(stringr)
#install_github("zdk123/SpiecEasi")
library(SpiecEasi)
library(igraph)
library(Matrix)
library(Hmisc)
library(corrplot)
library(mixOmics)
library(dplyr)

#significance threshold for what feature correlations are shown
sig_threshold <- 0.01

#should the network only show features that are significantly correlated across datasets?
only_use_sig_cross_dataset_correlation_features <- TRUE

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
col_order_df <- col_order_df[, c("Original_Compound_Name", "Class")]
col_order_df <- col_order_df[order(col_order_df$Class), ]
col_order_df[which(is.na(col_order_df$Class)), "Class"] <- "None"
metabolite_relative_abundance_subset_nodup_ct <- metabolite_relative_abundance_subset_nodup_ct[, col_order_df$"Original_Compound_Name"]

#abbreviate long metabolites by replacing with class identifier
if (length(which(str_length(colnames(metabolite_relative_abundance_subset_nodup_ct)) > 0))) { #check if any metabolite names are over 50 characters
    long_metabolite_names <- which(str_length(colnames(metabolite_relative_abundance_subset_nodup_ct)) >= 50)
    output_abbreviation_key <- colnames(metabolite_relative_abundance_subset_nodup_ct)[long_metabolite_names] #maps original annotation to updated annotation
    for (p in 1:length(long_metabolite_names)) {
      index <- which(metabolite_classes$Original_Compound_Name == output_abbreviation_key[p])
      if (length(index > 0)) {
        if (is.na(metabolite_classes$Class[index])) { #long-named metabolite was not given a class name
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
omics_combined <- cbind(otu_relative_abundance_subset_nodup, metabolite_relative_abundance_subset_nodup_ct)

#filter to those features that were significantly correlated across datasets
features_list <- read.csv('~/Desktop/Spiromics_BAL_Microbiome_Metabolome_V2/Spiromics_Colorado_Stratified/Spiromics_No_Severe_No_Healthy/Output/Heatmap/DIABLO_Correlation_Heatmap/BALF_Heatmap_Group_Selected_Features_OTU_vs_Met_No_Unknowns_Significant.csv', stringsAsFactors = FALSE, check.names = FALSE)
features <- colnames(features_list)
features <- c(features, features_list[, 1])
features <- features[-1]
features <- data.frame(features)
colnames(features)[1] <- "Feature"
#if (only_use_sig_cross_dataset_correlation_features == TRUE) {
#  omics_combined <- omics_combined[, which(colnames(omics_combined) %in% features$Feature)]
#}

#Generate pearson correlation matrix
cormat <- cor(omics_combined, method = "pearson")

#Get p-values for correlations and adjust to create p-value matrix
res1 <- cor.mtest(as.matrix(omics_combined), conf.level = .95) #return a list of p-values from significance test for these correlations
pAdj <- p.adjust(c(res1[[1]]), method = "BH") #correct p-values using Benjamini-Hochberg correction
resAdj <- matrix(pAdj, ncol = dim(res1[[1]])[1]) #p-value matrix

#filter cormat to only those features in the correlation heatmap
if (only_use_sig_cross_dataset_correlation_features == TRUE) {
  #omics_combined <- omics_combined[, which(colnames(omics_combined) %in% features$Feature)]
  sig_feature_indices <- which(rownames(cormat) %in% as.character(features$Feature))
  cormat <- cormat[sig_feature_indices, ]
  cormat <- cormat[, sig_feature_indices]
  resAdj <- resAdj[sig_feature_indices, ]
  resAdj <- resAdj[, sig_feature_indices]
}

#set insignificant correlations to 0
insig_correlation_indices <- which(resAdj >= sig_threshold)
cormat_mod <- cormat
cormat_mod[insig_correlation_indices] <- 0

#remove features without any significant correlations from rows and columns
insig_features <- which(rowSums(abs(cormat_mod)) == 1)
if (length(insig_features) > 0) {
  cormat_mod <- cormat_mod[-insig_features, , drop = FALSE] #features with no significant correlations would have only a correlation of 1 (with itself)
  cormat_mod <- cormat_mod[, -insig_features, drop = FALSE]
}

#create table to match vertex indices of graph to feature names
map_table <- data.frame(dimnames(cormat_mod)[2])
map_table$"Vertex Index" <- rownames(map_table)
colnames(map_table)[1] <- "Feature"

#remove vertex names from graph
dimnames(cormat_mod) <- NULL

#create igraph object
ig <- graph.adjacency(
  as.matrix(as.dist(abs(cormat_mod))),
  mode="undirected",
  weighted=TRUE,
  diag=FALSE
)

#remove isolated vertices from igraph object
ig <- delete.vertices(simplify(ig), degree(ig) == 0)

#plot raw igraph without isolated vertices
set.seed(1000)
#layout for igraph
l <- layout_with_fr(ig)
pdf('~/Desktop/Spiromics_BAL_Microbiome_Metabolome_V2/Spiromics_Colorado_Stratified/Spiromics_No_Severe_No_Healthy/Output/Heatmap/DIABLO_Pearson_Network/Pearson_Network.pdf')
plot(ig, layout = l)
dev.off()

#cluster vertices on plot using greedy clustering
clusters <- cluster_fast_greedy(ig)

#count number of clusters with more than 1 member and sort in ascending order
dupl_clusters <- sort(unique(clusters$membership[which(duplicated(clusters$membership) == TRUE)]))

#iterate through each cluster and find out which features belonged in that cluster
cluster_df <- NULL #data frame with clusters and features belonging to each cluster
for (i in dupl_clusters) {
  cluster_indices <- which(clusters$membership == i)
  cluster_features <- as.character(map_table$Feature[which(map_table$"Vertex Index" %in% cluster_indices)])
  i_df <- data.frame(cluster_indices, cluster_features, rep(i, length(cluster_features)))
  cluster_df <- rbind(cluster_df, i_df)
}

#save table of features and the cluster they are in
colnames(cluster_df) <- c("Vector Index", "Feature", "Cluster")
write.csv(cluster_df, file = '~/Desktop/Spiromics_BAL_Microbiome_Metabolome_V2/Spiromics_Colorado_Stratified/Spiromics_No_Severe_No_Healthy/Output/Heatmap/DIABLO_Pearson_Network/Features_by_Pearson_Cluster.csv', row.names = FALSE)

#plot igraph colored by greedy clustering algorithm (same plot as before, just colored by cluster)
V(ig)$community <- clusters$membership #assign points to cluster membership
rain <- rainbow(length(dupl_clusters), alpha=.5) #color scheme
V(ig)$color <- rain[V(ig)$community] #assign colors to vertices based on cluster membership

#identify which vertices are OTUs and which are metabolites
otu_vertices <- map_table$"Vertex Index"[which(str_sub(map_table$"Feature", start = 1, end = 3) == "Otu")]
met_vertices <- map_table$"Vertex Index"[which(str_sub(map_table$"Feature", start = 1, end = 3) != "Otu")]

#identify edges that cross from metabolite to OTU datasets or vice versa
cross_dataset_edges <- E(ig)[from(V(ig)[otu_vertices])]
cross_dataset_edges <- cross_dataset_edges[-which(head_of(ig, cross_dataset_edges) %in% otu_vertices)]

#color cross-dataset edges black and within-dataset edges green
E(ig)$color <- "gray" #by default all edges are colored gray
E(ig)$width <- 1
#E(ig)$lty <- "dashed"
E(ig)$color[cross_dataset_edges] <- "black"
E(ig)$width[cross_dataset_edges] <- 4

#plot final graph
pdf('~/Desktop/Spiromics_BAL_Microbiome_Metabolome_V2/Spiromics_Colorado_Stratified/Spiromics_No_Severe_No_Healthy/Output/Heatmap/DIABLO_Pearson_Network/Pearson_Clustered.pdf', width = 20, height = 20)
plot(ig, vertex.size = 7, vertex.label.cex = 2, vertex.shape = ifelse(V(ig) %in% otu_vertices, "square", "circle"), layout = l) #OTUs are squares, metabolites are circles
dev.off()

#*****************#
#Additional Graphs#
#*****************#

#collapse vertices into their greedy clusters
#c_g <- fastgreedy.community(ig)
#res_g <- simplify(contract(ig, membership(c_g)))

#set.seed(1234)
#pdf('~/Desktop/Spiromics_BAL_Microbiome_Metabolome_V2/Spiromics_Colorado_Stratified/Spiromics_No_Severe_No_Healthy/Output/Heatmap/DIABLO_Pearson_Network/Pearson_Clustered_Clusters_Only.pdf', width = 10, height = 10)
#plot(res_g) #OTUs are squares, metabolites are circles
#dev.off()

#layout for igraph
order_of_vertices <- c(unlist(clusters[1]), unlist(clusters[2]), unlist(clusters[3])) 

l <- layout_as_star(ig, center = V(ig)[8], order = order_of_vertices)
pdf('~/Desktop/Spiromics_BAL_Microbiome_Metabolome_V2/Spiromics_Colorado_Stratified/Spiromics_No_Severe_No_Healthy/Output/Heatmap/DIABLO_Pearson_Network/Pearson_Clustered_Circle.pdf')
plot(ig, vertex.size = 9, vertex.shape = ifelse(V(ig) %in% otu_vertices, "square", "circle"), layout = l) #OTUs are squares, metabolites are circles
dev.off()



