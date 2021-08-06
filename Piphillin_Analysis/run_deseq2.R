library("DESeq2")

#for DEseq2, counts have samples as columns and metadata has samples as rows
meta <- read.csv(file="~/Desktop/Spiromics_PiPhillin/Data/Master_Metadata_11_13_20_Updated.csv",row.names = 1)
counts <- read.table(file="~/Desktop/Spiromics_PiPhillin/Piphillin_Results/20201114011459__Siddharth_Madapoosi_keggPiphillin/ko_pathway_abund_table_unnorm.txt", sep = "\t", row.names = 1, check.names = FALSE, header = TRUE)
counts <- floor(counts) #floor abundances to the nearest integer

#line to ensure samples are in the same order
counts_reorder <- counts[,rownames(meta)]

#drop missing subjects
missing_subjects <- which(is.na(meta$PFV45_DERV_ctb))
if(length(missing_subjects) > 0) {
	meta <- meta[-missing_subjects, ]
	counts_reorder <- counts[, -missing_subjects]
}

#enter your R forumla on the right size of the ~
des = DESeqDataSetFromMatrix(counts_reorder, meta, ~ PFV45_DERV_ctb)

gm_mean2 = function(x, na.rm=TRUE){
  exp(sum(log(x[x > 0]), na.rm=na.rm) / length(x))
}
geoMeans = apply(counts(des), 1, gm_mean2)
des = estimateSizeFactors(des, geoMeans = geoMeans)
des = DESeq(des, fitType="local")

#adjust p-value threshold
alpha <- 1

res <- results(des)
res2 <- res[!is.na(res$padj),]


sigtab <- res2[(res2$padj < alpha), ]
write.table(file = "~/Desktop/Spiromics_Piphillin/DESeq2_Results/PFV45_DERV_ctb_DESeq2_Results.tsv", quote=FALSE, sep='\t', col.names = NA,sigtab)
