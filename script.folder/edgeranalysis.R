#careen Naitore
#7.12.2018
# This script is for EdgeR Analysis
rm(list=ls())
library(edgeR)
rm(list=ls())
## Set working directory.
setwd("/home/icipe/Documents/careenwork/objectivethree/novelmiRNA/Rawdata/")

data <- read.delim("/home/icipe/Documents/careenwork/objectivethree/novelmiRNA/Rawdata/Rawcounts_data.csv", sep = "," , row.names = 1)


#Adultmale_vs_Wildtypemale analysis

rm(list=ls())
## Set working directory.
setwd("/home/icipe/Documents/careenwork/objectivethree/novelmiRNA/Rawdata/")

data <- read.delim("/home/icipe/Documents/careenwork/objectivethree/novelmiRNA/Rawdata/Rawcounts_data.csv", sep = "," , row.names = 1)
Adultmale_vs_wildtype=(data [ ,c( 16, 17, 18, 22, 23, 24)])

rm("data")
write.csv(Adultmale_vs_wildtype, file = "/home/icipe/Documents/careenwork/objectivethree/novelmiRNA/Rawdata/Adultmale_vs_Wildtypemale.csv")
#creates a grouping factor which defines the sample types and replicates.
group <- factor(c("Adultmales", "Adultmales","Adultmales","wildtype", "wildtype", "wildtype"))

#Define the reference value as Control.
group <- relevel(group,ref = "Adultmales")

#Create a EdgeR DGE list object which stores the counts for each gene per sample. The grouping factor is used to define the samples.
y <- DGEList(counts=Adultmale_vs_wildtype[1:6],group=group)

#Generate a list with which to filter out genes with low expression values. To be included gene must have at least 1 count per million in at least 2 samples.
keep <- rowSums(cpm(y)>10) >= 2

#Use the keep list to remove genes with low expression values from the DGELIST object.
y <- y[keep, , keep.lib.sizes=FALSE]

#Calculate normalization factors for the DGELIST object.
y <- calcNormFactors(y)

#Print the samples, their respective library sizes and their nomalization factor.
y$samples

#Apply sample grouping to a data object.
design <- model.matrix(~group)

#Estimate dispersions for tags
y <- estimateDisp(y,design)

#Fit a generalized likelyhood model to the DGELIST using the Sample Grouping Data
fit <- glmFit(y,design)

#Perform likelyhood ratio test on to identify differentially expressed genes
lrt <- glmLRT(fit,coef = 2)

#Print the genes with the most significant differences.
topTags(lrt)
topGenes <- topTags(lrt, n=Inf)

#Output a summary of the top differentially expressed genes at 5% FDR
summary(de <- decideTestsDGE(lrt))

#Convert counts per million per gene to log counts per million for further downstream analyses.
logcpm <- cpm(y, prior.count = 10, log = TRUE)

#Create a data frame with final analysis data
Adultmale_vs_wildtype_DE_data <- topGenes$table
#create a data matrix with original count data 
counts <- fit$counts
#Convert data matrix to a data frame
countsdf <- as.data.frame(counts)
#Merge the analysis data with the original count data
Adultmale_vs_wildtype_data_combined <- merge(Adultmale_vs_wildtype_DE_data, countsdf, by="row.names")
#Add gene descriptions to data frame
#Extract gene subset used in the analysis
gene_names <- subset(Adultmale_vs_wildtype, keep)
#Append description data onto the data fram
Adultmale_vs_wildtype_data_combined$description <- gene_names$Description
Adultmale_vs_wildtype_data_combined$de <- de

#output the data frame to a csv file
write.csv(Adultmale_vs_wildtype_data_combined, file = "/home/icipe/Documents/careenwork/objectivethree/novelmiRNA/edgeranalysis/Adultmale_vs_wildtypemale_Analysis.csv")

#Extract the subset of upregulated genes from the gene names table
upregulated_genes <- subset(Adultmale_vs_wildtype_data_combined, PValue < 1E-3 & de == 1, select=c(Row.names, logCPM,logFC,PValue))
#Extract the subset of downregulated genes from the gene names table
downregulated_genes <- subset(Adultmale_vs_wildtype_data_combined, PValue < 1E-3 ,select=c(Row.names, logCPM,logFC, PValue))

#Plot log-fold change against log-counts per million, with differentially expressed genes highlighted
detags <- rownames(y)[as.logical(de)]
plotSmear(lrt, de.tags = detags, smooth.scatter = FALSE, pch= 16, cex= .7, lowess = FALSE, xlab = "Average log Counts per Million", ylab = "log Fold Change" , ylim=c(-10,10) , main = "log Abundance for Adultmale against Wildtypemale")
abline(h=c(-1,1),col="blue")
#Use the information in the up and downregulated genes tables to label the smear plot
#Label with Gene Identifiers
text(x=upregulated_genes$logCPM, y=upregulated_genes$logFC, labels = upregulated_genes$Row.names, cex = .5, pos=3)
text(x=downregulated_genes$logCPM, y=downregulated_genes$logFC, labels = downregulated_genes$Row.names, cex = .5, pos=1)

#Label with Gene Descriptions
text(x=upregulated_genes$logCPM, y=upregulated_genes$logFC, labels = upregulated_genes$Description, cex = .5, pos=3)
text(x=downregulated_genes$logCPM, y=downregulated_genes$logFC, labels = downregulated_genes$Description, cex = .5, pos=1)

#Larvae_vs_pupae analysis


rm(list=ls())
## Set working directory.
setwd("/home/icipe/Documents/careenwork/objectivethree/novelmiRNA/Rawdata/")

data <- read.delim("/home/icipe/Documents/careenwork/objectivethree/novelmiRNA/Rawdata/Rawcounts_data.csv", sep = "," , row.names = 1)

Larvae_vs_pupae=(data [ ,c( 1, 2, 4, 5, 6)])
rm("data")
write.csv(Larvae_vs_pupae, file = "/home/icipe/Documents/careenwork/objectivethree/novelmiRNA/Rawdata/Larvae_vs_pupae.csv")
#creates a grouping factor which defines the sample types and replicates.
group <- factor(c("Larvae", "Larvae" ,"Pupae", "Pupae", "Pupae"))

#Define the reference value as Control.
group <- relevel(group,ref = "Larvae")

#Create a EdgeR DGE list object which stores the counts for each gene per sample. The grouping factor is used to define the samples.
y <- DGEList(counts=Larvae_vs_pupae[1:5],group=group)

#Generate a list with which to filter out genes with low expression values. To be included gene must have at least 1 count per million in at least 2 samples.
keep <- rowSums(cpm(y)>10) >= 2

#Use the keep list to remove genes with low expression values from the DGELIST object.
y <- y[keep, , keep.lib.sizes=FALSE]

#Calculate normalization factors for the DGELIST object.
y <- calcNormFactors(y)

#Print the samples, their respective library sizes and their nomalization factor.
y$samples

#Apply sample grouping to a data object.
design <- model.matrix(~group)

#Estimate dispersions for tags
y <- estimateDisp(y,design)

#Fit a generalized likelyhood model to the DGELIST using the Sample Grouping Data
fit <- glmFit(y,design)

#Perform likelyhood ratio test on to identify differentially expressed genes
lrt <- glmLRT(fit,coef = 2)

#Print the genes with the most significant differences.
topTags(lrt)
topGenes <- topTags(lrt, n=Inf)

#Output a summary of the top differentially expressed genes at 5% FDR
summary(de <- decideTestsDGE(lrt))

#Convert counts per million per gene to log counts per million for further downstream analyses.
logcpm <- cpm(y, prior.count = 10, log = TRUE)

#Create a data frame with final analysis data
Larvae_vs_pupae_DE_data <- topGenes$table
#create a data matrix with original count data 
counts <- fit$counts
#Convert data matrix to a data frame
countsdf <- as.data.frame(counts)
#Merge the analysis data with the original count data
Larvae_vs_pupae_data_combined <- merge(Larvae_vs_pupae_DE_data, countsdf, by="row.names")
#Add gene descriptions to data frame
#Extract gene subset used in the analysis
gene_names <- subset(Larvae_vs_pupae, keep)
#Append description data onto the data fram
Larvae_vs_pupae_data_combined$description <- gene_names$Description
Larvae_vs_pupae_data_combined$de <- de

#output the data frame to a csv file
write.csv(Larvae_vs_pupae_data_combined, file = "/home/icipe/Documents/careenwork/objectivethree/novelmiRNA/edgeranalysis/Larvae_vs_pupae_Analysis.csv")

#Extract the subset of upregulated genes from the gene names table
upregulated_genes <- subset(Larvae_vs_pupae_data_combined, PValue < 1E-3 & de == 1 , select=c(Row.names, logCPM,logFC,PValue))
#Extract the subset of downregulated genes from the gene names table
downregulated_genes <- subset(Larvae_vs_pupae_data_combined, PValue < 1E-3  ,select=c(Row.names, logCPM,logFC, PValue))

#Plot log-fold change against log-counts per million, with differentially expressed genes highlighted
detags <- rownames(y)[as.logical(de)]
plotSmear(lrt, de.tags = detags, smooth.scatter = FALSE, pch= 16, cex= .7, lowess = FALSE, xlab = "Average log Counts per Million", ylab = "log Fold Change" , ylim=c(-10,10) , main = "log Abundance for Larvae against pupae")
abline(h=c(-1,1),col="blue")
#Use the information in the up and downregulated genes tables to label the smear plot
#Label with Gene Identifiers
text(x=upregulated_genes$logCPM, y=upregulated_genes$logFC, labels = upregulated_genes$Row.names, cex = .5, pos=3)
text(x=downregulated_genes$logCPM, y=downregulated_genes$logFC, labels = downregulated_genes$Row.names, cex = .5, pos=1)

#Label with Gene Descriptions
text(x=upregulated_genes$logCPM, y=upregulated_genes$logFC, labels = upregulated_genes$Description, cex = .5, pos=3)
text(x=downregulated_genes$logCPM, y=downregulated_genes$logFC, labels = downregulated_genes$Description, cex = .5, pos=1)

#Pupa_vs_teneralfemale analysis

setwd("/home/icipe/Documents/careenwork/objectivethree/novelmiRNA/Rawdata/")

data <- read.delim("/home/icipe/Documents/careenwork/objectivethree/novelmiRNA/Rawdata/Rawcounts_data.csv", sep = "," , row.names = 1)
Pupae_vs_teneralfemale=(data [ ,c( 4, 5, 6, 7, 8, 9)])

rm("data")
write.csv(Pupae_vs_teneralfemale, file = "/home/icipe/Documents/careenwork/objectivethree/novelmiRNA /Rawdata/ Pupae_vs_teneralfemale.csv")
#creates a grouping factor which defines the sample types and replicates.
group <- factor(c("Adultmales", "Adultmales","Adultmales","wildtype", "wildtype", "wildtype"))

#Define the reference value as Control.
group <- relevel(group,ref = "Adultmales")

#Create a EdgeR DGE list object which stores the counts for each gene per sample. The grouping factor is used to define the samples.
y <- DGEList(counts=Pupae_vs_teneralfemale[1:6],group=group)

#Generate a list with which to filter out genes with low expression values. To be included gene must have at least 1 count per million in at least 2 samples.
keep <- rowSums(cpm(y)>10) >= 2

#Use the keep list to remove genes with low expression values from the DGELIST object.
y <- y[keep, , keep.lib.sizes=FALSE]

#Calculate normalization factors for the DGELIST object.
y <- calcNormFactors(y)

#Print the samples, their respective library sizes and their nomalization factor.
y$samples

#Apply sample grouping to a data object.
design <- model.matrix(~group)

#Estimate dispersions for tags
y <- estimateDisp(y,design)

#Fit a generalized likelyhood model to the DGELIST using the Sample Grouping Data
fit <- glmFit(y,design)

#Perform likelyhood ratio test on to identify differentially expressed genes
lrt <- glmLRT(fit,coef = 2)

#Print the genes with the most significant differences.
topTags(lrt)
topGenes <- topTags(lrt, n=Inf)

#Output a summary of the top differentially expressed genes at 5% FDR
summary(de <- decideTestsDGE(lrt))

#Convert counts per million per gene to log counts per million for further downstream analyses.
logcpm <- cpm(y, prior.count = 10, log = TRUE)

#Create a data frame with final analysis data
Adultmale_vs_wildtype_DE_data <- topGenes$table
#create a data matrix with original count data 
counts <- fit$counts
#Convert data matrix to a data frame
countsdf <- as.data.frame(counts)
#Merge the analysis data with the original count data
Adultmale_vs_wildtype_data_combined <- merge(Adultmale_vs_wildtype_DE_data, countsdf, by="row.names")
#Add gene descriptions to data frame
#Extract gene subset used in the analysis
gene_names <- subset(Pupae_vs_teneralfemale, keep)
#Append description data onto the data fram
Adultmale_vs_wildtype_data_combined$description <- gene_names$Description
Adultmale_vs_wildtype_data_combined$de <- de

#output the data frame to a csv file
write.csv(Adultmale_vs_wildtype_data_combined, file = "/home/icipe/Documents/careenwork/objectivethree/novelmiRNA/edgeranalysis/Pupae_vs_teneralfemale_Analysis.csv")

#Extract the subset of upregulated genes from the gene names table
upregulated_genes <- subset(Adultmale_vs_wildtype_data_combined, PValue < 1E-3 & de == 1, select=c(Row.names, logCPM,logFC,PValue))
#Extract the subset of downregulated genes from the gene names table
downregulated_genes <- subset(Adultmale_vs_wildtype_data_combined, PValue < 1E-3 & de == -1,select=c(Row.names, logCPM,logFC, PValue))

#Plot log-fold change against log-counts per million, with differentially expressed genes highlighted
detags <- rownames(y)[as.logical(de)]
plotSmear(lrt, de.tags = detags, smooth.scatter = FALSE, pch= 16, cex= .7, lowess = FALSE, xlab = "Average log Counts per Million", ylab = "log Fold Change" , ylim=c(-10,10) , main = "log Abundance for Pupae against Teneralfemale")
abline(h=c(-1,1),col="blue")
#Use the information in the up and downregulated genes tables to label the smear plot
#Label with Gene Identifiers
text(x=upregulated_genes$logCPM, y=upregulated_genes$logFC, labels = upregulated_genes$Row.names, cex = .5, pos=3)
text(x=downregulated_genes$logCPM, y=downregulated_genes$logFC, labels = downregulated_genes$Row.names, cex = .5, pos=1)

#Label with Gene Descriptions
text(x=upregulated_genes$logCPM, y=upregulated_genes$logFC, labels = upregulated_genes$Description, cex = .5, pos=3)
text(x=downregulated_genes$logCPM, y=downregulated_genes$logFC, labels = downregulated_genes$Description, cex = .5, pos=1)



#pupae_vs_teneralmale analysis

rm(list=ls())
## Set working directory.
setwd("/home/icipe/Documents/careenwork/objectivethree/novelmiRNA/Rawdata/")

data <- read.delim("/home/icipe/Documents/careenwork/objectivethree/novelmiRNA/Rawdata/Rawcounts_data.csv", sep = "," , row.names = 1)
Pupae_vs_teneralmale=(data [ ,c( 4, 5, 6, 10, 11, 12)])

rm("data")
write.csv(Pupae_vs_teneralmale, file = "/home/icipe/Documents/careenwork/objectivethree/novelmiRNA /Rawdata/pupae_vs_teneralmale.csv")
#creates a grouping factor which defines the sample types and replicates.
group <- factor(c("Adultmales", "Adultmales","Adultmales","wildtype", "wildtype", "wildtype"))

#Define the reference value as Control.
group <- relevel(group,ref = "Adultmales")

#Create a EdgeR DGE list object which stores the counts for each gene per sample. The grouping factor is used to define the samples.
y <- DGEList(counts=Pupae_vs_teneralmale[1:6],group=group)

#Generate a list with which to filter out genes with low expression values. To be included gene must have at least 1 count per million in at least 2 samples.
keep <- rowSums(cpm(y)>10) >= 2

#Use the keep list to remove genes with low expression values from the DGELIST object.
y <- y[keep, , keep.lib.sizes=FALSE]

#Calculate normalization factors for the DGELIST object.
y <- calcNormFactors(y)

#Print the samples, their respective library sizes and their nomalization factor.
y$samples

#Apply sample grouping to a data object.
design <- model.matrix(~group)

#Estimate dispersions for tags
y <- estimateDisp(y,design)

#Fit a generalized likelyhood model to the DGELIST using the Sample Grouping Data
fit <- glmFit(y,design)

#Perform likelyhood ratio test on to identify differentially expressed genes
lrt <- glmLRT(fit,coef = 2)

#Print the genes with the most significant differences.
topTags(lrt)
topGenes <- topTags(lrt, n=Inf)

#Output a summary of the top differentially expressed genes at 5% FDR
summary(de <- decideTestsDGE(lrt))

#Convert counts per million per gene to log counts per million for further downstream analyses.
logcpm <- cpm(y, prior.count = 10, log = TRUE)

#Create a data frame with final analysis data
Adultmale_vs_wildtype_DE_data <- topGenes$table
#create a data matrix with original count data 
counts <- fit$counts
#Convert data matrix to a data frame
countsdf <- as.data.frame(counts)
#Merge the analysis data with the original count data
Adultmale_vs_wildtype_data_combined <- merge(Adultmale_vs_wildtype_DE_data, countsdf, by="row.names")
#Add gene descriptions to data frame
#Extract gene subset used in the analysis
gene_names <- subset(Pupae_vs_teneralmale, keep)
#Append description data onto the data fram
Adultmale_vs_wildtype_data_combined$description <- gene_names$Description
Adultmale_vs_wildtype_data_combined$de <- de

#output the data frame to a csv file
write.csv(Adultmale_vs_wildtype_data_combined, file = "/home/icipe/Documents/careenwork/objectivethree/novelmiRNA/edgeranalysis/Pupae_vs_teneralmale_Analysis.csv")

#Extract the subset of upregulated genes from the gene names table
upregulated_genes <- subset(Adultmale_vs_wildtype_data_combined, PValue < 1E-3 & de == 1, select=c(Row.names, logCPM,logFC,PValue))
#Extract the subset of downregulated genes from the gene names table
downregulated_genes <- subset(Adultmale_vs_wildtype_data_combined, PValue < 1E-3 & de == -1,select=c(Row.names, logCPM,logFC, PValue))

#Plot log-fold change against log-counts per million, with differentially expressed genes highlighted
detags <- rownames(y)[as.logical(de)]
plotSmear(lrt, de.tags = detags, smooth.scatter = FALSE, pch= 16, cex= .7, lowess = FALSE, xlab = "Average log Counts per Million", ylab = "log Fold Change" , ylim=c(-10,10) , main = "log Abundance for Pupae against teneralmale")
abline(h=c(-1,1),col="blue")
#Use the information in the up and downregulated genes tables to label the smear plot
#Label with Gene Identifiers
text(x=upregulated_genes$logCPM, y=upregulated_genes$logFC, labels = upregulated_genes$Row.names, cex = .5, pos=3)
text(x=downregulated_genes$logCPM, y=downregulated_genes$logFC, labels = downregulated_genes$Row.names, cex = .5, pos=1)

#Label with Gene Descriptions
text(x=upregulated_genes$logCPM, y=upregulated_genes$logFC, labels = upregulated_genes$Description, cex = .5, pos=3)
text(x=downregulated_genes$logCPM, y=downregulated_genes$logFC, labels = downregulated_genes$Description, cex = .5, pos=1)



#Adultfemale_vs_teneralfemale analysis

rm(list=ls())
## Set working directory.
setwd("/home/icipe/Documents/careenwork/objectivethree/novelmiRNA/Rawdata/")

data <- read.delim("/home/icipe/Documents/careenwork/objectivethree/novelmiRNA/Rawdata/Rawcounts_data.csv", sep = "," , row.names = 1)
teneralfemale_vs_Adultfemale=(data [ ,c( 7, 8, 9, 13,14, 15)])

rm("data")
write.csv(teneralfemale_vs_Adultfemale, file = "/home/icipe/Documents/careenwork/objectivethree/novelmiRNA/Rawdata/teneralfemale_vs_Adultfemale.csv")
#creates a grouping factor which defines the sample types and replicates.
group <- factor(c("Adultmales", "Adultmales","Adultmales","wildtype", "wildtype", "wildtype"))

#Define the reference value as Control.
group <- relevel(group,ref = "Adultmales")

#Create a EdgeR DGE list object which stores the counts for each gene per sample. The grouping factor is used to define the samples.
y <- DGEList(counts=teneralfemale_vs_Adultfemale[1:6],group=group)

#Generate a list with which to filter out genes with low expression values. To be included gene must have at least 1 count per million in at least 2 samples.
keep <- rowSums(cpm(y)>10) >= 2

#Use the keep list to remove genes with low expression values from the DGELIST object.
y <- y[keep, , keep.lib.sizes=FALSE]

#Calculate normalization factors for the DGELIST object.
y <- calcNormFactors(y)

#Print the samples, their respective library sizes and their nomalization factor.
y$samples

#Apply sample grouping to a data object.
design <- model.matrix(~group)

#Estimate dispersions for tags
y <- estimateDisp(y,design)

#Fit a generalized likelyhood model to the DGELIST using the Sample Grouping Data
fit <- glmFit(y,design)

#Perform likelyhood ratio test on to identify differentially expressed genes
lrt <- glmLRT(fit,coef = 2)

#Print the genes with the most significant differences.
topTags(lrt)
topGenes <- topTags(lrt, n=Inf)

#Output a summary of the top differentially expressed genes at 5% FDR
summary(de <- decideTestsDGE(lrt))

#Convert counts per million per gene to log counts per million for further downstream analyses.
logcpm <- cpm(y, prior.count = 10, log = TRUE)

#Create a data frame with final analysis data
Adultmale_vs_wildtype_DE_data <- topGenes$table
#create a data matrix with original count data 
counts <- fit$counts
#Convert data matrix to a data frame
countsdf <- as.data.frame(counts)
#Merge the analysis data with the original count data
Adultmale_vs_wildtype_data_combined <- merge(Adultmale_vs_wildtype_DE_data, countsdf, by="row.names")
#Add gene descriptions to data frame
#Extract gene subset used in the analysis
gene_names <- subset(teneralfemale_vs_Adultfemale, keep)
#Append description data onto the data fram
Adultmale_vs_wildtype_data_combined$description <- gene_names$Description
Adultmale_vs_wildtype_data_combined$de <- de

#output the data frame to a csv file
write.csv(Adultmale_vs_wildtype_data_combined, file = "/home/icipe/Documents/careenwork/objectivethree/novelmiRNA/edgeranalysis/teneralfemale_vs_Adultfemale_Analysis.csv")

#Extract the subset of upregulated genes from the gene names table
upregulated_genes <- subset(Adultmale_vs_wildtype_data_combined, PValue < 1E-3 & de == 1, select=c(Row.names, logCPM,logFC,PValue))
#Extract the subset of downregulated genes from the gene names table
downregulated_genes <- subset(Adultmale_vs_wildtype_data_combined, PValue < 1E-3 & de == -1,select=c(Row.names, logCPM,logFC, PValue))

#Plot log-fold change against log-counts per million, with differentially expressed genes highlighted
detags <- rownames(y)[as.logical(de)]
plotSmear(lrt, de.tags = detags, smooth.scatter = FALSE, pch= 16, cex= .7, lowess = FALSE, xlab = "Average log Counts per Million", ylab = "log Fold Change" , ylim=c(-10,10) , main = "log Abundance for Teneralfemale against Adultfemale")
abline(h=c(-1,1),col="blue")
#Use the information in the up and downregulated genes tables to label the smear plot
#Label with Gene Identifiers
text(x=upregulated_genes$logCPM, y=upregulated_genes$logFC, labels = upregulated_genes$Row.names, cex = .5, pos=3)
text(x=downregulated_genes$logCPM, y=downregulated_genes$logFC, labels = downregulated_genes$Row.names, cex = .5, pos=1)

#Label with Gene Descriptions
text(x=upregulated_genes$logCPM, y=upregulated_genes$logFC, labels = upregulated_genes$Description, cex = .5, pos=3)
text(x=downregulated_genes$logCPM, y=downregulated_genes$logFC, labels = downregulated_genes$Description, cex = .5, pos=1)



#Adultmale_vs_teneralmale analysis

rm(list=ls())
## Set working directory.
setwd("/home/icipe/Documents/careenwork/objectivethree/novelmiRNA/Rawdata/")

data <- read.delim("/home/icipe/Documents/careenwork/objectivethree/novelmiRNA/Rawdata/Rawcounts_data.csv", sep = "," , row.names = 1)
teneralmale_vs_Adultmale=(data [ ,c( 10, 11, 12, 16, 17, 18)])

rm("data")
write.csv(teneralmale_vs_Adultmale, file = "/home/icipe/Documents/careenwork/objectivethree/novelmiRNA/Rawdata/teneralmale_vs_Adultmale.csv")
#creates a grouping factor which defines the sample types and replicates.
group <- factor(c("Adultmales", "Adultmales","Adultmales","wildtype", "wildtype", "wildtype"))

#Define the reference value as Control.
group <- relevel(group,ref = "Adultmales")

#Create a EdgeR DGE list object which stores the counts for each gene per sample. The grouping factor is used to define the samples.
y <- DGEList(counts=teneralmale_vs_Adultmale[1:6],group=group)

#Generate a list with which to filter out genes with low expression values. To be included gene must have at least 1 count per million in at least 2 samples.
keep <- rowSums(cpm(y)>10) >= 2

#Use the keep list to remove genes with low expression values from the DGELIST object.
y <- y[keep, , keep.lib.sizes=FALSE]

#Calculate normalization factors for the DGELIST object.
y <- calcNormFactors(y)

#Print the samples, their respective library sizes and their nomalization factor.
y$samples

#Apply sample grouping to a data object.
design <- model.matrix(~group)

#Estimate dispersions for tags
y <- estimateDisp(y,design)

#Fit a generalized likelyhood model to the DGELIST using the Sample Grouping Data
fit <- glmFit(y,design)

#Perform likelyhood ratio test on to identify differentially expressed genes
lrt <- glmLRT(fit,coef = 2)

#Print the genes with the most significant differences.
topTags(lrt)
topGenes <- topTags(lrt, n=Inf)

#Output a summary of the top differentially expressed genes at 5% FDR
summary(de <- decideTestsDGE(lrt))

#Convert counts per million per gene to log counts per million for further downstream analyses.
logcpm <- cpm(y, prior.count = 10, log = TRUE)

#Create a data frame with final analysis data
Adultmale_vs_wildtype_DE_data <- topGenes$table
#create a data matrix with original count data 
counts <- fit$counts
#Convert data matrix to a data frame
countsdf <- as.data.frame(counts)
#Merge the analysis data with the original count data
Adultmale_vs_wildtype_data_combined <- merge(Adultmale_vs_wildtype_DE_data, countsdf, by="row.names")
#Add gene descriptions to data frame
#Extract gene subset used in the analysis
gene_names <- subset(teneralmale_vs_Adultmale, keep)
#Append description data onto the data fram
Adultmale_vs_wildtype_data_combined$description <- gene_names$Description
Adultmale_vs_wildtype_data_combined$de <- de

#output the data frame to a csv file
write.csv(Adultmale_vs_wildtype_data_combined, file = "/home/icipe/Documents/careenwork/objectivethree/novelmiRNA/edgeranalysis/teneralmale_vs_Adultmale_Analysis.csv")

#Extract the subset of upregulated genes from the gene names table
upregulated_genes <- subset(Adultmale_vs_wildtype_data_combined, PValue < 1E-3 & de == 1, select=c(Row.names, logCPM,logFC,PValue))
#Extract the subset of downregulated genes from the gene names table
downregulated_genes <- subset(Adultmale_vs_wildtype_data_combined, PValue < 1E-3 ,select=c(Row.names, logCPM,logFC, PValue))

#Plot log-fold change against log-counts per million, with differentially expressed genes highlighted
detags <- rownames(y)[as.logical(de)]
plotSmear(lrt, de.tags = detags, smooth.scatter = FALSE, pch= 16, cex= .7, lowess = FALSE, xlab = "Average log Counts per Million", ylab = "log Fold Change" , ylim=c(-10,10) , main = "log Abundance for Teneralmale against Adultmale")
abline(h=c(-1,1),col="blue")
#Use the information in the up and downregulated genes tables to label the smear plot
#Label with Gene Identifiers
text(x=upregulated_genes$logCPM, y=upregulated_genes$logFC, labels = upregulated_genes$Row.names, cex = .5, pos=3)
text(x=downregulated_genes$logCPM, y=downregulated_genes$logFC, labels = downregulated_genes$Row.names, cex = .5, pos=1)

#Label with Gene Descriptions
text(x=upregulated_genes$logCPM, y=upregulated_genes$logFC, labels = upregulated_genes$Description, cex = .5, pos=3)
text(x=downregulated_genes$logCPM, y=downregulated_genes$logFC, labels = downregulated_genes$Description, cex = .5, pos=1)



#teneralfemale_vs_teneralmale analysis

rm(list=ls())
## Set working directory.
setwd("/home/icipe/Documents/careenwork/objectivethree/novelmiRNA/Rawdata/")

data <- read.delim("/home/icipe/Documents/careenwork/objectivethree/novelmiRNA/Rawdata/Rawcounts_data.csv", sep = "," , row.names = 1)
teneralfemale_vs_teneralmale=(data [ ,c( 7, 8, 9, 10, 11, 12)])

rm("data")
write.csv(teneralfemale_vs_teneralmale, file = "/home/icipe/Documents/careenwork/objectivethree/novelmiRNA/Rawdata/teneralfemale_vs_teneralmale.csv")
#creates a grouping factor which defines the sample types and replicates.
group <- factor(c("Adultmales", "Adultmales","Adultmales","wildtype", "wildtype", "wildtype"))

#Define the reference value as Control.
group <- relevel(group,ref = "Adultmales")

#Create a EdgeR DGE list object which stores the counts for each gene per sample. The grouping factor is used to define the samples.
y <- DGEList(counts=teneralfemale_vs_teneralmale[1:6],group=group)

#Generate a list with which to filter out genes with low expression values. To be included gene must have at least 1 count per million in at least 2 samples.
keep <- rowSums(cpm(y)>10) >= 2

#Use the keep list to remove genes with low expression values from the DGELIST object.
y <- y[keep, , keep.lib.sizes=FALSE]

#Calculate normalization factors for the DGELIST object.
y <- calcNormFactors(y)

#Print the samples, their respective library sizes and their nomalization factor.
y$samples

#Apply sample grouping to a data object.
design <- model.matrix(~group)

#Estimate dispersions for tags
y <- estimateDisp(y,design)

#Fit a generalized likelyhood model to the DGELIST using the Sample Grouping Data
fit <- glmFit(y,design)

#Perform likelyhood ratio test on to identify differentially expressed genes
lrt <- glmLRT(fit,coef = 2)

#Print the genes with the most significant differences.
topTags(lrt)
topGenes <- topTags(lrt, n=Inf)

#Output a summary of the top differentially expressed genes at 5% FDR
summary(de <- decideTestsDGE(lrt))

#Convert counts per million per gene to log counts per million for further downstream analyses.
logcpm <- cpm(y, prior.count = 10, log = TRUE)

#Create a data frame with final analysis data
Adultmale_vs_wildtype_DE_data <- topGenes$table
#create a data matrix with original count data 
counts <- fit$counts
#Convert data matrix to a data frame
countsdf <- as.data.frame(counts)
#Merge the analysis data with the original count data
Adultmale_vs_wildtype_data_combined <- merge(Adultmale_vs_wildtype_DE_data, countsdf, by="row.names")
#Add gene descriptions to data frame
#Extract gene subset used in the analysis
gene_names <- subset(teneralfemale_vs_teneralmale, keep)
#Append description data onto the data fram
Adultmale_vs_wildtype_data_combined$description <- gene_names$Description
Adultmale_vs_wildtype_data_combined$de <- de

#output the data frame to a csv file
write.csv(Adultmale_vs_wildtype_data_combined, file = "/home/icipe/Documents/careenwork/objectivethree/novelmiRNA/edgeranalysis/teneralmale_vs_teneralfemale_Analysis.csv")

#Extract the subset of upregulated genes from the gene names table
upregulated_genes <- subset(Adultmale_vs_wildtype_data_combined, PValue < 1E-3 & de == 1, select=c(Row.names, logCPM,logFC,PValue))
#Extract the subset of downregulated genes from the gene names table
downregulated_genes <- subset(Adultmale_vs_wildtype_data_combined, PValue < 1E-3 ,select=c(Row.names, logCPM,logFC, PValue))

#Plot log-fold change against log-counts per million, with differentially expressed genes highlighted
detags <- rownames(y)[as.logical(de)]
plotSmear(lrt, de.tags = detags, smooth.scatter = FALSE, pch= 16, cex= .7, lowess = FALSE, xlab = "Average log Counts per Million", ylab = "log Fold Change" , ylim=c(-10,10) , main = "log Abundance for Teneralfemale against Teneralmale")
abline(h=c(-1,1),col="blue")
#Use the information in the up and downregulated genes tables to label the smear plot
#Label with Gene Identifiers
text(x=upregulated_genes$logCPM, y=upregulated_genes$logFC, labels = upregulated_genes$Row.names, cex = .5, pos=3)
text(x=downregulated_genes$logCPM, y=downregulated_genes$logFC, labels = downregulated_genes$Row.names, cex = .5, pos=1)

#Label with Gene Descriptions
text(x=upregulated_genes$logCPM, y=upregulated_genes$logFC, labels = upregulated_genes$Description, cex = .5, pos=3)
text(x=downregulated_genes$logCPM, y=downregulated_genes$logFC, labels = downregulated_genes$Description, cex = .5, pos=1)




#Gravidfemale_vs_Non-gravidfemale analysis

rm(list=ls())
## Set working directory.
setwd("/home/icipe/Documents/careenwork/objectivethree/novelmiRNA/Rawdata/")

data <- read.delim("/home/icipe/Documents/careenwork/objectivethree/novelmiRNA/Rawdata/Rawcounts_data.csv", sep = "," , row.names = 1)
nongravid_vs_gravid=(data [ ,c( 7, 8, 9, 19, 20, 21)])

rm("data")
write.csv(nongravid_vs_gravid, file = "/home/icipe/Documents/careenwork/objectivethree/novelmiRNA/Rawdata/nongravid_vs_gravid.csv")
#creates a grouping factor which defines the sample types and replicates.
group <- factor(c("Adultmales", "Adultmales","Adultmales","wildtype", "wildtype", "wildtype"))

#Define the reference value as Control.
group <- relevel(group,ref = "Adultmales")

#Create a EdgeR DGE list object which stores the counts for each gene per sample. The grouping factor is used to define the samples.
y <- DGEList(counts=nongravid_vs_gravid[1:6],group=group)

#Generate a list with which to filter out genes with low expression values. To be included gene must have at least 1 count per million in at least 2 samples.
keep <- rowSums(cpm(y)>10) >= 2

#Use the keep list to remove genes with low expression values from the DGELIST object.
y <- y[keep, , keep.lib.sizes=FALSE]

#Calculate normalization factors for the DGELIST object.
y <- calcNormFactors(y)

#Print the samples, their respective library sizes and their nomalization factor.
y$samples

#Apply sample grouping to a data object.
design <- model.matrix(~group)

#Estimate dispersions for tags
y <- estimateDisp(y,design)

#Fit a generalized likelyhood model to the DGELIST using the Sample Grouping Data
fit <- glmFit(y,design)

#Perform likelyhood ratio test on to identify differentially expressed genes
lrt <- glmLRT(fit,coef = 2)

#Print the genes with the most significant differences.
topTags(lrt)
topGenes <- topTags(lrt, n=Inf)

#Output a summary of the top differentially expressed genes at 5% FDR
summary(de <- decideTestsDGE(lrt))

#Convert counts per million per gene to log counts per million for further downstream analyses.
logcpm <- cpm(y, prior.count = 10, log = TRUE)

#Create a data frame with final analysis data
Adultmale_vs_wildtype_DE_data <- topGenes$table
#create a data matrix with original count data 
counts <- fit$counts
#Convert data matrix to a data frame
countsdf <- as.data.frame(counts)
#Merge the analysis data with the original count data
Adultmale_vs_wildtype_data_combined <- merge(Adultmale_vs_wildtype_DE_data, countsdf, by="row.names")
#Add gene descriptions to data frame
#Extract gene subset used in the analysis
gene_names <- subset(nongravid_vs_gravid, keep)
#Append description data onto the data fram
Adultmale_vs_wildtype_data_combined$description <- gene_names$Description
Adultmale_vs_wildtype_data_combined$de <- de

#output the data frame to a csv file
write.csv(Adultmale_vs_wildtype_data_combined, file = "/home/icipe/Documents/careenwork/objectivethree/novelmiRNA/edgeranalysis/nongravid_vs_gravid.csv")

#Extract the subset of upregulated genes from the gene names table
upregulated_genes <- subset(Adultmale_vs_wildtype_data_combined, PValue < 1E-3 & de == 1, select=c(Row.names, logCPM,logFC,PValue))
#Extract the subset of downregulated genes from the gene names table
downregulated_genes <- subset(Adultmale_vs_wildtype_data_combined, PValue < 1E-3 & de == -1,select=c(Row.names, logCPM,logFC, PValue))

#Plot log-fold change against log-counts per million, with differentially expressed genes highlighted
detags <- rownames(y)[as.logical(de)]
plotSmear(lrt, de.tags = detags, smooth.scatter = FALSE, pch= 16, cex= .7, lowess = FALSE, xlab = "Average log Counts per Million", ylab = "log Fold Change" , ylim=c(-10,10) , main = "log Abundance for Gravid against Non-gravid female")
abline(h=c(-1,1),col="blue")
#Use the information in the up and downregulated genes tables to label the smear plot
#Label with Gene Identifiers
text(x=upregulated_genes$logCPM, y=upregulated_genes$logFC, labels = upregulated_genes$Row.names, cex = .5, pos=3)
text(x=downregulated_genes$logCPM, y=downregulated_genes$logFC, labels = downregulated_genes$Row.names, cex = .5, pos=1)

#Label with Gene Descriptions
text(x=upregulated_genes$logCPM, y=upregulated_genes$logFC, labels = upregulated_genes$Description, cex = .5, pos=3)
text(x=downregulated_genes$logCPM, y=downregulated_genes$logFC, labels = downregulated_genes$Description, cex = .5, pos=1)





#Adultfemale_vs_Adultmale analysis

rm(list=ls())
## Set working directory.
setwd("/home/icipe/Documents/careenwork/objectivethree/novelmiRNA/Rawdata/")

data <- read.delim("/home/icipe/Documents/careenwork/objectivethree/novelmiRNA/Rawdata/Rawcounts_data.csv", sep = "," , row.names = 1)
Adultfemale_vs_Adultmale=(data [ ,c( 13, 14, 15, 16, 17, 18)])

rm("data")
write.csv(Adultfemale_vs_Adultmale, file = "/home/icipe/Documents/careenwork/objectivethree/novelmiRNA/Rawdata/Adultfemale_vs_Adultmale.csv")
#creates a grouping factor which defines the sample types and replicates.
group <- factor(c("Adultmales", "Adultmales","Adultmales","wildtype", "wildtype", "wildtype"))

#Define the reference value as Control.
group <- relevel(group,ref = "Adultmales")

#Create a EdgeR DGE list object which stores the counts for each gene per sample. The grouping factor is used to define the samples.
y <- DGEList(counts=Adultfemale_vs_Adultmale[1:6],group=group)

#Generate a list with which to filter out genes with low expression values. To be included gene must have at least 1 count per million in at least 2 samples.
keep <- rowSums(cpm(y)>10) >= 2

#Use the keep list to remove genes with low expression values from the DGELIST object.
y <- y[keep, , keep.lib.sizes=FALSE]

#Calculate normalization factors for the DGELIST object.
y <- calcNormFactors(y)

#Print the samples, their respective library sizes and their nomalization factor.
y$samples

#Apply sample grouping to a data object.
design <- model.matrix(~group)

#Estimate dispersions for tags
y <- estimateDisp(y,design)

#Fit a generalized likelyhood model to the DGELIST using the Sample Grouping Data
fit <- glmFit(y,design)

#Perform likelyhood ratio test on to identify differentially expressed genes
lrt <- glmLRT(fit,coef = 2)

#Print the genes with the most significant differences.
topTags(lrt)
topGenes <- topTags(lrt, n=Inf)

#Output a summary of the top differentially expressed genes at 5% FDR
summary(de <- decideTestsDGE(lrt))

#Convert counts per million per gene to log counts per million for further downstream analyses.
logcpm <- cpm(y, prior.count = 10, log = TRUE)

#Create a data frame with final analysis data
Adultmale_vs_wildtype_DE_data <- topGenes$table
#create a data matrix with original count data 
counts <- fit$counts
#Convert data matrix to a data frame
countsdf <- as.data.frame(counts)
#Merge the analysis data with the original count data
Adultmale_vs_wildtype_data_combined <- merge(Adultmale_vs_wildtype_DE_data, countsdf, by="row.names")
#Add gene descriptions to data frame
#Extract gene subset used in the analysis
gene_names <- subset(Adultfemale_vs_Adultmale, keep)
#Append description data onto the data fram
Adultmale_vs_wildtype_data_combined$description <- gene_names$Description
Adultmale_vs_wildtype_data_combined$de <- de

#output the data frame to a csv file
write.csv(Adultmale_vs_wildtype_data_combined, file = "/home/icipe/Documents/careenwork/objectivethree/novelmiRNA/edgeranalysis/Adultfemale_vs_Adultmale.csv")

#Extract the subset of upregulated genes from the gene names table
upregulated_genes <- subset(Adultmale_vs_wildtype_data_combined, PValue < 1E-3 & de == 1, select=c(Row.names, logCPM,logFC,PValue))
#Extract the subset of downregulated genes from the gene names table
downregulated_genes <- subset(Adultmale_vs_wildtype_data_combined, PValue < 1E-3 & de == -1,select=c(Row.names, logCPM,logFC, PValue))

#Plot log-fold change against log-counts per million, with differentially expressed genes highlighted
detags <- rownames(y)[as.logical(de)]
plotSmear(lrt, de.tags = detags, smooth.scatter = FALSE, pch= 16, cex= .7, lowess = FALSE, xlab = "Average log Counts per Million", ylab = "log Fold Change" , ylim=c(-10,10) , main = "log Abundance for Adultmale against Adult-female")
abline(h=c(-1,1),col="blue")
#Use the information in the up and downregulated genes tables to label the smear plot
#Label with Gene Identifiers
text(x=upregulated_genes$logCPM, y=upregulated_genes$logFC, labels = upregulated_genes$Row.names, cex = .5, pos=3)
text(x=downregulated_genes$logCPM, y=downregulated_genes$logFC, labels = downregulated_genes$Row.names, cex = .5, pos=1)

#Label with Gene Descriptions
text(x=upregulated_genes$logCPM, y=upregulated_genes$logFC, labels = upregulated_genes$Description, cex = .5, pos=3)
text(x=downregulated_genes$logCPM, y=downregulated_genes$logFC, labels = downregulated_genes$Description, cex = .5, pos=1)




