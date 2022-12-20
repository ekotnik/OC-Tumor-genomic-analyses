## RNA-Seq analysis - batch effect correction with SVA

###### QC AND NORMALIZATION #####
setwd() ##to make a new folder to save figures, need heatmapCol.R script and edgeR.data file
getwd() ##check which working directory you are in. 

options(stringsAsFactors = FALSE)

library(affycoretools)
library(limma)
library(edgeR)
library(rgl) 
library(rtracklayer)
library(sva)

config
head(config)

##check that config file looks similar to below example
#        sample group                                                                            path
# 1 ST-002-M1500406_1     1 /gscmnt/gc2799/fuh/tli/emilee/KALLISTO/H_SO-OV-SEQ-002-M1500406_1/abundance.tsv
# 2 ST-002-M1500488_1     1 /gscmnt/gc2799/fuh/tli/emilee/KALLISTO/H_SO-OV-SEQ-002-M1500488_1/abundance.tsv
# 3 ST-004-M1501201_1     1 /gscmnt/gc2799/fuh/tli/emilee/KALLISTO/H_SO-OV-SEQ-004-M1501201_1/abundance.tsv
# 4 ST-004-M1501761_1     1 /gscmnt/gc2799/fuh/tli/emilee/KALLISTO/H_SO-OV-SEQ-004-M1501761_1/abundance.tsv
# 5 ST-007-M1500407_1     1 /gscmnt/gc2799/fuh/tli/emilee/KALLISTO/H_SO-OV-SEQ-007-M1500407_1/abundance.tsv
# 6 ST-007-M1500409_1     1 /gscmnt/gc2799/fuh/tli/emilee/KALLISTO/H_SO-OV-SEQ-007-M1500409_1/abundance.tsv


####Sept2020: use rawcounts from edgeR data 
load("STvLT_all_Sept2020.RData")
all_rawcounts <- merge(rm1_outliers_rawcounts, rm1_rawcounts, by=0)
rownames(all_rawcounts) <- all_rawcounts$Row.names
all_rawcounts[,1] <-NULL

####Subset config file to include samples to analyze for DE
##config for STvLT
#config <- all_config[c(TRUE,TRUE,TRUE,TRUE, TRUE, TRUE,TRUE, TRUE, TRUE, FALSE, FALSE, TRUE, FALSE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, FALSE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE),]
##config for M_STvLT
#config <- all_config[c(TRUE, FALSE, TRUE, TRUE, TRUE, FALSE, TRUE, TRUE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, TRUE, FALSE, FALSE, TRUE, FALSE, TRUE, FALSE, FALSE, TRUE, TRUE, FALSE, FALSE, TRUE, FALSE, FALSE, FALSE, TRUE, FALSE, TRUE, FALSE, TRUE, TRUE, FALSE, FALSE, TRUE, TRUE, FALSE, FALSE, TRUE, FALSE, TRUE, TRUE, FALSE, TRUE, FALSE, TRUE, FALSE, TRUE, FALSE, TRUE, FALSE, TRUE, FALSE, TRUE, FALSE, TRUE, FALSE, TRUE, FALSE, TRUE, FALSE, TRUE, FALSE, TRUE, FALSE, TRUE, FALSE, TRUE, FALSE),]
##config for P_STvLT
config <- all_config[c(FALSE, TRUE, FALSE, FALSE, FALSE, TRUE, FALSE, FALSE, TRUE, FALSE, FALSE, TRUE, FALSE, TRUE, FALSE, TRUE, TRUE, FALSE, TRUE, FALSE, TRUE, TRUE, FALSE, FALSE, TRUE, TRUE, FALSE, TRUE, TRUE, TRUE, FALSE, TRUE, FALSE, TRUE, FALSE, FALSE, TRUE, FALSE, FALSE, FALSE, TRUE, TRUE, FALSE, TRUE, FALSE, FALSE, TRUE, FALSE, TRUE, FALSE, TRUE, FALSE, TRUE, FALSE, TRUE, FALSE, TRUE, FALSE, TRUE, FALSE,TRUE, FALSE, TRUE, FALSE, TRUE, FALSE, TRUE, FALSE, TRUE, FALSE, TRUE, FALSE, TRUE),]

####Subset rawcounts file to include samples to analyze for DE
##rawcounts for STvLT
#rawcounts2 <- all_rawcounts[,c(TRUE,TRUE,TRUE,TRUE, TRUE, TRUE,TRUE, TRUE, TRUE, FALSE, FALSE, TRUE, FALSE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, FALSE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE,TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE)]
##rawcounts for M_STvLT
#rawcounts2 <- all_rawcounts[,c(TRUE, FALSE, TRUE, TRUE, TRUE, FALSE, TRUE, TRUE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, TRUE, FALSE, FALSE, TRUE, FALSE, TRUE, FALSE, FALSE, TRUE, TRUE, FALSE, FALSE, TRUE, FALSE, FALSE, FALSE, TRUE, FALSE, TRUE, FALSE, TRUE, TRUE, FALSE, FALSE, TRUE, TRUE, FALSE, FALSE, TRUE, FALSE, TRUE, TRUE, FALSE, TRUE, FALSE, TRUE, FALSE, TRUE, FALSE, TRUE, FALSE, TRUE, FALSE, TRUE, FALSE, TRUE, FALSE, TRUE, FALSE, TRUE, FALSE, TRUE, FALSE, TRUE, FALSE, TRUE, FALSE, TRUE, FALSE)]
##for P_STvLT
rawcounts2 <- all_rawcounts[,c(FALSE, TRUE, FALSE, FALSE, FALSE, TRUE, FALSE, FALSE, TRUE, FALSE, FALSE, TRUE, FALSE, TRUE, FALSE, TRUE, TRUE, FALSE, TRUE, FALSE, TRUE, TRUE, FALSE, FALSE, TRUE, TRUE, FALSE, TRUE, TRUE, TRUE, FALSE, TRUE, FALSE, TRUE, FALSE, FALSE, TRUE, FALSE, FALSE, FALSE, TRUE, TRUE, FALSE, TRUE, FALSE, FALSE, TRUE, FALSE, TRUE, FALSE, TRUE, FALSE, TRUE, FALSE, TRUE, FALSE, TRUE, FALSE, TRUE, FALSE,TRUE, FALSE, TRUE, FALSE, TRUE, FALSE, TRUE, FALSE, TRUE, FALSE, TRUE, FALSE, TRUE)]

####Make DEGlist from raw counts
y <- DGEList(counts=rawcounts2);


#### Make new column with group name for graphing and modeling easier later:
config$GpF <- config$group
#config$group <- as.factor(c(rep("ST", 40), rep("LT", 28)))
#config$group <- as.factor(c(rep("P", 23), rep("M", 22))) ##not going to work when P and M samples are mixed up in the config list (at least I think thats how this list works)
#config$group <- as.factor(c("P","M","M","P","M","P","M","P","P","M","P","M","P","M","P","M","P","M","M","P","P","M","M","P","P","M","P","M","P","P","M","P","M","P","M","M","P","P","M","M","P","P","M","P","M"))
#config$group <- as.factor(c("P","M","M","P","M","P","M","P","P","M","P","M","P","M","P","M","P","M","M","P","P","M","M","P","P","M","P","M","P","P","M","P","M","P","M","M","P","P","M","M","P","P","M","P","M", "M", "P", "M", "P", "M", "P", "M", "P", "M", "P", "M", "P", "M", "P", "M", "P", "M", "P", "M", "P", "M", "P", "M", "P", "M", "P", "M", "P"))
#config$group <- as.factor(c(rep("ST", 18), rep("LT", 14)))
#config$group <- as.factor(c(rep("ST", 22), rep("LT", 14)))

#####prepare labeling for data
###matched FF PvM
#config$group <- as.factor(c("P","M", "M", "P", "M", "P", "P", "M", "P", "M", "P", "M", "M", "P", "P", "M"))
###rmoutliers July2020
#config$group <- as.factor(c(rep("ST", 39), rep("LT", 28)))
####rmoutliers STvLT Sept2020
#config$group <- as.factor(c(rep("ST", 41), rep("LT", 28)))
#### M STvLT Sept2020
#config$group <- as.factor(c(rep("ST", 20), rep("LT", 14)))
#### P STvLT Sept2020
config$group <- as.factor(c(rep("ST", 21), rep("LT", 14)))

config$group
config$GpF

#### Factor data can be easily transformed into numbers, which later can be interpreted as colors
#### for graphing purposes:

config$col <- as.numeric(config$group)
head(config)


#### Read in the count data using the function readDGE from edgeR 
###note that you can use any new config and load it here as .Rdata and make new config a y object. 
d = y

class(d)
# [1] "DGEList"
# attr(,"package")
# [1] "edgeR"

dim(d) #make sure this matches the number of genes and number of samples

####example outputs
########33817 69 #for rmoutliers Sept2020
########33817 34 #for M STvLT Sept2020
########33817 35 #for P STvLT Sept2020

names(d)
# "counts" "samples" 

d$samples   

####example output
#                       group lib.size norm.factors
# ST-002-M1500406_1         1 20588113            1
# ST-002-M1500488_1         1 31787689            1
# ST-004-M1501201_1         1 57687911            1
# ST-004-M1501761_1         1  3060289            1
# ST-007-M1500407_1         1 29522946            1
# ST-007-M1500409_1         1 32118895            1

### correct groups

d$samples$group = config$group
d$samples$GpF = config$GpF
d$samples$col = config$col

#### Look at the range of counts per sample:

apply(d$counts, 2, summary)

#### All 71 samples have some 0 counts
#### Check that the most highly expressed gene in each sample doesn't make up more than ~10% of the reads:

apply(d$counts, 2, max) / d$samples$lib.size
# all good



#### Examine the overall distributions of counts to see if any of the samples are different. 
#### Because of the extreme range in the data, converting to log2 scale helps in plotting; 
#### however, you can't take the log of 0, so we need to add a small constant before taking the logs 
#### and plotting the distributions. Since the smallest count value is 1, pick a constant smaller 
#### than that, such as 0.1

jpeg("Plots_QCandNorm/plotDensities_RawData.jpeg", width=7, height=7, units="in", res=300, quality=100)
plotDensities(log2(d$counts + 0.1), group = d$samples$GpF, col = 1:4, legend = "topright" )
dev.off()

# The shapes are roughly similar, except for some samples that have higher bumps. 
# Ideally you wouldn't be able to see group differences here. These are raw counts, 
# so they might be explained by differences in library sizes. 


#### To check the library size for each sample, sum each column and divide by 1e6
#### to get library size in millions

jpeg("Plots_QCandNorm/LibrarySizes_RawData.jpeg", width=7, height=7, units="in", res=300, quality=100)
barplot(colSums(d$counts) / 1e6, ylab = "Total Number of Reads (Million)", las = 2, 
        col = d$samples$col, main = "Library Sizes", cex.axis = 0.8 )
dev.off()


#### There often can be ~50% difference in library sizes, which can affect sampling 
#### of low expression genes. Check the relationship between number genes with 
#### 0 reads and total library size:

jpeg("Plots_QCandNorm/GenesWithZeroReadsVSLibZise.jpeg", width=7, height=7, units="in", res=300, quality=100)
plot(colSums(d$counts), colSums(d$counts == 0), pch = 16, col = d$samples$col)
dev.off()

# There is some correlation. If you see this, then FPKM might not be good to use 

#### Check how many genes in the data set have zero counts for all samples:

sum(rowSums(d$counts) == 0)

#975 genes for rmoutliers Sept2020
#1900 genes for M STvLT Sept2020
#1383 genes for P STvLT Sept2020

#### Do simple clustering on the raw counts. They need to have a small constant added and the log2 taken. 
#### Then, since we want to cluster the arrays, the data matrix must be transposed so that the 
#### rows X column are arrays X probes. Finally, we calculate a distance metric between the arrays 
#### and perform the hierarchical cluster. All this can be done in:

hc.raw <- hclust(dist(t(log2(d$counts + 0.1))), method = "average")

#### Then plot it:

jpeg("Plots_QCandNorm/HierCluster_RawData.jpeg", width=12, height=7, units="in", res=300, quality=100)
plot(hc.raw, hang = -1, main = "Hierarchical Cluster - Raw Values", sub = "", xlab = "", cex = 0.9)
dev.off()
# samples from the two groups appear to group together; ST-004 could be an outliet


#### PCA plot. This will do a principle components analysis, which compresses the
#### the info from all genes into just a few 'principle components', and plots the
#### results. The plot can either be a scree plot (which shows how much variation
#### each principle component explains) or a plot of 2 PCs.

#### First, check the screeplot to see how many of the PCs explain much variation.
#### By definition, PC1 > PC2 > PC3, etc.

jpeg("Plots_QCandNorm/Screeplot_RawData.jpeg", width=7, height=7, units="in", res=300, quality=100)
plotPCA(log2(d$counts + 0.1), screeplot = T)
dev.off()

#### Second, plot PC1 vs. PC2, coloring the samples by group and adding
#### on the sample names:

jpeg("Plots_QCandNorm/PCA_RawData.jpeg", width=15, height=15, units="in", res=300, quality=100)
plotPCA(log2(d$counts + 0.1), pch = 16, col = d$samples$col, groupnames = levels(
  d$samples$group), addtext = d$samples$Sample, main = "PCA on Raw Counts", legend=(0))
dev.off()

#### FILTERING ####

#### We need to filter out genes with 0 or few counts total across all samples. 
#### Look at the distribution of total counts per gene:

jpeg("Plots_QCandNorm/TotalCountsperGene_log2scale.jpeg", width=7, height=7, units="in", res=300, quality=100)
hist(log2(rowSums(d$counts + 0.1)), 1000, main = "Total Counts per Gene, log2 scale")
dev.off()

#### It's difficult to quantify what is "real" expression, but a threshold of 1 count 
#### per million (CPM) mapped reads is often used. When filtering on CPM values, include 
#### the TMM normalization factors:

d <- calcNormFactors(d)
d$samples


#### Make a barplot of the norm.factors to see if there are any group differences:

jpeg("Plots_QCandNorm/TMMNormFactors.jpeg", width=7, height=7, units="in", res=300, quality=100)
barplot(d$samples$norm.factors, col = d$samples$col, main = "TMM Normalization Factors",
        names.arg = d$samples$Label, las = 2)
dev.off()

#### To calculate "modified" CPM values with the TMM normalizations, do:

cpm.values <- cpm(d$counts)
head(cpm.values)

#### To ask how many samples had a value of at least 1 cpm, do:

above1cpm <- rowSums(cpm.values  >= 1)
table(above1cpm)

#### Most genes do have 1 cpm in all samples. Can also do a histogram:

jpeg("Plots_QCandNorm/Above1CPM.jpeg", width=7, height=7, units="in", res=300, quality=100)
hist(above1cpm, xlab = "Number of Samples with > 1 CPM", ylab = "Number of Genes")
dev.off()

#### These "present" histograms are usually bi-modal - most genes either present in 
#### all samples or absent in all samples.
#### Require at least 1cpm in half of the samples ~36 samples)

sum(above1cpm >= 17) #change to half the half n(#samples)
#########15838 for rmoutliers Sept2020, changed number to 34 (69/2)
#########15770 for M STvLT Sept2020, changed number to 17 (34/2)
#########16118 for P STvLT Sept2020, changed number to 17 (35/2)

sum(above1cpm >= 18) #22.5 is half the sample size in ST_PvM cohort
#####16036 (>=18) for P_STvLT
#######15854 (>=33) for rm outliers july2020
########15963 (>=33) for rm outliers Sept2020

sum(above1cpm >= 17) / length(above1cpm)
# 0.4683443 for rm outliers Sept2020
# 0.4718633 for M STvLT Sept 2020
# 0.4766242 for P STvLT Sept 2020

sum(above1cpm >= 33) / length(above1cpm)



#### Check on the total counts per gene for the kept genes vs. filtered genes:

jpeg("Plots_QCandNorm/AllvsKeptvsFilteredGenes.jpeg", width=7, height=7, units="in", res=300, quality=100)
temp <- log2(rowSums(d$counts + 0.1))
layout(matrix( 1:3 , 3, 1 ))
hist(temp, 1000, main = "All Genes", xlim = c(0, 20))
hist(temp[above1cpm >= 12], 1000, main = "Kept Genes", xlim = c(0, 20))
hist(temp[above1cpm < 12], 1000, main = "Filtered Genes", xlim = c(0, 20))
dev.off()

#### Now filter the raw data; you can subset directly on DGEList objects using rows
#### and/or columns, but there is an extra argument keep.lib.sizes that will re-calculate
#### the d$samples$lib.size for you after cutting out genes:

d.filt <- d[above1cpm >= 17, , keep.lib.sizes = FALSE]
#d.filt <- d[above1cpm >= 18, , keep.lib.sizes = FALSE]


#### Check what percentage of the total counts and genes were kept after filtering:

d.filt$samples$lib.size / d$samples$lib.size


nrow(d.filt) / nrow(d)
# About 99.99% of counts were kept, even though only ~46% of genes were kept after filtering out 
# no- and low-information genes.

####make sure to save progress as you go, can change these names for each analysis 
save.image("P_STvLT_Sept2020.RData")
savehistory("P_STvLT_Sept2020.Rhistory")

#### NORMALIZATION USING edgeR ####

#### Since we filtered out some genes, we should probably also re-do the TMM norm factors,
#### although in practice this makes very little difference.

d.filt <- calcNormFactors(d.filt)

#### Remember the TMM factors are multipliers for the library size to get an "effective"
#### library size to use for the adjustment of read depth:

layout(1)

jpeg("Plots_QCandNorm/ActualLibSizes.jpeg", width=7, height=7, units="in", res=300, quality=100)
barplot( d.filt$samples$lib.size / 1e6, main = "Actual Library Sizes", 
         las = 2, names.arg = d$samples$Label, col = d$samples$col, ylab = "Millions of Reads", 
         cex.axis = 0.8 )
dev.off()

jpeg("Plots_QCandNorm/EffectiveLibSizes.jpeg", width=7, height=7, units="in", res=300, quality=100)
barplot( d.filt$samples$lib.size * d.filt$samples$norm.factors / 1e6, main = "Effective Library Sizes", 
         las = 2, names.arg = d$samples$Label, col = d$samples$col, ylab = "Millions of Reads", 
         cex.axis = 0.8 )
dev.off()

#### Cluster the samples based on the normalized values...
#### edgeR has a function to cluster the samples based on multidimensional scaling,
#### which is sort of like PCA, but has a distance measure appropriate for count data:

jpeg("Plots_QCandNorm/MDSplot.jpeg", width=7, height=7, units="in", res=300, quality=100)
plotMDS(d.filt, main = "MDS plot", col = d.filt$samples$col, top = 1000 )
dev.off()

#### The edgeR vignette suggests using the cpm() function to get modified CPM values
#### like before, except on the log2 scale for clustering or heatmaps. Since you can't
#### take a log of 0, a small value needs to be added. The default is to add an average 
#### of 0.25 counts to each gene, proportional to library size:

logCPM <- cpm(d.filt, log = T)
head(logCPM)

#### Examine the logCPM values; first, the distribution of values:

jpeg("Plots_QCandNorm/plotDensities_NormData.jpeg", width=7, height=7, units="in", res=300, quality=100)
plotDensities(logCPM, group = d$samples$GpF, col = 1:4)
dev.off()

#qq bump still seems to be higher than the others.

#### Second, basic hierarchical clustering:

hc.edgeR <- hclust(dist(t(logCPM)), method = "average")

jpeg("Plots_QCandNorm/HierCluster_edgeRlogCPMValues.jpeg", width=12, height=7, units="in", res=300, quality=100)
plot(hc.edgeR, hang = -1, main = "Hierarchical Cluster - edgeR logCPM Values", sub = "", xlab = 
       "", cex = 0.9)
dev.off()

#### Third, can do both PCA plots. First, check the screeplot to see how much variation 
#### each PC explains:

jpeg("Plots_QCandNorm/Screeplot_NormData.jpeg", width=7, height=7, units="in", res=300, quality=100)
plotPCA(logCPM,screeplot=T)
dev.off()

#### Plot PC1 vs. PC2, coloring the samples by group and adding on the sample names

jpeg("Plots_QCandNorm/PCA_NormData.jpeg", width=15, height=15, units="in", res=300, quality=100)
plotPCA(logCPM, pch = 16, col = d.filt$samples$col, groupnames = levels(d.filt$samples$group), 
        addtext = d.filt$samples$Sample, main = "PCA - edgeR logCPM", legend=(0))
dev.off()


####======== End of the QC and normalization using the edgeR package.========#### 

####save progress, can rename these files
save.image("P_STvLT_Sept2020.RData")
savehistory("P_STvLT_Sept2020.Rhistory")


#### STATISTICAL ANALYSIS USING EDGER #####

#### If starting new session from here, first run:

options(stringsAsFactors = FALSE)

#had to install xcode on my Mac to download WGCNA
library(WGCNA)
library(affycoretools)
library(limma)
library(edgeR)

####load data if started from new session
###load("RNASeq_Analysis2.RData")

class(d.filt)
# [1] "DGElist"
# attr(,"package")
# "edgeR"

#### To analyze the RNA-Seq data using either the edgeR package or the limma package using voom values, 
#### it is necessary to construct design and contrast matrices, which we can make using the limma package.
#### Create the design matrix with one coefficient (column) for each of the 2 groups (qq and QQ):
#d.filt
design <- model.matrix(~0 + d.filt$samples$group)
#design
colnames(design) <- levels(d.filt$samples$group)
rownames(design) <- d.filt$samples$Label
design
####example output
#      LT ST
# [1,]  0  1
# [2,]  0  1
# [3,]  0  1
# [4,]  0  1
# ...
# ...
# attr(,"assign")
# [1] 1 1
# attr(,"contrasts")
# attr(,"contrasts")$`d.filt$samples$group`
# [1] "contr.treatment"

#### Specify contrasts of interest. Here we will a logical pairwise
#### comparisons between the two groups: 

cont.matrix <- makeContrasts(STvsLT = ST - LT, levels=design)
#cont.matrix <- makeContrasts(PvsM = P - M, levels=design)

cont.matrix
#     Contrasts
# Levels STvsLT
# LT     -1
# ST      1

#########HERE IS WHERE TO START DOING EDGER WITHOUT SVA, TO SEE RESULTS ##################

####======== PROBLEMS OBSERVED: ========####

#### Results from normalization were not very good.Samples are not separating as expecting in PCA plots.
#### Probably due to the presence of outliers or other unknown influencing covariates.
#### In order to improve such results, we can use the sva package to estimate unknown batch effects.

#### Surrogate Variables Analysis (SVA) package for removing batch effects and other unwanted 
#### variation in high-throughput experiments.  It contains  functions  for  identifying  and  
#### building surrogate variables for high-dimensional data sets.  Surrogate variables are covariates 
#### constructed directly from  high-dimensional  data  (like  gene  expression/RNA  
#### sequencing/methylation/brain  imaging  data) that can be used in subsequent analyses to adjust 
#### for unknown, unmodeled, or latent sources of noise. The sva package can be used to remove 
#### artifacts in two ways:  (1) identifying and estimating surrogate variables  for  unknown  
#### sources  of  variation  in  high-throughput  experiments  and  (2)  directly  removing known 
#### batch effects using ComBat. 
#### The sva package adjusts for unkown potential covariates that may be influencing the results.

#### SVA + edgeR ####

library(sva)


##to remove samples, do it with the d.filt object before this step. 

mod <- model.matrix(~d.filt$samples$group) # Full model
colnames(mod)[-1] <- paste0(levels(d.filt$samples$group)[-1], "vs", levels(d.filt$samples$group)[1])
head(mod)
####example output
#   (Intercept) STvsLT
# 1           1      1
# 2           1      1
# 3           1      1
# 4           1      1
# 5           1      1
# 6           1      1


mod0 <- model.matrix(~1, data=d.filt$samples$group) ## Null model
head(mod0)
####example output
#   (Intercept)
# 1           1
# 2           1
# 3           1
# 4           1
# 5           1
# 6           1

cont.mod <- makeContrasts(STvsLT = STvsLT, levels = mod)
#cont.mod <- makeContrasts(PvsM = PvsM, levels = mod)

cont.mod
# Contrasts
# Levels      STvsLT
# Intercept      0
# STvsLT         1


n.sv.eR <- num.sv(logCPM, mod)
n.sv.eR
#####7 variables for STvLT rmoutliers Sept2020
#####4 variables for M STvLT Sept2020
#####5 variables for P STvLT Sept2020

svobj.eR <- sva(logCPM, mod,mod0, n.sv=n.sv.eR)
##example output
# Number of significant surrogate variables is:  8
# Iteration (out of 5 ):1  2  3  4  5  

head(svobj.eR$sv)
####example output
#              [,1]        [,2]        [,3]        [,4]          [,5]        [,6]        [,7]        [,8]
# [1,] -0.483419970 -0.28944819  0.26245247 -0.70145826  0.0254851025 -0.02305705  0.09836440 -0.04928099
# [2,] -0.422633351 -0.03671694 -0.06458532 -0.01077011 -0.0006042416 -0.01516047 -0.28589129  0.05901358
# [3,] -0.009081903  0.13667050  0.18628191  0.04576136 -0.0468786462  0.07200384 -0.13284860  0.07308750
# [4,]  0.054236880 -0.61685953  0.53247975  0.43290479  0.1668004614  0.06064549  0.04912447 -0.05592462
# [5,] -0.058385899  0.11349504  0.03469986  0.05547762  0.0749636349 -0.20533021  0.08124823 -0.09401865
# [6,] -0.108794094  0.09074725  0.01999002  0.12055146 -0.1149266250 -0.32755260  0.19360126 -0.15798359

##ploting batch effects from sva
for (i in 1:ncol(svobj.eR$sv)){
  jpeg(filename = paste("Plots_SVA/SV", i, ".jpeg"), width=7, height=7, units="in", res=300, quality=100)
  barplot(svobj.eR$sv[,i], col = d.filt$samples$col)
  dev.off()
}


## Add all 8 surrogate variables (SV1 - SV8) to new design and contrast matrices, depending on number of variables for each analysis 
## and proceed with edgeR:

modSV <- model.matrix(~0 + d.filt$samples$group + svobj.eR$sv)
colnames(modSV)[1:2] <- c("ST", "LT")
colnames(modSV)[3:7] <- paste0("sv", 1:5) #adjust these numbers based on the number of variables (for STvLT analysis, the numbers were 3:10....1:8)
####3:9...1:7 #rmoutliers Sept2020
####3:6...1:4 #M STvLT Sept2020
####3:7...1:5 #P STvLT Sept2020

rownames(modSV) <- d.filt$samples$Label
head(modSV)


cont.modSV <- makeContrasts(STvsLT = ST - LT, levels=modSV)
###interpretation should be that if its high/low in ST, and opposite in LT
#cont.modSV <- makeContrasts(PvsM = P - M, levels=modSV)

cont.modSV
#     Contrasts
# Levels STvsLT
# LT      -1
# ST       1
# sv1      0
# sv2      0
# sv3      0
# sv4      0
# sv5      0
# sv6      0
# sv7      0
# sv8      0

#### Now do DE testing using edgeR.

d.filt.sv <- estimateGLMCommonDisp(d.filt, modSV, verbose = TRUE)
#Disp = 0.38715 , BCV = 0.6222 for STvLT, rmoutliers Sept2020
#Disp = 0.40458 , BCV = 0.6361 for M STvLT Sept2020
#Disp = 0.3817 , BCV = 0.6178 for P STvLT Sept2020

d.filt.sv <- estimateGLMTrendedDisp(d.filt.sv, modSV)
d.filt.sv <- estimateGLMTagwiseDisp(d.filt.sv, modSV)

jpeg("Plots_SVA/plotBCV.jpeg", width=7, height=7, units="in", res=300, quality=100)
plotBCV(d.filt.sv)
dev.off()

#### Step 1: Estimate model coefficients from count data and design matrix:

##fitting model
fit.edgeR.sv <- glmFit(d.filt.sv, modSV)

#### Steps 2 and 3: Specify contrasts of interest and do empirical Bayes "shrinkage" of 
#### variances and calculate test statistics. In edgeR, both of these are done in the same
#### function, glmLRT. Also, unlike regular limma modeling, must fit each
#### column of the design matrix separately:

##calculates DE between STvLT, and between each variable from the batch effects, change the er.SV numbers depending on number of variables for each analysis 
eR.STvsLT.sv <- glmLRT(fit.edgeR.sv, contrast = cont.modSV[ , 1]) #CHANGED eR.STvsLT.SV
eR.SV1.sv <- glmLRT(fit.edgeR.sv, coef = "sv1")
eR.SV2.sv <- glmLRT(fit.edgeR.sv, coef = "sv2")
eR.SV3.sv <- glmLRT(fit.edgeR.sv, coef = "sv3")
eR.SV4.sv <- glmLRT(fit.edgeR.sv, coef = "sv4")
eR.SV5.sv <- glmLRT(fit.edgeR.sv, coef = "sv5")
#eR.SV6.sv <- glmLRT(fit.edgeR.sv, coef = "sv6")
#eR.SV7.sv <- glmLRT(fit.edgeR.sv, coef = "sv7")
#eR.SV8.sv <- glmLRT(fit.edgeR.sv, coef = "sv8")
#eR.SV9.sv <- glmLRT(fit.edgeR.sv, coef = "sv9")

#RUN THIS PART TO EXCLUDE SVA
#eR.STvsLT.sv <- glmLRT(fit.edgeR, contrast = cont.mod[ , 1])

#### Step 4: Correct for multiple tests and extract relevant data;
##RUN THIS PART TO EXCLUDE SVA
#PvsM = decideTestsDGE(eR.PvsM)

##combines results from variables into a frame
edgeR.coded.0.05FDR.sv <- new("TestResults", cbind(STvsLT = decideTestsDGE(eR.STvsLT.sv)[ , 1],
                                                   sv1 = decideTestsDGE(eR.SV1.sv)[ ,1],
                                                   sv2 = decideTestsDGE(eR.SV2.sv)[ ,1],
                                                   sv3 = decideTestsDGE(eR.SV3.sv)[ ,1],
                                                   sv4 = decideTestsDGE(eR.SV4.sv)[ ,1],
                                                   sv5 = decideTestsDGE(eR.SV5.sv)[ ,1]))

rownames(edgeR.coded.0.05FDR.sv) <- rownames(eR.STvsLT.sv$table)
#rownames(edgeR.coded.0.05FDR.sv) <- rownames(eR.PvsM.sv$table)
summary(edgeR.coded.0.05FDR.sv, p.value = 0.05)

summary(decideTestsDGE(eR.STvsLT.sv, p.value = 0.01))


summary(decideTestsDGE(eR.STvsLT.sv, p.value = 0.1))


###what if the pvalue is changed to 0.1 (updated June 2020), cahnge change to see different cutoff levels 
summary(decideTestsDGE(eR.STvsLT.sv, p.value = 0.1))



####***WGCNA for network analysis for looking at network and pathways since this is a lot of DE genes

eR.STvsLT.detailed.sv <- topTags(eR.STvsLT.sv, n = Inf, sort.by = "none")$table  
###this object shows the log fold change
eR.STvsLT.detailed.sv$FC <- 2 ^ abs(eR.STvsLT.detailed.sv$logFC) * 
  sign(eR.STvsLT.detailed.sv$logFC)

hist(eR.STvsLT.sv$table$PValue, 1000)
hist(eR.STvsLT.detailed.sv$FDR, 1000)
hist(eR.STvsLT.detailed.sv$FDR, 1000, xlim = c(0,0.4), ylim=c(0,100))

edgeR.coded.0.01FDR.sv <- new("TestResults", cbind(STvsLT = decideTestsDGE(eR.STvsLT.sv, p.value = 0.01)[ , 1]))
rownames(edgeR.coded.0.01FDR.sv) <- rownames(eR.STvsLT.sv$table)
summary(edgeR.coded.0.01FDR.sv)

##the following was added June 2020, to see DEGs of FDR <0.1. 
edgeR.coded.0.1FDR.sv <- new("TestResults", cbind(STvsLT = decideTestsDGE(eR.STvsLT.sv, p.value = 0.1)[ , 1]))
rownames(edgeR.coded.0.1FDR.sv) <- rownames(eR.STvsLT.sv$table)
summary(edgeR.coded.0.1FDR.sv)



sva.rm.eR <- removeBatchEffect(logCPM, design = design, covariates = modSV[,3:7])  
##changed to modSV[,3:9] for rmoutliers Sept2020
##changed to modSV[,3:6] for M STvLT Sept2020
##changed to modSV[,3:7] for P STvLT Sept2020
#sva.rm.eR <- removeBatchEffect(logCPM, design = design, covariates = modSV[,3:5]) #for M_STvLT

jpeg("Plots_SVA/PCA_edgeRlogCPM_noSV.jpeg", width=15, height=15, units="in", res=300, quality=100)
plotPCA(sva.rm.eR, pch = 16, col = d.filt.sv$samples$col, groupnames = levels(d.filt.sv$samples$group), 
        addtext = d.filt.sv$samples$Sample, main = "PCA - PvsM")
dev.off()

####change this file name to what you want the DEG results to be saved as 
write.csv(eR.STvsLT.detailed.sv, file = "RNASeq_STvsLT_FDR_batchCorrected_results.csv")

####save progress
save.image("P_STvLT_Sept2020.RData")
savehistory("P_STvLT_Sept2020.Rhistory")

####check environment from session
sessionInfo()



##### GENE LISTS #####

options(stringsAsFactors = FALSE)

library(WGCNA)
library(affycoretools)
library(limma)
library(edgeR)

####load R data if starting from new session
#load("STvsLT.RData")

#### The d.filt.sv$samples object has the phenotypic information of the samples:

head(d.filt.sv$samples)
####example output
#                   group lib.size norm.factors GpF col            Sample
# ST-002-M1500406_1    ST 20332459    1.1437124   1   2 ST-002-M1500406_1
# ST-002-M1500488_1    ST 31519919    0.8771427   1   2 ST-002-M1500488_1
# ST-004-M1501201_1    ST 57204395    0.9084580   1   2 ST-004-M1501201_1
# ST-004-M1501761_1    ST  2947210    0.5831815   1   2 ST-004-M1501761_1
# ST-007-M1500407_1    ST 29169958    1.1435391   1   2 ST-007-M1500407_1
# ST-007-M1500409_1    ST 31784180    1.1158211   1   2 ST-007-M1500409_1

#### This data was analyzed in edgeR using the following design matrix:

head(modSV)
####example output
#      LT ST          sv1         sv2         sv3         sv4           sv5         sv6         sv7         sv8
# [1,]  0  1 -0.483419970 -0.28944819  0.26245247 -0.70145826  0.0254851025 -0.02305705  0.09836440 -0.04928099
# [2,]  0  1 -0.422633351 -0.03671694 -0.06458532 -0.01077011 -0.0006042416 -0.01516047 -0.28589129  0.05901358
# [3,]  0  1 -0.009081903  0.13667050  0.18628191  0.04576136 -0.0468786462  0.07200384 -0.13284860  0.07308750
# [4,]  0  1  0.054236880 -0.61685953  0.53247975  0.43290479  0.1668004614  0.06064549  0.04912447 -0.05592462
# [5,]  0  1 -0.058385899  0.11349504  0.03469986  0.05547762  0.0749636349 -0.20533021  0.08124823 -0.09401865
# [6,]  0  1 -0.108794094  0.09074725  0.01999002  0.12055146 -0.1149266250 -0.32755260  0.19360126 -0.15798359

#### And the following contrast matrix:

cont.modSV
#     Contrasts
# Levels STvsLT
# LT      -1
# ST       1
# sv1      0
# sv2      0
# sv3      0
# sv4      0
# sv5      0
# sv6      0
# sv7      0
# sv8      0

#### Here is the summary of "significant" (FDR < 0.01) genes for QQvsqq:

summary(edgeR.coded.0.01FDR.sv)


#### Here is the summary of "significant" (FDR < 0.05) genes for QQvsqq:

summary(edgeR.coded.0.05FDR.sv)   

summary(edgeR.coded.0.1FDR.sv) #added June 2020


#### Recall that the question in which we are interested is:
# 1. Does expression differ between ST and LT ?

# At FDR = 0.01

sum(eR.STvsLT.detailed.sv$FDR < 0.01)


# At FDR = 0.05

sum(eR.STvsLT.detailed.sv$FDR < 0.05)


# At FDR = 0.1

sum(eR.STvsLT.detailed.sv$FDR < 0.1) 

topTags(eR.STvsLT.sv)
####example output
# Coefficient:  -1*LT 1*ST 
#                 logFC    logCPM       LR        PValue           FDR
# ZNF738       6.303626 3.3829658 479.9282 2.215017e-106 3.500834e-102
# TATDN2P2     6.238937 2.2479866 462.6549 1.271021e-102  1.004424e-98
# C8orf49      7.165724 0.2066449 427.5103  5.659444e-95  2.451539e-91
# RP1-130L23.1 7.050164 0.9992780 427.3268  6.204465e-95  2.451539e-91
# KCNJ6        6.600935 4.6101964 423.3667  4.514823e-94  1.427135e-90
# ZNF252P      8.109938 6.6297308 422.5798  6.697814e-94  1.764316e-90
# HMGB3P9      6.797267 0.9029699 421.8731  9.544266e-94  2.154959e-90
# RP11-693N9.2 6.634157 1.5403811 404.1492  6.881952e-90  1.359616e-86
# SZRD1        5.262078 5.6602128 398.7594  1.025635e-88  1.801129e-85
# AC023271.1   9.331270 2.3603870 396.7762  2.771500e-88  4.380356e-85



##--- Gene lists for STvsLT---##

G.all.0.01FDR <- rownames(edgeR.coded.0.01FDR.sv)[edgeR.coded.0.01FDR.sv[,1] !=0]
length(G.all.0.01FDR)


G.all.0.05FDR <- rownames(edgeR.coded.0.05FDR.sv)[edgeR.coded.0.05FDR.sv[,1] !=0]
length(G.all.0.05FDR)


#added June 2020
G.all.0.1FDR <- rownames(edgeR.coded.0.1FDR.sv)[edgeR.coded.0.1FDR.sv[,1] !=0]
length(G.all.0.1FDR)

####save progress
save.image("P_STvsLT_Sept2020.RData")
savehistory("P_STvsLT_Sept2020.Rhistory")



##### HEATMAPS #####

options(stringsAsFactors = FALSE)

library(affycoretools)
library(limma)
library(edgeR)
library(Heatplus) #can install with BiocManager::install("Heatplus")
library(gplots)
library(WGCNA) #can install with BiocManager::install("WGCNA"); say yes to install from source; only need this package if you want the heatmaps with gene clusters

#### To make a heatmap, we first need to select a subset of genes and pull out their
#### logCPM expression values. Use the lists of genes created previously:

head(G.all.0.01FDR)
 

length(G.all.0.01FDR)


head(G.all.0.05FDR)
 

length(G.all.0.05FDR)


#added June 2020
head(G.all.0.1FDR)
length(G.all.0.1FDR)


#### 0.01FDR:
heatdata.0.01FDR <- sva.rm.eR[G.all.0.01FDR, ]
dim(heatdata.0.01FDR)


#### 0.05FDR:
heatdata.0.05FDR <- sva.rm.eR[G.all.0.05FDR, ]
dim(heatdata.0.05FDR)


### 0.1 FDR #added June 2020
heatdata.0.1FDR <- sva.rm.eR[G.all.0.1FDR, ]
dim(heatdata.0.1FDR)


#### Check the distribution of expression values in these objects
#dir.create("Plots_Heatmaps")
#Warning message:
#In dir.create("Plots_Heatmaps") : 'Plots_Heatmaps' already exists
jpeg("Plots_Heatmaps/heatdata0.01.jpeg", width=7, height=7, units="in", res=300, quality=100)
hist(heatdata.0.01FDR)
dev.off()

jpeg("Plots_Heatmaps/heatdata0.05.jpeg", width=7, height=7, units="in", res=300, quality=100)
hist(heatdata.0.05FDR)
dev.off()

#added June 2020
jpeg("Plots_Heatmaps/heatdata0.1.jpeg", width=7, height=7, units="in", res=300, quality=100)
hist(heatdata.0.1FDR)
dev.off()

# These values are not centered on zero, as they are modified CPM values (log2 scale).
# For the heatmap, we need values centered on zero

#### Scale the data to get data centered on zero with same standard
#### deviation. Instead of using the auto-scale in the heatmap.2 function, we will
#### do it ourselves in order to compute an appropriate color scale.
#### Note that the scale function works on columns, but we want to scale each row 
#### (gene). Therefore, we must transpose the data before and after scaling.

heatdata.0.01FDR.scaled <- t(scale(t(heatdata.0.01FDR)))

jpeg("Plots_Heatmaps/heatdata0.01scaled.jpeg", width=7, height=7, units="in", res=300, quality=100)
hist(heatdata.0.01FDR.scaled, 100)
dev.off()
##*** pick lim=2 for color --> smaller numbers will give you sharper colors; you basically want to pick the value that marks the lower and upper 25th percentile of your histogram

heatdata.0.05FDR.scaled <- t(scale(t(heatdata.0.05FDR)))

jpeg("Plots_Heatmaps/heatdata0.05scaled.jpeg", width=7, height=7, units="in", res=300, quality=100)
hist(heatdata.0.05FDR.scaled, 100)
dev.off()
## pick lim=2 for color

##added June 2020
heatdata.0.1FDR.scaled <- t(scale(t(heatdata.0.1FDR)))

jpeg("Plots_Heatmaps/heatdata0.1scaled.jpeg", width=7, height=7, units="in", res=300, quality=100)
hist(heatdata.0.1FDR.scaled, 100)
dev.off()

source("heatmapCol.R")

color.0.01FDR.scaled <- heatmapCol(data = heatdata.0.01FDR.scaled, lim = 2, 
                                   col =colorRampPalette(c("blue","white","red"))(128))

jpeg("Plots_Heatmaps/Heatmap_0.01FDR.jpeg", width=25, height=30, units="in", res=300, quality=100)
heatmap.2(heatdata.0.01FDR.scaled, col = color.0.01FDR.scaled, scale = "none", 
          labRow = F, trace = "none", density.info = "none", margins = c(10,4), 
          symbreaks= FALSE, key.xlab = "SD from mean", main = "Scaled Heatmap = PvsM (FDR < 0.01)") 

dev.off()


color.0.05FDR.scaled <- heatmapCol(data = heatdata.0.05FDR.scaled, lim = 2, 
                                   col =colorRampPalette(c("blue","white","red"))(128))

jpeg("Plots_Heatmaps/Heatmap_0.05FDR.jpeg", width=25, height=30, units="in", res=300, quality=100)
heatmap.2(heatdata.0.05FDR.scaled, col = color.0.05FDR.scaled, scale = "none", 
          labRow = F, trace = "none", density.info = "none", margins = c(10,4), 
          symbreaks= FALSE, key.xlab = "SD from mean", main = "Scaled Heatmap - PvsM (FDR < 0.05)")

dev.off()

#added June 2020
color.0.1FDR.scaled <- heatmapCol(data = heatdata.0.1FDR.scaled, lim = 2, 
                                   col =colorRampPalette(c("blue","white","red"))(128))

jpeg("Plots_Heatmaps/Heatmap_0.1FDR.jpeg", width=25, height=30, units="in", res=300, quality=100)
heatmap.2(heatdata.0.1FDR.scaled, col = color.0.1FDR.scaled, scale = "none", 
          labRow = F, trace = "none", density.info = "none", margins = c(10,4), 
          symbreaks= FALSE, key.xlab = "SD from mean", main = "Scaled Heatmap = PvsM (FDR < 0.1)") 

dev.off()

# Calculate your own gene dendrogram and split into clusters using a hybrid
# computation method:

cluster.0.01FDR<- hclust(dist(heatdata.0.01FDR.scaled))

cluster.0.05FDR <- hclust(dist(heatdata.0.05FDR.scaled))

#added June 2020
cluster.0.1FDR <- hclust(dist(heatdata.0.1FDR.scaled))

cutree.clusters.0.01FDR <- cutreeDynamic(dendro = cluster.0.01FDR, method = "hybrid", 
                                         minClusterSize = 40, distM = as.matrix(dist(heatdata.0.01FDR.scaled),
                                                                                deepSplit = 0))


cutree.clusters.0.05FDR <- cutreeDynamic(dendro = cluster.0.05FDR, method = "hybrid", 
                                         minClusterSize = 40, distM = as.matrix(dist(heatdata.0.05FDR.scaled),
                                                                                deepSplit = 0))

#added June 2020
cutree.clusters.0.1FDR <- cutreeDynamic(dendro = cluster.0.1FDR, method = "hybrid", 
                                         minClusterSize = 40, distM = as.matrix(dist(heatdata.0.1FDR.scaled),
                                                                                deepSplit = 0))



# Write out a high-res heatmap to the current working directory:

pdf("Plots_Heatmaps/Heatmap_0.01FDR_clusters.pdf", width=180, height=180, pointsize = 180)
out.heat.0.01 <- heatmap.2(heatdata.0.01FDR.scaled, col = color.0.01FDR.scaled, scale = "none", 
                           labRow = F, trace = "none", density.info = "none", margins = c(15,6), 
                           symbreaks= FALSE, key.xlab = "SD from mean",
                           Colv = T, dendrogram = "row", Rowv = as.dendrogram(cluster.0.01FDR),
                           keysize = 0.5, main = "Scaled Heatmap with Gene Clusters - STvsLT (FDR < 0.01)",
                           RowSideColors=labels2colors(cutree.clusters.0.01FDR))
dev.off()

pdf("Plots_Heatmaps/Heatmap_0.05FDR_clusters.pdf", width=180, height=180, pointsize = 180)
out.heat.0.05 <- heatmap.2(heatdata.0.05FDR.scaled, col = color.0.05FDR.scaled, scale = "none", 
                           labRow = F, trace = "none", density.info = "none", margins = c(15,6), 
                           symbreaks= FALSE, key.xlab = "SD from mean",
                           Colv = T, dendrogram = "row", Rowv = as.dendrogram(cluster.0.05FDR),
                           keysize = 0.5, main = "Scaled Heatmap with Gene Clusters - STvsLT (FDR < 0.05)",
                           RowSideColors=labels2colors(cutree.clusters.0.05FDR))
dev.off()

#added June 2020
pdf("Plots_Heatmaps/Heatmap_0.1FDR_clusters.pdf", width=180, height=180, pointsize = 180)
out.heat.0.1 <- heatmap.2(heatdata.0.1FDR.scaled, col = color.0.1FDR.scaled, scale = "none", 
                           labRow = F, trace = "none", density.info = "none", margins = c(15,6), 
                           symbreaks= FALSE, key.xlab = "SD from mean",
                           Colv = T, dendrogram = "row", Rowv = as.dendrogram(cluster.0.1FDR),
                           keysize = 0.5, main = "Scaled Heatmap with Gene Clusters - STvsLT (FDR < 0.1)",
                           RowSideColors=labels2colors(cutree.clusters.0.1FDR))
dev.off()

# Now output the values plotted in the heatmap, along with their module colors;

heatVal.0.01FDR.out <- data.frame(ID=rownames(heatdata.0.01FDR.scaled), 
                                  heatdata.0.01FDR.scaled[,out.heat.0.01$colInd],
                                  ModuleColor=labels2colors(cutree.clusters.0.01FDR))
heatVal.0.01FDR.out <- heatVal.0.01FDR.out[rev(out.heat.0.01$rowInd),]
heatVal.0.01FDR.out$TopToBottom <- 1:nrow(heatVal.0.01FDR.out)
write.csv(heatVal.0.01FDR.out,file="heatVal_0.01FDR_out.csv",row.names=F)

heatVal.0.05FDR.out <- data.frame(ID=rownames(heatdata.0.05FDR.scaled), 
                                  heatdata.0.05FDR.scaled[,out.heat.0.05$colInd],
                                  ModuleColor=labels2colors(cutree.clusters.0.05FDR))
heatVal.0.05FDR.out <- heatVal.0.05FDR.out[rev(out.heat.0.05$rowInd),]
heatVal.0.05FDR.out$TopToBottom <- 1:nrow(heatVal.0.05FDR.out)
write.csv(heatVal.0.05FDR.out,file="heatVal_0.05FDR_out.csv",row.names=F)

#added June 2020
heatVal.0.1FDR.out <- data.frame(ID=rownames(heatdata.0.1FDR.scaled), 
                                  heatdata.0.1FDR.scaled[,out.heat.0.1$colInd],
                                  ModuleColor=labels2colors(cutree.clusters.0.1FDR))
heatVal.0.1FDR.out <- heatVal.0.1FDR.out[rev(out.heat.0.1$rowInd),]
heatVal.0.1FDR.out$TopToBottom <- 1:nrow(heatVal.0.1FDR.out)
write.csv(heatVal.0.1FDR.out,file="heatVal_0.1FDR_out.csv",row.names=F)


####save progress, can change these file names
save.image("P_STvLT_Sept2020.RData")
savehistory("P_STvLT_Sept2020.Rhistory")


####how to generate the DEG list of gene symbols and ensembl IDs key in R studio:
####upload excel file of DE results from SVA into Rstudio

####change file names to paths, examples written here 
library(readxl)
####example line of code and path saves
STvLT_1FDR_out <- read_excel("/Users/emileekotnik/Documents/DE_Sept2020/STvLT_sept2020_0.1FDR_out.xlsx")
STvLT_01FDR_out <- read_excel("/Users/emileekotnik/Documents/DE_Sept2020/STvLT_sept2020_0.01FDR_out.xlsx")
STvLT_05FDR_out <- read_excel("/Users/emileekotnik/Documents/DE_Sept2020/STvLT_sept2020_0.05FDR_out.xlsx")



View(heatVal_0_01FDR_out)

##make a list of gene symbols from first column of table ("ID" is the heading of the column we need)
symbols1 <- STvLT_1FDR_out$ID
symbols2 <- STvLT_01FDR_out$ID
symbols3 <- STvLT_05FDR_out$ID


##load ensembl package
library("org.Hs.eg.db")

##convert gene symbols to ensembl IDs
ensembls1 <- mapIds(org.Hs.eg.db, keys = symbols1 , keytype="SYMBOL", column="ENSEMBL")
ensembls2 <- mapIds(org.Hs.eg.db, keys = symbols2 , keytype="SYMBOL", column="ENSEMBL")
ensembls3 <- mapIds(org.Hs.eg.db, keys = symbols3 , keytype="SYMBOL", column="ENSEMBL")


##combine the gene sysmbols and ensembl ID characters into a matrix
gene_key1 <- cbind(symbols1, ensembls1)
gene_key2 <- cbind(symbols2, ensembls2)
gene_key3 <- cbind(symbols3, ensembls3)



## write matrix to txt file
write.table(gene_key1, file="/Users/emileekotnik/Documents/DE_Sept2020/STvLT_1FDR_gene_ensembl_key.txt", row.names=FALSE, col.names=FALSE)
write.table(gene_key2, file="/Users/emileekotnik/Documents/DE_Sept2020/STvLT_01FDR_gene_ensembl_key.txt", row.names=FALSE, col.names=FALSE)
write.table(gene_key3, file="/Users/emileekotnik/Documents/DE_Sept2020/STvLT_05FDR_gene_ensembl_key.txt", row.names=FALSE, col.names=FALSE)




####session info from September 2020 analysis shown below
sessionInfo()


# R version 3.6.1 (2019-07-05)
# Platform: x86_64-apple-darwin15.6.0 (64-bit)
# Running under: macOS Mojave 10.14.6
# 
# Matrix products: default
# BLAS:   /System/Library/Frameworks/Accelerate.framework/Versions/A/Frameworks/vecLib.framework/Versions/A/libBLAS.dylib
# LAPACK: /Library/Frameworks/R.framework/Versions/3.6/Resources/lib/libRlapack.dylib
# 
# Random number generation:
#   RNG:     Mersenne-Twister 
# Normal:  Inversion 
# Sample:  Rounding 
# 
# locale:
#   [1] en_US.UTF-8/en_US.UTF-8/en_US.UTF-8/C/en_US.UTF-8/en_US.UTF-8
# 
# attached base packages:
#   [1] parallel  stats     graphics  grDevices utils     datasets  methods   base     
# 
# other attached packages:
#   [1] WGCNA_1.68            fastcluster_1.1.25    dynamicTreeCut_1.63-1 gplots_3.0.1.1        Heatplus_2.32.0       affycoretools_1.58.3  Biobase_2.46.0        BiocGenerics_0.32.0   edgeR_3.28.0         
# [10] limma_3.42.0         
# 
# loaded via a namespace (and not attached):
#   [1] backports_1.1.5             GOstats_2.52.0              Hmisc_4.3-0                 BiocFileCache_1.10.2        plyr_1.8.4                  lazyeval_0.2.2              GSEABase_1.48.0            
# [8] splines_3.6.1               BiocParallel_1.20.0         GenomeInfoDb_1.22.0         ggplot2_3.2.1               robust_0.4-18.1             digest_0.6.23               foreach_1.4.7              
# [15] ensembldb_2.10.2            htmltools_0.4.0             GO.db_3.10.0                gdata_2.18.0                magrittr_1.5                checkmate_1.9.4             memoise_1.1.0              
# [22] BSgenome_1.54.0             fit.models_0.5-14           doParallel_1.0.15           cluster_2.1.0               gcrma_2.58.0                Biostrings_2.54.0           annotate_1.64.0            
# [29] matrixStats_0.55.0          R.utils_2.9.0               ggbio_1.34.0                askpass_1.1                 prettyunits_1.0.2           colorspace_1.4-1            rrcov_1.4-9                
# [36] blob_1.2.0                  rappdirs_0.3.1              xfun_0.11                   dplyr_0.8.3                 crayon_1.3.4                RCurl_1.95-4.12             graph_1.64.0               
# [43] genefilter_1.68.0           impute_1.60.0               zeallot_0.1.0               survival_3.1-8              VariantAnnotation_1.32.0    iterators_1.0.12            glue_1.3.1                 
# [50] gtable_0.3.0                zlibbioc_1.32.0             XVector_0.26.0              DelayedArray_0.12.0         Rgraphviz_2.30.0            DEoptimR_1.0-8              scales_1.1.0               
# [57] mvtnorm_1.0-11              DBI_1.0.0                   GGally_1.4.0                Rcpp_1.0.3                  xtable_1.8-4                progress_1.2.2              htmlTable_1.13.2           
# [64] foreign_0.8-72              bit_1.1-14                  OrganismDbi_1.28.0          preprocessCore_1.48.0       Formula_1.2-3               stats4_3.6.1                AnnotationForge_1.28.0     
# [71] htmlwidgets_1.5.1           httr_1.4.1                  RColorBrewer_1.1-2          acepack_1.4.1               ff_2.2-14                   pkgconfig_2.0.3             reshape_0.8.8              
# [78] XML_3.98-1.20               R.methodsS3_1.7.1           nnet_7.3-12                 dbplyr_1.4.2                locfit_1.5-9.1              tidyselect_0.2.5            rlang_0.4.2                
# [85] reshape2_1.4.3              AnnotationDbi_1.48.0        munsell_0.5.0               tools_3.6.1                 RSQLite_2.1.3               stringr_1.4.0               knitr_1.26                 
# [92] bit64_0.9-7                 robustbase_0.93-5           oligoClasses_1.48.0         caTools_1.17.1.3            purrr_0.3.3                 AnnotationFilter_1.10.0     RBGL_1.62.1                
# [99] R.oo_1.23.0                 biomaRt_2.42.0              compiler_3.6.1              rstudioapi_0.10             curl_4.3                    affyio_1.56.0               PFAM.db_3.10.0             
# [106] tibble_2.1.3                geneplotter_1.64.0          pcaPP_1.9-73                stringi_1.4.3               GenomicFeatures_1.38.0      lattice_0.20-38             ProtGenerics_1.18.0        
# [113] Matrix_1.2-18               vctrs_0.2.0                 pillar_1.4.2                lifecycle_0.1.0             BiocManager_1.30.10         data.table_1.12.6           bitops_1.0-6               
# [120] rtracklayer_1.46.0          GenomicRanges_1.38.0        R6_2.4.1                    latticeExtra_0.6-28         affy_1.64.0                 hwriter_1.3.2               KernSmooth_2.23-16         
# [127] gridExtra_2.3               IRanges_2.20.1              codetools_0.2-16            dichromat_2.0-0             MASS_7.3-51.4               gtools_3.8.1                assertthat_0.2.1           
# [134] SummarizedExperiment_1.16.0 openssl_1.4.1               DESeq2_1.26.0               Category_2.52.1             ReportingTools_2.26.0       GenomicAlignments_1.22.1    Rsamtools_2.2.1            
# [141] S4Vectors_0.24.1            GenomeInfoDbData_1.2.2      hms_0.5.2                   grid_3.6.1                  rpart_4.1-15                biovizBase_1.34.0           base64enc_0.1-3            