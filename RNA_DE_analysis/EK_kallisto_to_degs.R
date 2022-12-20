source("http://bioconductor.org/biocLite.R")
#biocLite("rhdf5")
#biocLite("edgeR", suppressUpdates=TRUE)
#if (!requireNamespace("BiocManager", quietly = TRUE))
#    install.packages("BiocManager")

#BiocManager::install("edgeR")

library(edgeR);
#install.packages("gplots")
library(gplots);
#install.packages("RColorBrewer")
library(RColorBrewer);
#biocLite("tximport", suppressUpdates=TRUE)
library(tximport);

install.packages("statmod")
library(statmod)

###if (!requireNamespace("BiocManager", quietly = TRUE))
###    install.packages("BiocManager")
biocLite("GO.db", suppressUpdates=TRUE)
###BiocManager::install("GO.db")
library(GO.db)

###if (!requireNamespace("BiocManager", quietly = TRUE))
###    install.packages("BiocManager")
###BiocManager::install("org.Hs.eg.db")
biocLite("org.Hs.eg.db", suppressUpdates=TRUE)
library( org.Hs.eg.db) 
# takes three arguments - config file, transcript to gene table, and output directory

# config file specifies the samples to import, groupings, and paths to abundance.tsv files from kallisto
# groups should be either 0 or 1
# header: sample \t group \t /path/to/abundance.tsv


##---------------------------------------------------------------------------------
## Functions

## make a qc plot
mdsPlot <- function(y,file="mdsplot.pdf"){
  pdf(file=file,width=6,height=6,useDingbats=FALSE);
  plotMDS(y, top=1000, method="bcv", cex=0.5);
  ##legend("bottomleft",as.character(unique(y$samples$group)),col=1:3,pch=20);
  dev <- dev.off();
}

## hierarchical clustering/heatmap
doClustering <- function(data, type="top", num=1000, outviz="heatmap.pdf", outfile="clustered.tsv"){
  rmean = rowMeans(data); # row mean
  rsd = apply(data,1,sd); # row SD
  
  ## do row-normalization:
  data.z = (data-rmean)/rsd;
  data = data.z;
  
  measurements = length(data[1,]);

  ## turn data into numerical matrix:
  data.mat <- data.matrix(data);

  ## create heatmap
  col.sep=seq(0,length(data.mat[1,]),1);
  row.sep=seq(0,length(data.mat[,1]),1);

  ##blue:
  blues = colorRampPalette(c("white","blue"))(n=49);
  blue.breaks = c(seq(0,500,length=50));
  
  ##red-green:
  colors = unique(c(seq(-3,-0.5,length=27),seq(-0.5,0.5,length=25),seq(0.5,3,length=26)));

  if(type=="all"){
    ## huge heatmap of everything  - mem intensive
    pdf(outviz,width=8,height=10,useDingbats=FALSE);
    #alternately: distfun=function(x) as.dist(1-abs(cor(t(x))))
    #or           distfun=as.dist(1-cor(t(x)))
    heat.out <- heatmap.2(data.mat, lhei=c(1,4), lwid=c(4,16),
              key=TRUE, key.xlab="Count",trace="none",density.info="none",dendrogram="col",Rowv=TRUE,
              Colv=TRUE,col=bluered(75),breaks=colors,symbreaks=FALSE,labCol=names(data.mat),labRow=names(data.mat),
              cexCol=0.7,scale="none", margins=c(12,2));
    dev.off();
    ordered.data <- data.mat[heat.out$rowInd,heat.out$colInd]; # z-transformed data as ordered in the heatmap
    write.table(file=outfile, ordered.data, sep="\t", quote=FALSE,col.names=NA);
    return()
;
  }

  if(type=="top"){
    ##heatmap of top N differential genes
    rsd.topN = sort(rsd, decreasing=TRUE)[num];
    data.topN = data.z[which(rsd >= rsd.topN),]
    data.mat <- data.matrix(data.topN);
  
    pdf(outviz,width=8,height=10,useDingbats=FALSE);

    heat.out <- heatmap.2(data.mat, lhei=c(1,4), lwid=c(4,16),
                          key=TRUE, key.xlab="Count",trace="none",density.info="none",dendrogram="col",
                          Rowv=TRUE,Colv=TRUE,col=bluered(75), breaks=colors, symbreaks=FALSE,
                          labCol=names(data.mat),labRow=names(data.mat), cexCol=0.7, cexRow=0.1, scale="none", margins=c(12,2));
    ##                    was sepwidth=0.01,0.01, margins=10,10 for genewise plot
    dev.off();
    ordered.data <- data.mat[heat.out$rowInd,heat.out$colInd]; # z-transformed data as ordered in the heatmap
    write.table(file=outfile, ordered.data, sep="\t", quote=FALSE,col.names=NA);
    return();
  }

## ##do the same thing with euclidean distance
## pdf("heatmap.top1000.euclidian.pdf",width=8,height=24,useDingbats=FALSE);
## heat.out <- heatmap.2(data.mat, lhei=c(0.5,9.5), lwid=c(3,10), key=TRUE, key.xlab="Count",trace="none",
##                       density.info="none", dendrogram="col", Rowv=TRUE, Colv=TRUE, col=greenred(75),
##                       breaks=colors, symbreaks=FALSE, labCol=names(data.mat), labRow="", cexCol=0.7,
##                       scale="none", margins=c(15,5));# was sepwidth=0.01,0.01, margins=10,10 for genewise plot

## ordered.data <- data.mat[heat.out$rowInd,heat.out$colInd]; # z-transformed data as ordered in the heatmap
## write.table(file="out.lengthScaledTPM.counts.normalized.ordered1.tsv", ordered.data, sep="\t", quote=FALSE);
## dev.off();

  stop(paste0("unknown type ",type," for heatmap"))
}

##identify diff expr genes
edgeRclassic <- function(comparison, y, counts.cpm, outdir, tumor_group, type_group) { # “comparison” is just a label/name given to the comparison, y and counts.cpm are defined above
  ###print(counts.cpm[,1])
  print("below warning about 'Design matrix not provided' is expected")
  print("should not have design matrix warning when accounting for batch effects with an additive linear model") ###EK
  print("Let's just see what following the edgeR guide looks like first...")
  ###design <- model.matrix(~type_group+tumor_group)
  design <- model.matrix(~tumor_group+type_group)
  print("display the design matrix")
  print(design)
  logFC <- predFC(y,design,prior.count=1,dispersion=0.05)
  ###print("print the logFC")
  ###print(logFC)
  rownames(design) <-colnames(y)
  print("design after changing colnames")
  print(design)
  ###y.classic <- estimateDisp(y); # this performs dispersion (common and tagwise (ie gene wise)) estimation in the classic version of EdgeR
  #####y.classic <- estimateDisp(y, design, robust=TRUE); ###EK
  y.classic <- estimateDisp(y, design)
  y.classic_common <- estimateGLMCommonDisp(y, design)
  ###y.classic <- estimateDisp(y)
  ###print(counts.cpm)
  #####y.classic <- removeBatchEffect(y, batch=type_group)
  ###print(y.classic) 
  print("display common dispersion")
  print(y.classic$common.dispersion)
  print("first common dispersion done. next...")
  print(y.classic_common)
  print("done displaying glm common dispersion")
  fit <- glmQLFit(y.classic, design)
  et <- exactTest(y.classic); # recall that y.classic already has the group information
  qlf <- glmQLFTest(fit)
  print("done making the glmqlf test")
  print(topTags(qlf))
  print("done printing the top tags")
  top <-rownames(topTags(qlf))
  print(cpm(y.classic)[top,])
  print("done printing the cpm")
  print(summary(decideTests(qlf)))
  print("finished the toptags part")
  ###summary(decideTests(qlf)) ###EK addition
  #########total = length(qlf@.Data[[1]]$PValue);
  print("print the qlf total part")
  total_qlf<-length(qlf@.Data[[1]]);
  print(total_qlf)
  print("print the et total part")
  total = length(et@.Data[[1]]$PValue);
  print(total)
  print("print the total et part without the PValue")
  total_et = length(et@.Data[[1]]);
  print(total_et)
  #####results <- topTags(et, adjust.method="BH", sort.by="PValue", n=total); # table of the top differntially-expressed "tags" ie "genes"
  results <- topTags(qlf, adjust.method="BH", sort.by="PValue", n=total_qlf); # table of the top differntially-expressed "tags" ie "genes"
  ### GO terms - EK ###
  ###print("print the lrt")
  ###lrt <- glmLRT(fit)
  ###topTags(lrt)
  ###geneids <- counts.cpm[,1]
  ###go <- goana.DGELRT(lrt, de=results)
  ###go <- goana(lrt, geneid=geneids)
  ###print("GO Up")
  ###topGo(go, ont="BP", sort="Up", n=25, truncate=25)
  ###print("GO down")
  ###topGo(go, ont="BP", sort="Down", n=25, truncate=25)

  results.table <- results@.Data[[1]];
  results.pos <- results.table[which(results.table$logFC >= 0 & results.table$FDR <= 0.05),];
  results.neg <- results.table[which(results.table$logFC < 0 & results.table$FDR<= 0.05),];
  data.pos <- counts.cpm[rownames(results.pos),] # note that this uses x.0.wt.cpm
  data.neg <- counts.cpm[rownames(results.neg),]
  write.table(results@.Data[[1]],file=sprintf("%s/DEGs.%s.tsv", outdir, comparison), quote=FALSE,sep='\t',row.names=TRUE, col.names=NA);
  write.table(results.pos,file=sprintf("%s/DEGs.%s.pos.tsv", outdir, comparison), quote=FALSE,sep='\t',row.names=TRUE, col.names=NA);
  write.table(results.neg,file=sprintf("%s/DEGs.%s.neg.tsv", outdir, comparison), quote=FALSE,sep='\t',row.names=TRUE, col.names=NA);
  write.table(data.pos,file=sprintf("%s/DEGs.%s.pos.CPM.tsv", outdir, comparison), quote=FALSE,sep='\t',row.names=TRUE, col.names=NA);
  write.table(data.neg,file=sprintf("%s/DEGs.%s.neg.CPM.tsv", outdir, comparison), quote=FALSE,sep='\t',row.names=TRUE, col.names=NA);
  

  ###EK addition
  pdf(file="All_Dispersion_BCV.pdf") ###EK attempt to write BCV plot to pdf file
  plotBCV(y.classic) ###EK added plot to display dispersion output
  plotMD(qlf) ###EK plot log-fold change against log-counts per million, with DE gene highlighted
  ###abline(h=c(-1,1), col="blue")
  plotQLDisp(fit)
  dev.off(); ###EK end of pdf file writing?

  result = list(table=results, pos=data.pos, neg=data.neg);
  return(result);
}

#make groups vector
getGroups <- function(n,config){
  groups=c();
  for(i in n){
    g = config[config$sample==i,]$group
    if(is.na(g) | g==""){
      stop(paste0("sample name ",i," not found in config file"))
    } else {
      groups=c(groups,g)
    }
  }
  return(groups)
}
 
###make second group  vector (for FF vs FFPE samples)
getGroups2 <- function(n,config){
  groups=c();
  for(i in n){
    g = config[config$sample==i,]$type
    if(is.na(g) | g==""){
      stop(paste0("sample name ",i," not found in config file"))
    } else {
      groups=c(groups,g)
    }
  }
  return(groups)
}


#create matrix from the given files
createMatrix <- function(config, conv.table){
  ## list of kallisto output files to be read as input
  k=config$path
  names(k) = config$sample

  # transcript-gene conversions for ENSEMBL annotations
  tx2gene=read.table(conv.table, sep='\t',header=TRUE,quote='',check.names=FALSE, stringsAsFactors=F);  

  writeTables <- function(txi,suffix){
    write.table(txi$abundance, file=paste0(outdir,"/matrix.abundance.tsv"), sep='\t', quote=FALSE);
    write.table(txi$counts, file=paste0(outdir,"/matrix.counts.tsv"), sep='\t', quote=FALSE);
    write.table(txi$length, file=paste0(outdir,"/matrix.length.tsv"), sep='\t', quote=FALSE);
  }

  # length-scaled genewise counts from abundance:
  txi <- tximport(k, type="kallisto", tx2gene=tx2gene, countsFromAbundance="lengthScaledTPM",ignoreTxVersion=T);
  print("txi created")
  writeTables(txi,"lengthScaledTPM")

  ## # scaled genewise counts from abundance:
  ## txi <- tximport(files, type="kallisto", tx2gene=tx2gene, countsFromAbundance="scaledTPM");
  ## writeTables(txi,"scaledTPM")
  ## # transcript-wise counts from counts:
  ## txi <- tximport(files, type="kallisto", tx2gene=tx2gene, countsFromAbundance="no", txOut=TRUE); # output compiled transcript counts from kallisto?
  ## writeTables(txi,"txCount")
  ## # gene-wise counts from counts:
  ## txi <- tximport(files, type="kallisto", tx2gene=tx2gene, countsFromAbundance="no"); # output gene-level counts from kallisto without correction/normalization
  ## writeTables(txi,"geneCount")
  return(paste0(outdir,"/matrix.counts.tsv"))
}

##------------------------------------------------------
args <- commandArgs(trailingOnly = TRUE)
config=read.table(args[1], header=F, stringsAsFactors=F, sep="\t")
###names(config) = c("sample","group","path")
names(config) = c("sample", "group","type","path")

#create the output directory
outdir=args[3];
dir.create(outdir)

mat = createMatrix(config,args[2])
###print("This is the mat line")
###mat
###print("This is args 1")
###args[1]
###print("This is args 2")
###args[2]
###print("commented out the rest of the code below this line")


#read in the counts
print("reading count data")
rawcounts <- read.table(mat, sep='\t', quote="", check.names=FALSE,header=TRUE,row.names=1)
#subset to only those in the config file
y <- DGEList(counts=rawcounts);
###print("print the DGEList")
###y

###print("everything below this line is commented out")
###comment <-function() {
print("filtering/normalizing...")
## filter here for low expression, requiring counts-per-million (CPM) be at least 1 in at least half of the samples, and recalculate library sizes:
min.thresh=ceiling(length(names(rawcounts))/2);
keep <- rowSums(cpm(y)>1) >= min.thresh;
y <- y[keep, , keep.lib.sizes=FALSE]; # forces recalculation of library sizes

#normalize the data
y <- calcNormFactors(y);
counts.cpm <- cpm(y, normalized.lib.sizes=TRUE)
write.table(counts.cpm,file=paste0(outdir,"/counts.normalized.tsv"),quote=FALSE, sep='\t',col.names=NA);
###print("printing the normalized y value")
###y
###print("print the cpm counts")
###counts.cpm


###print("everything below this line is commented out")
###comment <-function() {
print("creating mdsplot...")
##make mds plot for QC
mdsPlot(y,file=paste0(outdir,"/mdsplot.pdf"));

###print("skip clustering for now...") ###EK added and uncommented below
print("doing hierarchical clustering/creating heatmap...")
doClustering(counts.cpm, outviz=paste0(outdir,"/heatmap.pdf"), outfile=paste0(outdir,"/heatmap.clustered.top1000.tsv"), num=1000, type="top");
print("clustering done")

save.image("EdgeR_STvLT.Rdata")


###_______________left off here, was able to grab sample types like tumor group but for FF and FFPE _______________________
###print("everything below this line is commented out")
###comment <-function() {
##now, before we calculate diff expr, let's subset to just the samples we actually want to compare
print("subsetting data")
counts.cpm.comp = counts.cpm[,colnames(counts.cpm) %in% config$sample]
###print("printing counts.cpm.comp")
###counts.cpm.comp
y2 <- DGEList(counts=counts.cpm.comp, group=getGroups(colnames(counts.cpm.comp),config));
mdsPlot(y2,file=paste0(outdir,"/mdsplot_2.pdf"));
tumor_group <- getGroups(colnames(counts.cpm.comp),config)
###print("print out the tumor group")
###tumor_group
###print("factoring test")
tumors <- factor(tumor_group)
tumors
types <-getGroups2(colnames(counts.cpm.comp),config)
###print("print attempt to get sample type FF or FFPE")
###types
types_f <- factor(types)
types_f

###print("print y2 before DE, the subset samples")
###y2

###print("before we run the edgeRclassic function, can I make a design matrix outside the function....")
###design <- model.matrix(~types_f+types_f:tumors)
###design
###design2 <-model.matrix(~types_f+tumors)
###print("print design as it should be within the function")
###design2
###rownames(design2) <- colnames(y2)
###design2
###logFC <- predFC(y,design,prior.count=1,dispersion=0.05)
###logFC
###cor(logFC[,2:4])

###design2



###y.classic <- estimateDisp(y2, design2, robust=TRUE); ###EK
###print("display common dispersion")
###y.classic$common.dispersion
###et <- exactTest(y.classic); # recall that y.classic already has the group information
###summary(decideTests(et))
###pdf(file="Dispersion_BCV.pdf") ###EK attempt to write BCV plot to pdf file
###plotBCV(y.classic) ###EK added plot to display dispersion output
###plotMD(et) ###EK plot log-fold change against log-counts per million, with DE gene highl
###dev.off()
###print("just seeing what the glmQLFTest from edgeR example looks like...")
###fit <- glmQLFit(y.classic, design, robust=TRUE)
###qlf <-glmQLFTest(fit, coef=1:2)
###topTags(qlf)
###FDR <- p.adjust(qlf$table$PValue, method="BH")
###sum(FDR <0.05)
###qlf2 <-glmQLFTest(fit)
###topTags(qlf2)
###top <- rownames(topTags(qlf2))
###cpm(y.classic)[top,]
###summary(decideTests(qlf))

###pdf(file="glm_Dispersion_BCV.pdf")	
###plotBCV(y.classic) ###EK added plot to display dispersion output
###plotMD(qlf2) ###EK plot log-fold change against log-counts per million, with DE gene highlighted
  ###abline(h=c(-1,1), col="blue")
###plotQLDisp(fit)
###dev.off();

###comment <-function() {
print("calculating diff expression...")
results = edgeRclassic("comp1", y2, counts.cpm.comp, outdir, tumors, types_f) 

###print("everything below this line is commented out")
###comment <-function() {

#make a heatmap of just the differentially expressed genes across all samples
#grab the DE genes
results.table=results$table@.Data[[1]];
results.sig <- results.table[which(results.table$FDR <= 0.05),];
counts.cpm.de <- counts.cpm[which(rownames(counts.cpm) %in% rownames(results.sig)),] # note that this uses x.0.wt.cpm
write.table(counts.cpm.de,file=paste0(outdir,"/DEG.cpm.tsv"),quote=FALSE,sep='\t',row.names=TRUE, col.names=NA);

z = merge(counts.cpm.de,results.sig,by=0)
write.table(z,file=paste0(outdir,"/DEG.cpm.pvals.tsv"),quote=FALSE,sep='\t',row.names=TRUE, col.names=NA);

print("doing hierarchical clustering/creating heatmap...")
doClustering(counts.cpm.de, outviz=paste0(outdir,"/heatmap.degenes.pdf"), outfile=paste0(outdir,"/heatmap.clustered.de.tsv"), type="all");


