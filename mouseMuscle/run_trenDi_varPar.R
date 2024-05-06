#### Loading data
metaData <- read.delim("/fs/cbcb-lab/rob/students/noor/Uncertainity/real_datasets/GSE100505_EDL_MAST/SRR_Acc_List.txt")
quantDir <- "/fs/cbcb-lab/rob/students/noor/Uncertainity/real_datasets/GSE100505_EDL_MAST/sal_out/mode_gcbias=True/posttype=gibbs_npost=100_tf=100"
samples <- metaData$RUN
files <- file.path(quantDir, samples, "quant.sf")
colData <- cbind(data.frame(files = files, names = samples), condition = as.factor(metaData$TissueName))

saveDir <- "../environment/mouseMuscle"
clustFile <- "/fs/cbcb-lab/rob/students/noor/Uncertainity/real_datasets/GSE100505_EDL_MAST/term_out/mode_gcbias=True/posttype=gibbs_npost=100_tf=100/no_threshold0/cluster_nwk.txt"

tseCons <- beaveR::buildTSE(treeTermFile = clustFile, coldata = colData)
treeCons <- TreeSummarizedExperiment::rowTree(tseCons)
l <- length(treeCons$tip)

yAll <- beaveR::computeSizeFactors(tseCons)
yAll <- beaveR::scInfReps(yAll)
yAll <- fishpond::labelKeep(yAll)
metadata(yAll)[["infRepsScaled"]] <- TRUE

#### Running trenDi
set.seed(10)
yTxps <- fishpond::swish(yAll[1:l,], x = "condition")
print(sum(mcols(yTxps)[["qvalue"]] <= 0.1, na.rm = T))

set.seed(10)
yInn <- fishpond::swish(yAll[(l+1):nrow(yAll),], x = "condition")

pvals <- c(mcols(yTxps)[["pvalue"]], mcols(yInn)[["pvalue"]])


alphas <- c(0.01, 0.05, 0.1)
minPs <- c(0.60, 0.65, seq(0.75, 0.90, 0.05))
mirvThresh <- c(seq(0.20, 0.35, 0.05), seq(0.50, 1, 0.1))
trenDiResPar <- vector(mode="list", length(minPs) + length(mirvThresh))
names(trenDiResPar) <- c(paste("minP", minPs,sep="="), paste("mirvThresh", mirvThresh,sep="="))

for(i in seq_along(minPs)) {
	    trenDiResPar[[i]] <- lapply(alphas, function(alpha) trenDi::trenDi(yAll, x="condition", pvalues=pvals, 
									               minP=minPs[i], alpha=alpha, cores=4))
    
}
save(trenDiResPar, file=file.path(saveDir,"trenDiResPar.RData"))
ll <-length(minPs)
for(i in seq_along(mirvThresh)) {
	    j <- ll+i
    trenDiResPar[[j]] <- lapply(alphas, function(alpha) trenDi::trenDi(yAll, x="condition", pvalues=pvals, 
								               minP=0.7, alpha=alpha, mIRVThresh=mirvThresh[i], cores=4))
        
}
save(trenDiResPar, file=file.path(saveDir,"trenDiResPar.RData"))
