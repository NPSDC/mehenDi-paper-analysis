quantDir <- "/fs/cbcb-scratch/jfan03/treeterm-paper/output/seed=1_fc=2:6/post_type=gibbs_nrep=100_tf=100/salmon_quants"
saveDir <- "../../environment/brain_sim_normal"
samples <- as.vector(outer(c(1:6), c(1,2), function(x,y) paste(x,y,sep="_")))
files <- file.path(quantDir, samples, "quant.sf")
coldata <- data.frame(files = files, names = samples, condition = as.factor(rep(c(1,2),each=6)))
clustFile <- "/fs/cbcb-scratch/jfan03/treeterm-paper/output/seed=1_fc=2:6/post_type=gibbs_nrep=100_tf=100/terminus/no_threshold0/cluster_nwk.txt"
tseCons <- beaveR::buildTSE(treeTermFile = clustFile, coldata = coldata)
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
mirvThresh <- c(seq(0, 0.35, 0.05), seq(0.50, 1, 0.1))
trenDiResPar <- vector(mode="list", length(minPs) + length(mirvThresh))
names(trenDiResPar) <- c(paste("minP", minPs,sep="="), paste("mirvThresh", mirvThresh,sep="="))

for(i in seq_along(minPs)) {
	    trenDiResPar[[i]] <- lapply(alphas, function(alpha) trenDi::trenDi(yAll, x="condition", pvalues=pvals, 
									               minP=minPs[i], alpha=alpha))
    
}
ll <-length(minPs)
for(i in seq_along(mirvThresh)) {
	    j <- ll+i
    trenDiResPar[[j]] <- lapply(alphas, function(alpha) trenDi::trenDi(yAll, x="condition", pvalues=pvals, 
								               minP=0.7, alpha=alpha, mIRVThresh=mirvThresh[i]))
        
}
save(trenDiResPar, file=file.path(saveDir,"trenDiResPar.RData"))
