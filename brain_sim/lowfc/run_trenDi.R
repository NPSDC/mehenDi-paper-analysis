#### Loading data
suppressPackageStartupMessages(library(pracma))
quantDir <- "/fs/cbcb-scratch/jfan03/treeterm-paper/output/seed=1_fc=1.4:2.8/post_type=gibbs_nrep=100_tf=100/salmon_quants"
saveDir <- "../../environment/brain_sim_lowfc"
samples <- as.vector(outer(c(1:6), c(1,2), function(x,y) paste(x,y,sep="_")))
files <- file.path(quantDir, samples, "quant.sf")
coldata <- data.frame(files = files, names = samples, condition = as.factor(rep(c(1,2),each=6)))
clustFile <- "/fs/cbcb-scratch/jfan03/treeterm-paper/output/seed=1_fc=1.4:2.8/post_type=gibbs_nrep=100_tf=100/terminus/no_threshold0/cluster_nwk.txt"
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
for(i in seq_along(alphas)) {
    tic()
    print(system.time(trenDiRes <- trenDi::trenDi(yAll, x="condition", pvalues=pvals, 
        minP=0.7, alpha=alphas[i])))
    toc()
}
save(trenDiRes, file=file.path(saveDir,"trenDiRes.RData"))

alphas <- c(0.01, 0.05, 0.1)
for(i in seq_along(alphas)) {
    tic()
    print(system.time(ss <- trenDi::trenDi(yAll, x="condition", pvalues=pvals, 
        minP=0.7, alpha=alphas[i], cores = 4)))
    toc()
}
