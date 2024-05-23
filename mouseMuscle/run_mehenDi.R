#### Loading data
suppressPackageStartupMessages(library(pracma))

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

#### Running mehenDi
set.seed(10)
yTxps <- fishpond::swish(yAll[1:l,], x = "condition")
print(sum(mcols(yTxps)[["qvalue"]] <= 0.1, na.rm = T))

set.seed(10)
yInn <- fishpond::swish(yAll[(l+1):nrow(yAll),], x = "condition")

pvals <- c(mcols(yTxps)[["pvalue"]], mcols(yInn)[["pvalue"]])

mehenDiRes <- list()
alphas <- c(0.01, 0.05, 0.1)
for(i in seq_along(alphas)) {
    tic()
    print(system.time(mehenDiRes[[i]] <- mehenDi::mehenDi(yAll, x="condition", pvalues=pvals, 
        minP=0.7, alpha=alphas[i])))
    toc()
}
save(mehenDiRes, file=file.path(saveDir,"mehenDiRes.RData"))

alphas <- c(0.01, 0.05, 0.1)
for(i in seq_along(alphas)) {
    tic()
    print(system.time(ss <- mehenDi::mehenDi(yAll, x="condition", pvalues=pvals, 
        minP=0.7, alpha=alphas[i], cores = 4)))
    toc()
}
