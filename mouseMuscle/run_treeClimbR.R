#### Loading TSE experiment
suppressPackageStartupMessages(library(treeclimbR))
metaData <- read.delim("/fs/cbcb-lab/rob/students/noor/Uncertainity/real_datasets/GSE100505_EDL_MAST/SRR_Acc_List.txt")
quantDir <- "/fs/cbcb-lab/rob/students/noor/Uncertainity/real_datasets/GSE100505_EDL_MAST/sal_out/mode_gcbias=True/posttype=gibbs_npost=100_tf=100"
samples <- metaData$RUN
files <- file.path(quantDir, samples, "quant.sf")
colData <- cbind(data.frame(files = files, names = samples), condition = as.factor(metaData$TissueName))

saveDir <- "../environment/mouseMuscle"
clustFile <- "/fs/cbcb-lab/rob/students/noor/Uncertainity/real_datasets/GSE100505_EDL_MAST/term_out/mode_gcbias=True/posttype=gibbs_npost=100_tf=100/no_threshold0/cluster_nwk.txt"

tseCons <- beaveR::buildTSE(treeTermFile = clustFile, coldata = colData)
treeCons <- rowTree(tseCons)
l <- length(treeCons$tip)


#### Swish on the entire dataset 
yAll <- beaveR::computeSizeFactors(tseCons)
yAll <- beaveR::scInfReps(yAll)
yAll <- labelKeep(yAll)
metadata(yAll)[["infRepsScaled"]] <- TRUE
set.seed(10)
yAll <- swish(yAll, x = "condition")

set.seed(10)
yTxps <- swish(yAll[1:l,], x = "condition")
print(sum(mcols(yTxps)[["qvalue"]] <= 0.1, na.rm = T))

set.seed(10)
yInn <- swish(yAll[(l+1):nrow(yAll),], x = "condition")

pvals <- c(mcols(yTxps)[["pvalue"]], mcols(yInn)[["pvalue"]])
signs <- computeSign(yAll, "condition")
mIRV <- mcols(yAll)[["meanInfRV"]]

### Running treeclimbR
mcols(yAll)[["pvalue"]] <- c(mcols(yTxps)[["pvalue"]], mcols(yInn)[["pvalue"]])
swishRes <- data.frame(mcols(yAll))
swishRes$node <- seq(1:nrow(yAll))
swishRes <- swishRes[,!colnames(swishRes) %in% c("keep")]

alphas <-  c(0.01,0.05,0.1)

print(system.time(cSwish <- getCand(treeCons,
                 score_data = swishRes,
                 node_column = "node",
                 p_column = "pvalue",
                 sign_column = "log2FC"
)))
save(cSwish, file=file.path(saveDir, "cSwishCons.RData"))
gc()

bSwish <- list()
for(i in seq_along(alphas)) {
print(system.time(bSwish[[i]] <- evalCand(tree = treeCons, levels = cSwish$candidate_list,
                      score_data = swishRes, node_column = "node",
                      p_column = "pvalue", sign_column = "log2FC", limit_rej=alphas[[i]]
)))
}
save(bSwish, file=file.path(saveDir, "bSwishCons.RData"))

