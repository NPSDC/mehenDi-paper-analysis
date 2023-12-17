#### Loading TSE experiment
suppressPackageStartupMessages(library(treeclimbR))
quantDir <- "/fs/cbcb-scratch/jfan03/treeterm-paper/output/seed=1_fc=2:6/post_type=gibbs_nrep=100_tf=100/salmon_quants"
saveDir <- "../../environment/brain_sim_normal"
samples <- as.vector(outer(c(1:6), c(1,2), function(x,y) paste(x,y,sep="_")))
files <- file.path(quantDir, samples, "quant.sf")
coldata <- data.frame(files = files, names = samples, condition = as.factor(rep(c(1,2),each=6)))
clustFile <- "/fs/cbcb-scratch/jfan03/treeterm-paper/output/seed=1_fc=2:6/post_type=gibbs_nrep=100_tf=100/terminus/no_threshold0/cluster_nwk.txt"

tseCons <- beaveR::buildTSE(treeTermFile = clustFile, coldata = coldata)
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

