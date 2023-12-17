#### Loading data
suppressPackageStartupMessages(library(BOUTH))
suppressPackageStartupMessages(library(SummarizedExperiment))
suppressPackageStartupMessages(source("../../tree_helper.R"))

quantDir <- "/fs/cbcb-scratch/jfan03/treeterm-paper/output/seed=1_fc=2:6/post_type=gibbs_nrep=100_tf=100/salmon_quants"
saveDir <- "../../environment/brain_sim_normal"
samples <- as.vector(outer(c(1:6), c(1,2), function(x,y) paste(x,y,sep="_")))
files <- file.path(quantDir, samples, "quant.sf")
coldata <- data.frame(files = files, names = samples, condition = as.factor(rep(c(1,2),each=6)))

tseCons <- beaveR::buildTSE(treeTermFile = clustFile, coldata = coldata)
treeCons <- rowTree(tseCons)
l <- length(treeCons$tip)

y <- fishpond::scaleInfReps(tseCons[1:l,])
y <- fishpond::labelKeep(y)
metadata(y)[["infRepsScaled"]] <- TRUE
mcols(y)[["keep"]] <- T ## BOUTH doesn't handle NA p-values
set.seed(10)
y <- fishpond::swish(y, x="condition")

bouthTree <- getTreeDf(treeCons)
# save(bouthTree, file=file.path(saveDir, "bouthTreeCons.RData"))
pvalues <- mcols(y)[["pvalue"]]
bouthBrain <- list()
print(system.time(bouthBrain[['0.01']] <- bouth(anno.table = bouthTree, pvalue.leaves = pvalues, na.symbol = "unknown", far = 0.01, is.weighted = TRUE)))
save(bouthBrain, file = file.path(saveDir, "bouthBrain.RData"))
print(system.time(bouthBrain[['0.05']] <- bouth(anno.table = bouthTree, pvalue.leaves = pvalues, na.symbol = "unknown", far = 0.05, is.weighted = TRUE)))
save(bouthBrain, file = file.path(saveDir, "bouthBrain.RData"))
print(system.time(bouthBrain[['0.1']] <- bouth(anno.table = bouthTree, pvalue.leaves = pvalues, na.symbol = "unknown", far = 0.1, is.weighted = TRUE)))
save(bouthBrain, file = file.path(saveDir, "bouthBrain.RData"))