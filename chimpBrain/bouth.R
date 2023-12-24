#### Loading data
suppressPackageStartupMessages(library(BOUTH))
suppressPackageStartupMessages(library(TreeSummarizedExperiment))
suppressPackageStartupMessages(source("../tree_helper.R"))

saveDir <- "../environment/chimpBrain"
load(file.path(saveDir, "y.RData"))

treeCons <- rowTree(yAll)
bouthTree <- getTreeDf(treeCons)
save(bouthTree, file=file.path(saveDir, "bouthTreeCons.RData"))
pvalues <- mcols(y)[["pvalue"]]
bouthBrain <- list()
print(system.time(bouthBrain[['0.01']] <- bouth(anno.table = bouthTree, pvalue.leaves = pvalues, na.symbol = "unknown", far = 0.01, is.weighted = TRUE)))
save(bouthBrain, file = file.path(saveDir, "bouthBrain.RData"))
print(system.time(bouthBrain[['0.05']] <- bouth(anno.table = bouthTree, pvalue.leaves = pvalues, na.symbol = "unknown", far = 0.05, is.weighted = TRUE)))
save(bouthBrain, file = file.path(saveDir, "bouthBrain.RData"))
print(system.time(bouthBrain[['0.1']] <- bouth(anno.table = bouthTree, pvalue.leaves = pvalues, na.symbol = "unknown", far = 0.1, is.weighted = TRUE)))
save(bouthBrain, file = file.path(saveDir, "bouthBrain.RData"))

which(bouthBrain[[1]][["results.by.node"]][["is.driver"]])
which(bouthBrain[[2]][["results.by.node"]][["is.driver"]])
which(bouthBrain[[3]][["results.by.node"]][["is.driver"]])