#### Loading data
suppressPackageStartupMessages(library(BOUTH))
suppressPackageStartupMessages(library(SummarizedExperiment))
suppressPackageStartupMessages(source("../tree_helper.R"))

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

y <- fishpond::scaleInfReps(tseCons[1:l,])
y <- fishpond::labelKeep(y)
metadata(y)[["infRepsScaled"]] <- TRUE
mcols(y)[["keep"]] <- T ## BOUTH doesn't handle NA p-values
set.seed(10)
y <- fishpond::swish(y, x="condition")

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