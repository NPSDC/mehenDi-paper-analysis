#### Loading data
suppressPackageStartupMessages(library(TreeSummarizedExperiment))
metaData <- read.delim("/fs/cbcb-lab/rob/students/noor/Uncertainity/real_datasets/GSE100505_EDL_MAST/SRR_Acc_List.txt")
quantDir <- "/fs/cbcb-lab/rob/students/noor/Uncertainity/real_datasets/GSE100505_EDL_MAST/sal_out/mode_gcbias=True/posttype=gibbs_npost=100_tf=100"
samples <- metaData$RUN
files <- file.path(quantDir, samples, "quant.sf")
colData <- cbind(data.frame(files = files, names = samples), condition = as.factor(metaData$TissueName))
saveDir <- "../environment/mouseMuscle"
clustFile <- "/fs/cbcb-lab/rob/students/noor/Uncertainity/real_datasets/GSE100505_EDL_MAST/term_out/mode_gcbias=True/posttype=gibbs_npost=100_tf=100/no_threshold0/cluster_nwk.txt"
termFile <- "/fs/cbcb-lab/rob/students/noor/Uncertainity/real_datasets/GSE100505_EDL_MAST/term_out/mode_gcbias=True/posttype=gibbs_npost=100_tf=100/old/SRR5758630/clusters.txt"

tseCons <- beaveR::buildTSE(treeTermFile = clustFile, coldata = colData)
treeCons <- rowTree(tseCons)
l <- length(treeCons$tip)

seMuscle <- tximeta::tximeta(colData)

save(treeCons, file=file.path(saveDir, "treeCons.RData"))

#### Swish on transcripts
y <- fishpond::scaleInfReps(tseCons[1:l,])
y <- fishpond::labelKeep(y)
metadata(y)[["infRepsScaled"]] <- TRUE
set.seed(10)
y <- fishpond::swish(y, x="condition")
                           
#### Swish on genes
gse <- tximeta::summarizeToGene(seMuscle)
gy <- fishpond::scaleInfReps(gse)
set.seed(10)
gy <- fishpond::labelKeep(gy)
gy <- fishpond::swish(gy, x="condition")
                           
#### Swish on Terminus
suppressPackageStartupMessages(source("../terminus_helper.R"))
yTermThrNS <- tseCons[1:l,]
groupsClust <- parseClustFile(termFile, yTermThrNS)
mInds <- seq(nrow(yTermThrNS) + length(groupsClust))
yAggTermThrNS <- prepTerm(yTermThrNS, mInds, groupsClust)
yTerm <- yAggTermThrNS[-unlist(groupsClust),]
yTerm <- fishpond::scaleInfReps(yTerm)
yTerm <- fishpond::labelKeep(yTerm)
set.seed(10)
yTerm <- fishpond::swish(yTerm, x="condition")
                           
#### Loading features from tree based methods
load(file.path(saveDir, "trenDiRes.RData"))
load(file.path(saveDir, "bouthBrain.RData"))
load(file.path(saveDir, "bSwishCons.RData"))

#### Storing features
detNodes <- list()
detNodes[["Txps"]] <- lapply(c(0.01, 0.05, 0.1), function(x) which(mcols(y)[,"qvalue"] <= x ))
detNodes[["Genes"]] <- lapply(c(0.01, 0.05, 0.1), function(x) which(mcols(gy)[,"qvalue"] <= x ))
detNodes[["Terminus"]] <- lapply(c(0.01, 0.05, 0.1), function(x) {
        nodes <- rownames(yTerm)[which(mcols(yTerm)[,"qvalue"] <= x)]
        match(nodes, names(yAggTermThrNS))
    })
detNodes[["trenDi"]] <- lapply(trenDiRes, function(ta) ta[["candNodes"]])
detNodes[["treeClimbR(N)"]] <- lapply(bSwish, function(sw) sw$output[sw$output$signal.node,][["node"]])
detNodes[["treeClimbR(L)"]] <- lapply(bSwish, function(sw) unlist(phangorn::Descendants(treeCons,sw$output[sw$output$signal.node,][["node"]])))
save(detNodes, file=file.path(saveDir,"detNodes.RData"))

negNodes <- list()
negNodes[["Txps"]] <- lapply(detNodes[["Txps"]], function(nodes) setdiff(seq(l), nodes))
negNodes[["Genes"]] <- lapply(detNodes[["Genes"]], function(nodes) setdiff(seq(nrow(gy)), nodes))
negNodes[["Terminus"]] <- lapply(c(0.01, 0.05, 0.1), function(x) {
        nodes <- rownames(yTerm)[setdiff(seq(nrow(yTerm)),
                                                 which(mcols(yTerm)[,"qvalue"] <= x))]
        match(nodes, rownames(yAggTermThrNS))
    })
negNodes[["trenDi"]] <- lapply(detNodes[["trenDi"]], function(nodes) setdiff(seq(nrow(y)), unlist(phangorn::Descendants(treeCons, nodes))))
negNodes[["treeClimbR(N)"]] <- lapply(detNodes[["treeClimbR(N)"]], function(det) setdiff(seq(nrow(y)), unlist(phangorn::Descendants(treeCons, det))))
negNodes[["treeClimbR(L)"]] <- lapply(detNodes[["treeClimbR(L)"]], function(det) setdiff(seq(nrow(y)), det))
save(negNodes, file=file.path(saveDir,"negNodes.RData"))

which(bouthBrain[[1]][["results.by.node"]][["is.driver"]])
which(bouthBrain[[2]][["results.by.node"]][["is.driver"]])
which(bouthBrain[[3]][["results.by.node"]][["is.driver"]])
