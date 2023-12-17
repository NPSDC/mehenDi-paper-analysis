#### Loading data
suppressPackageStartupMessages(library(TreeSummarizedExperiment))
quantDir <- "/fs/cbcb-scratch/jfan03/treeterm-paper/output/seed=1_fc=2:6/post_type=gibbs_nrep=100_tf=100/salmon_quants"
saveDir <- "../../environment/brain_sim_normal"
samples <- as.vector(outer(c(1:6), c(1,2), function(x,y) paste(x,y,sep="_")))
                           
files <- file.path(quantDir, samples, "quant.sf")
coldata <- data.frame(files = files, names = samples, condition = as.factor(rep(c(1,2),each=6)))
clustFile <- "/fs/cbcb-scratch/jfan03/treeterm-paper/output/seed=1_fc=2:6/post_type=gibbs_nrep=100_tf=100/terminus/no_threshold0/cluster_nwk.txt"
tseCons <- beaveR::buildTSE(treeTermFile = clustFile, coldata = coldata)
treeCons <- rowTree(tseCons)
l <- length(treeCons$tip)
termFile <- "/fs/cbcb-scratch/jfan03/treeterm-paper/output/seed=1_fc=2:6/post_type=gibbs_nrep=100_tf=100/terminus/old/1_1/clusters.txt"
save(treeCons, file=file.path(saveDir, "treeCons.RData"))

#### Swish on transcripts
y <- fishpond::scaleInfReps(tseCons[1:l,])
y <- fishpond::labelKeep(y)
metadata(y)[["infRepsScaled"]] <- TRUE
set.seed(10)
y <- fishpond::swish(y, x="condition")
                           
#### Swish on genes
gse <- tximeta::summarizeToGene(tseCons[1:l,])
gy <- fishpond::scaleInfReps(gse)
set.seed(10)
gy <- fishpond::labelKeep(gy)
gy <- fishpond::swish(gy, x="condition")
                           
#### Swish on Terminus
suppressPackageStartupMessages(source("../../terminus_helper.R"))
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
