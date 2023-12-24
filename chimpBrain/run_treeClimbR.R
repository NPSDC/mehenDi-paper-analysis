#### Loading TSE experiment
suppressPackageStartupMessages(library(treeclimbR))
suppressPackageStartupMessages(library(treeSummarizedExperiment))
saveDir <- "../environment/chimpBrain"
load(file.path(saveDir, "yAll.RData"))

treeCons <- rowTree(yAll)
l <- length(treeCons$tip)

pvals <- mcols(yAll)[["pvalue"]]

### Running treeclimbR
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

#### Saving features
load(file.path(saveDir,"detNodes.RData"))
load(file.path(saveDir,"negNodes.RData"))

detNodes[["treeClimbR(N)"]] <- lapply(bSwish, function(sw) sw$output[sw$output$signal.node,][["node"]])
detNodes[["treeClimbR(L)"]] <- lapply(bSwish, function(sw) unlist(phangorn::Descendants(treeCons,sw$output[sw$output$signal.node,][["node"]])))
negNodes[["treeClimbR(N)"]] <- lapply(detNodes[["treeClimbR(N)"]], function(det) setdiff(seq(l), unlist(phangorn::Descendants(treeCons, det))))
negNodes[["treeClimbR(L)"]] <- lapply(detNodes[["treeClimbR(L)"]], function(det) setdiff(seq(l), det))

save(file.path(saveDir,"detNodes.RData"))
save(file.path(saveDir,"negNodes.RData"))