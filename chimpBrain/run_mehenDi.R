#### Loading data
suppressPackageStartupMessages(library(pracma))
suppressPackageStartupMessages(library(treeSummarizedExperiment))

saveDir <- "../environment/chimpBrain"
clustFile <- "/fs/cbcb-lab/rob/students/noor/Uncertainity/real_datasets/GSE100505_EDL_MAST/term_out/mode_gcbias=True/posttype=gibbs_npost=100_tf=100/no_threshold0/cluster_nwk.txt"
load(file.path(saveDir, "yAll.RData"))

treeCons <- rowTree(yAll)
l <- length(treeCons$tip)

#### Running mehenDi
alphas <- c(0.01, 0.05, 0.1)
mehenDiRes <- list()
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
load(file.path(saveDir,"detNodes.RData"))
load(file.path(saveDir,"negNodes.RData"))

detNodes[["mehenDi"]] <- lapply(mehenDiRes, function(ta) ta[["candNodes"]])
negNodes[["mehenDi"]] <- lapply(detNodes[["mehenDi"]], function(nodes) setdiff(seq(nrow(y)), unlist(phangorn::Descendants(treeCons, nodes))))
                               
save(file.path(saveDir,"detNodes.RData"))
save(file.path(saveDir,"negNodes.RData"))