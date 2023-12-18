getTreeDf <- function(tree) {
    nleaves <- length(tree$tip.label)
    maxDepth <- max(node.depth(tree, 2))
    
    df <- data.frame(matrix("unknown",nrow = nleaves, ncol = maxDepth, dimnames=list(c(tree$tip.label), paste("Group", c(1:maxDepth), sep="_"))))
    df[,maxDepth] <- tree$tip.label
    colnames(df)[maxDepth] <- "LEAF"
    des <- lapply(seq(1:nleaves), function(i) phangorn::Ancestors(tree, i))
    df[,1] <- rep("Root", nrow(df))
    for(i in seq_along(des)) {
        l <- length(des[[i]])
        if(l > 1)
            df[i,2:l] <- rev(des[[i]])[2:l]
    }
    df
}
                  
computeAggNodes <- function(tree, nodeID, se_counts, group_inds = NULL) {
    performRowAgg <- function(counts, col_inds = NULL) {
        if(is.null(dim(counts)))
            stop("counts has to be matrix/dataframe")
        if(is.null(rownames(counts)))
            stop("rows must be named")
        
        if(is.null(col_inds))
            return(counts)
        
        else
        {
            
            df <- matrix(0, nrow = nrow(counts), ncol = length(col_inds), dimnames = list(rownames(counts)))
            for(i in seq_along(col_inds))
                df[,i] = rowMeans(counts[,col_inds[[i]]])
            return(df)
        }
    }
    performColAgg <- function(counts, row_inds = NULL)
    {
        # if(is.null(dim(counts)))
        #   stop("counts has to be matrix/dataframe")
        if(is.null(row_inds)) {
          if(is.null(dim(counts)))
            return(matrix(counts, nrow=length(counts), ncol=1))
          return(counts)
        }
            
        if(is.null(names(row_inds)))
            stop("row indexes must be named")
        
        sInd <- F ###if not matrix then column sums
        
        if(dim(counts)[2] == 1)
          sInd <- T
        df <- matrix(0, nrow = length(row_inds), ncol = ncol(counts))
        for(i in seq_along(row_inds)){
     #     print(sInd)
          if(!sInd)
            df[i,] <- colSums(counts[row_inds[[i]],])
          else
            df[i,] <- sum(counts[row_inds[[i]]])
        }
          
        rownames(df) <- names(row_inds)
        colnames(df) <- colnames(counts)
        return(df)
    }
    
    if(!is.numeric(nodeID))
    {
        nodeID <- as.numeric(nodeID)
        if(sum(is.na(nodeID)) > 0)
            stop("Node ids contain a non numeric")
    }
    
    mat <- matrix(0, nrow=0, ncol=ncol(se_counts))
    
    leaves <- which(nodeID <= nrow(se_counts))
    innNodes <- which(nodeID > nrow(se_counts))
    
    lInds <- phangorn::Descendants(tree, nodeID[innNodes], type = "tips")
    names(lInds) <- as.character(nodeID[innNodes])
    ls <- sapply(lInds, length)
    #print(dim(se_counts[nodeID[leaves],]))
    
    if(length(leaves) > 0) {
      mat <- rbind(mat, performColAgg(se_counts[nodeID[leaves],]))
    }
    
    if(length(innNodes) > 0) {
      mat <- rbind(mat, performColAgg(se_counts, lInds))
      #print("ss")
    }
        
    mat <- performRowAgg(mat, group_inds)
  
    return(mat)
}

computeTPFP <- function(tSig, sSig, y, logFC = NULL, tree = NULL, type = "all", pTab = F)
{
  # nodeDf <- bouth_ob$tree@node
  # nodeDf$id[nodeDf$id=="Root"] <- as.character(nrow(y)+1)
  # inds <- match(nodeDf$id, as.character(seq_along(logFC)))
  # pvalues <- bouth_ob$tree@test$pvalue[inds]
    if(type == "all")
    {
      truth <- getTruth(length(logFC), tSig)
      simTruth <- getTruth(length(logFC), sSig)    
    }
    else 
    {
      if(is.null(tree))
        stop("Tree cannot be null")
      tSig <- unique(unlist(Descendants(tree, tSig, type = "tips")))
      sSig <- unique(unlist(Descendants(tree, sSig, type = "tips")))
      truth <- getTruth(nrow(y), tSig)
      simTruth <- getTruth(nrow(y), sSig)
    }
    
    tab <- table(simTruth, truth)
    tpr <- tab[2,2]/colSums(tab)[2]
    fdr <- tab[2,1]/rowSums(tab)[2]
    if(pTab)
      print(tab)
    list(fdr=fdr, tpr=tpr)
}
                  
computeMetOut <- function(detNodes, logFCTrue, negNodes = NULL, nodesLooked = NULL, tree = NULL, onlyDet = F, lfcThresh = 0.067, res = T)
{
  if(onlyDet)
  {
    if(is.null(tree))
      stop("Tree cannot be null")
    if(!is.null(negNodes))
      stop("neg nodes should be null")
    descLeaves <- unlist(Descendants(tree, detNodes, "tips"))
    nodesLooked <- sort(c(setdiff(seq_along(tree$tip.label), descLeaves), detNodes))
  }
  if(!is.null(negNodes))
    nodesLooked <- sort(c(negNodes, detNodes))
  
  logFCTrue <- logFCTrue[nodesLooked]
  allNodes <- seq_along(nodesLooked)
  names(allNodes) <- as.character(nodesLooked)
  sNodes <- allNodes[as.character(detNodes)]
  tpNodes <- which(abs(logFCTrue) >= lfcThresh)
  tp = length(intersect(sNodes, tpNodes))/length(sNodes)
  fp = length(setdiff(sNodes, tpNodes))/length(sNodes)
  
  met <- computeTPFP(tpNodes, sNodes, y = NULL, logFC = logFCTrue, type = "all")
  if(res)
    return(met)
  
  fps <- as.numeric(names(allNodes)[setdiff(sNodes, tpNodes)])
  fns <- as.numeric(names(allNodes)[setdiff(tpNodes, sNodes)])
  tps <- as.numeric(names(allNodes)[tpNodes])
  tns <- as.numeric(names(allNodes)[setdiff(allNodes, tpNodes)])
  return(list(fps = fps, fns = fns, tps = tps, tns = tns))
}
                  
getTruth <- function(n, trueInds)
{
  truth <- rep(0,n)
  truth[trueInds] <- 1
  as.factor(truth)
}
                  
compFdratEffSizes <- function(nodes) {
    effSizes <- sort(abs(mcols(yA)[nodes, "log2FC"]), decreasing = T)
    fdrs <- c()
    ls <- c()
    for(i in seq(effSizes)){
        rNodes <- nodes[abs(mcols(yA)[nodes, "log2FC"]) >= effSizes[i]]
        ls <- c(ls, length(rNodes))
        fdrs <- c(fdrs,sum(abs(logFCNodes[rNodes]) < rootFC)/length(rNodes))
    }
    return(list("fdrs"=fdrs, "lengths"=ls, "effSizes"=effSizes, nodes=rNodes))
}
plotScatPlot <- function(vals, size=20) {
    df <- data.frame(vals)
    pFDR <- ggplot(df, aes(x=effSizes,y=fdrs)) +
        geom_point() +
        xlab("Absolute LFC") +
        ylab("FDR") +
        theme_bw() +
        theme(text=element_text(size=size))
#         ggtitle(paste("Unique nodes obtained at FDR", nFDR, "for", type))
    
    pL <- ggplot(df, aes(x=effSizes,y=lengths)) +
        geom_point() +
        xlab("Absolute LFC") +
        ylab("Number of nodes") +
        theme_bw() +
        theme(text=element_text(size=size))
#         ggtitle(paste("Unique nodes obtained at FDR", nFDR, "for", type))
    
    pp <- ggarrange(pFDR, pL, common.legend = T)
    pp
}