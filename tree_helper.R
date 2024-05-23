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
plotScatPlot <- function(vals,ps=2, size=20,diff=4) {
    df <- data.frame(vals)
    pFDR <- ggplot(df, aes(x=effSizes,y=fdrs)) +
        geom_point(size=ps) +
        xlab("Absolute LFC") +
        ylab("FDR") +
        theme_bw() +
        theme(title =element_text(size=size-diff, margin=unit(c(0, -10, 0, 0), "mm")),
             axis.title = element_text(size = size-diff),
             axis.text = element_text(size = size, , vjust=1))
#         ggtitle(paste("Unique nodes obtained at FDR", nFDR, "for", type))
    
    pL <- ggplot(df, aes(x=effSizes,y=lengths)) +
        geom_point(size=ps) +
        xlab("Absolute LFC") +
        ylab("Number of nodes") +
        theme_bw() +
        theme(title =element_text(size=size-diff, margin=unit(c(0, -10, 0, 0), "mm")),
             axis.title = element_text(size = size-diff),
             axis.text = element_text(size = size, vjust=1))
#         ggtitle(paste("Unique nodes obtained at FDR", nFDR, "for", type))
    
    pp <- ggarrange(pFDR, pL, common.legend = T)
#     pp <- pp + ggtitle("yo")
    pp
    
}

### This function returns the useful elements needed for plot construction, if things can be automated, otherwise we shall start from scratch
#' @ tse - TreeSummarizedExperiment
#' @ tL - tree list
#' @ indlist - for the index set
#' @ gM - for the index set
#' @ i - FDR ind
#' @ j - node index in the index set
extPreInf <- function(tse, y, tL, indList, txpM, gM, i, j) {
    treeCons <- rowTree(tse)
    iNode <- tL[[i]][indList[[i]][[j]]]
    txShow <- rownames(tse)[Descendants(treeCons,iNode)[[1]]]
    print(txShow)
    anc <- Ancestors(treeCons, iNode)
    anc <- ifelse(length(anc)==1, iNode, anc[length(anc)-1])
    treeSub <- tidytree::tree_subset(treeCons, anc, levels_back = 0)
    print(treeSub)
    gs <- txpM %>% 
        filter(tx_name %in% treeSub$tip) %>%
        tibble::as_tibble() %>%
        dplyr::select(ensID) %>%
        unlist %>%
        unique
    print(paste("Genes", gs))
    
    gTxps <- txpM %>% 
        filter(ensID == gs[1]) %>%
        tibble::as_tibble() %>%
        dplyr::select(tx_name) %>%
        unlist
    
    g <- gM %>% 
        filter(ensID == gs[1]) %>%
        tibble::as_tibble()

    print(gTxps)
    print(treeSub$tip)
    print(all(treeSub$tip %in% gTxps))
    print(all(gTxps %in% treeSub$tip))
    txpMin <- gTxps[which.min(mcols(y)[gTxps, "pvalue"])]
    
    if(length(txpMin)==0) {
        txpMin <- gTxps[which.max(abs(mcols(y)[gTxps, "log10mean"]))]
    }
    minTSInd <- match(txpMin, treeSub$tip)
    return(list(tSub=treeSub, txpMin=txpMin, txShow=txShow, minTInd=minTSInd, iNode=iNode, g=g))
}

### This function plots the infRep plot
plotIReps <- function(y, txpMin, iNode, lp="right",x="condition", cex=1.6) {
    legPos = ifelse(lp=="right", "topright", "topleft")
    pTxp <- as.grob(function() plotInfReps(y, txpMin, x = x, legend=TRUE,
                              main=txpMin, legendTitle=TRUE, legendCex=cex,
                             legendPos = legPos, mainCex=cex, axisCex=cex,labCex=cex))
    pInn <- as.grob(function() plotInfReps(y, iNode, x = x, legend=TRUE,
                              main="mehenDi Selected Node", legendTitle=TRUE, legendCex=cex,
                             legendPos = legPos, mainCex=cex, axisCex=cex,labCex=cex))
    return(list(pTxp, pInn))    
}

plotTree <- function(treeSub, iNode, txNode, of=20.5, ofex=4, xlim=NA, size=20,tip.size=7) {
    suppressPackageStartupMessages(library(ggtree))
    suppressPackageStartupMessages(library(ggplot2))
    xx <- ifelse(is.na(xlim), 80, xlim)
    pTree <- ggtree(treeSub) + #ggtitle("Tree representing the transcripts covered by gene ENSMUSG00000070509") +
        xlim(NA, xx) +
        geom_tiplab(size=tip.size, hjust=-0.1) +
        geom_point2(aes(subset=(node==iNode), color="red"), 
                size=5, fill='red', show.legend=T) +
        geom_point2(aes(subset=(node==txNode), color = "black"), 
                size=5, fill='black', show.legend=T) +
        geom_cladelab(node = iNode, label = "", textcolour="red", barsize=2,
                      barcolour="red",  fontsize=5, offset = of) + 
        geom_cladelab(node = txNode, label = "", textcolour="black",
                      fontsize=5, offset = of+ofex) +
        theme(legend.position = "bottom", legend.text=element_text(size=size),
             plot.title=element_text(size=size, face="bold", hjust=0)) +
        scale_color_manual(name = "", labels=c("Transcript with the lowest p-value", "Selected Node"),
               values=c("black", "red"))
    return(pTree)
}

parF <- function(g, txShow, treeSub, chromSt=100, chromEnd=200, fs=14, assemb="mm10") {
    par <- plotgardener::pgParams(
          chrom = as.character(g[["seqnames"]]), 
          chromstart = g[["start"]]-chromSt, chromend = g[["end"]]+chromEnd,
          assembly = assemb, just = c("left", "bottom"), fontsize = fs,
        default.units = "inches"
        )
    
    hilite <- data.frame(transcript=c(txShow, setdiff(treeSub$tip, txShow)), 
                     color=c(rep("red", length(txShow)),
                             rep("blue", length(setdiff(treeSub$tip, txShow)))))
    
    return(list(par, hilite))
}

plotInfReps <- function(y, idx, x, cov=NULL,
                        colsDrk=c("dodgerblue","goldenrod4","royalblue4",
                                  "red3","purple4","darkgreen"),
                        colsLgt=c("lightblue1","goldenrod1","royalblue1",
                                  "salmon1","orchid1","limegreen"),
                        xaxis, xlab, ylim,
                        main, mainCol,
                        legend=FALSE,
                        legendPos="topleft",
                        legendTitle=FALSE,
                        legendCex=1,
                        useMean=TRUE,
                        applySF=FALSE,
                        mainCex=1, 
                        axisCex=1, 
                        labCex=1,
                        reorder,
                        thin) {
  
  hasInfReps <- any(grepl("infRep", assayNames(y)))
  # logical switch for if variance assay is present
  # if not, just show the count or mean
  showVar <- TRUE
  if (!hasInfReps) {
    showVar <- "variance" %in% assayNames(y)
    if (useMean & !("mean" %in% assayNames(y))) {
      message("using 'counts' assay, as 'mean' is missing, see argument 'useMean'")
      useMean <- FALSE
    }
  }
  stopifnot(x %in% names(colData(y)))
  condition <- colData(y)[[x]]
  # whether x is a factor variable or numeric (e.g. pseudotime)
  xfac <- is(condition, "factor")
  if (!xfac) stopifnot(is(condition, "numeric"))
  stopifnot(length(colsDrk) == length(colsLgt))
  if (xfac) {
    ncond <- nlevels(condition)
    stopifnot(ncond <= length(colsDrk))
    colsDrk <- colsDrk[seq_len(ncond)]
    colsLgt <- colsLgt[seq_len(ncond)]
  } else {
    if (hasInfReps) {
      message("boxplot of inf. reps over numeric 'x' not supported yet")
      hasInfReps <- FALSE
      stopifnot("variance" %in% assayNames(y))
    }
  }
  if (missing(xaxis)) {
    if (xfac) {
      xaxis <- ncol(y) < 30
    } else {
      xaxis <- TRUE
    }
  }
  if (missing(thin)) {
    thin <- if (ncol(y) >= 400) 2 else if (ncol(y) >= 150) 1 else 0
  } else {
    stopifnot(thin >= 0 & thin <= 2)
  }
  if (!is.null(cov)) {
    stopifnot(cov %in% names(colData(y)))
    covariate <- factor(colData(y)[[cov]])
    stopifnot(is(covariate, "factor"))
    ngrp <- nlevels(covariate)
  }
  infRepsScaled <- FALSE
  if (!is.null(metadata(y)$infRepsScaled)) {
    infRepsScaled <- metadata(y)$infRepsScaled
  }
  # single cell?
  sc <- FALSE 
  if (!is.null(metadata(y)$tximetaInfo$type)) {
    if (metadata(y)$tximetaInfo$type == "alevin") {
      sc <- TRUE
    }
  }
  if (missing(xlab)) {
    if (xfac) {
      xlab <- if (sc) "cells" else "samples"
    } else {
      xlab <- x
    }
  }
  if (missing(reorder)) {
    if (xfac) {
      reorder <- sc
    } else {
      reorder <- FALSE
    }
  } else {
    if (!xfac & reorder) stop("reorder not used when 'x' is numeric")
  }
  ylab <- if (infRepsScaled) "scaled counts" else "counts"
  # this is a dummy variable used when making the plot()
  # if we don't put x-axis ticks, then we will move
  # the label up higher using title()
  xlabel <- if (xaxis) xlab else ""
  if (missing(main)) {
    if (missing(mainCol)) {
      if (is.character(idx)) {
        main <- idx
      } else {
        main <- rownames(y)[idx]
      }
    } else {
      stopifnot(mainCol %in% names(mcols(y)))
      main <- mcols(y)[idx,mainCol]
    }
  }
  if (reorder) {
    if (hasInfReps) {
      infReps <- assays(y[idx,])[grep("infRep",assayNames(y))]
      value <- colMeans(unlist(infReps))
    } else {
      which.assay <- if (useMean) "mean" else "counts"
      value <- assays(y)[[which.assay]][idx,]
      if (applySF & !is.null(y$sizeFactor)) {
        value <- value/y$sizeFactor
      }
    }
    if (is.null(cov)) {
      o <- order(condition, value)
    } else {
      o <- order(covariate, condition, value)
    }
    # not reordering:
  } else {
    if (xfac) {
      if (is.null(cov)) {
        o <- order(condition)
      } else {
        o <- order(covariate, condition)
      }
      # numeric 'x', don't reorder
    } else {
      o <- seq_along(condition)
    }
  }
  # col - dark color for boxplot border, points
  # col.hglt - highlight color for inside boxplot, lines when n >= 400
  if (xfac) {
    if (is.null(cov)) {
      samp.nums <- unlist(lapply(table(condition), seq_len))
      col <- rep(colsDrk, table(condition))
      col.hglt <- rep(colsLgt, table(condition))
    } else {
      vec.tab <- as.vector(table(condition, covariate))
      samp.nums <- unlist(lapply(vec.tab, seq_len))
      col <- rep(rep(colsDrk, ngrp), vec.tab)
      col.hglt <- rep(rep(colsLgt, ngrp), vec.tab)
    }
  } else {
    if (is.null(cov)) {
      col <- colsDrk[1]
      col.hglt <- colsLgt[1]
    } else {
      col <- colsDrk[covariate]
      col.hglt <- colsLgt[covariate]
    }
  }
  ### boxplot ###
  if (hasInfReps) {
    infReps <- assays(y[idx,])[grep("infRep",assayNames(y))]
    cts <- unlist(infReps)[,o]
    ymax <- max(cts)
    ymin <- if (is.null(cov)) 0 else -0.02 * ymax
    if (missing(ylim)) {
      ylim <- c(ymin,ymax)
    } else {
      stopifnot(length(ylim) == 2)
    }
    boxplot2(cts, col=col, col.hglt=col.hglt, ylim=ylim,
             xlab=xlabel, ylab=ylab, main=main, 
             mainCex=mainCex, axisCex=axisCex, labCex=labCex)
    ### point and line plot by 'x' levels or numeric 'x' ###
  } else {
    which.assay <- if (useMean) "mean" else "counts"
    cts <- assays(y)[[which.assay]][idx,o]
    if (showVar) {
      sds <- sqrt(assays(y)[["variance"]][idx,o])
      Q <- qnorm(.975)
    }
    if (applySF & !is.null(y$sizeFactor)) {
      cts <- cts / y$sizeFactor[o]
      if (showVar) {
        sds <- sds / y$sizeFactor[o]
      }
      ylab <- "scaled counts"
    }
    ymax <- if (showVar) max(cts + Q*sds) else max(cts)
    ymin <- 0
    if (xfac & !is.null(cov)) {
      ymin <- -0.02 * ymax
    }
    if (missing(ylim)) {
      ylim <- c(ymin, ymax)
    } else {
      stopifnot(length(ylim) == 2)
    }
    if (xfac) {
      plot(cts, type="n", main=main,
           xaxt="n", ylim=ylim,
           xlab=xlabel, ylab=ylab)
    } else {
      plot(condition, cts, type="n", main=main,
           xaxt="n", ylim=ylim,
           xlab=xlabel, ylab=ylab)
    }
    seg.lwd <- if (thin == 0) 2 else if (thin == 1) 1 else 3
    seg.col <- if (thin < 2) col else col.hglt
    if (showVar) {
      if (xfac) {
        segments(seq_along(cts), pmax(cts - Q*sds, 0),
                 seq_along(cts), cts + Q*sds,
                 col=seg.col, lwd=seg.lwd)
      } else {
        segments(condition, pmax(cts - Q*sds, 0),
                 condition, cts + Q*sds,
                 col=seg.col, lwd=seg.lwd)
      }
    }
    pts.pch <- if (thin == 0) 22 else 15
    pts.lwd <- if (thin == 0) 1. else 1
    pts.cex <- if (thin == 0) 1 else 0.5
    if (xfac) {
      points(cts,
             col=col, pch=pts.pch, bg=col.hglt,
             cex=pts.cex, lwd=pts.lwd)
    } else {
      points(condition, cts,
             col=col, pch=pts.pch, bg=col.hglt,
             cex=pts.cex, lwd=pts.lwd)
    }
  }
  if (xaxis) {
    if (xfac) {
      axis(1, seq_along(condition), samp.nums, cex.axis=axisCex)
    } else {
      axis(1, cex.axis=axisCex)
    }
  }
  if (!xaxis) {
    title(xlab=xlab, mgp=c(1,1,0))
  }
  if (xfac & !is.null(cov)) {
    cuts <- cumsum(table(covariate))
    segments(c(1,cuts[-ngrp]+1),ymin,cuts,ymin,lwd=3,
             col=rep(c("black","grey60"),length=ngrp))
  }
  if (legend & (xfac | !is.null(cov))) {
    group <- if (xfac) condition else covariate
    group.name <- if (xfac) x else cov
    ltitle <- if (legendTitle) group.name else NULL
    legend(legendPos, legend=levels(group), title=ltitle,
           col=colsDrk, pt.bg=colsLgt, pch=22,
           cex=legendCex, bg="white", box.col=NA, inset=.01)
  }
}

boxplot2 <- function(x, w=.4, ylim, col, col.hglt, mainCex=1, axisCex=1, labCex=1, xlab="", ylab="", main="") {
  qs <- matrixStats::rowQuantiles(t(x), probs=0:4/4)
  if (missing(ylim)) {
    ylim <- c(min(x),max(x))
  }
  par(mgp=c(3,1,0)) 
  plot(qs[,3], type="n", xlim=c(0.5,ncol(x)+.5), xaxt="n",
       xlab=xlab, ylab=ylab, main=main, ylim=ylim,
      cex.main=mainCex, cex.axis=axisCex, cex.lab=labCex)
  s <- seq_len(ncol(x))
  rect(s-w,qs[,2],s+w,qs[,4], col=col.hglt, border=col)
  segments(s-w, qs[,3], s+w, qs[,3], col=col, lwd=3, lend=1)
  segments(s, qs[,2], s, qs[,1], col=col, lend=1)
  segments(s, qs[,4], s, qs[,5], col=col, lend=1)
  segments(s-w/2, qs[,1], s+w/2, qs[,1], col=col)
  segments(s-w/2, qs[,5], s+w/2, qs[,5], col=col)
}