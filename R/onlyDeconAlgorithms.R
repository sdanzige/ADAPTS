#' Hierarchical Deconvolution
#' @description Deconvolve cell types based on clusters detected by an n-pass spillover matrix
#'
#' @param sigMatrix  The deconvolution matrix, e.g. LM22 or MGSM27
#' @param geneExpr  The source gene expression matrix used to calculate sigMatrix
#' @param toPred  The gene expression to ultimately deconvolve
#' @param hierarchData  The results of hierarchicalSplit OR hierarchicalSplit.sc (DEFAULT: NULL, ie hierarchicalSplit)
#' @param pdfDir  A fold to write the pdf file to (DEFAULT: tempdir())
#' @param oneCore Set to TRUE to disable parallelization (DEFAULT: FALSE)
#' @param nPasses  The maximum number of iterations for spillToConvergence (DEFAULT: 100)
#' @param remZinf Set to TRUE to remove any ratio with zero or infinity when generating gList (DEFAULT: FALSE)
#' @param method  One of 'DCQ', 'SVMDECON', 'DeconRNASeq', 'proportionsInAdmixture', 'nnls' (DEFAULT: DCQ)
#' @param useRF  Set to TRUE to use ranger random forests to build the seed matrix (DEFAULT: TRUE)
#' @param incNonCluster  Set to TRUE to include a 'nonCluster' in each of the sub matrices (DEFAULT: TRUE)
#' @export
#' @return a matrix of cell counts
#' @examples
#' #This toy example 
#' library(ADAPTS)
#' fullLM22 <- ADAPTS::LM22[1:30, 1:4]
#' smallLM22 <- fullLM22[1:25,] 
#' 
#' cellCounts <- hierarchicalClassify(sigMatrix=smallLM22, geneExpr=fullLM22, toPred=fullLM22, 
#'     oneCore=TRUE, nPasses=10, method='DCQ')
hierarchicalClassify <- function(sigMatrix, geneExpr, toPred, hierarchData=NULL, pdfDir=tempdir(), oneCore=FALSE, nPasses=100, remZinf=TRUE, method='DCQ', useRF=TRUE, incNonCluster=TRUE) {
  if(is.null(hierarchData)) {
    hierarchData <- hierarchicalSplit(sigMatrix, geneExpr, oneCore=oneCore, nPasses=nPasses, remZinf=remZinf, useRF=useRF, incNonCluster=incNonCluster)
  }
  
  #Step 1: Baseline deconvolution
  toPred.sub <- toPred[rownames(toPred) %in% rownames(sigMatrix),,drop=FALSE]
  colnames(toPred.sub) <- make.names(colnames(toPred.sub), unique=TRUE)
  toPred.sub.imp <- missForest.par(toPred.sub)
  
  initDecon <- estCellPercent(refExpr = sigMatrix, geneExpr=toPred.sub.imp, method=method)
  
  #Step 2: Build the clustered Deconvolution
  clusters <- NULL
  for (i in 1:length(hierarchData$allClusters)) {
    clusters <- rbind(data.frame(cell=hierarchData$allClusters[[i]], clust=i), clusters)
  }
  clustIDs <- clusters$clust
  names(clustIDs) <- clusters$cell
  clustIDs <- clustIDs[rownames(initDecon)]
  initDecon.clust <- apply(initDecon, 2, function(x){tapply(x,clustIDs, sum)})
  
  #Step 3: Split the clusters based on the smaller groups
  curDecon.break.list <- list()
  for (i in as.numeric(rownames(initDecon.clust))) {
    curCellTypes <- hierarchData$allClusters[[i]]
    if(length(curCellTypes) > 1) {
      #curSigMat <- hierarchData$sigMatList[[i]][,curCellTypes,drop=FALSE]
      curSigMat <- hierarchData$sigMatList[[i]][,,drop=FALSE]
      toPred.sub <- toPred[rownames(toPred) %in% rownames(curSigMat),,drop=FALSE] #Filter genes
      colnames(toPred.sub) <- make.names(colnames(toPred.sub), unique=TRUE)
      toPred.sub.imp <- missForest.par(toPred.sub)
      
      curDecon <- estCellPercent(refExpr = curSigMat, geneExpr=toPred.sub.imp, method=method)
      curDecon <- curDecon[curCellTypes,,drop=FALSE]
      curDecon.frac <- apply(curDecon, 2, function(x){x/sum(x)})
      curDecon.frac <- curDecon.frac[rownames(curDecon.frac) != 'others',,drop=FALSE]
      nas <- apply(curDecon.frac, 2, function(x){ any(is.na(x)) })
      curDecon.frac[,nas] <- rep(1/nrow(curDecon.frac), nrow(curDecon.frac))
      #curDecon.frac.  Each row is a component of the current cluster.  
      #  Each column is a particular sample 
      
      # initDecon.clust[as.character(i),]  The fraction of cells in the current sample
      
      curDecon.break <- apply(curDecon.frac, 1, function(x){x * initDecon.clust[as.character(i),]}) #initDecon.clust[as.character(i),,drop=FALSE] * curDecon.frac
      
      if(ncol(toPred) > 1) {
        curDecon.break <- t(curDecon.break)
      } else {
        curDecon.break <- data.frame(first=curDecon.break)
        colnames(curDecon.break) <- colnames(toPred)
      }
      
      #curDecon.break should have broken up the amount in each column initDecon.clust[as.character(i),] scaled by each column of curDecon.frac
    } else {
      #Note, potentially for each single cluster, we could rescale this to self vs other, and rescale accordingly.
      #  This could then be fixed in the later renormalization.  However this does get a little wierd.
      
      curDecon.break <- initDecon.clust[as.character(i),,drop=FALSE]
      rownames(curDecon.break) <- curCellTypes
    } #if(length(curCellTypes) > 1) {
    curDecon.break.list[[as.character(i)]] <- curDecon.break
  } #for (i in as.numeric(rownames(initDecon.clust))) {
  
  curDecon.new <- do.call(rbind, curDecon.break.list)
  rownames(curDecon.new) <- unlist(sapply(curDecon.break.list, rownames))
  curDecon.new <- rbind(curDecon.new, initDecon['others',,drop=FALSE])
  curDecon.new.rescale <- apply(curDecon.new, 2, function(x){100*x/sum(x)}) #Fix rounding errors.
  outFile <- paste('hierarchicalSplit', Sys.Date(),'pdf', sep='.')
  outFile <- file.path(pdfDir, outFile)
  grDevices::pdf(outFile)
  res <- try(pheatmap::pheatmap(t(curDecon.new.rescale), fontsize=6, fontsize_row = 4, cluster_cols=FALSE, cluster_rows = TRUE), silent=TRUE)
  if(inherits(res, 'try-error')) { try(pheatmap::pheatmap(t(curDecon.new.rescale), fontsize=6, fontsize_row = 4, cluster_cols=FALSE, cluster_rows = FALSE), silent=TRUE) }
  grDevices::dev.off()
  
  return(curDecon.new.rescale)
}

#' Build hierarchical cell clusters.
#' @description Attempt to deconvolve cell types by building a hierarchy of cell types using
#'   spillToConvergence to determine cell types that are not signficantly different.
#'   First deconvolve those clusters of cell types.
#'   Deconvolution matrices are then built to separate the cell types that formerly could
#'   not be separated.
#'
#' @param sigMatrix  The deconvolution matrix, e.g. LM22 or MGSM27
#' @param geneExpr  The source gene expression matrix used to calculate sigMatrix
#' @param oneCore Set to TRUE to disable parallelization (DEFAULT: FALSE)
#' @param nPasses  The maximum number of iterations for spillToConvergence (DEFAULT: 100)
#' @param deconMatrices  Optional pre-computed results from spillToConvergence (DEFAULT: NULL)
#' @param remZinf Set to TRUE to remove any ratio with zero or infinity when generating gList (DEFAULT: FALSE)
#' @param method  One of 'DCQ', 'SVMDECON', 'DeconRNASeq', 'proportionsInAdmixture', 'nnls' (DEFAULT: DCQ)
#' @param useRF  Set to TRUE to use ranger random forests to build the seed matrix (DEFAULT: TRUE)
#' @param incNonCluster  Set to TRUE to include a 'nonCluster' in each of the sub matrices (DEFAULT: TRUE)
#' @export
#' @return A list of clusters and a list of signature matrices for breaking those clusters
#' @examples
#' #This toy example 
#' library(ADAPTS)
#' fullLM22 <- ADAPTS::LM22[1:30, 1:4]
#' smallLM22 <- fullLM22[1:25,] 
#' 
#' clusters <- hierarchicalSplit(sigMatrix=smallLM22, geneExpr=fullLM22, oneCore=TRUE, nPasses=10,
#'     deconMatrices=NULL, remZinf=TRUE, method='DCQ', useRF=TRUE, incNonCluster=TRUE)
hierarchicalSplit <- function(sigMatrix, geneExpr, oneCore=FALSE, nPasses=100, deconMatrices=NULL, remZinf=TRUE, method='DCQ', useRF=TRUE, incNonCluster=TRUE) {
  allClusters.rv <- clustWspillOver(sigMatrix, geneExpr, nPasses=nPasses, deconMatrices=deconMatrices, method=method)
  allClusters <- allClusters.rv$allClusters
  deconMatrices <- allClusters.rv$deconMatrices
  
  #Step 1: Do the level 1 deconvolution
  cNames <- sub('\\.+[0-9]+$', '', colnames(geneExpr))
  if(!all(cNames == colnames(geneExpr))) {
    message('Stripping .[0-9]+ from the end of gene expression column names.')
    colnames(geneExpr) <- cNames
  }
  
  #Make new signature matrices for each split.  How to determine # of genes for only 2 cell types?
  sigMatList <- list()
  for (i in 1:length(allClusters)) {
    curLen <- length(allClusters[[i]])
    if (curLen == 1) {
      sigMatList[[i]] <- as.matrix(sigMatrix)
      next;
    }
    
    curGeneExpr <- geneExpr[,colnames(geneExpr) %in% allClusters[[i]]]
    
    if(incNonCluster == TRUE) {
      #Add the other cell types
      curGeneExpr.other <- geneExpr[,!(colnames(geneExpr) %in% allClusters[[i]])]
      colnames(curGeneExpr.other) <- rep('nonCluster', length=ncol(curGeneExpr.other))
      curGeneExpr <- cbind(curGeneExpr, curGeneExpr.other)
    } #if(incNonCluster == TRUE) {
    
    naBool <- apply(curGeneExpr, 1, function(x){ any(is.na(x)) })
    if(any(naBool))  {
      message(paste('Removing', sum(naBool), 'genes due to NAs'))
      curGeneExpr <- curGeneExpr[!naBool,]
      #Note:  Why is it imputing later if I've removed all of these??? I need to fix the NA problem better.
    }
    
    colnames(curGeneExpr) <- sub('\\.[0-9]+$', '', colnames(curGeneExpr))
    if(useRF==TRUE) {
      gList <- gListFromRF(trainSet = curGeneExpr, oneCore=oneCore)
    } else {
      gList <- rankByT(geneExpr = curGeneExpr, qCut=0.3, oneCore=oneCore, remZinf=remZinf)
    } # if(useRF==TRUE) {
    if(length(gList) == 1) {
      otherCellType <- allClusters[[i]][!allClusters[[i]] %in% names(gList)]
      gList[[otherCellType]] <- gList[[1]]
    }
    origMatrix <- sigMatrix[,allClusters[[i]]]
    origMatrix.sm <- origMatrix[names(utils::tail(sort(apply(origMatrix,1,stats::var)),ceiling(nrow(sigMatrix)/10))),]
    newMatData <- try(AugmentSigMatrix(origMatrix=origMatrix.sm, fullData=curGeneExpr, newData=curGeneExpr, gList=gList,
                                       nGenes=1:100, plotToPDF=TRUE, imputeMissing=TRUE, condTol=1.01, postNorm=FALSE,
                                       minSumToRem=NA, addTitle=paste(allClusters[[i]],collapse='_'), autoDetectMin=TRUE,
                                       calcSpillOver=FALSE), silent = TRUE)
    if(inherits(newMatData, 'try-error')) {
      sigMatList[[i]] <- as.matrix(origMatrix.sm)
    } else {
      #Revision 08-20-18 -  Add in all of cell types, but with just genes in newMatData
      
      sigMatList[[i]] <- newMatData
    }
  } #for (i in 1:length(allClusters)) {
  
  return(list(allClusters=allClusters, sigMatList=sigMatList, deconMatrices=deconMatrices))
}

#' Cluster with spillover
#' @description Build clusters based on n-pass spillover matrix
#'
#' @param sigMatrix  The deconvolution matrix, e.g. LM22 or MGSM27
#' @param geneExpr  The source gene expression matrix used to calculate sigMatrix.
#' @param nPasses  The maximum number of iterations for spillToConvergence (DEFAULT: 100)
#' @param deconMatrices  Optional pre-computed results from spillToConvergence (DEFAULT: NULL)
#' @param method  One of 'DCQ', 'SVMDECON', 'DeconRNASeq', 'proportionsInAdmixture', 'nnls' (DEFAULT: DCQ)
#' @export
#' @return  Cell types grouped by cluster
#' @examples
#' #This toy example 
#' library(ADAPTS)
#' fullLM22 <- ADAPTS::LM22[1:30, 1:4]
#' smallLM22 <- fullLM22[1:25,] 
#' 
#' clusters <- clustWspillOver(sigMatrix=smallLM22, geneExpr=fullLM22, nPasses=10)
clustWspillOver <- function(sigMatrix, geneExpr, nPasses=100, deconMatrices=NULL, method='DCQ') {
  
  if(is.null(deconMatrices)) {
    curGeneExpr <- geneExpr
    naBool <- apply(curGeneExpr, 1, function(x){ any(is.na(x)) })
    if(any(naBool))  {
      message(paste('clustWspillOver: Removing', sum(naBool), 'genes due to NAs'))
      curGeneExpr <- curGeneExpr[!naBool,,drop=FALSE]
      
      keepBool <- rownames(sigMatrix) %in% rownames(geneExpr)
      message(paste('clustWspillOver: Trimming', sum(!keepBool), '/', length(keepBool), 'genes from sigMatrix due to missingness'))
      sigMatrix <- sigMatrix[keepBool,,drop=FALSE]
      #Note:  Why is it imputing later if I've removed all of these??? I need to fix the NA problem better.
    } #if(is.null(deconMatrices)) {
    
    deconMatrices <- spillToConvergence(sigMatrix=sigMatrix, geneExpr=curGeneExpr,plotIt=FALSE, nPasses=nPasses, method=method)
  }
  curExpr <- estCellCounts.nPass(geneExpr=sigMatrix, deconMatrices=deconMatrices, method=method)
  
  #Any two identical columns belong in a cluster.
  curCor <- stats::cor(curExpr)
  allClusters <- list()
  while(nrow(curCor) > 0) {
    curLabel <- rownames(curCor)[1]
    clustIdx <- (round(curCor[,curLabel],2) - 1) == 0
    curRows <- rownames(curCor)[clustIdx]
    allClusters[[length(allClusters)+1]] <- curRows
    curCor <- curCor[!clustIdx,,drop=FALSE]
  } #while(ncol(curCor) > 1) {
  
  return(list(allClusters=allClusters, deconMatrices=deconMatrices))
}

#' Deconvolve with an n-pass spillover matrix
#' @description curExpr <- estCellCounts.nPass(sigMatrix, deconMatrices)
#'
#' @param geneExpr  The gene expression matrix
#' @param deconMatrices  The results from spillToConvergence()
#' @param method  One of 'DCQ', 'SVMDECON', 'DeconRNASeq', 'proportionsInAdmixture', 'nnls' (DEFAULT: DCQ)
#' @export
#' @return An estimate of cell counts
#' @examples
#' #This toy example 
#' library(ADAPTS)
#' fullLM22 <- ADAPTS::LM22[1:30, 1:4]
#' smallLM22 <- fullLM22[1:25,] 
#' 
#' deconMatrices <- spillToConvergence(sigMatrix=smallLM22, geneExpr=fullLM22, nPasses=10)
#' cellCounts <- estCellCounts.nPass(geneExpr=fullLM22, deconMatrices=deconMatrices, method='DCQ')
estCellCounts.nPass <- function(geneExpr, deconMatrices, method='DCQ') {
  curExpr <- geneExpr
  for (curDecon in deconMatrices) {
    curExpr <- estCellPercent(refExpr = curDecon, geneExpr = curExpr, method=method)
  }
  return(curExpr)
}

#' Spillover to convergence
#' @description Build an n-pass spillover matrix, continuing until the results converge into clusters of cell types
#'
#' deconMatrices <- spillToConvergence(sigMatrix, geneExpr, 100, FALSE, TRUE)
#'
#' @param sigMatrix  The deconvolution matrix, e.g. LM22 or MGSM27
#' @param geneExpr  The source gene expression matrix used to calculate sigMatrix
#' @param nPasses  The maximum number of iterations (DEFAULT: 100)
#' @param plotIt  Set to TRUE to plot it (DEFAULT: FALSE)
#' @param imputNAs  Set to TRUE to imput genes with missing values & cache the imputed.  FALSE will just remove them (DEFAULT: FALSE)
#' @param method  One of 'DCQ', 'SVMDECON', 'DeconRNASeq', 'proportionsInAdmixture', 'nnls' (DEFAULT: DCQ)
#' @export
#' @return  A list of signature matrices
#' @examples
#' #This toy example 
#' library(ADAPTS)
#' fullLM22 <- ADAPTS::LM22[1:30, 1:4]
#' smallLM22 <- fullLM22[1:25,] 
#' 
#' deconMatrices <- spillToConvergence(sigMatrix=smallLM22, geneExpr=fullLM22, nPasses=10, plotIt=TRUE)
spillToConvergence <- function(sigMatrix, geneExpr, nPasses=100, plotIt=FALSE, imputNAs=FALSE, method='DCQ') {
  keepBool <- sub('\\.[0-9]+$', '', colnames(geneExpr)) %in% colnames(sigMatrix)
  if (!all(keepBool)) {
    message(paste('Removing', sum(!keepBool), '/', length(keepBool), 'from geneExpr that are not in sigMatrix'))
    geneExpr <- geneExpr[,keepBool,drop=FALSE]
  }
  
  keepBool <- sub('\\.[0-9]+$', '', colnames(sigMatrix)) %in% colnames(geneExpr)
  if (!all(keepBool)) {
    message(paste('Removing', sum(!keepBool), '/', length(keepBool), 'from sigMatrix that are not in geneExpr'))
    sigMatrix <- sigMatrix[,keepBool,drop=FALSE]
  }
  
  
  missingGeneBool <- !rownames(sigMatrix) %in% rownames(geneExpr)
  if (sum(missingGeneBool) > 0) {
    message(paste('Removing', sum(missingGeneBool), 'that are in sigMatrix but not geneExpr'))
    sigMatrix <- sigMatrix[!missingGeneBool,]
  }
  
  geneExpr.sub <- geneExpr[rownames(sigMatrix),]
  naGeneBool <- apply(geneExpr.sub, 1, function(x) {any(is.na(x))})
  if(any(naGeneBool)) {
    if (imputNAs==TRUE) {
      message(paste('Imputing for', sum(naGeneBool), '/', length(naGeneBool), 'genes with missing values'))
      saveFile <- paste('spillToConvergence.geneExpr.sub.imp',nrow(sigMatrix),ncol(sigMatrix),nrow(geneExpr),ncol(geneExpr),sum(naGeneBool),length(naGeneBool),'RData',sep='.')
      saveFile <- file.path(tempdir(), saveFile)
      newMatrix <- NULL
      if(file.exists(saveFile)) {
        message(paste('Loading pre-imputed', saveFile))
        newMatrix <- get(load(saveFile)[1])
        if(any(rownames(newMatrix) != rownames(geneExpr.sub)) | any(colnames(newMatrix) != colnames(geneExpr.sub))) {
          message(paste('Loaded matrix does not match genes/experiments in input matrix.'))
          newMatrix <- NULL
        }
      } #if(file.exists(saveFile)) {
      
      if(is.null(newMatrix)) {
        message('Imputing')
        newMatrix <- missForest.par(dataMat = geneExpr.sub, parallelize = "variables")
        save(newMatrix, file=saveFile)
      }
      geneExpr.sub <- newMatrix
    } else {
      message(paste('Removing', sum(naGeneBool), '/', length(naGeneBool), 'genes with missing values'))
      keepGenes <- names(naGeneBool)[!naGeneBool]
      geneExpr.sub <- geneExpr.sub[keepGenes,]
      sigMatrix <- sigMatrix[keepGenes,]
    } #if (imputNAs==TRUE) {
  }
  
  #E_0
  
  cellEst <- estCellPercent(refExpr=sigMatrix,  geneExpr=geneExpr.sub, method=method)
  if(is.null(cellEst)) {message('spillToConvergence deconvolution failed'); return(NULL)}
  #S_1
  newSig <- t(apply(cellEst, 1, function(x) { tapply(x, sub('\\.[0-9]+$', '', colnames(cellEst)), mean, na.rm=TRUE)}))
  #estToPlot <- newSig[c(sort(colnames(newSig)),'others'),sort(colnames(newSig))]
  if(plotIt==TRUE) {pheatmap::pheatmap(t(cellEst), main='First Decon Results\n y=Purified, x=Decon As', fontsize = 4)}#, cluster_rows = FALSE, cluster_cols = FALSE)
  cellEst.last <- cellEst
  cellEst.last2 <- cellEst
  
  addPassList <- list()
  addPassList[[1]] <- sigMatrix
  addPassList[[2]] <- newSig
  for (curPass in 3:nPasses) {
    #E_1, etc
    cellEst.next <- estCellPercent(refExpr=addPassList[[curPass-1]],  geneExpr=cellEst.last, method=method)
    #S_2
    newSig.next <- t(apply(cellEst.next, 1, function(x) { tapply(x, sub('\\.[0-9]+$', '', colnames(cellEst.next)), mean, na.rm=TRUE)}))
    rnames <- c(sort(rownames(cellEst.next)[rownames(cellEst.next) %in% colnames(cellEst.next)]), sort(rownames(cellEst.next)[!rownames(cellEst.next) %in% colnames(cellEst.next)]))
    cnames <- c(sort(colnames(cellEst.next)[colnames(cellEst.next) %in% rownames(cellEst.next)]), sort(colnames(cellEst.next)[!colnames(cellEst.next) %in% rownames(cellEst.next)]))
    cellEst.next <- cellEst.next[rnames,colnames(cellEst.next) %in% cnames]
    
    pairedTypes <- sort(rownames(cellEst.next)[rownames(cellEst.next) %in% colnames(cellEst.next)])
    estWself3 <- sapply(pairedTypes, function(i) {cellEst.next[i,i]})
    if(plotIt==TRUE) {graphics::barplot(estWself3, col=grDevices::rainbow(length(estWself3)), main=paste('Self Identification in Purified Samples\nPass',curPass), ylim=c(0,100), cex.names=0.66, las=2)}
    
    titleStr3 <- paste0(curPass, 'x-Re-Deconvolving results with spillover matrix\nMean Self % = ', round(mean(estWself3)))
    if(plotIt==TRUE) {pheatmap::pheatmap(t(newSig.next),cluster_rows = FALSE, cluster_cols = FALSE, main=titleStr3, xlab='DCQ', ylab='Reference', fontsize = 4)} #xlab and ylab don't work
    
    #addPassList[[curPass]] <- list(res=res3, cellEst=cellEst.third, estWself=estWself3, titleStr=titleStr3)
    addPassList[[curPass]] <- newSig.next
    
    #Test to see if there was any change from the last pass.  If not, then stop
    diffs <- sapply(colnames(cellEst.next), function(x){ sqrt(mean((cellEst.last[,x]-cellEst.next[,x])^2))})
    if(all(diffs == 0)) { break;}
    
    #Test to see if the signature matrices are oscilating
    diffs <- sapply(colnames(cellEst.next), function(x){ sqrt(mean((cellEst.last2[,x]-cellEst.next[,x])^2))})
    if(all(diffs == 0)) { break;}
    
    cellEst.last2 <- cellEst.last
    cellEst.last <- cellEst.next
  } #for (curPass in 3:nPasses) {
  
  if(plotIt==TRUE) {pheatmap::pheatmap(t(newSig.next),cluster_rows = TRUE, cluster_cols = TRUE, main=titleStr3, xlab='DCQ', ylab='Reference', fontsize = 4)} #xlab and ylab don't work
  
  return(addPassList)
}


#' Build a spillover matrix
#' @description Build a spillover matrix, i.e. what do purified samples deconvolve as?
#'
#' spillExpr <- buildSpilloverMat(refExpr, geneExpr, method='DCQ')
#'
#' @param refExpr  The deconvolution matrix, e.g. LM22 or MGSM27
#' @param geneExpr  The full gene expression for purified cell types.  Multiple columns (examples) for each column in the reference expr.
#' @param method  One of 'DCQ', 'SVMDECON', 'DeconRNASeq', 'proportionsInAdmixture', 'nnls' (DEFAULT: DCQ)
#' @export
#' @return A spillover matrix showing how purified cell types deconvolve
#' @examples
#' #This toy example 
#' library(ADAPTS)
#' fullLM22 <- ADAPTS::LM22[1:30, 1:4]
#' smallLM22 <- fullLM22[1:25,] 
#' 
#' spillover <- buildSpilloverMat(refExpr=smallLM22, geneExpr=fullLM22, method='DCQ')
buildSpilloverMat <- function(refExpr, geneExpr, method='DCQ') {
  if(any(grepl('\\.[0-9]+$', unique(colnames(geneExpr))))) {
    print('###Note: Some of the column names in geneExpr end in a . followed by numbers###')
    print('Different samples with the same cell types should have exactly the same column names')
    print('')
  }
  
  olGenes <- rownames(refExpr)[rownames(refExpr) %in% rownames(geneExpr)]
  refExpr <- refExpr[olGenes,]
  failedGene1 <- apply(refExpr, 1, function(x){any(is.na(x))})
  geneExpr <- geneExpr[olGenes,]
  failedGene2 <- apply(geneExpr, 1, function(x){any(is.na(x))})
  keepGenes <- !failedGene1 & !failedGene2
  
  cellEst <- estCellPercent(refExpr=refExpr[keepGenes,],  geneExpr=geneExpr[keepGenes,], method=method)
  
  res <- res.bk <- apply(cellEst, 1, function(x) { tapply(x, colnames(cellEst), mean, na.rm=TRUE)})
  res <- res[,colnames(res)!='others']
  res <- res[order(toupper(rownames(res))),order(toupper(colnames(res)))]
  res <- cbind(res, data.frame(others = res.bk[,'others']))
  
  return(res)
}

#' Estimate cell percentage from spillover
#' @description Use a spillover matrix to deconvolve a samples
#'
#' @param spillExpr  A spill over matrix, as calculated by buildSpilloverMat(). (e.g. LM22.spillover.csv.gz)
#' @param refExpr  a data frame representing immune cell expression profiles. Each row represents an expression of a gene, and each column represents a different immune cell type. colnames contains the name of each immune cell type and the rownames includes the genes' symbol. The names of each immune cell type and the symbol of each gene should be unique. Any gene with missing expression values must be excluded.
#' @param geneExpr  a data frame representing RNA-seq or microarray gene-expression profiles of a given complex tissue. Each row represents an expression of a gene, and each column represents a different experimental sample. colnames contain the name of each sample and rownames includes the genes' symbol. The name of each individual sample and the symbol of each gene should be unique. Any gene with missing expression values should be excluded.
#' @param method  One of 'DCQ', 'SVMDECON', 'DeconRNASeq', 'proportionsInAdmixture', 'nnls' (DEFAULT: DCQ)
#' @param ...  Parameters for estCellPercent.X (e.g. number_of_repeats for .DCQ)
#' @export
#' @return  a matrix of estimate cell type percentages in samples
#' @examples
#' #This toy example 
#' library(ADAPTS)
#' fullLM22 <- ADAPTS::LM22[1:30, 1:4]
#' smallLM22 <- fullLM22[1:25,] 
#' 
#' spillover <- buildSpilloverMat(refExpr=smallLM22, geneExpr=fullLM22) 
#' cellEst <- estCellPercent.spillOver(spillExpr=spillover, refExpr=smallLM22, geneExpr=fullLM22)
estCellPercent.spillOver <- function(spillExpr, refExpr,  geneExpr, method='DCQ', ...) {
  if(method == 'DCQ') {
    estCellPercent.X <- estCellPercent.DCQ
  } else if (method == 'SVMDECON') {
    estCellPercent.X <- estCellPercent.svmdecon
  } else if (method == 'DeconRNASeq')  {
    estCellPercent.X <- estCellPercent.DeconRNASeq
  } else if (method == 'proportionsInAdmixture') {
    estCellPercent.X <- estCellPercent.proportionsInAdmixture
  } else if (method == 'nnls') {
    estCellPercent.X <- estCellPercent.nnls
  }
  
  cellEst <- estCellPercent.X(refExpr=refExpr,  geneExpr=geneExpr, ...)
  cellEst2 <- estCellPercent.X(refExpr=t(spillExpr),  geneExpr=cellEst, ...)
  
  return(cellEst2)
}

#' DCQ Deconvolution
#' @description Use DCQ to estimate the cell count percentage
#' Requires installation of package 'ComICS'
#'   To Do: Also report the standard deviation as a confidence metric
#'
#' @param refExpr  a data frame representing immune cell expression profiles. Each row represents an expression of a gene, and each column represents a different immune cell type. colnames contains the name of each immune cell type and the rownames includes the genes' symbol. The names of each immune cell type and the symbol of each gene should be unique. Any gene with missing expression values must be excluded.
#' @param geneExpr  a data frame representing RNA-seq or microarray gene-expression profiles of a given complex tissue. Each row represents an expression of a gene, and each column represents a different experimental sample. colnames contain the name of each sample and rownames includes the genes' symbol. The name of each individual sample and the symbol of each gene should be unique. Any gene with missing expression values should be excluded.
#' @param marker_set  data frames of one column, that includes a preselected list of genes that likely discriminate well between the immune-cell types given in the reference data. (DEFAULT: NULL, i.e. one for each gene in the refExpr)
#' @param number_of_repeats  using one repeat will generate only one output model. Using many repeats, DCQ calculates a collection of models, and outputs the average and standard deviation for each predicted relative cell quantity. (DEFAULT: 1)
#' @param alpha  The elasticnet mixing parameter, with 0 <= alpha <= 1. alpha=1 is the lasso penalty, and alpha=0 the ridge penalty. (DEFAULT: 0.05)
#' @param lambda A minimum value for the elastic net lambda parameter (DEFAULT: 0.2)
#' @export
#' @return A matrix with cell type estimates for each samples
#' @examples
#' #This toy example 
#' library(ADAPTS)
#' fullLM22 <- ADAPTS::LM22[1:30, 1:4]
#' smallLM22 <- fullLM22[1:25,] 
#' 
#' cellEst <- estCellPercent.DCQ(refExpr=smallLM22, geneExpr=fullLM22)
estCellPercent.DCQ <- function(refExpr,  geneExpr, marker_set=NULL, number_of_repeats=10, alpha=0.05, lambda=0.2) {
  
  if(any(is.na(geneExpr))) {
    message('There are some NAs in geneExpr, please impute or remove NAs')
    return(NULL)
  }
  
  if(is.null(marker_set)) {marker_set <- data.frame(marker_set=rownames(refExpr))}
  if(ncol(geneExpr)==1) {geneExpr <- cbind(geneExpr, geneExpr); fixOneCol<-TRUE} else {fixOneCol<-FALSE}
  suppressWarnings(sink(""))
  cellCounts <- try(ComICS::dcq(reference_data = refExpr, mix_data = geneExpr, marker_set = marker_set, number_of_repeats=number_of_repeats,lambda_min = lambda,alpha_used = alpha))
  sink()
  if(inherits(cellCounts, 'try-error')) {return(cellCounts)}
  #DCQ returns a list with two matrices: average quantities for each cell type, stdev over all repeats for each cell type
  #Question for MM-010 data, how many cells?? Frank says 1-10 million
  cellCoefVar <- cellCoefVar.m <- t(cellCounts$stdev / cellCounts$average)
  cellCountsPercent <- cellCountsPercent.m <- t(cellCounts$average*100)
  
  cellCountsPercent.m[cellCountsPercent.m > 0] <- 0
  cellCoefVar.m[cellCoefVar.m > 0] <- 0
  others <- abs(colSums(cellCountsPercent.m))
  others.m <- abs(colMeans(cellCoefVar.m, na.rm=TRUE))
  
  cellCountsPercent[cellCountsPercent < 0] <- 0
  cellCountsPercent <- rbind(cellCountsPercent, others)
  cellCountsPercent <- round(apply(cellCountsPercent, 2, function(x){100*x/sum(x)}),2)
  
  cellCoefVar <- rbind(cellCoefVar, others.m)
  cellCoefVar[is.na(cellCoefVar)] <- 0
  cellCountsPercent.std <- cellCoefVar * cellCountsPercent
  
  #Note: cellCountsPercent.std is very small.  This DCQ error is not good enough for the estimate
  if(fixOneCol==TRUE) {
    cellCountsPercent <- cellCountsPercent[,1,drop=FALSE]
  }
  
  return (cellCountsPercent)
}

#' SVMDECON deconvolution
#' @description Use SVMDECON to estimate the cell count percentage
#' Performs considerably worse in deconvolution than DCQ
#'
#' cellEst <- estCellPercent.svmdecon(refExpr,  geneExpr)
#'
#' @param refExpr  a data frame representing immune cell expression profiles. Each row represents an expression of a gene, and each column represents a different immune cell type. colnames contains the name of each immune cell type and the rownames includes the genes' symbol. The names of each immune cell type and the symbol of each gene should be unique. Any gene with missing expression values must be excluded.
#' @param geneExpr  a data frame representing RNA-seq or microarray gene-expression profiles of a given complex tissue. Each row represents an expression of a gene, and each column represents a different experimental sample. colnames contain the name of each sample and rownames includes the genes' symbol. The name of each individual sample and the symbol of each gene should be unique. Any gene with missing expression values should be excluded.
#' @param marker_set  data frames of one column, that includes a preselected list of genes that likely discriminate well between the immune-cell types given in the reference data. (DEFAULT: NULL, i.e. one for each gene in the refExpr)
#' @param useOldVersion  Set the TRUE to 2^ the data (DEFAULT: FALSE)
#' @param progressBar  Set to TRUE to show a progress bar  (DEFAULT: TRUE)
#' @export
#' @return A matrix with cell type estimates for each samples
#' #This toy example 
#' library(ADAPTS)
#' fullLM22 <- ADAPTS::LM22[1:30, 1:4]
#' smallLM22 <- fullLM22[1:25,] 
#' 
#' cellEst <- estCellPercent.svmdecon(refExpr=smallLM22, geneExpr=fullLM22)
estCellPercent.svmdecon <- function(refExpr,  geneExpr, marker_set=NULL, useOldVersion=F,progressBar = T) {
  if(is.null(marker_set)) {marker_set <- data.frame(marker_set=rownames(refExpr))}
  #if(ncol(geneExpr)==1) {geneExpr <- cbind(geneExpr, geneExpr)}
  
  refExprMatrix <- as.matrix(refExpr)
  
  refGenes <- rownames(refExprMatrix)[rownames(refExprMatrix) %in% rownames(geneExpr)]
  refExprMatrix <- refExprMatrix[refGenes,]
  geneExpr <- geneExpr[refGenes,,drop=FALSE]
  
  if(useOldVersion == F){
    geneExpr <- 2^geneExpr
    refExprMatrix <- 2^refExprMatrix
  }
  
  proportions <- matrix(nrow=ncol(refExprMatrix), ncol=ncol(geneExpr))
  for (column in 1:ncol(geneExpr)) {
    # SVMDECON returns the estimated proportions of each cell-type in the given sample
    dataCol <- as.matrix(geneExpr[,column,drop=FALSE])
    proportions[, column] <- t(SVMDECON(dataCol, refExprMatrix))
  }
  
  others <- numeric(ncol(proportions))
  proportions <- rbind(proportions, others)
  
  colnames(proportions) <- colnames(geneExpr)
  rownames(proportions) <- c(colnames(refExprMatrix), "others")
  
  cellCountsPercent <- round(proportions * 100, 2)
  return(cellCountsPercent)
}

#' DeconRNASeq deconvolution
#' @description Use DeconRNASeq to estimate the cell count percentage
#' Performs with similar effectiveness as DCQ, but identifies different proportions of cell-types
#' Requires installation of package 'DeconRNASeq':
#'    source("https://bioconductor.org/biocLite.R")
#'    biocLite("DeconRNASeq")
#'
#'    <joseph.szustakowski@novartis.com> TGJDS (2013). DeconRNASeq: Deconvolution of Heterogeneous Tissue Samples for mRNA-Seq data. R package version 1.18.0.
#'
#' cellEst <- estCellPercent.DeconRNASeq(refExpr,  geneExpr, marker_set=NULL)
#'
#' @param refExpr  a data frame representing immune cell expression profiles. Each row represents an expression of a gene, and each column represents a different immune cell type. colnames contains the name of each immune cell type and the rownames includes the genes' symbol. The names of each immune cell type and the symbol of each gene should be unique. Any gene with missing expression values must be excluded.
#' @param geneExpr  a data frame representing RNA-seq or microarray gene-expression profiles of a given complex tissue. Each row represents an expression of a gene, and each column represents a different experimental sample. colnames contain the name of each sample and rownames includes the genes' symbol. The name of each individual sample and the symbol of each gene should be unique. Any gene with missing expression values should be excluded.
#' @param marker_set  data frames of one column, that includes a preselected list of genes that likely discriminate well between the immune-cell types given in the reference data. (DEFAULT: NULL, i.e. one for each gene in the refExpr)
#' @export
#' @return A matrix with cell type estimates for each samples
#' @examples
#' \donttest{
#' #This toy example, donttest due to performance issues in windows development build 
#' library(ADAPTS)
#' fullLM22 <- ADAPTS::LM22[1:30, 1:4]
#' smallLM22 <- fullLM22[1:25,] 
#' 
#' cellEst <- estCellPercent.DeconRNASeq(refExpr=smallLM22, geneExpr=fullLM22)
#' }
estCellPercent.DeconRNASeq <- function(refExpr,  geneExpr, marker_set=NULL) {
  if(!'DeconRNASeq' %in% rownames(utils::installed.packages()) | !exists('DeconRNASeq')) {
    message('estCellPercent.DeconRNASeq requires DeconRNASeq')
    message('https://www.bioconductor.org/packages/release/bioc/html/DeconRNASeq.html')
    message('Run library(DeconRNASeq) after installation')
    return(NULL)
  }
  
  if(!exists('DeconRNASeq')) {
    DeconRNASeq <- function(datasets = NULL, signatures=NULL) { stop('DeconRNASeq not loaded') }
  }
  
  if(is.null(marker_set)) {marker_set <- data.frame(marker_set=rownames(refExpr))}
  if(ncol(geneExpr)==1) {geneExpr <- cbind(geneExpr, geneExpr)}
  
  pca <- pcaMethods::pca  #Something has clobbered PCA
  curDecon <- try(DeconRNASeq(datasets=as.data.frame(geneExpr), signatures=refExpr))
  if(inherits(curDecon, 'try-error')) {
    message('Please update all packages called by DeconRNASeq')
    return(NULL)
  }
  cellProportions <- t(curDecon$out.all)
  
  cellCountsPercent <- round(cellProportions * 100, 2)
  others <- numeric(ncol(cellCountsPercent))
  cellCountsPercent <- rbind(cellCountsPercent, others)
  colnames(cellCountsPercent) <- colnames(geneExpr)
  rownames(cellCountsPercent) <- c(colnames(refExpr), "others")
  
  return (cellCountsPercent)
}

#' WGCNA::proportionsInAdmixture deconvolution
#' @description Use R function proportionsInAdmixture to estimate the cell count percentage
#' Uses the 'WGCNA' package
#'
#' cellEst <- estCellPercent.proportionsInAdmixture(refExpr)
#'
#' @param refExpr  a data frame representing immune cell expression profiles. Each row represents an expression of a gene, and each column represents a different immune cell type. colnames contains the name of each immune cell type and the rownames includes the genes' symbol. The names of each immune cell type and the symbol of each gene should be unique. Any gene with missing expression values must be excluded.
#' @param geneExpr  a data frame representing RNA-seq or microarray gene-expression profiles of a given complex tissue. Each row represents an expression of a gene, and each column represents a different experimental sample. colnames contain the name of each sample and rownames includes the genes' symbol. The name of each individual sample and the symbol of each gene should be unique. Any gene with missing expression values should be excluded.
#' @param marker_set  data frames of one column, that includes a preselected list of genes that likely discriminate well between the immune-cell types given in the reference data. (DEFAULT: NULL, i.e. one for each gene in the refExpr)
#' @export
#' @return A matrix with cell type estimates for each samples
#' @examples
#' \donttest{
#' #This toy example 
#' library(ADAPTS)
#' fullLM22 <- ADAPTS::LM22[1:30, 1:4]
#' smallLM22 <- fullLM22[1:25,] 
#' 
#' cellEst <- estCellPercent.proportionsInAdmixture(refExpr=smallLM22, geneExpr=fullLM22)
#' }
estCellPercent.proportionsInAdmixture <- function(refExpr,  geneExpr, marker_set=NULL) {
  if(!'WGCNA' %in% rownames(utils::installed.packages()) | !exists('proportionsInAdmixture')) {
    message('WGCNA required for proportionsInAdmixture deconvolution')
    message('https://cran.r-project.org/web/packages/WGCNA/')
    message('Run library(WGCNA) after installation')
    return(NULL)
  }
  
  if(!exists('proportionsInAdmixture')) {
    proportionsInAdmixture <- function(MarkerMeansPure = NULL, datE.Admixture=NULL) { stop('proportionsInAdmixture not loaded') }
  }
  
  if(is.null(marker_set)) {marker_set <- data.frame(marker_set=rownames(refExpr))}
  if(ncol(geneExpr)==1) {geneExpr <- cbind(geneExpr, geneExpr)}
  
  # Filter the gene expression matrix and reference expression matrix to only having the same genes
  refExprMatrix <- as.matrix(refExpr)
  refGenes <- rownames(refExprMatrix)[rownames(refExprMatrix) %in% rownames(geneExpr)]
  refExprMatrix <- refExprMatrix[refGenes,]
  geneExprMatrix <- geneExpr[refGenes,]
  
  # Changing the input matrices into relevant formats (data frames of corresponding sizes)
  refExprDF <- as.data.frame(refExprMatrix)
  refExprDF <- cbind(rownames(refExprDF), refExprDF)
  geneExprDF <- as.data.frame(t(geneExprMatrix))
  
  # function proportionsInAdmixture takes in a data frame whose first column reports the gene names and remaining columns report
  # the gene expression in specific cell types, and a data frame whose rows represent each sample and columns represent the gene expression.
  # proportionInAdmixture returns a list where rows are samples and columns are cell types
  proportionList <- try(proportionsInAdmixture(MarkerMeansPure = refExprDF, datE.Admixture = geneExprDF))
  if(inherits(proportionList, 'try-error')) { message('proportionsInAdmixture failed'); return(NULL); }
  proportionList <- proportionList["PredictedProportions"]
  proportionMatrix <- t(matrix(unlist(proportionList), ncol = ncol(refExprMatrix), byrow = FALSE))
  
  cellCountsPercent <- round(proportionMatrix * 100, 2)
  others <- numeric(ncol(cellCountsPercent))
  cellCountsPercent <- rbind(cellCountsPercent, others)
  colnames(cellCountsPercent) <- colnames(geneExprMatrix)
  rownames(cellCountsPercent) <- c(colnames(refExprMatrix), "others")
  
  return (cellCountsPercent)
}

#' Non-negative least squares deconvolution
#' @description Use non-negative least squares regression to deconvolve a sample
#'     This is going to be to simple to be useful
#'     This might be more interesting if I used non-positive least squares to detect 'other'
#'
#' @param refExpr  a data frame representing immune cell expression profiles. Each row represents an expression of a gene, and each column represents a different immune cell type. colnames contains the name of each immune cell type and the rownames includes the genes' symbol. The names of each immune cell type and the symbol of each gene should be unique. Any gene with missing expression values must be excluded.
#' @param geneExpr  a data frame representing RNA-seq or microarray gene-expression profiles of a given complex tissue. Each row represents an expression of a gene, and each column represents a different experimental sample. colnames contain the name of each sample and rownames includes the genes' symbol. The name of each individual sample and the symbol of each gene should be unique. Any gene with missing expression values should be excluded.
#' @export
#' @return A matrix with cell type estimates for each samples
#' @examples
#' #This toy example 
#' library(ADAPTS)
#' fullLM22 <- ADAPTS::LM22[1:30, 1:4]
#' smallLM22 <- fullLM22[1:25,] 
#' 
#' cellEst <- estCellPercent.nnls(refExpr=smallLM22, geneExpr=fullLM22)
estCellPercent.nnls <- function(refExpr,  geneExpr) {
  marker_set <- data.frame(marker_set=rownames(refExpr))
  #if(ncol(geneExpr)==1) {geneExpr <- cbind(geneExpr, geneExpr)}
  
  refExprMatrix <- as.matrix(refExpr)
  
  refGenes <- rownames(refExprMatrix)[rownames(refExprMatrix) %in% rownames(geneExpr)]
  refExprMatrix <- refExprMatrix[refGenes,,drop=FALSE]
  geneExpr <- geneExpr[refGenes,,drop=FALSE]
  
  proportions <- matrix(nrow=ncol(refExprMatrix), ncol=ncol(geneExpr))
  for (column in 1:ncol(geneExpr)) {
    # SVMDECON returns the estimated proportions of each cell-type in the given sample
    dataCol <- as.matrix(geneExpr[,column,drop=FALSE])
    reg <- nnls::nnls(refExprMatrix, dataCol)
    curDec <- reg$x
    names(curDec) <- colnames(refExprMatrix)
    proportions[, column] <- curDec
  }
  
  others <- numeric(ncol(proportions))
  proportions <- rbind(proportions, others)
  
  colnames(proportions) <- colnames(geneExpr)
  rownames(proportions) <- c(colnames(refExprMatrix), "others")
  
  cellCountsPercent <- round(proportions * 100, 2)
  return(cellCountsPercent)
}

#' SVMDECONV helper function
#' @description Use weightNorm to normalize the SVM weights.  Used for SVMDECONV
#'
#' w1 <- weightNorm(w)
#'
#' @param w  The weight vector from fitting an SVM, something like something like t(fit1$coefs) \%*\% fit1$SV, where fit comes from  <- svm(m~B, nu=0.25, kernel="linear"))
#' @return a weight vector
weightNorm <- function(w) {
  w[w<0] <- 0
  return(w/sum(w))
}


#' Support vector machine deconvolution
#' @description Use SVMDECONV to estimate the cell count percentage
#' David L Gibbs, dgibbs@systemsbiology.org
#' June 9, 2017
#'
#' v-SVR is applied with a linear kernel to solve for f,
#' and the best result from three values of v = {0.25, 0.5, 0.75}
#' is saved, where ‘best’ is defined as the lowest root mean squared error
#' between m and the deconvolution result, f x B.
#'
#' Our current implementation executes v-SVR using the
#' ‘svm’ function in the R package, ‘e1071’.
#'
#' w2 <- SVMDECON(m, B)
#'
#' @param m  a matrix represenging the mixture (genes X 1 sample)
#' @param B  a matrix representing the references (genes X cells), m should be subset to match B
#'
#' @return A matrix with cell type estimates for each samples
SVMDECON <- function(m,B) {
  
  # three models are fit with different values of nu
  fit1 <- e1071::svm(m~B, nu=0.25, kernel="linear", scale=T, type="nu-regression")
  fit2 <- e1071::svm(m~B, nu=0.50, kernel="linear", scale=T, type="nu-regression")
  fit3 <- e1071::svm(m~B, nu=0.75, kernel="linear", scale=T, type="nu-regression")
  # these w's are the cell fractions
  w1 <- weightNorm(t(fit1$coefs) %*% fit1$SV)
  w2 <- weightNorm(t(fit2$coefs) %*% fit2$SV)
  w3 <- weightNorm(t(fit3$coefs) %*% fit3$SV)
  # return the model with the smallest mean sq error
  err1 <- sqrt( sum( (m - B %*% t(w1))^2 )/nrow(m) )
  err2 <- sqrt( sum( (m - B %*% t(w2))^2 )/nrow(m) )
  err3 <- sqrt( sum( (m - B %*% t(w3))^2 )/nrow(m) )
  resIdx <- which(c(err1,err2,err3) == min(c(err1,err2,err3)))[1]
  return(list(w1,w2,w3)[[resIdx]])
}


#' Collapse cell types
#' @description Collapse the cell types (in rows) to super-classes
#' Including MGSM36 cell types
#'
#' @param cellCounts A matrix with cell counts
#' @param method The method for combining cell types ('Default: 'Pheno2')
#'               Pheno1: Original cell-type based combinations
#'               Pheno2: Original cell-type based combinations, omitting Macrophages
#'               Pheno3: Alt Phenotype definitions based on WMB deconvolution correlations
#'               Pheno4: Consensus cell types
#'               Pheno5: Consensus cell types, combined myeloma & plasma
#'               Spillover1: Empirical combinations based on compToLM22source
#'               Spillover2: More agressive combination based on empirical combinations based on compToLM22source
#'               Spillover3: Combinations determined by spillToConvergence on 36 cell types
#' @return NULL
#' @export
#' @return a cell estimate matrix with the names changed
#' @examples
#' #This toy example 
#' library(ADAPTS)
#' fullLM22 <- ADAPTS::LM22[1:30, 1:4]
#' smallLM22 <- fullLM22[1:25,] 
#' 
#' cellEst <- estCellPercent.DCQ(refExpr=smallLM22, geneExpr=fullLM22)
#' collapseCounts <- collapseCellTypes(cellCounts=cellEst)
collapseCellTypes <- function(cellCounts, method='Pheno4') {
  
  if(method == "Pheno1" | method == 'Pheno2') {
    combList <- list(MastCells=c('Mast.cells.activated', 'Mast.cells.resting'),
                     NKcells=c('NK.cells.activated', 'NK.cells.resting'),
                     HealthyBcells=c('B.cells.memory', 'B.cells.naive'),
                     CD8s=c('T.cells.CD8'),
                     CD4s=c('T.cells.CD4.memory.resting', 'T.cells.CD4.memory.activated', 'T.cells.CD4.naive',
                            'T.cells.follicular.helper', 'T.cells.regulatory..Tregs.'),
                     DendriticCells=c('Dendritic.cells.resting', 'Dendritic.cells.activated'),
                     Myeloma=c('CustomMM', 'PlasmaMemory', 'Plasma.cells', 'MM.plasma.cell'),
                     BMFibroblast=c('BMFibroblast','uamsBMFibroblast')
    )
  }
  if(method == "Pheno1" | method == 'Pheno2') {
    newList <- list(Macrophages=c('Macrophages.M0', 'Macrophages.M1', 'Macrophages.M2'))
    combList <- c(combList, newList)
  }
  
  if(method == 'Pheno3') {
    combList <- list(MastCells=c('Mast.cells.activated', 'Mast.cells.resting'),
                     NKcells=c('NK.cells.activated', 'NK.cells.resting'),
                     HealthyBcells=c('B.cells.memory', 'B.cells.naive'),
                     CD8s=c('T.cells.CD8'),
                     CD4s=c('T.cells.CD4.memory.resting', 'T.cells.CD4.naive',
                            'T.cells.follicular.helper', 'T.cells.regulatory..Tregs.'), #Keep 'T.cells.CD4.memory.activated' separate
                     DendriticCells=c('Dendritic.cells.resting', 'Dendritic.cells.activated'),
                     PlasmaCellMemory=c('PlasmaMemory', 'Plasma.cells'),
                     Myeloma=c('CustomMM', 'MM.plasma.cell'),
                     MonocyteNeutrophil=c('Monocytes', 'Neutrophils'),
                     BMFibroblast=c('BMFibroblast','uamsBMFibroblast')
    )
  }
  
  if(method == 'Pheno4') {
    combList <- list(MastCells=c('Mast.cells.activated', 'Mast.cells.resting'),
                     NKcells=c('NK.cells.activated', 'NK.cells.resting'),
                     Bcells=c('B.cells.memory', 'B.cells.naive'),
                     CD8s=c('T.cells.CD8'),
                     CD4s=c('T.cells.CD4.memory.resting', 'T.cells.CD4.naive', 'T.cells.CD4.memory.activated',
                            'T.cells.follicular.helper', 'T.cells.regulatory..Tregs.'),
                     DendriticCells=c('Dendritic.cells.resting', 'Dendritic.cells.activated'),
                     HealthyPlasmaCell=c('PlasmaMemory', 'Plasma.cells'),
                     Myeloma=c('CustomMM', 'MM.plasma.cell'),
                     BMFibroblast=c('BMFibroblast','uamsBMFibroblast')
    )
  }
  
  if(method == 'Pheno5') {
    combList <- list(MastCells=c('Mast.cells.activated', 'Mast.cells.resting'),
                     NKcells=c('NK.cells.activated', 'NK.cells.resting'),
                     Bcells=c('B.cells.memory', 'B.cells.naive'),
                     CD8s=c('T.cells.CD8'),
                     CD4s=c('T.cells.CD4.memory.resting', 'T.cells.CD4.naive', 'T.cells.CD4.memory.activated',
                            'T.cells.follicular.helper', 'T.cells.regulatory..Tregs.'),
                     DendriticCells=c('Dendritic.cells.resting', 'Dendritic.cells.activated'),
                     Myeloma=c('CustomMM', 'MM.plasma.cell', 'PlasmaMemory', 'Plasma.cells'),
                     BMFibroblast=c('BMFibroblast','uamsBMFibroblast')
    )
  }
  
  if(method == 'Spillover1') {
    #Comb B-cells
    #Comb DC.act & M1
    #Comb DC.rest & M0 & M2
    #Comb Mast
    #Comb Mono/Neutro
    #Comb CD4.resting, CD4.naive, FHs, Tregs
    #Comb NKs
    #Plasma Cells alone
    #CD4.activated alone
    #Leave Eo alone
    #Leave CD8s alone
    #Leave gdTs alone
    combList <- list(MastCells=c('Mast.cells.activated', 'Mast.cells.resting'),
                     NKcells=c('NK.cells.activated', 'NK.cells.resting'),
                     HealthyBcells=c('B.cells.memory', 'B.cells.naive'),
                     MonocyteNeutrophil=c('Monocytes', 'Neutrophils'),
                     APC.activated=c('Dendritic.cells.activated', 'Macrophages.M1'),
                     APC.inactive=c('Dendritic.cells.resting', 'Macrophages.M0', 'Macrophages.M2'),
                     CD8s=c('T.cells.CD8'),
                     CD4s.MemoryActivated= c('T.cells.CD4.memory.activated'),
                     CD4s.others=c('T.cells.CD4.memory.resting', 'T.cells.CD4.naive',
                                   'T.cells.follicular.helper', 'T.cells.regulatory..Tregs.'), #Keep separate
                     PlasmaCellMemory=c('PlasmaMemory', 'Plasma.cells'),
                     Myeloma=c('CustomMM', 'MM.plasma.cell'),
                     BMFibroblast=c('BMFibroblast','uamsBMFibroblast')
    )
  }
  
  if(method == 'Spillover2') {
    # Eosinphils alone
    # Monocytes+Neutrophils
    # Mast.cells.activated/resting
    # M1 macrophages & activated dendritic cells
    # Adipocyte / Osteoblast
    # B.cells memory/naive
    # cd138p/PlasmaMemory/Plasma.cell/MM.plasma.cell
    # T.cells.gamma.delta
    # NK.cells act/resting
    # osteoclast/M2/resting dendritic/M0
    # CD8s
    # Treg / TcelFH / memory resting / naive
    combList <- list(Eosinophils=c('Eosinophils'),
                     MonocyteNeutrophil=c('Monocytes', 'Neutrophils'),
                     MastCells=c('Mast.cells.activated', 'Mast.cells.resting'),
                     APC.activated=c('Dendritic.cells.activated', 'Macrophages.M1'),
                     AdipoOsteo=c('adipocyte', 'osteoblast'),
                     HealthyBcells=c('B.cells.memory', 'B.cells.naive'),
                     Myeloma=c('CustomMM', 'MM.plasma.cell', 'PlasmaMemory', 'Plasma.cells'),
                     GammaDeltaT=c('T.cells.gamma.delta'),
                     NKcells=c('NK.cells.activated', 'NK.cells.resting'),
                     APC.inactive=c('osteoblast', 'Dendritic.cells.resting', 'Macrophages.M0', 'Macrophages.M2'),
                     CD8s=c('T.cells.CD8'),
                     CD4s.MemoryActivated= c('T.cells.CD4.memory.activated'),
                     CD4s.others=c('T.cells.CD4.memory.resting', 'T.cells.CD4.naive',
                                   'T.cells.follicular.helper', 'T.cells.regulatory..Tregs.'), #Keep separate
                     BMFibroblast=c('BMFibroblast','uamsBMFibroblast')
    ) #combList
  }
  
  if(method == 'Spillover3') {
    #MDSCs and CAFs should be left alone
    combList <- list(EosinoMegakaryoEryth=c('Eosinophils','Megakaryocyte','Erythroid'),
                     MonocyteNeutrophil=c('Monocytes', 'Neutrophils'),
                     MastCells=c('Mast.cells.activated', 'Mast.cells.resting'),
                     APC.activated=c('Dendritic.cells.activated', 'Macrophages.M1'),
                     APC.resting=c('Dendritic.cells.resting', 'Macrophages.M0', 'Macrophages.M2', 'osteoclast'),
                     AdipoLym=c('adipocyte', 'LymEndothelial'),
                     FibroOsteoBlast = c('uamsBMFibroblast','BMFibroblast','osteoblast'),
                     HealthyBcells=c('B.cells.memory', 'B.cells.naive', 'Plasma.cells'),
                     Myeloma=c('CustomMM', 'MM.plasma.cell', 'PlasmaMemory', 'Haemetopoetic'),
                     GammaDeltaT.CD4naive.FH=c('T.cells.gamma.delta','T.cells.CD4.naive', 'T.cells.follicular.helper'),
                     NKcells=c('NK.cells.activated', 'NK.cells.resting'),
                     CD8s.CD4sOthers=c('T.cells.CD8', 'T.cells.CD4.memory.resting', 'T.cells.CD4.memory.activated', 'T.cells.regulatory..Tregs.')
    ) #combList
  }
  
  countMatrix <- cellCounts[!grepl('_', rownames(cellCounts)),]
  
  #Add tumor & make it optional (DEFAULT to collapse to these types)
  countMatrix.comT <- countMatrix
  for(combCell in names(combList)) {
    curTypes <- combList[[combCell]]
    curTypes <- curTypes[curTypes %in% rownames(countMatrix)]
    if(length(curTypes) == 0) { next; }
    
    combCounts <- colSums(countMatrix[curTypes,,drop=FALSE])
    newDF <- data.frame(newData=combCounts)
    colnames(newDF) <- combCell
    countMatrix.comT <- rbind(countMatrix.comT[!(rownames(countMatrix.comT) %in% curTypes),], t(newDF))
  }
  
  countMatrix.comT <- countMatrix.comT[order(toupper(rownames(countMatrix.comT))),]
  countMatrix.comT <- rbind(countMatrix.comT, cellCounts[grepl('_', rownames(cellCounts)),])
  
  return(countMatrix.comT)
}
#' Wrapper for deconvolution methods
#' @description A wrapper function to call any of the estCellPercent functions
#'    Modified on June 16th 2021 to quantile normalize the geneExpr data to match refExpr
#'    Set preNormalize to FALSE for previous behavior.
#'
#' @param refExpr  a data frame representing immune cell expression profiles. Each row represents an expression of a gene, and each column represents a different immune cell type. colnames contains the name of each immune cell type and the rownames includes the genes' symbol. The names of each immune cell type and the symbol of each gene should be unique. Any gene with missing expression values must be excluded.
#' @param geneExpr  a data frame representing RNA-seq or microarray gene-expression profiles of a given complex tissue. Each row represents an expression of a gene, and each column represents a different experimental sample. colnames contain the name of each sample and rownames includes the genes' symbol. The name of each individual sample and the symbol of each gene should be unique. Any gene with missing expression values should be excluded.
#' @param method  One of 'DCQ', 'SVMDECON', 'DeconRNASeq', 'proportionsInAdmixture', 'nnls' (DEFAULT: DCQ)
#' @param preNormalize  Set to TRUE to quantile normalize geneExpr to match refExpr (DEFAULT: TRUE)
#' @param verbose  Set to TRUE to echo the results of parameters (DEFAULT: TRUE)
#' @param ...  Parameters for estCellPercent.X (e.g. number_of_repeats for .DCQ)
#' @export
#' @return A matrix with cell type estimates for each samples
#' @examples
#' #This toy example 
#' library(ADAPTS)
#' fullLM22 <- ADAPTS::LM22[1:30, 1:4]
#' smallLM22 <- fullLM22[1:25,] 
#' 
#' cellEst <- estCellPercent(refExpr=smallLM22, geneExpr=fullLM22, preNormalize=FALSE, verbose=TRUE)
#' 
estCellPercent <- function(refExpr, geneExpr, preNormalize=TRUE, verbose=TRUE, method='DCQ', ...) {
  if (preNormalize == TRUE) {
    if(verbose==TRUE) {message('Quantile Normalizing geneExpr to match refExpr')}
    rns <- rownames(geneExpr)
    cns <- colnames(geneExpr)
    geneExpr <- preprocessCore::normalize.quantiles.use.target(as.matrix(geneExpr), as.numeric(as.matrix(refExpr)))
    rownames(geneExpr) <- rns
    colnames(geneExpr) <- cns
  }
  
  if(verbose==TRUE) {message(paste('Setting method to:', method))}
  if(method == 'DCQ') {
    cellEst <- estCellPercent.DCQ(refExpr = refExpr, geneExpr=geneExpr, ...)
  } else if (method == 'SVMDECON') {
    cellEst <- estCellPercent.svmdecon(refExpr = refExpr, geneExpr=geneExpr, ...)
  } else if (method == 'DeconRNASeq')  {
    cellEst <- estCellPercent.DeconRNASeq(refExpr = refExpr, geneExpr=geneExpr, marker_set=NULL, ...)
  } else if (method == 'proportionsInAdmixture') {
    cellEst <- estCellPercent.proportionsInAdmixture(refExpr = refExpr, geneExpr=geneExpr, marker_set=NULL, ...)
  } else if (method == 'nnls') {
    cellEst <- estCellPercent.nnls(refExpr = refExpr, geneExpr=geneExpr, ...)
  }
  return(cellEst) 
}
