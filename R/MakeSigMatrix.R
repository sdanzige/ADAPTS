# Ugly workaround to make foreach pass CRAN syntax check
#http://r.789695.n4.nabble.com/R-CMD-check-and-foreach-code-td4687660.html
globalVariables(c('fe_cType', 'fe_curGene'))

#' Use parallel missForest to impute missing values.
#' @description  This wrapper is helpful because missForest crashes if you have more cores than variables.
#'   This will default to no parellelization for Windows
#'
#'   newMatrix <- missForest.par(dataMat)
#'
#' @param dataMat  Columns are features, Rows examples. The data with NA values.  'xmis' in missForest
#' @param parallelize  split on 'forests' or 'variables' (DEFAULT: 'variables')
#'
#' @export
#' @return a matrix including imputed values
#' @examples
#' library(ADAPTS)
#' LM22 <- ADAPTS::LM22
#' LM22[2,3] <- as.numeric(NA) #Make some missing data to impute
#' LM22.imp <- missForest.par(LM22)
missForest.par <- function(dataMat, parallelize = "variables") {

  if (.Platform$OS.type == 'windows') {
    parallelize <- 'no'
    fixCores <- FALSE
  } else {
    fixCores <- ncol(dataMat) < getDoParWorkers()
  
    if (fixCores) {
      oldCores <- getDoParWorkers()
      if(ncol(dataMat) == 1) {
        parallelize <- 'no'
      } else {
        options(mc.cores = ncol(dataMat))
        options(cores = ncol(dataMat))
        doParallel::registerDoParallel(cores = ncol(dataMat))
      } #if(ncol(dataMat) == 1 ) {
    } #if (fixCores) {
  } #if (.Platform$OS.type == 'windows') {

  newMatrix <- try(missForest::missForest(dataMat, parallelize = parallelize)$ximp)

  if(inherits(x=newMatrix,'try-error')) {
    message('missForest error')
  }

  if(fixCores) {
    options(mc.cores = oldCores)
    options(cores = oldCores)
    doParallel::registerDoParallel(cores = oldCores);
  }

  return(newMatrix)
}

#' Make an Augmented Signature Matrix
#' @description With the ADAPTSdata packge, it will use the full LM22 data matrix and add a few 
#' additional genes to cover osteoblasts, osteoclasts, Plasma.memory, MM.  In many ways this is 
#' just a convenient wrapper for AugmentSigMatrix that calculates and caches a gList.
#'
#'
#' @param exprData  The gene express data to use to augment LM22, e.g. ADAPTSdata::addMGSM27
#' @param fullLM22  LM22 data with all genes.  Available in ADAPTSdata2::fullLM22
#' @param smallLM22  The small LM22 matrix, if it includes new cell types in exprData those will not be overwritten (DEFAULT: NULL, i.e. buildLM22plus(useLM22genes = TRUE)
#' @param plotToPDF  TRUE: pdf, FALSE: standard display (DEFAULT: TRUE)
#' @param condTol  The tolerance in the reconstruction algorithm.  1.0 = no tolerance, 1.05 = 5\% tolerance (DEFAULT: 1.01)
#' @param postNorm  Set to TRUE to normalize new signatures to match old signatures.  To Do: Redo Kappa curve? (DEFAULT: TRUE)
#' @param autoDetectMin Set to true to automatically detect the first local minima. GOOD PRELIMINARY RESULTS (DEAFULT: FALSE)
#' @param pdfDir  A fold to write the pdf file to if plotToPDF=TRUE (DEFAULT: tempdir())
#' @param oneCore Set to TRUE to disable parallelization (DEFAULT: FALSE)
#' @export
#' @return a cell type signature matrix
#' @examples
#' #This toy example treats the LM22 deconvolution matrix as if it were all of the data
#' #  For a real example, look at the vignette or comments in exprData, fullLM22, small LM22
#' library(ADAPTS)
#' fullLM22 <- ADAPTS::LM22[1:200, 1:8]
#' #Make a fake signature matrix out of 100 genes and the first 8 cell types
#' smallLM22 <- fullLM22[1:100, 1:8] 
#' 
#' #Make fake data representing two replicates of purified Mast.cells 
#' exprData <- ADAPTS::LM22[1:200, c("Mast.cells.resting","Mast.cells.activated")]
#' colnames(exprData) <- c("Mast.cells", "Mast.cells")
#' newSig <- remakeLM22p(exprData=exprData, fullLM22=fullLM22, smallLM22=smallLM22, 
#'     plotToPDF=FALSE, oneCore=TRUE)
remakeLM22p <- function(exprData, fullLM22, smallLM22=NULL, plotToPDF=TRUE, condTol = 1.01, postNorm=TRUE, autoDetectMin = FALSE, pdfDir=tempdir(), oneCore=FALSE) {
  exprData <- as.data.frame(exprData)

  if (is.null(smallLM22)) {
    smallLM22 <- ADAPTS::LM22
  }

  #Combine additonalMM data and the full LM22 dataset
  colnames(exprData) <- sub('.[0-9]+$', '', colnames(exprData))
  cNames <- c(colnames(fullLM22), colnames(exprData))

  #Problem 04-19-17 - If I just use the original LM22 data, there's too many NA, I end up adding 950 genes
  #  and 913 of them have NA values in the original data.  I will artifically limit the dataset to genes that I
  #  have data for
  ##rNames <- unique(c(rownames(fullLM22), rownames(exprData)))
  rNames.1 <- rownames(fullLM22)[apply(fullLM22, 1, function(x) {mean(is.na(x)) < 0.25})]
  rNames.2 <- rownames(exprData)[apply(exprData, 1, function(x) {mean(is.na(x)) < 0.25})]
  rNames <- rNames.1[rNames.1 %in% rNames.2]

  geneExpr <- cbind(fullLM22[rNames,], as.data.frame(exprData)[rNames,])
  colnames(geneExpr) <- sub('.[0-9]+$', '', colnames(geneExpr))

  fName <- paste('gList', paste(rev(unique(colnames(geneExpr))), collapse="_"), 'RData', sep='.')
  if(nchar(fName) > 240) { print('Truncating name list.  File name may not be unique') }
  fName <- paste0(strtrim(fName, 240),'.RData')  #Avoid too long filesnames, but introduct possible bug where two specs can generate the same file
  fName <- file.path(tempdir(), fName)
  
  if(!file.exists(fName)) {
    gList <- rankByT(geneExpr = geneExpr, qCut=0.3, oneCore=oneCore)
    save(gList, file=fName)
  } else {
    gList <- get(load(fName)[1])
  }

  #Normalize the new expression data against the full data?
  newMatData <- AugmentSigMatrix(origMatrix=smallLM22, fullData=fullLM22[rNames,], newData=exprData[rNames,], gList=gList, plotToPDF = plotToPDF, condTol = condTol, postNorm=postNorm, autoDetectMin=autoDetectMin, pdfDir=pdfDir)

  return(newMatData=as.data.frame(newMatData))
}

#' Make an augmented signature matrix
#' @description Build an augmented signature matrix from an initial signature matrix, source data, and a list of 
#' differentially expressed genes (gList).  The user might want to modify gList to make certain that particular 
#' genes are included in the matrix.  The algorithm will be to add one additional gene from each new cell type
#' Record the condition number, and plot those.  Will only consider adding rows shared by fullData and newData
#'
#'  newMatData <- AugmentSigMatrix(origMatrix, fullData, newData, gList)
#'
#' @param origMatrix  The original signature matrix
#' @param fullData  The full data for the signature matrix
#' @param newData  The new data to add signatures from
#' @param gList  The ordered list of genes from running rankByT() on newData. NOTE: best genes at the bottom!!
#' @param nGenes  The number of additional genes to consider (DEFAULT: 1:100)
#' @param plotToPDF  Plot the output condition numbers to a pdf file. (DEFAULT: TRUE)
#' @param imputeMissing  Set to TRUE to impute missing values. NOTE: adds stoachasiticity (DEFAULT: TRUE)
#' @param condTol  Setting higher tolerances will result in smaller numbers extra genes. 1.00 minimizes compliment number (DEFAULT: 1.00)
#' @param postNorm  Set to TRUE to normalize new signatures to match old signatures.  (DEFAULT: FALSE)
#' @param minSumToRem  Set to non-NA to remove any row with the sum(abs(row)) < minSumToRem (DEFAULT: NA)
#' @param addTitle  An optional string to add to the plot and savefile (DEFAULT: NULL)
#' @param autoDetectMin Set to true to automatically detect the first local minima. GOOD PRELIMINARY RESULTS (DEAFULT: FALSE)
#' @param calcSpillOver Use the training data to calculate a spillover matrix (DEFAULT: FALSE)
#' @param pdfDir  A fold to write the pdf file to if plotToPDF=TRUE (DEFAULT: tempdir())
#' @export
#' @return an augmented cell type signature matrix
#' @examples
#' #This toy example treats the LM22 deconvolution matrix as if it were all of the data
#' #  For a real example, look at the vignette or comments in exprData, fullLM22, small LM22
#' library(ADAPTS)
#' fullLM22 <- ADAPTS::LM22[1:200, 1:8]
#' #Make a fake signature matrix out of 100 genes and the first 8 cell types
#' smallLM22 <- fullLM22[1:100, 1:8] 
#' 
#' #Make fake data representing two replicates of purified Mast.cells 
#' exprData <- ADAPTS::LM22[1:200, c("Mast.cells.resting","Mast.cells.activated")]
#' colnames(exprData) <- c("Mast.cells", "Mast.cells")
#' 
#' #Fake source data with replicates for all purified cell types.
#' #  Note in this fake data set, many cell types have exactly one replicate
#' fakeAllData <- cbind(fullLM22, as.data.frame(exprData)) 
#' gList <- rankByT(geneExpr = fakeAllData, qCut=0.3, oneCore=TRUE)
#' 
#' newSig <- AugmentSigMatrix(origMatrix=smallLM22, fullData=fullLM22, newData=exprData, 
#'     gList=gList, plotToPDF=FALSE)
AugmentSigMatrix <- function(origMatrix, fullData, newData, gList, nGenes=1:100, plotToPDF=TRUE, imputeMissing=TRUE, condTol=1.01, postNorm=FALSE, minSumToRem=NA, addTitle=NULL, autoDetectMin=FALSE, calcSpillOver=FALSE, pdfDir=tempdir()) {
  origMatrix <- as.data.frame(origMatrix)
  if(autoDetectMin == TRUE) {
    if(!is.null(addTitle)) {
      addTitle <- paste('Auto', addTitle, sep='.')
    } else {
      addTitle <- 'Auto'
    }
  } #if(autoDetectMin == TRUE) {

  if(any(!(colnames(origMatrix) %in% c(colnames(fullData), colnames(newData))))) {
    missingData <- colnames(origMatrix)[!(colnames(origMatrix) %in% c(colnames(fullData), colnames(newData)))]
    print(paste(missingData, 'in origMatrix but not fullData or newData'))
    return(NULL)
  }

  #Make sure that all columns in the newData already exist in the origMatrix

  missingCols <- unique(colnames(newData)[!(colnames(newData) %in% colnames(origMatrix))])
  if(length(missingCols) > 0) {
    augData <- matrix(as.numeric(NA), ncol=length(missingCols), nrow=nrow(origMatrix),
                  dimnames=list(rownames(origMatrix), missingCols))
    olGenes <- rownames(origMatrix)[rownames(origMatrix) %in% rownames(newData)]

    for(missingCol in missingCols) {
      newExp <- apply(newData[olGenes, colnames(newData) == missingCol, drop=FALSE],1, stats::median, na.rm=TRUE)
      augData[names(newExp),missingCol] <- newExp
    } #for(missingCol in missingCols) {
    origMatrix <- cbind(origMatrix, augData)
  }

  unAug <- origMatrix[colnames(origMatrix)[colnames(origMatrix) %in% colnames(fullData)]]
  cNums.orig <- kappa(unAug)

  #Baseline condition number.
  newCtypes <- unique(colnames(newData))
  if(any(!(newCtypes %in% names(gList)))) {
    outStr <- paste(paste(newCtypes[!(newCtypes %in% names(gList))], sep=', '), 'missing')
    print(outStr)
    return(NULL)
  }
  gList <- gList[newCtypes]


  allGenes <- rownames(fullData)[rownames(fullData) %in% rownames(newData)]
  fullData <- fullData[allGenes,]
  newData <- newData[allGenes,]

  if(any(is.na(origMatrix))) {
    if(imputeMissing == TRUE) {
      origMatrix.imp <- t(missForest.par(t(origMatrix)))
    } else {
      remBool <- apply(origMatrix, 1, function(x){any(is.na(x))})
      origMatrix.imp <- origMatrix[!remBool,]
    }
  } else {
    origMatrix.imp <- origMatrix
  }
  cNums.new <- kappa(origMatrix.imp)

  selGenes <- list()
  newMatrix <- origMatrix

  for(gNum in nGenes) {
    newGenes <- NULL
    for (cType in newCtypes) {
      gNames <- rev(rownames(gList[[cType]]))
      #New matrix will be iteratively built and code will make sure that new genes are not in matrix
      newGenes <- c(newGenes, gNames[which(!(gNames %in% rownames(newMatrix)))[1]])
    } #for (cType in newCtypes) {
    newGenes <- unique(newGenes)

    if(all(is.na(newGenes))) { next; }
    newGenes <- newGenes[!is.na(newGenes)]
    
    augData.new <- cbind(fullData[newGenes,,drop=FALSE], newData[newGenes,,drop=FALSE])
    augData <- apply(augData.new, 1, function(x) {
      tapply(x, colnames(augData.new), stats::median, na.rm=TRUE)
    })
    augData <- t(augData)

    newMatrix <- rbind(newMatrix, augData[, colnames(newMatrix), drop=FALSE])
    selGenes[[as.character(gNum)]] <- newGenes
  } #for(gNum in nGenes) {

  #Impute the full matrix and back-calculate the kappa
  if(any(is.na(newMatrix))) {
    if(imputeMissing == TRUE) {
      impMatrix <- t(missForest.par(t(newMatrix)))
    } else {
      remBool <- apply(newMatrix, 1, function(x){any(is.na(x))})
      impMatrix <- newMatrix[!remBool,]
    }
    
  } else {
    impMatrix <- newMatrix
  }
  impMatrix[impMatrix > max(origMatrix, na.rm=TRUE)] <- max(origMatrix, na.rm=TRUE)
  impMatrix[impMatrix < min(origMatrix, na.rm=TRUE)] <- min(origMatrix, na.rm=TRUE)

  cNums <- cNums.new
  nGenes <- nrow(origMatrix)
  for (i in 1:length(selGenes)) {
    newGenes <- unlist(selGenes[1:i])
    newGenes <- newGenes[!is.na(newGenes)]
    curMatGenes <- c(rownames(origMatrix), newGenes)
    curMat <- impMatrix[curMatGenes[curMatGenes %in% rownames(impMatrix)],]
    cNums <- c(cNums, kappa(curMat))
    nGenes <- c(nGenes, nrow(curMat))
  }

  #Plot the results

  smData <- stats::smooth(cNums)
  #smData <- ges(cNums)

  #Modification: use standard smoothing and then pick point with smallest # features that is no more than 1% higher

  if(autoDetectMin == TRUE) {
    smData2 <- stats::predict(stats::smooth.spline(smData))$y
    mins <- quantmod::findPeaks(-smData2)
    bestMin <- mins[1]
    mVal <- smData[bestMin]
    if(is.na(mVal)) {
      message('autoDetectMin failed, reverting to absolute min')
    }
  } 
  if(autoDetectMin == FALSE || is.na(mVal)) {
    #mVal will be NA if autoDetectMin fails.
    mVal <- min(smData)
    bestMin <- which(smData == mVal)[1]
  }

  smallMin <- which(smData <= condTol*mVal)[1]

  #newGenes <- unlist(selGenes[1:which.min(smData)])
  #newGenes <- unique(unlist(selGenes[1:smallMin]))
  newGenes <- unique(unlist(selGenes[1:(smallMin-1)]))
  if(imputeMissing == TRUE) {
    sigMatrix <- impMatrix[rownames(impMatrix) %in% c(rownames(origMatrix), newGenes),]
  } else {
    sigMatrix <- newMatrix[rownames(newMatrix) %in% c(rownames(origMatrix), newGenes),]
    sigMatrix <- sigMatrix[apply(sigMatrix, 1, function(x){!any(is.na(x))}),]
  } #if(imputeMissing == TRUE) {

  if (postNorm==TRUE) {
    #message(paste('Pre-normalization Kappa:', kappa(sigMatrix)))
    newPartBool <- colnames(sigMatrix) %in% colnames(newData)
    if(sum(!newPartBool)==0) {
      renormNewPart <- preprocessCore::normalize.quantiles(x=as.matrix(sigMatrix[,newPartBool]))
    } else {
      renormNewPart <- preprocessCore::normalize.quantiles.use.target(x=as.matrix(sigMatrix[,newPartBool]), target=as.vector(as.matrix(sigMatrix[,!newPartBool])))
    } #if(sum(!newPartBool)) {
    sigMatrix[,newPartBool] <- renormNewPart
    #message(paste('Post-normalization Kappa:', kappa(sigMatrix)))
  }

  titleStr <- paste('Augmenting Signature Matrix ( tol =', condTol, ')\n# Cell-types:', ncol(unAug), '->', ncol(newMatrix),
                  '| # Genes:', nrow(unAug), '->', nrow(sigMatrix))
  if(!is.null(addTitle)) { titleStr <- paste(addTitle, titleStr) }

  if(plotToPDF == TRUE) {
    pdfFile <- paste('AugmentSigMatrix', condTol, Sys.Date(), 'pdf', sep='.')
    pdfFile <- file.path(pdfDir, pdfFile)
    if(!is.null(addTitle)) { pdfFile <- sub('.pdf$', paste0('.', addTitle, '.pdf'), pdfFile) }
    if(imputeMissing) { pdfFile <- sub('.pdf', '.impute.pdf', pdfFile) }
    if(postNorm) { pdfFile <- sub('.pdf', '.postNorm.pdf', pdfFile) }
    grDevices::pdf(pdfFile)
  }

  legText <- c('Unagumented Signature Matrix', 'Minimum Smoothed Condition Number', 'Best Augmented Signature Matrix')
  pchs <- c('o', 'x', 'x')
  cols <- c('red', 'purple', 'blue')
  ylims <- c(min(cNums.orig,cNums,kappa(sigMatrix))*0.95, max(cNums.orig, cNums)*1.05)
  graphics::plot(x=nGenes, y=cNums,
     xlab='Number of Genes', ylab='Condition Number (lower is more stable)',
     main=titleStr, ylim=ylims)
  graphics::points(x=nrow(unAug), y=cNums.orig, col='red', pch='o', cex=1.5 )
  graphics::lines(x=nGenes, y=smData, col='green')
  graphics::points(x=nGenes[bestMin], y=cNums[bestMin], col='purple', pch='x', cex=1.5)
  graphics::points(x=nGenes[smallMin], y=cNums[smallMin], col='blue', pch='x', cex=1.5)
  if (postNorm==TRUE) {
    legText <- c(legText, 'Post-Normalized')
    pchs <- c(pchs, 'x')
    cols <- c(cols, 'violetred')
    graphics::points(x=nGenes[smallMin], y=kappa(sigMatrix), col='violetred', pch='x', cex=1.5)
  }
  graphics::legend('topright', legend=legText, pch=pchs, col=cols)

  origTitle <- 'Original Matrix'
  if(!is.null(addTitle)) { origTitle <- paste(origTitle, titleStr) }

  pheatmap::pheatmap(origMatrix, main=origTitle, fontsize_row = 4)
  pheatmap::pheatmap(sigMatrix, main=titleStr, fontsize_row = 4)

  if(!is.na(minSumToRem)) {
    sigMatrix <- sigMatrix[rowSums(abs(sigMatrix)) > minSumToRem,]
    titleStr <- paste('Augmenting Signature Matrix ( tol =', condTol, 'filter =', minSumToRem, ')\n# Cell-types:',
                      ncol(unAug), '->', ncol(sigMatrix), '| # Genes:', nrow(unAug), '->', nrow(sigMatrix))
    if(!is.null(addTitle)) { titleStr <- paste(addTitle, titleStr) }
    pheatmap::pheatmap(sigMatrix, main=titleStr, fontsize_row = 4)
  }

  if(plotToPDF == TRUE) {
    grDevices::dev.off()
  }

  #Remove genes one at a time until the condition number stops removing
  #note, this is a bad idea, it pretty much monotonically shrinks.
  postShrink <- FALSE
  if(postShrink==TRUE) {
    curComp <- kappa(sigMatrix)
    names(curComp) <- ""

    newSigMatrix <- sigMatrix
    for(i in 1:(nrow(sigMatrix)-1)) {
      #Remove genes one at a time and determine which will minimize the condition number
      newComps <- foreach(fe_curGene = rownames(newSigMatrix), .combine=c) %dopar% {
        kappa(newSigMatrix[rownames(newSigMatrix)!=fe_curGene,])
      }
      names(newComps) <- rownames(newSigMatrix)

      remGene <- names(newComps)[which.min(newComps)]
      curComp <- c(curComp, newComps[remGene])

      newSigMatrix <- newSigMatrix[rownames(newSigMatrix)!=remGene,]
    } #for(i in 1:(nrow(sigMatrix)-1)) {
  }

  #Optional, calculate a spillover matrix
  if(calcSpillOver == TRUE) {
    curAllDat <- cbind(fullData, newData)
    res <- buildSpilloverMat(sigMatrix, geneExpr=curAllDat)
    pheatmap::pheatmap(res, cluster_rows = FALSE, cluster_cols = FALSE, main=paste0('Spillover Matrix: ', nrow(res), ' Cell Types'))
    rv <- list(sigMatrix=sigMatrix, spillOver=res)
  } else {
    rv <- sigMatrix
  }

  return(rv)
}

#' Rank genes for each cell type
#' @description Use a t-test to rank to features for each cell type
#'
#' gList <- rankByT(geneExpr, qCut=0.3)
#'
#' @param geneExpr  The gene expression data
#' @param qCut  (DEFAULT: 0.3)
#' @param oneCore Set to TRUE to disable paralellization (DEFAULT: FALSE)
#' @param secondPval Set to TRUE to use p-Values as a second sort criteria (DEFAULT: TRUE)
#' @param remZinf Set to TRUE to remove any ratio with zero or infinity.  Good for scRNAseq. (DEFAULT: FALSE)
#' @param reqRatGT1 Set to TRUE to remove any gene with a ratio with less than 1.  Good for scRNAseq. (DEFAULT: FALSE)
#' @export
#' @return a list of cell types with data frames ranking genes
#' @examples
#' #This toy example treats the LM22 deconvolution matrix as if it were all of the data
#' #  For a real example, look at the vignette or comments in exprData, fullLM22, small LM22
#' library(ADAPTS)
#' fullLM22 <- ADAPTS::LM22[1:200, 1:8]
#' #Make a fake signature matrix out of 100 genes and the first 8 cell types
#' smallLM22 <- fullLM22[1:100, 1:8] 
#' 
#' #Make fake data representing two replicates of purified Mast.cells 
#' exprData <- ADAPTS::LM22[1:200, c("Mast.cells.resting","Mast.cells.activated")]
#' colnames(exprData) <- c("Mast.cells", "Mast.cells")
#' 
#' #Fake source data with replicates for all purified cell types.
#' #  Note in this fake data set, many cell types have exactly one replicate
#' fakeAllData <- cbind(fullLM22, as.data.frame(exprData)) 
#' gList <- rankByT(geneExpr = fakeAllData, qCut=0.3, oneCore=TRUE, reqRatGT1=FALSE)
rankByT <- function(geneExpr, qCut=0.3, oneCore=FALSE, secondPval=TRUE, remZinf=FALSE, reqRatGT1=FALSE) {
  colnames(geneExpr) <- sub(".[0-9]+$", '', colnames(geneExpr)) #Strip any trailing numbers added by make.names()
  cTypes <- unique(colnames(geneExpr))

  if(length(cTypes) > 2 & oneCore==FALSE) {
    gList <- foreach (fe_cType = cTypes) %dopar% {
      print(fe_cType)
      isType <- colnames(geneExpr) == fe_cType
      geneExpr.cur <- geneExpr
      
      tRes <- lapply(rownames(geneExpr.cur), function(x) {
        rv <- try(stats::t.test(geneExpr.cur[x,isType], geneExpr.cur[x,!isType], na.action=stats::na.omit), silent=TRUE)
        if(inherits(rv, 'try-error')) {rv <- list(estimate=c(1,1), statistic=0, p.value=1)}
        return(rv)
      }) #tRes <- mclapply(gNames, function(x) {

      geneDF <- do.call(rbind, lapply(tRes, function(x) {data.frame(rat=x$estimate[1]/x$estimate[2], t=x$statistic, pVal=x$p.value)}))
      rownames(geneDF) <- rownames(geneExpr.cur)
      geneDF$qVal <- stats::p.adjust(geneDF$pVal, method = 'fdr')

      if(secondPval==TRUE) {
        geneDF <- geneDF[order(abs(log2(geneDF$rat)), -1*log(geneDF$pVal)),]
      } else {
        geneDF <- geneDF[order(abs(log2(geneDF$rat))),]
      }
      geneDF <- geneDF[geneDF$qVal <= qCut,]

      geneDF <- geneDF[!is.na(geneDF$rat), ]
      #gList[[cType]] <- geneDF
      
      if(reqRatGT1==TRUE) { geneDF <- geneDF[geneDF$rat>1, ,drop=FALSE] }
      if (remZinf==TRUE) { 
        isZ <- geneDF$rat == 0
        isInf <- is.infinite(geneDF$rat)
        geneDF <- geneDF[!(isZ | isInf), ,drop=FALSE]
      }
      
      return(geneDF)
    } #for (cType in unique(colnames(geneExpr))) {
    names(gList) <- cTypes
  } else {
    if(length(cTypes) <= 2) {cTypes <- cTypes[2]}
    gList <- lapply(cTypes, function(fe_cType) {
      print(fe_cType)
      isType <- colnames(geneExpr) == fe_cType
      
      #if (remZinf) {
      #  isZ <- rowSums(geneExpr[,isType,drop=FALSE]) == 0
      #  notZ <- rowSums(geneExpr[,!isType,drop=FALSE]) == 0
      #  remBool <- isZ | notZ
      #  geneExpr.cur <- geneExpr[!remBool,]
      #} else {
        geneExpr.cur <- geneExpr
      #}
      
      if(oneCore==TRUE) {
        tRes <- lapply(rownames(geneExpr.cur), function(x) {
          rv <- try(stats::t.test(geneExpr.cur[x,isType], geneExpr.cur[x,!isType], na.action=na.omit), silent=TRUE)
          if(inherits(rv, 'try-error')) {rv <- list(estimate=c(1,1), statistic=0, p.value=1)}
          return(rv)
        }) #tRes <- lapply(gNames, function(x) {
      } else {
        tRes <- parallel::mclapply(rownames(geneExpr.cur), function(x) {
          rv <- try(stats::t.test(geneExpr.cur[x,isType], geneExpr.cur[x,!isType], na.action=na.omit), silent=TRUE)
          if(inherits(rv, 'try-error')) {rv <- list(estimate=c(1,1), statistic=0, p.value=1)}
          return(rv)
        }) #tRes <- mclapply(gNames, function(x) {
      } #if(oneCore==TRUE) {
        
      geneDF <- do.call(rbind, lapply(tRes, function(x) {data.frame(rat=x$estimate[1]/x$estimate[2], t=x$statistic, pVal=x$p.value)}))
      rownames(geneDF) <- rownames(geneExpr.cur)
      geneDF$qVal <- stats::p.adjust(geneDF$pVal, method = 'fdr')

      if(secondPval==TRUE) {
        geneDF <- geneDF[order(abs(log2(geneDF$rat)), -1*log(geneDF$pVal)),]
      } else {
        geneDF <- geneDF[order(abs(log2(geneDF$rat))),]
      }
      geneDF <- geneDF[geneDF$qVal <= qCut,]

      geneDF <- geneDF[!is.na(geneDF$rat), ]
      #gList[[cType]] <- geneDF
      
      if(reqRatGT1==TRUE) { geneDF <- geneDF[geneDF$rat>1, ,drop=FALSE] }
      if (remZinf==TRUE) { 
        isZ <- geneDF$rat == 0
        isInf <- is.infinite(geneDF$rat)
        geneDF <- geneDF[!(isZ | isInf), ,drop=FALSE]
      }
      
      return(geneDF)
    }) #for (cType in unique(colnames(geneExpr))) {
    names(gList) <- cTypes
  } #if(length(cTypes) > 1) {

  return(gList)
}

#' Load MGSM27
#' @description Load the MGSM27 signature matrix
#'
#' @export
#' @return  The MGSM27 signature matrix from Identifying a High-risk Cellular Signature in the Multiple Myeloma Bone Marrow Microenvironment
#' @examples
#' MGSM27 <- loadMGSM27()
loadMGSM27 <- function() {
  MGSM27 <- ADAPTS::MGSM27
  return(MGSM27)
}

#' LM22 look up table
#' @description Load a map of cell type names
#'
#' @export
#' @return a map of cell types names
#' @examples
#' cellMap <- getLM22cells()
getLM22cells <- function() {
  mapTypes <- list('naive B-cells',       # B.cells.naive
                   'Memory B-cells',      # B.cells.memory
                   'Plasma cells',        # Plasma.cells
                   c('CD8+ T-cells', 'CD8+ Tcm', 'CD8+ Tem', 'CD8+ naive T-cells'),        # T.cells.CD8
                   'CD8+ T-cells',
                   'CD4+ naive T-cells',  # T.cells.CD4.naive
                   'CD4+ memory T-cells', # T.cells.CD4.memory.resting   NOTE divide by 2
                   'CD4+ memory T-cells', # T.cells.CD4.memory.activated NOTE divide by 2
                   c('CD4+ T-cells', 'CD4+ Tem', 'CD4+ Tcm', 'Th1 cells', 'Th2 cells'),                    # T.cells.follicular.helper
                   'Tregs',               # T.cells.regulatory..Tregs
                   'Tgd cells',           # T.cells.gamma.delta
                   'NK cells',            # NK.cells.resting             NOTE divide by 2
                   'NK cells',            # NK.cells.activated           NOTE divide by 2
                   'Monocytes',           # Monocytes
                   'Macrophages',         # Macrophages.M0
                   'Macrophages M1',      # Macrophages.M1
                   'Macrophages M2',      # Macrophages.M2
                   c('cDC','pDC','DC'),   # Dendritic.cells.resting
                   'aDC',                 # Dendritic.cells.activated
                   'Mast cells',          # Mast.cells.resting           NOTE divide by 2
                   'Mast cells',          # Mast.cells.activated         NOTE divide by 2
                   'Eosinophils',         # Eosinophils
                   'Neutrophils',         # Neutrophils
                   'Plasma cells',                    # MM.plasma.cell
                   'Macrophages',                    # osteoclast
                   'Plasma cells',                    # PlasmaMemory
                   '')                    # Other, need something special here
  mapTypes <- unique(unlist(mapTypes))
  return(mapTypes)
}

#' LM22 to xCell LUT
#' @description Load the LM22 xCell map
#'
#' @export
#' @return A map between xCell cell type names and LM22 cell type names
#' @examples
#' xcellMap <- loadModMap()
loadModMap <- function() {
  modMap <- rbind(c('naive B-cells', 'B.cells.naive'),
                  c('Memory B-cells','B.cells.memory'),
                  c('Plasma cells','Plasma.cells'),
                  c('CD8+ T-cells','T.cells.CD8'),
                  c('CD4+ naive T-cells','T.cells.CD4.naive'),
                  c('CD4+ memory T-cells', 'T.cells.CD4.memory.resting'),
                  c('CD4+ memory T-cells', 'T.cells.CD4.memory.activated'),
                  c('CD4+ Tcm','T.cells.follicular.helper'),
                  c('Tregs','T.cells.regulatory..Tregs.'),
                  c('Tgd cells','T.cells.gamma.delta'),
                  c('NK cells',  'NK.cells.resting'),
                  c('NK cells',  'NK.cells.activated'),
                  c('Monocytes','Monocytes'),
                  c('Macrophages', 'Macrophages.M0'),
                  c('Macrophages M1','Macrophages.M1'),
                  c('Macrophages M2', 'Macrophages.M2'),
                  c('iDC','Dendritic.cells.resting'),
                  c('aDC', 'Dendritic.cells.activated'),
                  c('Mast cells', 'Mast.cells.resting'),
                  c('Mast cells', 'Mast.cells.activated'),
                  c('Eosinophils', 'Eosinophils'),
                  c('Neutrophils', 'Neutrophils'),
                  c('Adipocytes',  'adipocyte'),
                  c('Plasma cells', 'MM.plasma.cell'),
                  c('Osteoblast',  'osteoblast'),
                  c('Plasma cells', 'PlasmaMemory'))
  return(modMap)
}

#  Function removed to pass R CMD check --as-cran because xCell is on GitHub not CRAN. 
#
# #' Use xCellSignifcanceBetaDist in xCell to estimate the probability that a cell type is in a sample.
# #'
# #' cellPvals <- estxCellSig(geneExpr.pbmc, rnaseq = FALSE)
# #'
# #' @param geneExpr.pbmc  The gene expression data
# #' @param rnaseq  Set to TRUE if the data is RNAseq data (DEFAULT: FALSE)
# #' @export
# #' @return cell type p-values
# estxCellSig <- function(geneExpr.pbmc, rnaseq = FALSE) {

#  if(!require('xCell')) {
#    message('This function requires the xCell package')
#    message('It can be found here: https://github.com/dviraran/xCell')
#    return(NULL)
#  }
#  xCells <- xCellAnalysis(geneExpr.pbmc,rnaseq = rnaseq)
#  xCellSigs <- xCellSignifcanceBetaDist(xCells, rnaseq = rnaseq)

  #Map to LM22 where possible.
#  xCellMap <- loadModMap()
#  xCellMap <- xCellMap[xCellMap[,1] %in% rownames(xCellSigs),]

#  mappable <- xCellSigs[xCellMap[,1],]
#  rownames(mappable) <- xCellMap[,2]

#  notmappable <- xCellSigs[!(rownames(xCellSigs) %in% xCellMap[,1]),]
#  cellPvals <- rbind(mappable, notmappable)

#  return(cellPvals)
#}

#  Function removed to pass R CMD check --as-cran because xCell is on GitHub not CRAN
#
# #' Deconvolve PBMCs and see what you get.
# #'   This requires the ADAPTSdata packages from CRAN
# #'   https://github.com/sdanzige/ADAPTSdata
# #'   https://github.com/sdanzige/ADAPTSdata2
# #'
# #'
# #' matData <- deconvolvePBMC(geneExpr.pbmc=ADAPTSdata::PBMC)
# #'
# #' @param geneExpr.pbmc  The PBMC data to look at, for example: ADAPTSdata::PBMC
# #' @param refExpr  The reference matrix (DEFAULT: NULL, ie load LM22)
# #' @param deconName  The deconvolved matrix name (DEFAULT: 'LM22')
# #' @param dsName  The sample name (DEFAULT: 'PBMC')
# #' @param decons  Set to include already deconvolved values (DEFAULT: NULL)
# #' @param incXcell  Set to TRUE to include xCell estimates of the fraction of samples with this cell type (DEFAULT: TRUE)
# #' @param plotIT  Set to TRUE to plot it  (DEFAULT: FALSE)
# #' @export
# #' @return cell type estimates for each sample
# deconvolvePBMC <- function(geneExpr.pbmc, refExpr=NULL, deconName="LM22", dsName='PBMC', decons=NULL, incXcell=TRUE, plotIT=FALSE) {
#   if(is.null(decons)) {
#     if(is.null(refExpr)) {
#       refExpr <- utils::read.csv(gzfile('/GitHub/Deconvolution/LM22.csv.gz'), row.names = 1)
#     }

#     decons <- estCellPercent.DCQ(refExpr, geneExpr.pbmc, marker_set = NULL, number_of_repeats = 10)
#   } #if(!is.null(decons)) {

  #Plot it
#   means <- apply(decons, 1, base::mean, na.rm=TRUE)
#   sds <- apply(decons, 1, stats::sd, na.rm=TRUE)
#   medians <- apply(decons, 1, stats::median, na.rm=TRUE)

  #Also test with xCell
#   if(incXcell == TRUE) {
#     cellPvals <- estxCellSig(geneExpr.pbmc)
#     fracSig <- apply(cellPvals, 1, mean)
#     fracSig.inc <- fracSig[rownames(decons)[rownames(decons) %in% names(fracSig)]]
#     inBoth <- rownames(decons)[rownames(decons) %in% names(fracSig.inc)]
#     for (curBoth in inBoth) {
#       i <- which(rownames(decons) == curBoth)
#       rownames(decons)[i] <- paste0(rownames(decons)[i], ' (', round(100*fracSig.inc[curBoth]), '%)')
#     }
#   } #if(incXcell == TRUE) {

#   for(iter in 1:2) {
#     if(iter==2 && plotIT==TRUE) { grDevices::pdf(paste('deconvolvePBMC', deconName, dsName, Sys.Date(), 'pdf', sep='.')) }

#     titleStr <- paste('DCQ', deconName, 'plot of', dsName)
#     newMax <- max(means) + max(sds) * 1.1

#     graphics::par(mar=c(12,4.1,4.1,2.1))
#     graphics::plot(x=means, ylim=c(0, newMax), axes=FALSE, xlab=NA, main=titleStr)
#     graphics::points(medians, col='red',pch='x')
#     graphics::par(las=2)
#     graphics::axis(1, at=1:nrow(decons), labels=rownames(decons), cex.axis=0.75)
#     graphics::axis(2, labels=TRUE) #default way
#     graphics::box()

#     graphics::segments(x0=1:nrow(decons), y0=means-sds, y1=means+sds)
#     epsilon = 0.2
#     graphics::segments(x0=(1:nrow(decons))-epsilon, means-sds ,(1:nrow(decons))+epsilon, means-sds)
#     graphics::segments(x0=(1:nrow(decons))-epsilon, means+sds ,(1:nrow(decons))+epsilon, means+sds)

#     if(incXcell == TRUE) {
#       graphics::legend('topleft', legend=c('median','() xCell %'), pch=c('x', ''), col=c('red','black'))
#     } else {
#       graphics::legend('topleft', legend=c('median'), pch=c('x'), col=c('red'))
#     }
#     if(iter==2 && plotIT==TRUE) { grDevices::dev.off() }
#   } #for(iter = 1:2) {

#   invisible(decons)
# }
