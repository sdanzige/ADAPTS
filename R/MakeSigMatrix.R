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
    fixCores <- ncol(dataMat) < foreach::getDoParWorkers()
    
    if (fixCores) {
      oldCores <- foreach::getDoParWorkers()
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
#' @param cache_gList Set to TRUE to cache slow gList calculations (DEFAULT: TRUE)
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
#' #Make fake data representing two replicates of purified Mast.cells types 
#' exprData <- ADAPTS::LM22[1:200, c("Mast.cells.resting","Mast.cells.activated")]
#' colnames(exprData) <- c("Mast.cells", "Mast.cells")
#' newSig <- remakeLM22p(exprData=exprData, fullLM22=fullLM22, smallLM22=smallLM22, 
#'     plotToPDF=FALSE, oneCore=TRUE, cache_gList=FALSE)
remakeLM22p <- function(exprData, fullLM22, smallLM22=NULL, plotToPDF=TRUE, condTol = 1.01, postNorm=TRUE, autoDetectMin = FALSE, pdfDir=tempdir(), oneCore=FALSE, cache_gList=TRUE) {
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
  
  if(file.exists(fName) & cache_gList == TRUE) {
    gList <- get(load(fName)[1])
  } else {
    gList <- rankByT(geneExpr = geneExpr, qCut=0.3, oneCore=oneCore)
    if(cache_gList == TRUE) { save(gList, file=fName, compress = TRUE) }
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
#' @param plotIt  Set to FALSE to suppress non-PDF plotting (DEFAULT: TRUE)
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
AugmentSigMatrix <- function(origMatrix, fullData, newData, gList, nGenes=1:100, plotToPDF=TRUE, imputeMissing=TRUE, condTol=1.01, postNorm=FALSE, minSumToRem=NA, addTitle=NULL, autoDetectMin=FALSE, calcSpillOver=FALSE, pdfDir=tempdir(), plotIt=TRUE) {
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
  cNums.orig <- kappa(as.matrix(unAug))
  
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
  cNums.new <- kappa(as.matrix(origMatrix.imp))
  
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
    cNums <- c(cNums, kappa(as.matrix(curMat)))
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
  
  if(plotIt == TRUE) {
    legText <- c('Unagumented Signature Matrix', 'Minimum Smoothed Condition Number', 'Best Augmented Signature Matrix')
    pchs <- c('o', 'x', 'x')
    cols <- c('red', 'purple', 'blue')
    ylims <- c(min(cNums.orig,cNums,kappa(as.matrix(sigMatrix)))*0.95, max(cNums.orig, cNums)*1.05)
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
      graphics::points(x=nGenes[smallMin], y=kappa(as.matrix(sigMatrix)), col='violetred', pch='x', cex=1.5)
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
  } #if(plotIt = TRUE) {
  
  if(plotToPDF == TRUE) {
    grDevices::dev.off()
  }
  
  #Remove genes one at a time until the condition number stops removing
  #note, this is a bad idea, it pretty much monotonically shrinks.
  postShrink <- FALSE
  if(postShrink==TRUE) {
    curComp <- kappa(as.matrix(sigMatrix))
    names(curComp) <- ""
    
    newSigMatrix <- sigMatrix
    for(i in 1:(nrow(sigMatrix)-1)) {
      #Remove genes one at a time and determine which will minimize the condition number
      newComps <- foreach(fe_curGene = rownames(newSigMatrix), .combine=c) %dopar% {
        kappa(as.matrix(newSigMatrix[rownames(newSigMatrix)!=fe_curGene,]))
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
  colnames(geneExpr) <- sub("\\.[0-9]+$", '', colnames(geneExpr)) #Strip any trailing numbers added by make.names()
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

#New functions to work into ADAPTS
#' Plot condition numbers
#' @description  Plot the condition numbers during the growing and shrinking of signature matrices.
#'
#'    bonusPoints <- data.frame(legText = c('Unagumented Signature Matrix', 'Minimum Smoothed Condition Number', 'Best Augmented Signature Matrix'), 
#'                                pchs = c('o', 'x', 'x'), 
#'                                cols = c('red', 'purple', 'blue'), 
#'                                kappa = c(10, 15, 20), 
#'                                nGene = c(5, 10, 15))
#'
#' @param kappas  The condition numbers to plot
#' @param nGenes  The number of genes associated with each kapp
#' @param smData  Smoothed data to plot as a green line (DEFAULT: NULL)
#' @param titleStr  The title of the plot (DEFAULT: 'Shrink Signature Matrix')
#' @param bonusPoints  Set to plot additional points on the plot, see description (DEFAULT: NULL)
#' @param maxCond  Cap the condition number to maxCond (DEFAULT: 100)
#'
#' @export
#' @return a matrix including imputed values
#' @examples
#' nGenes <- 1:300
#' kappas <- log(abs(nGenes-250))
#' kappas[is.infinite(kappas)] <- 0
#' kappas <- kappas+runif(300, 0, 1)
#' smData <- stats::smooth(kappas)
#' bonusPoints <- data.frame(legText = 'Minimum Smoothed ', pchs='x', cols='purple', 
#' kappa=min(smData), nGenes=nGenes[which.min(smData)])
#' plotKappas(kappas=kappas, nGenes=nGenes, smData=smData, bonusPoints=bonusPoints, maxCond=100)
#' 
plotKappas <- function(kappas, nGenes, smData=NULL, titleStr='Shrink Signature Matrix', bonusPoints=NULL, maxCond=100) {
  
  testKappa <- kappas
  if(!is.null(bonusPoints)) {testKappa <- c(testKappa, bonusPoints$kappa)}
  
  if(any(kappas > maxCond)) {
    message(paste('Capping condition number to', maxCond))
    kappas[kappas > maxCond] <- maxCond
  }
  
  ylims <- c(min(testKappa)*0.95, max(testKappa)*1.05)
  graphics::plot(x=nGenes, y=kappas,
                 xlab='Number of Genes', ylab='Condition Number (lower is more stable)',
                 main=titleStr, ylim=ylims)
  
  if(!is.null(smData)) {graphics::lines(x=nGenes, y=smData, col='green')}
  if(!is.null(bonusPoints)) {
    graphics::points(x=bonusPoints$'nGene', y=bonusPoints$'kappa', 
                     col=as.character(bonusPoints$'cols'), pch=as.character(bonusPoints$'pchs'), cex=1.5)  
    graphics::legend('topright', legend=as.character(bonusPoints$'legText'), 
                     pch=as.character(bonusPoints$'pchs'), col=as.character(bonusPoints$'cols'))
  }
}

#' Calculate conditions numbers for signature subsets
#' @description  Remove genes by chunks by picking those the most improve the condition number.  
#' Will set any infinite condition numbers to max(kappas[!is.infinite(kappas)])+1
#' Return the condition numbers with their gene lists
#'
#' @param sigMatrix  The original signature matrix
#' @param numChunks  The number of groups of genes to remove (DEFAULT: NULL)
#' @param verbose  Print out the current chunk as is it's being calculated (DEFAULT: NULL)
#' @param plotIt  The title of the plot (DEFAULT: TRUE)
#' @param singleCore  Set to FALSE to use multiple cores to calculate condition numbers (DEFAULT: FALSE)
#' @param fastStop  Halt early when the condition number changes by less than 1 for 3 iterations (DEFAULT: FALSE)
#'
#' @export
#' @return A list with condition numbers and gene lists
#' @examples
#' library(ADAPTS)
#' LM22 <- ADAPTS::LM22
#' sigGenesList <- shrinkByKappa(sigMatrix=LM22[1:100,1:5], numChunks=4, 
#' verbose=FALSE, plotIt=FALSE, singleCore=TRUE, fastStop=TRUE)
#' 
shrinkByKappa <- function(sigMatrix, numChunks=NULL, verbose=TRUE, plotIt=TRUE, singleCore=FALSE, fastStop=TRUE) {
  curComp <- kappa(as.matrix(sigMatrix))
  names(curComp) <- ""
  
  newSigMatrix <- sigMatrix
  
  #It is much to slow for 4000+ genes in a sig matrix with one gene at a time
  if(is.null(numChunks)){numChunks <- nrow(sigMatrix)}
  stepSize <- max(floor(nrow(sigMatrix) / numChunks), 1)
  
  sigGenesList <- list()
  for(i in 1:numChunks) {
    if(verbose==TRUE){message(paste(i, '/', numChunks))}
    
    #Remove genes one at a time and determine which will minimize the condition number,
    #  It's really to slow to do recalcualte this every iteration.
    if(singleCore==TRUE) {
      newComps <- foreach(fe_curGene = rownames(newSigMatrix), .combine=c) %do% {
        kappa(as.matrix(newSigMatrix[rownames(newSigMatrix)!=fe_curGene,]))
      }
    } else {
      newComps <- foreach(fe_curGene = rownames(newSigMatrix), .combine=c) %dopar% {
        kappa(as.matrix(newSigMatrix[rownames(newSigMatrix)!=fe_curGene,]))
      }
    } #if(singleCore==TRUE) {
    names(newComps) <- rownames(newSigMatrix)
    
    #Alternative: Why not just rank them by variance?  Remove the lowest variance first??
    #  Really, we need soem clustering here to no remove say, B-cell specific genes
    #  That would require some prior information about clusters??
    
    #remGene <- names(newComps)[which.min(newComps)]
    remGene <- names(utils::head(sort(newComps), stepSize))
    curComp <- c(curComp, newComps[remGene])
    
    #Problem: if we remove all non-zero genes for a single sample, that will increase the kappa to Inf.
    #  This loop makes sure that this doesn't happen. 
    for(curRemGene in remGene) {
      tempSigMat <- newSigMatrix[rownames(newSigMatrix) != curRemGene,,drop=FALSE]
      if(nrow(tempSigMat)>0 && !is.infinite(kappa(as.matrix(tempSigMat)))) {newSigMatrix <- tempSigMat}
    }
    #newSigMatrix <- newSigMatrix[!rownames(newSigMatrix) %in% remGene,]
    
    condNum <- kappa(as.matrix(newSigMatrix))
    if( i == 1) {
      deltaKappa <- as.numeric(NA)
    } else {
      deltaKappa <- condNum - sigGenesList[[i-1]]$condNum
    } #if( i == 1) {
    
    if(nrow(newSigMatrix)>0){
      sigGenesList[[i]] <- list(sigGenes=rownames(newSigMatrix), condNum=condNum, deltaKappa=deltaKappa)
    }
    
    if(fastStop == TRUE) {
      if (i > 3) {
        kappas <- c(sigGenesList[[i-2]]$deltaKappa, sigGenesList[[i-1]]$deltaKappa, sigGenesList[[i]]$deltaKappa)
        if(all(abs(kappas) < 1)) { break; }
      } #if (i > 3) {
    } #if(fastStop == TRUE) {
  } #for(i in 1:(nrow(sigMatrix)-1)) {
  
  kappas <- sapply(sigGenesList, function(x){x$condNum})
  if (any(is.infinite(kappas))) {
    kappas[is.infinite(kappas)] <- max(kappas[!is.infinite(kappas)])+1
    for (i in 1:length(sigGenesList)) { sigGenesList[[i]]$condNum <- kappas[i]  }
  }
  
  #plot(y=kappas, x=nGenes, xlab='Number of Genes', ylab='Condition Number')
  if(plotIt==TRUE) {
    kappas <- sapply(sigGenesList, function(x){x$condNum})  #Redo to make sure Inf is gone
    nGenes <- sapply(sigGenesList, function(x){length(x$sigGenes)})
    
    plotKappas(kappas, nGenes=nGenes)
  }
  
  return(sigGenesList)
}

#' Shrink a signature matrix
#' @description  Use shrinkByKappa and automatic minima detection to reduce a signature matrix.
#' Select the new signature matrix with the minima and the maximum number of genes.  There is an
#' inherent difficult in that the condition number will tend to have a second peak at a relatively
#' small number of genes, and then crash so that smallest condition number has more or less one gene.
#' 
#' By default, the algorithm will tend to pick the detected minima with the largest nubmer of genes.
#' aggressiveMin=TRUE will try to find the minimum number of genes that has more genes than the 
#' maxima at the smallest number of genes   
#'
#' @param sigMatrix  The original signature matrix
#' @param numChunks  The number of groups of genes to remove. NULL is all genes (DEFAULT: 100)
#' @param verbose  Print out the current chunk as is it's being calculated (DEFAULT: NULL)
#' @param plotIt  Set to TRUE to plot (DEFAULT: FALSE)
#' @param aggressiveMin  Set to TRUE to aggresively seek the smallest number of genes (DEFAULT: TRUE)
#' @param sigGenesList  Set to use precomputed results from shrinkByKappa (DEFAULT: NULL)
#' @param singleCore  Set to FALSE to use multiple cores to calculate condition numbers (DEFAULT: FALSE)
#' @param fastStop  Halt early when the condition number changes by less than 1 for 3 iterations (DEFAULT: TRUE)
#'
#' @export
#' @return A list with condition numbers and gene lists
#' @examples
#' library(ADAPTS)
#' LM22 <- ADAPTS::LM22
#' newSigMat <- shrinkSigMatrix(sigMatrix=LM22[1:100,1:5], numChunks=4, verbose=FALSE, 
#' plotIt=FALSE, aggressiveMin=TRUE, sigGenesList=NULL, singleCore=TRUE, fastStop=FALSE)
#' 
shrinkSigMatrix <- function(sigMatrix, numChunks=100, verbose=FALSE, plotIt=FALSE, aggressiveMin=TRUE, sigGenesList=NULL, singleCore=FALSE, fastStop=TRUE) {
  if(fastStop == TRUE) {
    message('fastStop==TRUE overwriting aggressiveMin option.')
    aggressiveMin <- TRUE #The logic is getting convoluted.  Please refactor
  }
  else{
    if(is.null(numChunks) || numChunks > nrow(sigMatrix)-1 ){
      numChunks<-nrow(sigMatrix)-1
    }
  }
  
  if(is.null(sigGenesList)) {
    sigGenesList <- shrinkByKappa(sigMatrix=sigMatrix, numChunks=numChunks, verbose=verbose, plotIt=FALSE, singleCore=singleCore, fastStop=fastStop)
  }
  
  kappas <- sapply(sigGenesList, function(x){x$condNum})
  nGenes <- sapply(sigGenesList, function(x){length(x$sigGenes)})
  kappas <- kappas[order(nGenes)]
  nGenes <- nGenes[order(nGenes)]
  
  smData <- stats::smooth(kappas)
  smData2 <- try(stats::predict(stats::smooth.spline(smData))$y)
  if(inherits(smData2, 'try-error')) {
    smData2 <- stats::smooth(smData)
  }
  mins <- quantmod::findValleys(smData2)
  maxs <- quantmod::findPeaks(smData2)
  if(length(mins)==0) {mins <- which.min(kappas)}
  if(length(maxs)==0) {maxs <- which.max(kappas)}
  
  if (fastStop== TRUE) {
    legText <- 'min'
    pchs <- 2
    cols <- 'orange'
    bonusGenes <- min(nGenes) 
    bonusKappas <- kappas[which(nGenes==bonusGenes)[1]]
    bonusPoints <- data.frame(legText = legText, pchs=pchs, cols=cols, kappa=bonusKappas, nGenes=bonusGenes, stringsAsFactors = FALSE)
  } else {
    legText <- c(make.names(rep('min', length(mins)), unique=TRUE), make.names(rep('max', length(maxs)), unique=TRUE))
    pchs <- c(rep(2, length(mins)), rep(6, length(maxs)))
    cols <- c(rep('red', length(mins)), rep('blue', length(maxs)))
    bonusKappas <- kappas[c(mins, maxs)]
    bonusGenes <- nGenes[c(mins, maxs)]
    bonusPoints <- data.frame(legText = legText, pchs=pchs, cols=cols, kappa=bonusKappas, nGenes=bonusGenes, stringsAsFactors = FALSE)
    
    #There's two possible algorithms, 
    #  A) Find the min with the minimum condition number
    #  B) Find the min with maximum number of genes
    # First: use B
  } # if (fastStop== TRUE) {
  minIdx <-  sub('\\.+[0-9]$', '', bonusPoints$legText) == 'min' 
  bonus.min <- bonusPoints[minIdx, ]
  
  if (aggressiveMin == TRUE) {
    chosenPoint <- bonus.min$legText[which.max(bonus.min$nGenes)]
    bonusPoints[bonusPoints$legText==chosenPoint,'cols'] <- 'orange'
    bonusPoints[bonusPoints$legText==chosenPoint,'pchs'] <- 17
  } else {
    maxIdx <-  sub('\\.+[0-9]$', '', bonusPoints$legText) == 'max' 
    bonus.max <- bonusPoints[maxIdx, ]
    chosenMax <- bonus.max$legText[which.min(bonus.max$nGenes)]
    nGene.bool <- nGenes > bonus.max$nGenes[bonus.max$legText==chosenMax]
    if(all(nGene.bool == FALSE)) { nGene.bool <- nGenes >= bonus.max$nGenes[bonus.max$legText==chosenMax] }
    relMinKappa <- min(kappas[nGene.bool])
    chosenPointDF <- data.frame(legText = 'Chosen', pchs=17, cols='orange', 
                                kappa=relMinKappa, nGenes=nGenes[kappas==relMinKappa], 
                                stringsAsFactors = FALSE)
    bonusPoints <- rbind(bonusPoints, chosenPointDF)
    chosenPoint <- 'Chosen'
  }
  
  if (plotIt==TRUE) {
    plotKappas(kappas=kappas, nGenes=nGenes, smData = smData2, bonusPoints=bonusPoints)
  }
  
  numChosenGenes <- bonusPoints[bonusPoints$legText==chosenPoint,'nGenes']
  #Check the bestSig calculation
  bestSig <- sigGenesList[sapply(sigGenesList, function(x){length(x$sigGenes)}) == numChosenGenes][[1]]
  
  smallMatrix <- sigMatrix[bestSig$sigGenes,]
  return(smallMatrix)
}

#' Split a single cell dataset into multiple sets
#' @description  Take a matrix of single cell data with genes as rows and each column corresponding
#' to a single cells. Break it up into rougly equal subsets, taking care to make sure that each cell type is represented
#' in each set if possible
#'
#' @param RNAcounts  The single cell matrix
#' @param cellIDs  A vector will cell types for each column in scCountMatrix (DEFAULT: colnames(RNAcounts))
#' @param numSets  The number of sets to break it up into (DEFAULT: 3)
#' @param verbose  Set to TRUE to print cell counts as it goes (DEFAULT: TRUE)
#' @param randomize  Set to TRUE to randomize the sets (DEFAULT: TRUE)
#'
#' @export
#' @return a list with a multiple sets
#' @examples
#' RNAcounts <- matrix(0, nrow=10, ncol=30)
#' rownames(RNAcounts) <- make.names(rep('Gene', nrow(RNAcounts)), unique=TRUE)
#' colnames(RNAcounts) <- make.names(c('CellX', rep('CellY', 9), 
#' rep('CellZ', 10), rep('CellB', 10)), unique=TRUE)
#' RNAcounts[, grepl('CellY', colnames(RNAcounts))] <- 1
#' RNAcounts[, grepl('CellZ', colnames(RNAcounts))] <- 2
#' RNAcounts[, grepl('CellB', colnames(RNAcounts))] <- 3
#' splitSCdata(RNAcounts, numSets=3)
#' 
splitSCdata <- function(RNAcounts, cellIDs=colnames(RNAcounts), numSets=3, verbose=TRUE, randomize=TRUE) {
  cellTypes <- unique(sub('\\.[0-9]+$', '', cellIDs))
  seqList <- list()
  for (cellType in cellTypes) {
    matchStr <- sub("+",'\\+',paste0(cellType,'\\.*[0-9]*$'), fixed=TRUE)
    idxs <- grep(matchStr, cellIDs)
    
    maxSam <- ceiling(length(idxs)/numSets)
    
    #https://stackoverflow.com/questions/3318333/split-a-vector-into-chunks-in-r
    #> max <- 20
    #> x <- seq_along(d)
    #> d1 <- split(d, ceiling(x/max))
    x <- seq_along(idxs)
    if(randomize==TRUE) {x <- sample(x)}#Add option to not randomize?
    idxList <- split(idxs, ceiling(x/maxSam))
    message(c(paste(cellType, ':', length(idxs)), '; ', paste(sapply(idxList, length), collapse=', ')))
    for (y in names(idxList)) {seqList[[y]] <- c(seqList[[y]], idxList[[y]])}
  } #for (cellType in cellTypes) {
  
  #Reshape into lists for each group
  setList <- lapply(seqList, function(x){ RNAcounts[,x] })
  return(setList)
} #splitSCdata <- function(RNAcounts, cellIDs=colnames(RNAcounts), numSets=3) {

#' Build groupSize pools according to cellIDs
#' @description  This function is intended to collapse many single cells into 3 (groupsize) groups
#' with the average count across all cells in each of the groups.  These groups can then be used to perform a 
#' t-test (for example) between the 3 groups of CellX with 3 groups of CellY
#' 
#' @param RNAcounts  The single cell matrix
#' @param cellIDs  A vector will cell types for each column in scCountMatrix (DEFAULT: colnames(RNAcounts))
#' @param groupSize  The number of sets to break it up into (DEFAULT: 3)
#' @param randomize  Set to TRUE to randomize the sets (DEFAULT: TRUE)
#' @param mc.cores  The number of cores to use (DEFAULT: 1)
#'
#' @export
#' @return a list with a multiple sets
#' @examples
#' RNAcounts <- matrix(0, nrow=10, ncol=100)
#' rownames(RNAcounts) <- make.names(rep('Gene', nrow(RNAcounts)), unique=TRUE)
#' colnames(RNAcounts) <- make.names(c('CellX', rep('CellY', 39), 
#' rep('CellZ', 30), rep('CellB', 30)), unique=TRUE)
#' RNAcounts[, grepl('CellY', colnames(RNAcounts))] <- 1
#' RNAcounts[, grepl('CellZ', colnames(RNAcounts))] <- 2
#' RNAcounts[, grepl('CellB', colnames(RNAcounts))] <- 3
#' scSample(RNAcounts, groupSize=3)
#' 
scSample <- function(RNAcounts, cellIDs=colnames(RNAcounts), groupSize=3, randomize=TRUE, mc.cores=1) {
  cellTypes <- unique(sub('\\.[0-9]+$', '', as.character(cellIDs)))
  #combCellList <- list()
  #for (cellType in cellTypes) {
  combCellList <- parallel::mclapply(cellTypes, mc.cores=mc.cores, function(cellType) {
    matchStr <- sub("+",'\\+',paste0(cellType,'\\.*[0-9]*$'), fixed=TRUE)
    idxs <- grep(matchStr, cellIDs)
    if(randomize==TRUE) {idxs <- sample(x=idxs, length(idxs))}
    if(length(idxs) > groupSize) {
      maxSam <- ceiling(length(idxs)/groupSize)
      
      #https://stackoverflow.com/questions/3318333/split-a-vector-into-chunks-in-r
      #> max <- 20
      #> x <- seq_along(d)
      #> d1 <- split(d, ceiling(x/max))
      x <- seq_along(idxs)
      idxList <- split(idxs, ceiling(x/maxSam))
      
      curData.list <- lapply(idxList, function(x) {
        rowMeans(RNAcounts[,x,drop=FALSE], na.rm=TRUE)
      })
      curData <- do.call(cbind, curData.list)
      colnames(curData) <- rep(cellType, ncol(curData))
      #combCellList[[cellType]] <- curData
    } else {
      curData <- RNAcounts[,idxs,drop=FALSE]
      if(length(idxs) == 1) {
        colnames(curData) <- cellType
      } else {
        #colnames(curData) <- 1:length(idxs)#paste(cellType, 1:length(idxs), sep='.')
        colnames(curData) <-  rep(cellType, ncol(curData))
      }
      #combCellList[[cellType]] <- curData
    }
    curData
  } ) #for (cellType in cellTypes) {
  
  #Now convert it to a data.frame
  combDF <- do.call(cbind, combCellList)
  colnames(combDF) <- sub('\\..*$', '', colnames(combDF)) #Remove <cell>.<cell> name.  Need smarter RegEX
  return(combDF)
}

#' Calculate prediction accuracy
#' @description  Calculate correlation coeffifients, p-Values, MAE, RMSE for deconvolution predictions
#'
#' @param estimates  The estimated cell percentages
#' @param reference  The reference cell percentages
#'
#' @export
#' @return a list with a multiple sets
#' @examples
#' estimates <- sample(c(runif(8), 0 ,0))
#' estimates <- 100 * estimates / sum(estimates)
#' reference <- sample(c(runif(7), 0 , 0, 0))
#' reference <- 100 * reference / sum(reference)
#' calcAcc(estimates, reference)
#' 
calcAcc <- function(estimates, reference) {
  ct <- stats::cor.test(estimates, reference)
  ct.spear <- suppressWarnings(stats::cor.test(estimates, reference, method = 'spearman'))
  mae <- mean(abs(estimates - reference), na.rm=TRUE)
  rmse <- sqrt(mean((estimates - reference)^2, na.rm=TRUE))
  
  #Binarize to get at sensitivity and specificity
  F1mcc <- getF1mcc(estimate=estimates>0, reference=reference>0)
  Zs <- estimates == 0
  minForNon0 <- min(reference[!Zs & !is.na(reference)])
  
  out <- c(rho=ct$estimate, pVal=ct$p.value, spear=ct.spear$estimate, pVal.spear=ct.spear$p.value, 
           mae=mae, rmse=rmse, sens=F1mcc[['sensitivity']], minForNon0=minForNon0)
  return(out)
}

#' Get f1 / mcc 
#' @description Get f1 / mcc and other accuracy measurements for binary predictions.  
#' Provide either an estimate and reference vector
#' e.g. getF1mcc(estimate, reference)
#' Or TPs, FPs, etc.
#' e.g. getF1mcc(tps=3, fps=4, tns=7, fns=2)
#'
#' @param estimate  A binary vector of predictions
#' @param reference  a binary vector of actual values
#' @param tps  The number of TPs
#' @param fps  The number of FPs
#' @param tns  The number of TNs
#' @param fns  The number of FNs
#'
#' @export
#' @return A vector with sensitivity, specificity, fpr, fdr, f1, agreement, p.value, mcc, and mcc.p
#' @examples 
#' estimates <- sample(c(runif(8), 0 ,0))
#' reference <- sample(c(runif(7), 0 , 0, 0))
#' accuracyStats <- getF1mcc(estimate=estimates>0, reference=reference>0)
getF1mcc <- function(estimate=NULL, reference=NULL, tps=NULL, fps=NULL, tns=NULL, fns=NULL) {
  if(is.null(tps) | is.null(fps) | is.null(tns) | is.null(fns)) {
    if(is.null(estimate) | is.null(reference)) {
      stop('Must specific either estimates and reference or tps, fps, tns, fns.')
    } else {
      naIdx <- is.na(estimate) | is.na(reference) | is.null(estimate) | is.null(reference)
      estimate <- estimate[!naIdx]
      reference <- reference[!naIdx]
      tps <- sum(estimate==TRUE & reference==TRUE)
      fps <- sum(estimate==TRUE & reference==FALSE)
      tns <- sum(estimate==FALSE & reference==FALSE)
      fns <- sum(estimate==FALSE & reference==TRUE)
    }
  } #if(is.null(tps) | is.null(fps) | is.null(tns) | is.null(fns)) {
  
  sensitivity <- tps/(tps+fns)
  specificity <- tns/(fps+tns)
  fpr <- fps / (fps+tns)
  fdr <- fps / (fps+tps)
  f1 <- (2*tps) / (2*tps + fps + fns)
  N <- tns+tps+fns+fps
  S <- (tps+fns) / N
  P <- (tps+fps) / N
  
  mcc <- (tps/N - S * P) / sqrt(P*S*(1-S)*(1-P))
  mccData <- rbind( t(matrix(c(1,1), 2, max(tps,1))), t(matrix(c(0,0), 2, max(tns,1))),
                    t(matrix(c(1,0), 2, max(fps,1))), t(matrix(c(0,1), 2, max(fns,1))) )
  #Note: max assumes that you will get a correlation
  mcc.p <- stats::cor.test(mccData[,1], mccData[,2])$p.value
  if (is.na(mcc) | mcc==Inf) { mcc <- 0; mcc.p <- 1 }
  
  agreement <- (tps+tns)/(tps+fps+tns+fns)
  bgProb <- max(c( (tps+fns)/(tps+fps+tns+fns) , 1-(tps+fns)/(tps+fps+tns+fns) ))
  p.value <- stats::pbinom(q=tps+tns, size=tps+fps+tns+fns, prob=bgProb, lower.tail = FALSE)
  
  return(c(sensitivity=sensitivity, specificity=specificity, fpr=fpr, fdr=fdr, f1=f1, agreement=agreement, p.value=p.value, mcc=mcc, mcc.p=mcc.p))
}



#' Build a gList using random forest
#' @description Use ranger to select features and build a genesInSeed gene matrix
#'
#' @param trainSet  Each row is a gene, and each column is an example of a particular cell type, e.g. ADAPTS::scSample(trainSet, groupSize=30)
#' @param oneCore  SEt to TRUE to disable multicore (DEFAULT: FALSE)
#'
#' @export
#' @return A cell specific geneList for ADAPTS::AugmentSigMatrix()
#' @examples
#' library(ADAPTS)
#' ct1 <- runif(1000, 0, 100)
#' ct2 <- runif(1000, 0, 100)
#' dataMat <- cbind(ct1, ct1, ct1, ct1, ct1, ct1, ct2, ct2, ct2, ct2)
#' rownames(dataMat) <- make.names(rep('gene', nrow(dataMat)), unique=TRUE)
#' noise <- matrix(runif(nrow(dataMat)*ncol(dataMat), -2, 2), nrow = nrow(dataMat), byrow = TRUE)
#' dataMat <- dataMat + noise
#' gList <- gListFromRF(trainSet=dataMat, oneCore=TRUE)
#' 
gListFromRF <- function(trainSet, oneCore=FALSE) {
  clusterNames <- unique(sub('\\.[0-9]+$', '', colnames(trainSet)))
  trainSet.4reg <- t(trainSet)
  
  if(oneCore==TRUE) {
    gList.fromRF <- lapply (clusterNames, function(cn) {
      clusterBool <- colnames(trainSet) == cn
      
      rf1 <- ranger::ranger(x=trainSet.4reg, y=clusterBool, num.trees=1000, importance='impurity')
      imp <- ranger::importance(rf1)
      imp <- sort(imp[imp>0], decreasing = FALSE)
      curDF <- data.frame(rat=imp, t=0, pVal=0, qVal=0)
    })
  } else {
    gList.fromRF <- parallel::mclapply (clusterNames, function(cn) {
      clusterBool <- colnames(trainSet) == cn
      
      rf1 <- ranger::ranger(x=trainSet.4reg, y=clusterBool, num.trees=1000, importance='impurity')
      imp <- ranger::importance(rf1)
      imp <- sort(imp[imp>0], decreasing = FALSE)
      curDF <- data.frame(rat=imp, t=0, pVal=0, qVal=0)
    })
  } #if(oneCore==TRUE) {
  names(gList.fromRF) <- clusterNames
  return(gList.fromRF)
}

#' Make a GSVA genelist
#' @description Provide a gList and signature matrix with matched cell types to get signatures
#'   gene lists for GSVA and similar algorithms.
#'   gList=NULL select highest genes for each cell type, minimum of 3.
#'
#' @param sigMat  A signature matrix such as from ADAPTS::AugmentSigMatrix()
#' @param gList  A list of prioritized genes such as from ADAPTS::gListFromRF() (DEFAULT:NULL)
#'
#' @export
#' @return A list of genes for each cell types musually in sigMat and gList
#' @examples
#' library(ADAPTS)
#' ct1 <- runif(1000, 0, 100)
#' ct2 <- runif(1000, 0, 100)
#' dataMat <- cbind(ct1, ct1, ct1, ct1, ct1, ct1, ct2, ct2, ct2, ct2)
#' rownames(dataMat) <- make.names(rep('gene', nrow(dataMat)), unique=TRUE)
#' noise <- matrix(runif(nrow(dataMat)*ncol(dataMat), -2, 2), nrow = nrow(dataMat), byrow = TRUE)
#' dataMat <- dataMat + noise
#' gList <- ADAPTS::gListFromRF(trainSet=dataMat, oneCore=TRUE)
#' newSigMat <- ADAPTS::buildSeed(trainSet=dataMat, plotIt=FALSE)
#' geneLists <- matrixToGenelist(sigMat=newSigMat, gList=gList)
#' 
matrixToGenelist <- function(sigMat, gList=NULL) {
  geneLists <- list() #The output variable
  
  if(is.null(gList)) {
    bestCols <- apply(sigMat, 1, which.max)
    for (curName in colnames(sigMat)) {
      seedGenes <- names(tail(sort(sigMat[, curName]),3))
      bestGenes <- names(bestCols)[bestCols == which(colnames(sigMat)==curName)]
      allGenes <- unique(c(seedGenes, bestGenes))
      rats <- sapply(allGenes, function(x){ sigMat[x,curName] / sum(sigMat[x,colnames(sigMat) != curName])})
      geneLists[[curName]] <- names(sort(rats, decreasing = TRUE))
    }
  } else {
    olNames <- names(gList)[names(gList) %in% colnames(sigMat)]
    if(length(olNames) < length(gList) | length(olNames) < ncol(sigMat)) {
      message('Not all cell-type names match between gList and sigMat')
      message('  Only matched cell-types will be calculated')
    }
    for (curName in olNames) {
      olGenes <- rownames(gList[[curName]])[rownames(gList[[curName]]) %in% rownames(sigMat)]
      
      #Make sure that genes are high in current cell type.
      rats <- sapply(olGenes, function(x){ sigMat[x,curName] / sum(sigMat[x,colnames(sigMat) != curName])})
      rats <- rats[rats>1]
      geneLists[[curName]] <- names(sort(rats, decreasing = TRUE))
    }
  }
  
  return(geneLists)
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
