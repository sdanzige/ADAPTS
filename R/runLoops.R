
#' Build a deconvolution seed matrix, add the proportional option 
#' @description Use ranger to select features and build a genesInSeed gene matrix
#'
#' @param trainSet  Each row is a gene, and each column is an example of a particular cell type, ie from single cell data 
#' @param genesInSeed  The maximum number of genes in the returned seed matrix (DEFAULT: 200)
#' @param groupSize  The number of groups to break the trainSet into by ADAPTS::scSample (DEFAULT: 30)
#' @param randomize  Set to TRUE randomize the sets selected by ADAPTS::scSample (DEFAULT: TRUE)
#' @param num.trees  The number of trees to be used by ranger (DEFAULT: 1000)
#' @param plotIt  Set to TRUE to plot (DEFAULT: TRUE)
#' @param trainSet.3sam  Optional pre-calculated ADAPTS::scSample(trainSet, groupSize = 3) (DEFAULT: NULL)
#' @param trainSet.30sam  Optional pre-calculated ADAPTS::scSample(trainSet, groupSize=groupSize, randomize=randomize) (DEFAULT: NULL)
#' @param proportional  Set to true to make the training set cell type proportional.  Ignores group size (DEFAULT: FALSE)
#'
#' @export
#' @return A list with condition numbers and gene lists
#' @examples
#' library(ADAPTS)
#' ct1 <- runif(1000, 0, 100)
#' ct2 <- runif(1000, 0, 100)
#' dataMat <- cbind(ct1, ct1, ct1, ct1, ct1, ct1, ct2, ct2, ct2, ct2)
#' rownames(dataMat) <- make.names(rep('gene', nrow(dataMat)), unique=TRUE)
#' noise <- matrix(runif(nrow(dataMat)*ncol(dataMat), -2, 2), nrow = nrow(dataMat), byrow = TRUE)
#' dataMat <- dataMat + noise
#' newSigMat <- buildSeed(trainSet=dataMat)
#' 
buildSeed <- function(trainSet, genesInSeed=200, groupSize=30, randomize=TRUE, num.trees=1000, plotIt=TRUE, trainSet.3sam=NULL, trainSet.30sam=NULL, proportional=FALSE) {
  if(is.null(trainSet.3sam)) {trainSet.3sam <- ADAPTS::scSample(RNAcounts = trainSet, groupSize = 3, randomize = randomize)}
  if (proportional==TRUE) {
    #colnames(trainSet) <- sub('\\.[0-9]+$', '', colnames(trainSet))
    tsNames <- sub('\\.[0-9]+$', '', colnames(trainSet))
    cellProps <- table(tsNames)
    cellSampleCounts <- 3*ceiling(100*cellProps / sum(cellProps))
    trainList <- list()
    for (curCount in unique(cellSampleCounts)) {
      curClusts <- names(cellSampleCounts)[cellSampleCounts==curCount]
      trainList[[as.character(curCount)]] <- ADAPTS::scSample(RNAcounts = trainSet[, tsNames %in% curClusts], groupSize = curCount, randomize = randomize)
    }
    trainSet.30sam <- do.call(cbind, trainList)
  } else {
    if(is.null(trainSet.30sam)) {trainSet.30sam <- ADAPTS::scSample(RNAcounts = trainSet, groupSize = groupSize, randomize = randomize)}
  }
  
  clusterIDs <- factor(colnames(trainSet.30sam))
  trainSet.4reg <- t(trainSet.30sam)
  
  rf1 <- ranger::ranger(x=trainSet.4reg, y=clusterIDs, num.trees=num.trees, importance='impurity')
  imp <- ranger::importance(rf1)
  imp <- sort(imp[imp>0], decreasing = TRUE)
  
  
  topGenes <- names(imp)[1:min(genesInSeed, length(imp))]
  #topGenes[!topGenes %in% rownames(trainSet.3sam)]
  
  seedMat <- trainSet.3sam[rownames(trainSet.3sam) %in% topGenes,]
  cellTypes <- sub('\\.[0-9]+$', '', colnames(seedMat))
  seedMat <- t(apply(seedMat, 1, function(x){tapply(x, cellTypes, mean, na.rm=TRUE)}))
  if(plotIt==TRUE) {pheatmap:: pheatmap(seedMat, main=paste('Seed Matrix','\n# Cell Types:',ncol(seedMat),'| # Genes:',nrow(seedMat))) }
  return(seedMat)
}
  
  

   
   



#' Generate all the signature matrices one time with the option to leave out half of the data as a test set
#' @description  This wrapper is helpful for repetitively matrix generation. It generates seed matrix, all-gene matrix, augmented matrix, shrunk matrix,
#' and all the clustered matrices in one call.
#'
#' @param exprData The gene express data. Each row is a gene, and each column is an example of a particular cell type.
#' @param randomize Set to to TRUE randomize the sets selected by ADAPTS::scSample (DEFAULT: TRUE)
#' @param skipShrink Set to TRUE to skip shrinking the signatrure matrix (DEFAULT: TRUE)
#' @param proportional Set to true to make the training set cell type proportional.  Ignores group size (DEFAULT: FALSE)
#' @param handMetaCluster A List of pre-defined meta clusters.Set to NULL to automatically group indistinguishable cells 
#' into same cluster using clustWspillOver.(DEFAULT: NULL)
#' @param testOnHalf Set to TRUE to leave half the data as a test set
#' @param condTol  The tolerance in the reconstruction algorithm.  1.0 = no tolerance, 1.05 = 5\% tolerance (DEFAULT: 1.01)
#' @param numChunks  The number of groups of genes to remove while shrinking (DEFAULT: NULL, i.e. 1)
#' @param plotIt  Set to FALSE to suppress plots (DEFAULT: TRUE)
#' @param fastStop  Halt early when the condition number changes by less than 1 for 3 iterations (DEFAULT: TRUE)
#' @param singleCore  TRUE for a single core (DEFAULT: TRUE)
#'
#' @export
#' @return A list of results including prediction accuracy and cell enrichment
#'
#' @examples
#' ct1 <- runif(1000, 0, 100)
#' ct2 <- runif(1000, 0, 100)
#' ct3 <- runif(1000, 0, 100)
#' ct4 <- runif(1000, 0, 100)
#' dataMat <- cbind(ct1, ct1, ct1, ct1, ct1, ct1, ct2, ct2, ct2, ct2, ct3, ct3, ct3,ct3,ct4,ct4)
#' rownames(dataMat) <- make.names(rep('gene', nrow(dataMat)), unique=TRUE)
#' noise <- matrix(runif(nrow(dataMat)*ncol(dataMat), -2, 2), nrow = nrow(dataMat), byrow = TRUE)
#' dataMat <- dataMat + noise
#' metaList <- list()
#' colnames(dataMat) <- sub('\\..*','', colnames(dataMat))
#' metaList[[1]] <- c(unique(colnames(dataMat))[1])  #Cell Type 1
#' metaList[[2]] <- c(unique(colnames(dataMat))[2])  #Cell Type 2
#' metaList[[3]] <- c(unique(colnames(dataMat))[3])  #Cell Type 3
#' metaList[[4]] <- c(unique(colnames(dataMat))[4:length(unique(colnames(dataMat)))])  #Cell Type 4
#' #options(mc.cores=2)
#' #  This is a meta-function that calls other functions, 
#' #  The execution speed is too slow for the CRAN automated check
#' #testAllSigMatrices(exprData=dataMat, randomize = TRUE, skipShrink=FALSE, 
#' #    proportional=FALSE, handMetaCluster=metaList, testOnHalf=TRUE, numChunks=NULL)
    
testAllSigMatrices <- function(exprData, randomize = TRUE, skipShrink=FALSE, proportional=FALSE, handMetaCluster=NULL, testOnHalf=TRUE, condTol=1.01, numChunks=100, plotIt=TRUE, fastStop=TRUE, singleCore=TRUE) {
  
  if(randomize==TRUE) {set.seed(Sys.time())}
  resList <- list()
  
  if(testOnHalf == TRUE){
    trainTestSet <- ADAPTS::splitSCdata(exprData, numSets=2, randomize = randomize)
    trainSet <- trainTestSet[[1]]
    testSet <-trainTestSet[[2]]
  }
  else {
  trainSet <- exprData
  testSet <- exprData
}
  
  
  trainSet.30sam <- ADAPTS::scSample(RNAcounts = trainSet, groupSize = 30, randomize = randomize)
  trainSet.3sam <- ADAPTS::scSample(RNAcounts = trainSet, groupSize = 3, randomize = randomize)
  
  pseudobulk.test <- data.frame(test=rowSums(testSet))
  pseudobulk.test.counts<-table(sub('\\..*','',colnames(testSet)))
  actFrac.test <- 100 * pseudobulk.test.counts / sum(pseudobulk.test.counts)
  
  
  #seed
  genesInSeed<-100
  seedMat <- buildSeed(trainSet, genesInSeed=genesInSeed, groupSize=30, randomize=TRUE, num.trees=1000, plotIt=plotIt, trainSet.3sam=trainSet.3sam, trainSet.30sam=trainSet.30sam,proportional = proportional)
  resList[['matrix.seed']] <- seedMat
  
  estimates.onTest <- as.data.frame(ADAPTS::estCellPercent.DCQ(seedMat, pseudobulk.test))
  
  colnames(estimates.onTest) <- paste('Seed Matrix')
  estimates.onTest$actFrac.test <- round(actFrac.test[rownames(estimates.onTest)],2)
  resList[['estimates.onTest']] <- estimates.onTest
  
  resList[['testAcc.seed']] <- seed2TestAcc <- ADAPTS::calcAcc(estimates=estimates.onTest[,1], reference=estimates.onTest[,2])
  
  
  #All gene
  allGeneSig <- apply(trainSet.3sam, 1, function(x){tapply(x, colnames(trainSet.3sam), mean, na.rm=TRUE)})
  
  estimates.allGene <- as.data.frame(ADAPTS::estCellPercent.DCQ(t(allGeneSig), pseudobulk.test))
  colnames(estimates.allGene)<-'All Gene Matrix'
  
  estimates.onTest<-cbind(estimates.allGene,estimates.onTest)
  resList[['estimates.onTest']] <- estimates.onTest
  
  resList[['testAcc.all']] <- seed2TestAcc <-  ADAPTS::calcAcc(estimates=estimates.onTest[,1],reference=estimates.onTest[,ncol(estimates.onTest)])
  
  
  # Aug
  gList <- ADAPTS::gListFromRF(trainSet=trainSet.30sam)
  resList[['gList']] <- gList
  
  augTrain <- ADAPTS::AugmentSigMatrix(origMatrix = seedMat, fullData = trainSet.3sam, gList = gList, nGenes = 1:100, newData = trainSet.3sam, plotToPDF = FALSE, pdfDir = '.', condTol=condTol, plotIt=plotIt)
  
  resList[['matrix.aug']] <- augTrain
  resList[['matrix.aug.bestGenes']] <- ADAPTS::matrixToGenelist(augTrain, gList)
  
  
  estimates.augment <- as.data.frame(ADAPTS::estCellPercent.DCQ(augTrain, pseudobulk.test))
  colnames(estimates.augment) <- paste('Augmented Matrix')
  estimates.onTest <- cbind(estimates.augment, estimates.onTest)
  resList[['estimates.onTest']] <- estimates.onTest
  
  resList[['testAcc.aug']] <- seed2TestAcc <- ADAPTS::calcAcc(estimates=estimates.onTest[,1], reference=estimates.onTest[,ncol(estimates.onTest)])
  
  #shrink
  if(skipShrink == FALSE) {
    #augTrain.shrink <- ADAPTS::shrinkSigMatrix(sigMatrix=augTrain, numChunks=NULL, verbose=FALSE, plotIt = FALSE, aggressiveMin=TRUE,sigGenesList=NULL, fastStop=TRUE, singleCore=TRUE)
    augTrain.shrink <- ADAPTS::shrinkSigMatrix(sigMatrix=augTrain, verbose=FALSE, plotIt = plotIt, aggressiveMin=TRUE, fastStop=fastStop, singleCore=singleCore, numChunks=numChunks)
    #pheatmap(augTrain.shrink)
    resList[['matrix.shrink']] <- augTrain.shrink
    resList[['matrix.shrink.bestGenes']] <- ADAPTS::matrixToGenelist(augTrain.shrink, gList)
    
    
    estimates.shrink <- as.data.frame(ADAPTS::estCellPercent.DCQ(augTrain.shrink, pseudobulk.test))
    colnames(estimates.shrink) <- paste('Shrunk Matrix')
    resList[['estimates.onTest']] <- estimates.onTest <- cbind(estimates.shrink, estimates.onTest)
    
    resList[['testAcc.shrink']] <- seed2TestAcc<- ADAPTS::calcAcc(estimates=estimates.onTest[,1], reference=estimates.onTest[,ncol(estimates.onTest)])
  } else {
    augTrain.shrink <- augTrain  #Used for later clustering
  }
  
  
  #Clustering
  
  if(!is.null(handMetaCluster)) {
    resList[['allClusters']]  <- handMetaCluster
  } else {
    varClusts <- ADAPTS::clustWspillOver(sigMatrix = augTrain.shrink, geneExpr = trainSet.3sam)
    resList[['allClusters']]  <- varClusts$allClusters
  }
  
  resList[['allClusters']] 
  
  metaCluster.id <- list()
  for(i in 1:length(resList[['allClusters']])) {
    for (x in resList[['allClusters']][[i]]) {
      metaCluster.id[[x]] <- paste('Meta',i,sep='_')
    }
  }
  metaClust.LUT <- resList[['metaClust.LUT']]  <- unlist(metaCluster.id)
  names(resList[['allClusters']]) <- metaClust.LUT[sapply(resList[['allClusters']], function(x){x[1]})]
  
  #Update 05-19-20: More informative names
  metaNames <- sapply(unique(metaClust.LUT), function(x) { paste(names(metaClust.LUT)[metaClust.LUT == x], collapse='_')})
  metaClust.LUT.MI <- metaClust.LUT
  metaClust.LUT.MI <- metaNames[match(metaClust.LUT.MI, names(metaNames))]
  names(metaClust.LUT.MI) <- names(metaClust.LUT)
  
  metatrainSet<-trainSet
  metatestSet<-testSet
  
  colnames(metatrainSet) <- metaClust.LUT.MI[sub('\\..*','',colnames(metatrainSet))]
  colnames(metatestSet) <- metaClust.LUT.MI[sub('\\..*','',colnames(metatestSet))]
  
  metatrainSet.3sam <- ADAPTS::scSample(RNAcounts = metatrainSet, groupSize = 3, randomize = TRUE)
  metatrainSet.30sam <- ADAPTS::scSample(RNAcounts = metatrainSet, groupSize = 30, randomize = TRUE)
  
  metaclusterIDs <- factor(colnames(metatrainSet.30sam))
  metatrainSet.4reg <- t(metatrainSet.30sam)
  
  metapseudobulk.test <- data.frame(test=rowSums(metatestSet))
  metapseudobulk.test.counts <- table(sub('\\..*','',colnames(metatestSet)))
  meta.actFrac <- 100 * metapseudobulk.test.counts / sum(metapseudobulk.test.counts)
  
  #Metaseed
  genesInSeed<-100
  metaseedMat <-buildSeed(metatrainSet, genesInSeed=genesInSeed, groupSize=30, randomize=TRUE, num.trees=1000, plotIt=plotIt, trainSet.3sam=metatrainSet.3sam, trainSet.30sam=metatrainSet.30sam,proportional = proportional)
  
  resList[['matrix.metaSeed']] <- metaseedMat
  
  estimates.Meta.onTest <- as.data.frame(ADAPTS::estCellPercent.DCQ(metaseedMat, metapseudobulk.test))
  
  colnames(estimates.Meta.onTest) <- paste('Seed Meta')
  estimates.Meta.onTest$actFrac.test <- round(meta.actFrac[rownames(estimates.Meta.onTest)],2)
  
  resList[['estimates.onTest.meta']] <- estimates.Meta.onTest
  
  resList[['testAcc.metaSeed']] <- seed2TestAcc.meta <- ADAPTS::calcAcc(estimates=estimates.Meta.onTest[,1], reference=estimates.Meta.onTest[,2])
  
  
  #meta all gene
  metaallGeneSig <- apply(metatrainSet.3sam, 1, function(x){tapply(x, colnames(metatrainSet.3sam), mean, na.rm=TRUE)})
  metaestimates.allGene <- as.data.frame(ADAPTS::estCellPercent.DCQ(t(metaallGeneSig), metapseudobulk.test))
  colnames(metaestimates.allGene)<-paste('All Gene Meta')
  estimates.Meta.onTest <- cbind(metaestimates.allGene, estimates.Meta.onTest)
  resList[['estimates.onTest.meta']] <- estimates.Meta.onTest
  
  resList[['testAcc.metaAll']] <- seed2TestAcc <- ADAPTS::calcAcc(estimates=estimates.Meta.onTest[,1], reference=estimates.Meta.onTest[,ncol(estimates.Meta.onTest)])
  
  #meta aug
  metagList <- ADAPTS::gListFromRF(trainSet=metatrainSet.30sam)
  resList[['gList.meta']] <- metagList
  sapply(gList,dim)
  
  meta.augTrain <- ADAPTS::AugmentSigMatrix(origMatrix = metaseedMat, fullData = metatrainSet.3sam, gList = metagList, nGenes = 1:100, newData = metatrainSet.3sam, plotToPDF = FALSE, pdfDir = '.', condTol = condTol, plotIt=plotIt)
  
  resList[['matrix.metaAug']] <- meta.augTrain
  resList[['matrix.metaAug.bestGenes']] <- ADAPTS::matrixToGenelist(meta.augTrain, metagList)
  
  
  estimates.Meta.augment <- as.data.frame(ADAPTS::estCellPercent.DCQ(meta.augTrain, metapseudobulk.test))
  colnames(estimates.Meta.augment) <- paste('Augmented Meta')
  resList[['estimates.onTest.meta']] <- estimates.Meta.onTest <- cbind(estimates.Meta.augment, estimates.Meta.onTest)
  
  resList[['testAcc.metaAug']] <- seed2TestAcc.aug.meta <- ADAPTS::calcAcc(estimates=estimates.Meta.onTest[,1], reference=estimates.Meta.onTest[,ncol(estimates.Meta.onTest)])
  
  #meta shrink
  if(skipShrink == FALSE) {
    gc()
    #meta.augTrain.shrink <- ADAPTS::shrinkSigMatrix(sigMatrix=meta.augTrain, numChunks=NULL, verbose=FALSE, plotIt = FALSE, aggressiveMin=TRUE, sigGenesList=NULL,fastStop=TRUE, singleCore=TRUE)
    meta.augTrain.shrink <- ADAPTS::shrinkSigMatrix(sigMatrix=meta.augTrain, verbose=FALSE, plotIt = plotIt, aggressiveMin=TRUE, fastStop=fastStop, singleCore=singleCore, numChunks=numChunks)
    dim(meta.augTrain.shrink)
    #pheatmap(augTrain.shrink)
    resList[['matrix.metaAugShrink']] <- meta.augTrain.shrink
    resList[['matrix.metaAugShrink.bestGenes']] <- ADAPTS::matrixToGenelist(meta.augTrain.shrink, metagList)
    
    estimates.Meta.shrink <- as.data.frame(ADAPTS::estCellPercent.DCQ(meta.augTrain.shrink, metapseudobulk.test))
    colnames(estimates.Meta.shrink) <- paste('Shrunk Meta')
    resList[['estimates.onTest.meta']] <- estimates.Meta.onTest <- cbind(estimates.Meta.shrink, estimates.Meta.onTest)
    
    resList[['testAcc.metaAugShrink']] <- seed2TestAcc.shrink.meta <- ADAPTS::calcAcc(estimates=estimates.Meta.onTest[,1], reference=estimates.Meta.onTest[,ncol(estimates.Meta.onTest)])
  }
  
  return(resList)
}


#' Find out at which iteration the results converge, i.e. the mean results are stable.
#'
#' @param curSeq A sequence of results that generated from each iteration of the loop
#' @param changePer The maximum percentage of change allowed
#' @param winSize  The window size for mean calculation
#'
#' @return The minimum number of iterations needed for the results to converge

findConvergenceIter <- function(curSeq, changePer=1, winSize=5) {
  
  #Note, this will remove NAs.  Is that best?  They're caused by bad correlations
  runMean <- sapply(1:length(curSeq), function(x) {sum(curSeq[1:x], na.rm=TRUE)/x})
  #Criteria: running mean has changed less than 5% in the last 5? point
  convIter <- as.numeric(NA)
  if (length(runMean) > winSize) {
    winOffset <- winSize-1
    maxWinChange <- sapply(winSize:length(runMean), function(x) {
      win <- runMean[(x-winOffset):x]
      max(abs((win - mean(win))/mean(win)))
    }) #maxWinChange
    mwcBool <- maxWinChange < changePer/100
    if(any(mwcBool)) {convIter <- winOffset + which(mwcBool)[1]}
  }
  return (convIter)
}



#' A meta analysis for the results from multiple iterations
#' @description Calculate the mean and the standard deviation of the results from all the iterations, and also 
#' test for convergence by % of change with each additional iteration.
#' 
#' @param allResList A list of results generated from all the iterative calls of testAllSigMatrices
#' @param changePer  The maximum percentage of change allowed for convergence
#' @export
#'
#' @return The mean and standard deviation of all the results, along with the number of iterations needed for the results to converge.
#' A meta analysis for the results from multiple iterations
#' @description Calculate the mean and the standard deviation of the results from all the iterations, and also 
#' test for convergence by % of change with each additional iteration.
#' 
#' @param allResList A list of results generated from all the iterative calls of testAllSigMatrices
#' @param changePer  The maximum percentage of change allowed for convergence
#' @export
#'
#' @return The mean and standard deviation of all the results, along with the number of iterations needed for the results to converge.
meanResults <- function (allResList,changePer=1) {
  testNames <- unique(sub('^.*\\.', '', names(allResList[[1]])))
  testNames <- testNames[!testNames %in% c("onTest", "allClusters", "LUT", "meta", "gList")]
  compTypes <- names(allResList[[1]][[paste0('testAcc.', testNames[1])]])
  
  allResList <- allResList[!sapply(allResList, function(x){inherits(x,'try-error')})]
  
  compList <- list()
  for (curName in testNames) {
    compList[[curName]] <- list()
    for (curComp in compTypes) {
      compList[[curName]][[curComp]] <- sapply(allResList, function(x){ 
        y <- x[[paste0('testAcc.', curName)]]
        if(!inherits(y, 'try-error')) {return(y[[curComp]])} else {return(as.numeric(NA))}  
        #If only two clusters, set correlation to zero or NA?
      }) #compList[[curName]][[curComp]] <- sapply(allResList, function(x){ 
    } #for (curComp in compTypes) {
  } #for (curName in testNames) {
  
  curColN <- unlist(lapply(compTypes, function(x){c(x, paste0(x,'.sd'))}))
  resMat <- matrix(as.numeric(NA), nrow=length(testNames), ncol=length(curColN),
                   dimnames=list(testNames, curColN))
  convMat <- matrix(as.numeric(NA), nrow=length(testNames), ncol=length(compTypes),
                    dimnames=list(testNames, compTypes))
  for (curName in testNames) {
    if(is.null(compList[[curName]][[compTypes[1]]][[1]])) { message(paste('Omitting', curName, 'due to NULL results')) ; next; }
    
    for (curComp in compTypes) {
      curMean <- mean(compList[[curName]][[curComp]], na.rm = TRUE)
      curSD <- stats::sd(compList[[curName]][[curComp]], na.rm = TRUE)
      resMat[curName,curComp] <- curMean
      resMat[curName,paste0(curComp,'.sd')] <- curSD
      
      #Also test for convergence by % of change with each additional
      convMat[curName, curComp] <- findConvergenceIter(curSeq=compList[[curName]][[curComp]], changePer=changePer, winSize=5)
      
    } #for (curComp in compTypes) {
  } #for (curName in testNames) {
  colnames(convMat) <- paste0('convIt.',colnames(convMat))
  
  resMat <- as.data.frame(resMat)
  resMat$N <- length(allResList)
  resMat <- cbind(resMat, as.data.frame(convMat))
  
  remBool <- apply(resMat, 1, function(x){all(is.na(x))})
  resMat <- resMat[!remBool,]
  
  return(resMat)
}


#' Loop testAllSigMatrices until convergence
#' @description Iteratively call testAllSigMatrices numLoops times with the option to fast stop 
#'  if correlation, correlation spear, mae and rmse all converge
#' @param numLoops The number of iterations. Set to null to loop until results converge.
#' @param fastStop Set to TRUE to break the loop when correlation, correlation spear, mae and rmse all converge
#' @param exprData The single cell matrix
#' @param changePer The maximum percentage of change allowed for convergence 
#' @param handMetaCluster A List of pre-defined meta clusters. Set to NULL to automatically group indistinguishable 
#' cells into same cluster use clustWspillOver (DEFAULT: NULL)
#' @param testOnHalf Set to TRUE to leave half the data as a test set to validate all the matrices
#' @param condTol  The tolerance in the reconstruction algorithm.  1.0 = no tolerance, 1.05 = 5\% tolerance (DEFAULT: 1.01)
#'
#' @return  A list of results generated from all the iterative calls of testAllSigMatrices
#' @export
#'
#' @examples
#' ct1 <- runif(1000, 0, 100)
#' ct2 <- runif(1000, 0, 100)
#' ct3 <- runif(1000, 0, 100)
#' ct4 <- runif(1000, 0, 100)
#' dataMat <- cbind(ct1, ct1, ct1, ct1, ct1, ct1, ct2, ct2, ct2, ct2, ct3, ct3, ct3,ct3,ct4,ct4)
#' rownames(dataMat) <- make.names(rep('gene', nrow(dataMat)), unique=TRUE)
#' noise <- matrix(runif(nrow(dataMat)*ncol(dataMat), -2, 2), nrow = nrow(dataMat), byrow = TRUE)
#' dataMat <- dataMat + noise
#' #options(mc.cores=2)
#' #  This is a meta-function that calls other functions, 
#' #  The execution speed is too slow for the CRAN automated check
#' #loopTillConvergence(numLoops=10, fastStop=TRUE, exprData=dataMat, 
#' #    changePer=10,handMetaCluster=NULL, testOnHalf=TRUE)
loopTillConvergence <- function(numLoops, fastStop, exprData, changePer, handMetaCluster, testOnHalf, condTol=1.01){
  if(is.null(numLoops)){
    fastStop<-TRUE
    numLoops<-1000000
  }
  allResListOut <- list()
  for (i in 1:numLoops) {
    curName <- paste0('res', i)
    allResListOut[[curName]] <- try(testAllSigMatrices(exprData, randomize = TRUE, proportional=FALSE, handMetaCluster=handMetaCluster, testOnHalf=testOnHalf, condTol = condTol))
    if(fastStop==TRUE){
      covtmp<-meanResults(allResList=allResListOut,changePer)[ ,c("convIt.rho.cor", "convIt.spear.rho", "convIt.mae","convIt.rmse")]
      if(all(!is.na(covtmp))) break
    }}
  return(allResListOut)
}
