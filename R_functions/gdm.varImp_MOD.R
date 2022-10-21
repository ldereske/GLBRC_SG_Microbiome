#modified from 'gdm.varImp' from GDM package. 
#This code give you variable importance and p-values for models with 3 or less varibales
#This code does not do model selection
#Please use citation(package="gdm") for citation 

gdm.varImp_MOD<-function (spTable, geo, splines = NULL, knots = NULL, 
          nPerm = 50, pValue = 0.05, parallel = FALSE, cores = 2, 
          sampleSites = 1, sampleSitePairs = 1, outFile = NULL) 
{
  library(pbapply)
  permutateSitePair <- function(spTab, siteVarTab, indexTab, vNames){
    
    #################
    #spTab <- currSitePair    ##site-pair table
    #siteVarTab <- siteData   ##siteXvar table
    #indexTab <- indexTab     ##table of the index of sites
    #vNames <- varNames       ##variables names
    #vNames <- c("awcA", "phTotal", "shcA", "solumDepth", "bio5", "bio19")
    #################
    
    ##randomizes the row order of the given siteXvar table
    randVarTab <- siteVarTab[sample(nrow(siteVarTab), nrow(siteVarTab)), ]
    
    #site1x <- siteVarTab$xCoord[1]
    #site1y <- siteVarTab$yCoord[1]
    #checkingIn <- siteVarTab[siteVarTab$xCoord==site1x & siteVarTab$yCoord==site1y,]
    #checkX <- siteVarTab[siteVarTab$xCoord==site1x,]
    #checkingRand <- randVarTab[randVarTab$xCoord==site1x & randVarTab$yCoord==site1y,]
    
    ##sets up the coordinate values for the randomized site-pair table
    s1xCoord <- sapply(1:nrow(spTab), function(i){randVarTab[indexTab[i,1],1]})
    s1yCoord <- sapply(1:nrow(spTab), function(i){randVarTab[indexTab[i,1],2]})
    s2xCoord <- sapply(1:nrow(spTab), function(i){randVarTab[indexTab[i,2],1]})
    s2yCoord <- sapply(1:nrow(spTab), function(i){randVarTab[indexTab[i,2],2]})
    
    #print(vNames)
    ##extracts values of other variables
    varLists <- lapply(vNames, function(vn, rvTab, spt, inT){if(vn!="Geographic"){
      ###################
      #vn <- vNames[[2]]
      #rvTab=randVarTab
      #spt=spTab
      #inT=indexTab
      ###################
      ##identifies variable columns in randVarTab
      randCols <- grep(paste("^", vn, "$", sep=""), colnames(rvTab))
      #print(randCols)
      ##identifies variable columns in site-pair table
      spCols <- grep(vn, colnames(spt))
      
      s1var <- sapply(1:nrow(spt), function(i){rvTab[inT[i,1],randCols]})
      s2var <- sapply(1:nrow(spt), function(i){rvTab[inT[i,2],randCols]})
      
      return(list(s1var, s2var))
    }
    }, rvTab=randVarTab, spt=spTab, inT=indexTab)
    
    # unravels the varList into a data.frame of the variable portion of a site-pair table
    bySite <- lapply(1:2, function(i,vlist){sapply(vlist, function(vl,k){vl[[k]]},k=i)}, vlist=varLists)
    
    if(is(bySite[[1]], "list")){
      site1Vars <- do.call("cbind", bySite[[1]])
      site2Vars <- do.call("cbind", bySite[[2]])
    }else{
      site1Vars <- bySite[[1]]
      site2Vars <- bySite[[2]]
    }
    
    ##sets up new site-pair table
    newSP <- as.data.frame(cbind(spTab$distance, spTab$weights, s1xCoord, s1yCoord, s2xCoord, s2yCoord, site1Vars, site2Vars))
    colnames(newSP) <- colnames(spTab)
    class(newSP) <- c(class(spTab))
    
    #getCoords1 <- newSP[newSP$s1.xCoord==site1x,]
    #getCoords2 <- newSP[newSP$s2.xCoord==site1x,]
    
    return(newSP)
  }
  k <- NULL
  if (!is(spTable, "gdmData")) {
    warning("The spTable object is not of class 'gdmData'. See the formatsitepair function for help.")
  }
  if (!(is(spTable, "gdmData") | is(spTable, "matrix") | is(spTable, 
                                                            "data.frame"))) {
    stop("spTable argument needs to be of class 'gdmData', 'matrix', or 'data frame'")
  }
  if (ncol(spTable) < 6) {
    stop("spTable object requires at least 6 columns: distance, weights, s1.xCoord, s1.yCoord, s2.xCoord, s2.yCoord")
  }
  if (nrow(spTable) < 1) {
    stop("The spTable object contains zero rows of data.")
  }
  if (!(geo == TRUE | geo == FALSE)) {
    stop("The geo argument must be either TRUE or FALSE.")
  }
  if (is.null(splines) == FALSE & !is(splines, "numeric")) {
    stop("The splines argument needs to be a numeric data type.")
  }
  if (is.null(knots) == FALSE & !is(knots, "numeric")) {
    stop("The knots argument needs to be a numeric data type.")
  }
  if ((is.null(nPerm) == FALSE & is.numeric(nPerm) == FALSE) | 
      nPerm < 1) {
    stop("The nPerm argument needs to be a positive integer.")
  }
  if (!(parallel == TRUE | parallel == FALSE)) {
    stop("The parallel argument must be either TRUE or FALSE.")
  }
  if (parallel == TRUE & is.null(cores) == TRUE) {
    stop("If parallel==TRUE, the number of cores must be specified.")
  }
  if ((is.null(cores) == FALSE & is.numeric(cores) == FALSE) | 
      cores < 1) {
    stop("The cores argument needs to be a positive integer.")
  }
  if (is.numeric(sampleSites) == FALSE | sampleSites < 0 | 
      sampleSites > 1) {
    stop("The sampleSites argument needs to be a positive number between 0 and 1.")
  }
  if (is.numeric(sampleSitePairs) == FALSE | sampleSitePairs < 
      0 | sampleSitePairs > 1) {
    stop("The sampleSitePairs argument needs to be a positive number between 0 and 1.")
  }
  if (sampleSites == 0) {
    stop("A sampleSites value of 0 will remove all sites from the analysis.")
  }
  if (sampleSitePairs == 0) {
    stop("A sampleSitePairs value of 0 will remove all sites from the analysis.")
  }
  if (is.null(outFile) == FALSE) {
    if (is.character(outFile) == FALSE) {
      stop("The outFile argument needs to be a character string of the directory and file name you wish the tables to be written to")
    }
    outFileChar <- nchar(outFile)
    if (substr(outFile, outFileChar - 5, outFileChar) != 
        ".RData") {
      outFile <- paste(outFile, ".RData", sep = "")
    }
    if (length(strsplit(outFile, "/")[[1]]) > 1) {
      splitOutFile <- strsplit(outFile, "/")[[1]][-length(strsplit(outFile, 
                                                                   "/")[[1]])]
      dir.create(paste(splitOutFile, collapse = "/"))
    }
    else {
      outFile <- paste("./", outFile, sep = "")
    }
  }
  nPerm <- as.integer(nPerm)
  cores <- as.integer(cores)
  if (sampleSites < 1) {
    spTable <- subsample.sitepair(spTable, sampleSites = sampleSites)
    if (sampleSitePairs < 1) {
      warning("You have selected to randomly remove sites and/or site-pairs.")
    }
  }
  if (sampleSitePairs < 1) {
    numRm <- sample(1:nrow(spTable), round(nrow(spTable) * 
                                             (1 - sampleSitePairs)))
    spTable <- spTable[-c(numRm), ]
  }
  rtmp <- spTable[, 1]
  if (length(rtmp[rtmp < 0]) > 0) {
    stop("The spTable contains negative distance values. Must be between 0 - 1.")
  }
  if (length(rtmp[rtmp > 1]) > 0) {
    stop("The spTable contains distance values greater than 1. Must be between 0 - 1.")
  }
  nVars <- (ncol(spTable) - 6)/2
  varNames <- colnames(spTable[c(7:(6 + nVars))])
  varNames <- sapply(strsplit(varNames, "s1."), "[[", 2)
  if (geo == TRUE) {
    nVars <- nVars + 1
    varNames <- c("Geographic", varNames)
  }
  message(paste0("Fitting initial model with all ", nVars, 
                 " predictors..."))
  Sys.sleep(1)
  fullGDM <- gdm(spTable, geo = geo, splines = splines, knots = knots)
  thiscoeff <- 1
  thisquant <- 1
  sumCoeff <- NULL
  for (i in 1:length(fullGDM$predictors)) {
    numsplines <- fullGDM$splines[[i]]
    holdCoeff <- NULL
    for (j in 1:numsplines) {
      holdCoeff[j] <- fullGDM$coefficients[[thiscoeff]]
      thiscoeff <- thiscoeff + 1
    }
    sumCoeff[i] <- sum(holdCoeff)
  }
  zeroSum <- fullGDM$predictors[which(sumCoeff == 0)]
  if (length(zeroSum) > 0) {
    message(paste0("Sum of I-spline coefficients for predictors(s) ", 
                   zeroSum, " = 0"))
    Sys.sleep(1)
    message(paste0("Removing predictor(s) ", zeroSum, " and proceeding with permutation testing..."))
    Sys.sleep(1)
    for (z in zeroSum) {
      testVarCols1 <- grep(paste("^s1.", z, "$", sep = ""), 
                           colnames(spTable))
      testVarCols2 <- grep(paste("^s2.", z, "$", sep = ""), 
                           colnames(spTable))
      spTable <- spTable[, -c(testVarCols1, testVarCols2)]
    }
  }
  nVars <- (ncol(spTable) - 6)/2
  varNames <- colnames(spTable[c(7:(6 + nVars))])
  varNames <- sapply(strsplit(varNames, "s1."), "[[", 2)
  if (geo == TRUE) {
    nVars <- nVars + 1
    varNames <- c("Geographic", varNames)
  }
  if (cores > nVars) {
    cores <- nVars
  }
  splines <- rep(unique(splines), nVars)
  sortMatX <- sapply(1:nrow(spTable), function(i, spTab) {
    c(spTab[i, 3], spTab[i, 5])
  }, spTab = spTable)
  sortMatY <- sapply(1:nrow(spTable), function(i, spTab) {
    c(spTab[i, 4], spTab[i, 6])
  }, spTab = spTable)
  sortMatNum <- sapply(1:nrow(spTable), function(i) {
    c(1, 2)
  })
  sortMatRow <- sapply(1:nrow(spTable), function(i) {
    c(i, i)
  })
  fullSortMat <- cbind(as.vector(sortMatX), as.vector(sortMatY), 
                       as.vector(sortMatNum), as.vector(sortMatRow), rep(NA, 
                                                                         length(sortMatX)))
  siteByCoords <- as.data.frame(unique(fullSortMat[, 1:2]))
  numSites <- nrow(siteByCoords)
  for (i in 1:numSites) {
    fullSortMat[which(fullSortMat[, 1] == siteByCoords[i, 
                                                       1] & fullSortMat[, 2] == siteByCoords[i, 2]), 5] <- i
  }
  indexTab <- matrix(NA, nrow(spTable), 2)
  for (iRow in 1:nrow(fullSortMat)) {
    indexTab[fullSortMat[iRow, 4], fullSortMat[iRow, 3]] <- fullSortMat[iRow, 
                                                                        5]
  }
  rm(fullSortMat)
  rm(sortMatX)
  rm(sortMatY)
  rm(sortMatNum)
  rm(sortMatRow)
  rm(siteByCoords)
  exBySite <- lapply(1:numSites, function(i, index, tab) {
    rowSites <- which(index[, 1] %in% i)
    if (length(rowSites) < 1) {
      rowSites <- which(index[, 2] %in% i)
    }
    exSiteData <- tab[rowSites[1], ]
    return(exSiteData)
  }, index = indexTab, tab = spTable)
  outSite <- which(!(1:numSites %in% indexTab[, 1]))
  for (i in 1:length(exBySite)) {
    siteRow <- exBySite[[i]]
    if (i %in% outSite) {
      siteRow <- siteRow[grep("s2.", colnames(siteRow))]
      colnames(siteRow) <- sapply(strsplit(colnames(siteRow), 
                                           "s2."), "[[", 2)
    }
    else {
      siteRow <- siteRow[grep("s1.", colnames(siteRow))]
      colnames(siteRow) <- sapply(strsplit(colnames(siteRow), 
                                           "s1."), "[[", 2)
    }
    exBySite[[i]] <- siteRow
  }
  siteData <- do.call("rbind", exBySite)
  modelTestValues <- data.frame(matrix(NA, 4, 1, dimnames = list(c("Model deviance", 
                                                                       "Percent deviance explained", "Model p-value", "Fitted permutations"), 
                                                                     c("All predictors"))))
  varImpTable <- matrix(NA, nVars, 1)
  rownames(varImpTable) <- varNames
  colnames(varImpTable) <- c("All predictors")
  pValues <- nModsConverge <- varImpTable
  currSitePair <- spTable
  nullGDMFullFit <- 0
  message(paste0("Creating ", nPerm, " permuted site-pair tables..."))
  if (parallel == F | nPerm <= 50) {
    permSpt <- pbreplicate(nPerm, list(permutateSitePair(currSitePair, 
                                                         siteData, indexTab, varNames)))
  }
  if (parallel == T & nPerm > 50) {
    cl <- makeCluster(cores)
    registerDoParallel(cl)
    permSpt <- foreach(k = 1:nPerm, .verbose = F, .packages = c("gdm"), 
                       .export = c("permutateSitePair")) %dopar% permutateSitePair(currSitePair, 
                                                                                   siteData, indexTab, varNames)
    stopCluster(cl)
  }
  varNames.x <- varNames
  message("Starting model assessment...")
  for (v in 1:length(varNames)) {
    if (length(varNames.x) < 2) {
      break
    }
    if (is.numeric(splines)) {
      splines <- rep(unique(splines), length(varNames.x))
    }
    fullGDM <- gdm(currSitePair, geo = geo, splines = splines, 
                   knots = knots)
    message(paste0("Percent deviance explained by the full model =  ", 
                   round(fullGDM$explained, 3)))
    if (is.null(fullGDM) == TRUE) {
      warning(paste("The model did not converge when testing variable: ", 
                    varNames.x[v], ". Terminating analysis and returning output completed up to this point.", 
                    sep = ""))
      break
    }
    message("Fitting GDMs to the permuted site-pair tables...")
    permGDM <- lapply(permSpt, function(x) {
      gdm(x, geo = geo, splines = splines, knots = knots)
    })
    permModelDev <- sapply(permGDM, function(mod) {
      mod$gdmdeviance
    })
    modPerms <- length(which(sapply(permModelDev, is.null) == 
                               TRUE))
    if (modPerms > 0) {
      permModelDev <- unlist(permModelDev[-(which(sapply(permModelDev, 
                                                         is.null) == T))])
    }
    if (v > 1) {
      colnames(modelTestValues)[v] <- colnames(pValues)[v] <- colnames(varImpTable)[v] <- colnames(nModsConverge)[v]
    }
    modelTestValues[1, v] <- round(fullGDM$gdmdeviance, 
                                   3)
    modelTestValues[2, v] <- round(fullGDM$explained, 3)
    modelTestValues[3, v] <- round(sum(permModelDev <= fullGDM$gdmdeviance)/(nPerm - 
                                                                               modPerms), 3)
    modelTestValues[4, v] <- round(nPerm - modPerms, 0)
    if (parallel == TRUE) {
      if (length(varNames.x) < cores) {
        cores <- length(varNames.x)
      }
      cl <- makeCluster(cores)
      registerDoParallel(cl)
      permVarDev <- foreach(k = 1:length(varNames.x), 
                            .verbose = F, .packages = c("gdm"), .export = c("currSitePair")) %dopar% 
        {
          if (varNames.x[k] != "Geographic") {
            lll <- lapply(permSpt, function(x, spt = currSitePair) {
              idx <- grep(varNames.x[k], colnames(x))
              spt[, idx] <- x[, idx]
              return(spt)
            })
          }
          if (varNames.x[k] == "Geographic") {
            lll <- lapply(permSpt, function(x, spt = currSitePair) {
              s1 <- sample(1:nrow(spt), nrow(spt))
              s2 <- sample(1:nrow(spt), nrow(spt))
              s3 <- sample(1:nrow(spt), nrow(spt))
              s4 <- sample(1:nrow(spt), nrow(spt))
              spt[, 3] <- spt[s1, 3]
              spt[, 4] <- spt[s2, 4]
              spt[, 5] <- spt[s3, 5]
              spt[, 6] <- spt[s4, 6]
              return(spt)
            })
          }
          gdmPermVar <- lapply(lll, function(x) {
            try(gdm(x, geo = geo, splines = splines, 
                    knots = knots))
          })
          permModelDev <- sapply(gdmPermVar, function(mod) {
            mod$gdmdeviance
          })
          return(permModelDev)
        }
      stopCluster(cl)
    }
    if (parallel == FALSE) {
      permVarDev <- list()
      for (k in 1:length(varNames.x)) {
        if (varNames.x[k] != "Geographic") {
          message(paste0("Assessing importance of ", 
                         varNames.x[k], "..."))
          lll <- lapply(permSpt, function(x, spt = currSitePair) {
            idx <- grep(varNames.x[k], colnames(x))
            spt[, idx] <- x[, idx]
            return(spt)
          })
        }
        if (varNames.x[k] == "Geographic") {
          message("Assessing importance of geographic distance...")
          lll <- lapply(permSpt, function(x, spt = currSitePair) {
            s1 <- sample(1:nrow(spt), nrow(spt))
            s2 <- sample(1:nrow(spt), nrow(spt))
            s3 <- sample(1:nrow(spt), nrow(spt))
            s4 <- sample(1:nrow(spt), nrow(spt))
            spt[, 3] <- spt[s1, 3]
            spt[, 4] <- spt[s2, 4]
            spt[, 5] <- spt[s3, 5]
            spt[, 6] <- spt[s4, 6]
            return(spt)
          })
        }
        gdmPermVar <- lapply(lll, function(x) {
          try(gdm(x, geo = geo, splines = splines, knots = knots))
        })
        permVarDev[[k]] <- sapply(gdmPermVar, function(mod) {
          mod$gdmdeviance
        })
      }
    }
    names(permVarDev) <- varNames.x
    nullDev <- fullGDM$nulldeviance
    for (var in varNames.x) {
      grepper <- grep(var, names(permVarDev))
      varDevTab <- permVarDev[[grepper]]
      nConv <- length(which(is.null(varDevTab)))
      nModsConverge[which(rownames(varImpTable) == var), 
                    v] <- nPerm - nConv
      if (nConv > 0) {
        varDevTab <- unlist(varDevTab[-(which(sapply(is.null(varDevTab))))])
      }
      varDevExplained <- 100 * (nullDev - varDevTab)/nullDev
      varImpTable[which(rownames(varImpTable) == var), 
                  v] <- median(100 * abs((varDevExplained - fullGDM$explained)/fullGDM$explained))
      if (var != "Geographic") {
        testVarCols1 <- grep(paste("^s1.", var, "$", 
                                   sep = ""), colnames(currSitePair))
        testVarCols2 <- grep(paste("^s2.", var, "$", 
                                   sep = ""), colnames(currSitePair))
        testSitePair <- currSitePair[, -c(testVarCols1, 
                                          testVarCols2)]
        noVarGDM <- gdm(testSitePair, geo = geo, splines = splines[-1], 
                        knots = knots)
      }
      else {
        noVarGDM <- gdm(currSitePair, geo = F, splines = splines[-1], 
                        knots = knots)
      }
      permDevReduct <- noVarGDM$gdmdeviance - varDevTab
      pValues[which(rownames(pValues) == var), v] <- sum(permDevReduct >= 
                                                           (varDevTab - fullGDM$gdmdeviance))/(nPerm - 
                                                                                                 nConv)
    }
    if (max(na.omit(pValues[, v])) < pValue) {
      message("All remaining predictors are significant, ceasing assessment.")
      message(paste0("Percent deviance explained by final model = ", 
                     round(fullGDM$explained, 3)))
      message("Final set of predictors returned: ")
      for (vvv in 1:length(fullGDM$predictors)) {
        message(fullGDM$predictors[vvv])
      }
      break
    }
    else {
      geo <- F
    }
        if (v == 1) {
      message("One variable. Ceasing assessment.")
      message(paste0("Percent deviance explained by final model = ", 
                     round(fullGDM$explained, 3)))
      message("Final set of predictors returned: ")
      for (vvv in 1:length(fullGDM$predictors)) {
        message(fullGDM$predictors[vvv])
      }
      break
    }
  }
  if (v == 1) {
    modelTestVals <- data.frame(matrix(round(modelTestValues[, 
                                                             1], 3), ncol = 1))
    rownames(modelTestVals) <- rownames(modelTestValues)
    colnames(modelTestVals) <- "All predictors"
    varImpTab <- data.frame(matrix(round(varImpTable[, 1], 
                                         3), ncol = 1))
    rownames(varImpTab) <- rownames(varImpTable)
    colnames(varImpTab) <- "All predictors"
    pVals <- varImpTab
    pVals[, 1] <- round(pValues[, 1], 3)
    nModsConv <- varImpTab
    nModsConv[, 1] <- round(nModsConverge[, 1], 3)
    outObject <- list(modelTestVals, varImpTab, pVals, nModsConv)
    names(outObject) <- c("Model assessment", "Predictor Importance", 
                          "Predictor p-values", "Model Convergence")
  }
  else {
    outObject <- list(round(modelTestValues[, 1:v], 3), 
                      round(varImpTable[, 1:v], 3), round(pValues[, 1:v], 
                                                          3), nModsConverge[, 1:v])
    names(outObject) <- c("Model assessment", "Predictor Importance", 
                          "Predictor p-values", "Model Convergence")
  }
  if (is.null(outFile) == FALSE) {
    save(outObject, file = outFile)
  }
  return(outObject)
}
