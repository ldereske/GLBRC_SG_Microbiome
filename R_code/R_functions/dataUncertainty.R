#Modified from 'plotUncertainy' from GDM package. 
#This code gives you the 95% confidence intervals based off of bootstraps
#Please use citation(package="gdm") for citation 

dataUncertainty<-function (spTable, sampleSites, bsIters, geo = FALSE, splines = NULL, 
                           knots = NULL, parallel = FALSE, cores = 2) 
{
  if (!is(spTable, "gdmData")) {
    warning("The spTable object is not of class 'gdmData'. See the formatsitepair function for help.")
  }
  if (!(is(spTable, "gdmData") | is(spTable, "matrix") | is(spTable, 
                                                            "data.frame"))) {
    stop("The spTable object must be of class 'gdmData', 'matrix', or 'data.frame'.")
  }
  if (ncol(spTable) < 6) {
    stop("spTable object requires at least 6 columns: Observed, weights, s1.xCoord, s1.yCoord, s2.xCoord, s2.yCoord")
  }
  if (nrow(spTable) < 1) {
    stop("spTable object has < 1 rows.")
  }
  if (!(geo == TRUE | geo == FALSE)) {
    stop("geo argument must be either TRUE or FALSE")
  }
  if (is.null(splines) == FALSE & !is(splines, "numeric")) {
    stop("splines object must of of class = 'numeric'.")
  }
  if (is.null(knots) == FALSE & !is(knots, "numeric")) {
    stop("knots object must of of class = 'numeric'.")
  }
  if (!(parallel == TRUE | parallel == FALSE)) {
    stop("parallel argument must be either TRUE or FALSE")
  }
  if (parallel == TRUE & is.null(cores) == TRUE) {
    stop("If parallel==TRUE, the number of cores must be specified")
  }
  if ((is.null(cores) == FALSE & is.numeric(cores) == FALSE) | 
      cores < 1) {
    stop("argument cores needs to be a positive integer")
  }
  if ((is.null(bsIters) == FALSE & is.numeric(bsIters) == 
       FALSE) | bsIters < 1) {
    stop("argument bsIters needs to be a positive integer")
  }
  if (is.numeric(sampleSites) == FALSE) {
    stop("sampleSites must be a number between 0 and 1")
  }
  if (sampleSites < 0) {
    stop("sampleSites must be a number between 0 and 1")
  }
  if (sampleSites > 1) {
    stop("sampleSites must be a number between 0 and 1")
  }
  if (sampleSites == 0) {
    stop("a sampleSites value of 0 will remove all sites from the analysis (bad).")
  }
  cores <- as.integer(cores)
  bsIters <- as.integer(bsIters)
  k <- NULL
  pred_data <- NULL
  lstSP <- lapply(1:bsIters, function(i) {
    spTable
  })
  
  if (parallel == TRUE) {
    cl <- makeCluster(cores, outfile = "")
    registerDoParallel(cl)
    subSamps <- foreach(k = 1:length(lstSP), .verbose = F, 
                        .packages = c("gdm")) %dopar% subsample.sitepair(lstSP[[k]], 
                                                                         sampleSites = sampleSites)
    gdmMods <- foreach(k = 1:length(subSamps), .verbose = F, 
                       .packages = c("gdm")) %dopar% gdm(subSamps[[k]], 
                                                         geo = geo, splines = splines, knots = knots)
    stopCluster(cl)
  }
  else {
    subSamps <- lapply(lstSP, subsample.sitepair, sampleSites = sampleSites)
    gdmMods <- lapply(subSamps, gdm, geo = geo, splines = splines, 
                      knots = knots)
  }
  fullGDMmodel <- gdm(spTable, geo = geo, splines = splines, 
                      knots = knots)
  exUncertSplines <- lapply(gdmMods, isplineExtract)
  fullGDMsplines <- isplineExtract(fullGDMmodel)
  predVars <- colnames(exUncertSplines[[1]][[1]])
  totalYmin <- Inf
  totalYmax <- -Inf
  for (p in 1:length(predVars)) {
    predV <- predVars[p]
    for (nm in 1:length(exUncertSplines)) {
      selPlot <- exUncertSplines[[nm]]
      spYmax <- max(selPlot[[2]][, predV])
      spYmin <- min(selPlot[[2]][, predV])
      totalYmax <- max(c(totalYmax, spYmax))
      totalYmin <- min(c(totalYmin, spYmin))
    }
  }
  for (p in 1:length(predVars)) {
    predV <- predVars[p]
    totalXmin <- Inf
    totalXmax <- -Inf
    for (nm in 1:length(exUncertSplines)) {
      selPlot <- exUncertSplines[[nm]]
      spXmax <- max(selPlot[[1]][, predV])
      spXmin <- min(selPlot[[1]][, predV])
      if (spXmax > totalXmax) {
        totalXmax = spXmax
      }
      if (spXmin < totalXmin) {
        totalXmin = spXmin
      }
    }
    if (totalYmax != 0) {
      plotX <- NULL
      plotY <- NULL
      byVarMatX <- NULL
      byVarMatY <- NULL
      for (nn in 1:length(exUncertSplines)) {
        plotX[[nn]] <- exUncertSplines[[nn]][[1]]
        plotY[[nn]] <- exUncertSplines[[nn]][[2]]
        byVarMatY <- cbind(byVarMatY, plotY[[nn]][, 
                                                  predV])
        byVarMatX <- cbind(byVarMatX, plotX[[nn]][, 
                                                  predV])
      }
      fullPlotX <- fullGDMsplines[[1]]
      fullPlotX <- fullPlotX[, predV]
      fullPlotY <- fullGDMsplines[[2]]
      fullPlotY <- fullPlotY[, predV]
      sdX <- apply(as.matrix(byVarMatX), 1, sd)
      sdY <- apply(as.matrix(byVarMatY), 1, sd)
      highBoundX <- fullPlotX + sdX
      lowBoundX <- fullPlotX - sdX
      highBoundY <- fullPlotY + sdY
      lowBoundY <- fullPlotY - sdY
      
      
      pred_data_1<-data.frame(cbind(fullPlotX,fullPlotY,sdX,sdY,highBoundX,lowBoundX,highBoundY,lowBoundY),"factor"=rep(predV))
      pred_data=rbind(pred_data,pred_data_1)
    }
    
  }
  return(pred_data)
}
