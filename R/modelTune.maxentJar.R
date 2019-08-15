#################################################
#########	MODEL TUNE for maxent.jar #############
#################################################

modelTune.maxentJar <- function(pres, bg, env, nk, group.data, args.i, userArgs, 
                                rasterPreds, clamp, categoricals, path, outputformat) {
  
  # set up data: x is coordinates of occs and bg, 
  # p is vector of 0's and 1's designating occs and bg
  x <- rbind(pres, bg)
  p <- c(rep(1, nrow(pres)), rep(0, nrow(bg)))
  
  # build the full model from all the data
  full.mod <- dismo::maxent(x, p, args = c(args.i, userArgs), removeDuplicates = TRUE,
                            factors = categoricals, path = paste(path,args.i))  
  pred.args <- c(paste0("outputformat=",outputformat), ifelse(clamp==TRUE, "doclamp=true", "doclamp=false"))
  
  # if rasters selected, predict for the full model
  if (rasterPreds == TRUE) {
    predictive.map <- predict(full.mod, env, args = pred.args)  
  } else {
    predictive.map <- stack()
  }
  
  # set up empty vectors for stats
  AUC.TEST <- double()
  AUC.DIFF <- double()
  OR10 <- double()
  ORmin <- double()
  KAPPA <- double()
  MAX.F1 <- double()
  
  # cross-validation on partitions
  for (k in 1:nk) {
    # set up training and testing data groups
    train.val <- pres[group.data$occ.grp != k,, drop = FALSE]
    test.val <- pres[group.data$occ.grp == k,, drop = FALSE]
    bg.val <- bg[group.data$bg.grp != k,, drop = FALSE]
    # redefine x and p for partition groups
    x <- rbind(train.val, bg.val)
    p <- c(rep(1, nrow(train.val)), rep(0, nrow(bg.val)))
    
    # run the current test model
    mod <- dismo::maxent(x, p, args = c(args.i, userArgs), factors = categoricals)  
    eval <- dismo::evaluate(test.val, bg, mod)
    
    recall <- eval@TPR
    precision <- eval@confusion[,1]/(eval@confusion[,1]+eval@confusion[,2])
    
    MAX.F1[k] <- max(2*((precision*recall)/(precision+recall)), na.rm = TRUE)
    AUC.TEST[k] <- eval@auc
    AUC.DIFF[k] <- max(0, eval@auc - AUC.TEST[k])
    KAPPA[k] <- threshold(eval)$kappa
    
    # predict values for training and testing data
    p.train <- predict(mod, train.val, args = pred.args)
    p.test <- predict(mod, test.val, args = pred.args)  
    
    # figure out 90% of total no. of training records
    if (nrow(train.val) < 10) {
      n90 <- floor(nrow(train.val) * 0.9)
    } else {
      n90 <- ceiling(nrow(train.val) * 0.9)
    }
    train.thr.10 <- rev(sort(p.train))[n90]
    OR10[k] <- mean(p.test < train.thr.10)
    train.thr.min <- min(p.train)
    ORmin[k] <- mean(p.test < train.thr.min)
  }
  stats <- c(AUC.DIFF, AUC.TEST, OR10, ORmin, KAPPA, MAX.F1)
  out.i <- list(full.mod, stats, predictive.map)
  return(out.i)
}

