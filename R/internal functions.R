# I suggest creating another internal function that does all the necessary calculations and returns a list as described in the attached PDF file. 
# For this we just need to move a part of the code from fit.cox and fit.nonpar to separate internal functions - I put those pieces in the attached R file - you'll need to add the necessary input arguments in these functions - I think they will be similar to what we tried to pass yesterday. 
# Do you see where I moved these pieces of code from?

# The next step will be simple aggregation without any calculations, and you said you know how to do that, right?  

# Some remarks:
# -  fit.cox and fit.nonpar will call  these internal functions both for the original data and for bootstrap samples. 
# - the time vector from the original data should be passed to all the bootstrap replications. 
# - we need to figure out how to use seed in order to make our results reproducible.  

# Let me know what you think and if we need to connect on Zoom!



.nonpar.run <- function()
{
  if (wtype!="unadj")
  {
    ps.fit <- get.weights(df, A, C, wtype)
    est.0 <- .estimate.nonpar(X=X[trt==0], E=E[trt==0],
                              case.w =ps.fit$w[trt==0], cens=cens, time=time, E.set)
    est.1 <- .estimate.nonpar(X=X[trt==1], E=E[trt==1],
                              case.w=ps.fit$w[trt==1], cens=cens, time=time, E.set)
  }
  else
  {
    un.w <- rep(1, nobs)
    est.0 <- .estimate.nonpar(X=X[trt==0], E=E[trt==0],
                              case.w =un.w[trt==0], cens=cens, time=time, E.set)
    est.1 <- .estimate.nonpar(X=X[trt==1], E=E[trt==1],
                              case.w=un.w[trt==1], cens=cens, time=time, E.set)
  }

  res$trt.0 <- est.0 # save point estimates in the trt world 0
  res$trt.1 <- est.1 # save point estimates in the trt world 1

  # calculate and save trt effect measures:
  for (k in E.set)
  {
    # calculate log(CumHaz1/CumHaz0)
    b <- log(est.1[[paste("Ev=", k, sep="")]]$CumHaz) -
      log(est.0[[paste("Ev=", k, sep="")]]$CumHaz) # under PH, this is beta=log(HR_TRT)
    # without PH, it is log(HR_1(t)/HR_0(t))

    CIF.1 <- est.1[[paste("Ev=", k, sep="")]]$CIF
    CIF.0 <- est.0[[paste("Ev=", k, sep="")]]$CIF
    RMT.1 <- est.1[[paste("Ev=", k, sep="")]]$RMT
    RMT.0 <- est.0[[paste("Ev=", k, sep="")]]$RMT

    res$trt.eff[[paste("Ev=", k, sep="")]] <- list(log.CumHazR=b,
                                                   RD=CIF.1-CIF.0, RR=CIF.1/CIF.0,
                                                   ATE.RMT=RMT.1-RMT.0)
  } # res is a list with 4 fields: time, trt.0, trt.0, trt.eff
  res
}

.cox.run <- function(df, T, E, A, C, wtype, cens, E.set,time,trt,nobs,X,case.w)
{
  if (wtype!="unadj")
  {
    ps.fit <- get.weights(df, A, C, wtype,case.w = case.w)
    est <- .estimate.cox(X=X, E=E, trt=trt, case.w = case.w*ps.fit$w, cens=cens, time=time, E.set=E.set)
  }
  else
    est <- .estimate.cox(X=X, E=E, trt=trt, case.w=case.w, cens=cens, time=time, E.set=E.set)

  for (k in E.set)
  {
    # CumHaz:
    time.k <- est[[paste("Ev=", k, sep="")]]$bh$time
    res$trt.0[[paste("Ev=", k, sep="")]][["CumHaz"]] <-
      sapply(time, .get.CIF, time.k, est[[paste("Ev=", k, sep="")]]$bh$haz)
    res$trt.1[[paste("Ev=", k, sep="")]][["CumHaz"]] <-
      sapply(time, .get.CIF, time.k, est[[paste("Ev=", k, sep="")]]$bh$haz * exp(est[[paste("Ev=", k, sep="")]]$logHR))

    # CIF:
    bh.k <- .base.haz.std(est[[paste("Ev=", k, sep="")]]$bh)
    cif.0.k <- cumsum(sapply(bh.k$time, .get.S, est$OS$time, est$OS$OS.trt.0) * bh.k$haz)
    cif.1.k <- cumsum(sapply(bh.k$time, .get.S, est$OS$time, est$OS$OS.trt.1) * bh.k$haz * exp(est[[paste("Ev=", k, sep="")]]$logHR))
    res$trt.0[[paste("Ev=", k, sep="")]][["CIF"]] <- sapply(time, .get.CIF, time=bh.k$time, cuminc=cif.0.k)
    res$trt.1[[paste("Ev=", k, sep="")]][["CIF"]] <- sapply(time, .get.CIF, time=bh.k$time, cuminc=cif.1.k)
    # RMT:
    res$trt.0[[paste("Ev=", k, sep="")]][["RMT"]] <- sapply(time, .get.RMT, time=bh.k$time, cuminc=cif.0.k)
    res$trt.1[[paste("Ev=", k, sep="")]][["RMT"]] <- sapply(time, .get.RMT, time=bh.k$time, cuminc=cif.1.k)

    RD <- res$trt.1[[paste("Ev=", k, sep="")]]$CIF - res$trt.0[[paste("Ev=", k, sep="")]]$CIF
    RR <- res$trt.1[[paste("Ev=", k, sep="")]]$CIF / res$trt.0[[paste("Ev=", k, sep="")]]$CIF
    ATE.RMT <- res$trt.1[[paste("Ev=", k, sep="")]]$RMT - res$trt.0[[paste("Ev=", k, sep="")]]$RMT

    res$trt.eff[[paste("Ev=", k, sep="")]] <- list(log.CumHazR=est[[paste("Ev=", k, sep="")]]$logHR,
                                                   RD=RD, RR=RR,
                                                   ATE.RMT=ATE.RMT)
  } # res is a list with 4 fields: time, trt.0, trt.0, trt.eff
  res
}


my_list <- list(time = rep(1,10,1))
my_function <- function(object){
  object$trt <- seq(100,110,1);
  return(object)
}
my_list = my_function(my_list)