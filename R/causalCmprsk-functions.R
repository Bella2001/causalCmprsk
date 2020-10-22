.get.S <- function(t, time, surv)	# assumes that t is a scalar
{
  n <- length(time)
  if (t<time[1] && abs(t-time[1])>100*.Machine$double.eps) return(1)
  if (t>time[n] && abs(t-time[n])>100*.Machine$double.eps) return(surv[n])
  return( surv[sum(time<=t | abs(time-t)<100*.Machine$double.eps)] )
}

# obtain CIF or CumHaz for arbitrary t:
.get.CIF <- function(t, time, cuminc)	# assumes that t is a scalar
{
  n <- length(time)
  if (t<time[1] && abs(t-time[1])>100*.Machine$double.eps) return(0)
  if (t>time[n] && abs(t-time[n])>100*.Machine$double.eps) return(cuminc[n])
  return( cuminc[sum(time<=t | abs(time-t)<100*.Machine$double.eps)] )
}

# RMT is an integral of CIF:
.get.RMT <- function(t, time, cuminc)	# assumes that t is a scalar
{
  n <- length(time)
  if (t<time[1]) return(0)
  else
  {
    ind <- (1:n)[time <=t | abs(time-t)<100*.Machine$double.eps]
    time.dif <- diff(c(time[ind], t))
    return(sum(cuminc[ind]*time.dif))
  }
}

.base.haz.std <- function(bh)  # returns (non-cumulative) hazard
{
  n = nrow(bh)
  # bh$hazard is a cumulative hazard
  hazdiff = c(bh$hazard[1], apply(cbind(bh$hazard[1:(n-1)],bh$hazard[2:n]),1, diff))
  data.frame(time=bh$time, haz=hazdiff)
}

#' Fitting a logistic regression model for propensity scores and estimating weights
#'
#' Fits a propensity scores model by logistic regression and returns both estimated propensity scores
#' and requested weights. The estimated propensity scores can be used
#' for further diagnostics, e.g. for testing a positivity assumption and covariate balance.
#'
#' @param df a data frame that includes a treatment indicator \code{A}  and covariates \code{C}.
#' @param A a character specifying the name of the treatment/exposure variable.
#' It is assumed that \code{A} is a numeric binary indicator with 0/1 values, where \code{A}=1
#' is assumed a treatment group, and \code{A}=0 a control group.
#' @param C a vector of character strings with variable names (potential confounders)
#' in the logistic regression model for Propensity Scores, i.e. P(A=1|C=c).
#' The default value of \code{C} is NULL corresponding to \code{wtype}="unadj"
#' that will estimate treatment effects in the raw (observed) data.
#' @param wtype a character string variable indicating the type of weights that will define the target
#' population for which the ATE will be estimated.
#' The default is "stab.ATE" defined as P(A=a)/P(A=a|C=c) - see Hernan et al. (2000).
#' Other possible values are "ATE", "ATT", "ATC", and "overlap".
#' See Table 1 from Li, Morgan, and Zaslavsky (2018).
#' @param case.w a vector of case weights.
#'
#' @return  A list with the following fields:
#' \itemize{
#' \item{\code{wtype}} a character string indicating the type of the estimated weights
#' \item{\code{ps}} a vector of estimated propensity scores P(A=a|C=c)
#' \item{\code{w}} a vector of estimated weights
#' \item{\code{summary.glm}} a summary of the logistic regression fit which is done
#'  using \code{stats::glm}} function
#'
#' @references F. Li, K.L. Morgan, and A.M. Zaslavsky. 2018. Balancing Covariates via Propensity Score Weighting. Journal of the American Statistical Association 113 (521): 390–400.
#' @references M. Hernan, K.L. Morgan, and A.M. Zaslavsky. 2000. Balancing Covariates via Propensity Score Weighting. Journal of the American Statistical Association 113 (521): 390–400.
#'
#' @seealso \code{\link{fit.nonpar}}, \code{\link{fit.cox}}, \code{\link{causalCmprsk}}
#' @examples
#' # please see our package vignette for practical examples
#'
#' @export
get.weights <- function(df, A, C, wtype="stab.ATE", case.w=NULL)
{
  form.txt <- paste(A, " ~ ", sep="")
  p <- length(C) - 1
  for (i in 1:p)
    form.txt <- paste(form.txt, C[i], "+", sep="")
  form.txt <- paste(form.txt, C[p+1], sep="")
  form <- as.formula(form.txt)

  wei <- NULL
  trt <- df[[A]]

  if (is.null(case.w))
    case.w=rep(1, length(trt))

  suppressWarnings(p_logistic.fit <- glm(form, family = binomial(link = "logit"), data = df, weights=case.w) )
  # computing the propensity scores:
  pscore <- predict(p_logistic.fit, type = "response")

  pr.1 <- sum(trt*case.w)/sum(case.w)
  if (wtype=="stab.ATE")
    wei <- pr.1*trt/pscore + (1-trt)*(1-pr.1)/(1-pscore)
  if (wtype=="ATE")
    wei <- trt/pscore + (1-trt)/(1-pscore)
  if (wtype=="ATT")
    wei <- trt + (1-trt)*pscore/(1-pscore)
  if (wtype=="ATC")
    wei <- trt*(1-pscore)/pscore + (1-trt)
  if (wtype=="overlap")
    wei <- trt*(1-pscore) + (1-trt)*pscore

  list(wtype=wtype, ps=pscore, w=wei, summary.glm=summary(p_logistic.fit))
}


#' Number-at-risk in raw and weighted data
#'
#'
#' Obtaining time-varying number-at-risk statistic in both raw and weighted data
#'
#' @param df a data frame that includes time-to-event \code{T}, type of event \code{E},
#' a treatment indicator \code{A}  and covariates \code{C}.
#' @param T a character string specifying the name of the time-to-event variable in \code{df}.
#' @param E a character string specifying the name of the "event type" variable in \code{df}.
#' @param A a character specifying the name of the treatment/exposure variable.
#' It is assumed that \code{A} is a numeric binary indicator with 0/1 values, where \code{A}=1
#' is assumed a treatment group, and \code{A}=0 a control group.
#' @param C a vector of character strings with variable names (potential confounders)
#' in the logistic regression model for Propensity Scores, i.e. P(A=1|C=c).
#' The default value of \code{C} is NULL corresponding to \code{wtype}="unadj"
#' that will estimate treatment effects in the raw (observed) data.
#' @param wtype a character string variable indicating the type of weights that will define the target
#' population for which the ATE will be estimated.
#' The default is "unadj" - this will not adjust for possible
#' treatment selection bias and will not use propensity scores weighting. It can be used, for example,
#' in data from a randized controlled trial (RCT) where there is no need for emulation of baseline randomization.
#' Other possible values are "stab.ATE", "ATE", "ATT", "ATC" and "overlap".
#' See Table 1 from Li, Morgan, and Zaslavsky (2018).
#' "stab.ATE" is defined as P(A=a)/P(A=a|C=c) - see Hernan et al. (2000).
#' @param cens an integer value in \code{E} that corresponds to censoring times recorded in \code{T}.
#' By default \code{fit.nonpar} assumes \code{cens}=0
#'
#' @return  A list with two fields:
#' \itemize{
#' \item{\code{trt.0}} a matrix with three columns, \code{time}, \code{num} and \code{sample}
#' corresponding to the treatment arm with \code{A}=0.
#' The results for both weighted and unadjusted number-at-risk are returnd in  a long-format matrix.
#' The column \code{time} is a vector of time points at which we calculate the number-at-risk.
#' The column \code{num} is the number-at-risk.
#' The column \code{sample} is a factor variable that gets one of two values, "Weighted" or "Unadjusted".
#' The estimated number-at-risk in the weighted sample corresponds to the rows with \code{sample="Weighted"}.
#' \item{\code{trt.1}} a matrix with three columns, \code{time}, \code{num} and \code{sample}
#' corresponding to the treatment arm with \code{A}=1.
#' The results for both weighted and unadjusted number-at-risk are returnd in  a long-format matrix.
#' The column \code{time} is a vector of time points at which we calculate the number-at-risk.
#' The column \code{num} is the number-at-risk.
#' The column \code{sample} is a factor variable that gets one of two values, "Weighted" or "Unadjusted".
#' The estimated number-at-risk in the weighted sample corresponds to the rows with \code{sample="Weighted"}.}
#'
#' @seealso \code{\link{get.weights}}, \code{\link{get.pointEst}}, \code{\link{causalCmprsk}}
#'
#'
#' @export
get.numAtRisk <- function(df, T, E, A, C=NULL, wtype="unadj", cens=0)
{
  X <- df[[T]]
  E <- df[[E]]
  nobs <- length(X)
  trt <- df[[A]]

  res <- list()
  ps.fit <- get.weights(df=df, A=A, C=C, wtype=wtype)
  fit.w.1 <- survfit(Surv(time = X[trt==1], ev = as.numeric(E[trt==1]!=0) )~1, weights=ps.fit$w[trt==1])
  fit.unw.1 <- survfit(Surv(time = X[trt==1], ev = as.numeric(E[trt==1]!=0) )~1)
  fit.w.0 <- survfit(Surv(time = X[trt==0], ev = as.numeric(E[trt==0]!=0) )~1, weights=ps.fit$w[trt==0])
  fit.unw.0 <- survfit(Surv(time = X[trt==0], ev = as.numeric(E[trt==0]!=0) )~1)

  # trt.1
  num.at.risk.1 <- rbind(
    data.frame(time=fit.w.1$time, num=fit.w.1$n.risk, sample=1),
    data.frame(time=fit.unw.1$time, num=fit.unw.1$n.risk, sample=2)
  )

  num.at.risk.1$sample <- factor(num.at.risk.1$sample)
  levels(num.at.risk.1$sample) <- c("Weighted", "Unadjusted")

  # trt.0
  num.at.risk.0 <- rbind(
    data.frame(time=fit.w.0$time, num=fit.w.0$n.risk, sample=1),
    data.frame(time=fit.unw.0$time, num=fit.unw.0$n.risk, sample=2)
  )

  num.at.risk.0$sample <- factor(num.at.risk.0$sample)
  levels(num.at.risk.0$sample) <- c("Weighted", "Unadjusted")

  res$trt.1 <- num.at.risk.1
  res$trt.0 <-  num.at.risk.0

  return(res)
}


# time is a vector of arbitrary time points
.estimate.nonpar <- function(X, E, case.w, cens, time, E.set) # within one counterfactual world
{
  # overall survival:
  all.ev <- as.numeric(E!=cens) # an indicator of any event
  fit.all <- coxph(Surv(time=X, event=all.ev)~1, weights=case.w)
  bh.all <- basehaz(fit.all)
  OS <- exp(-bh.all$haz)

  res.a <- list() # a list of results; time is common for all the events
  # event-specific quantities:
  for (k in E.set)
  {
    Ek <- as.numeric(E==k)

    fit.k <- coxph(Surv(time = X, ev = Ek)~1, weights=case.w)
    cum.bh.k <- basehaz(fit.k, centered=FALSE)
    bh.k <- .base.haz.std(cum.bh.k)
    cif.k <- cumsum(bh.k$haz * sapply(bh.k$time, .get.S, bh.all$time, OS))

    CumHaz <- sapply(time, .get.CIF, cum.bh.k$time, cum.bh.k$haz)
    CIF.k <- sapply(time, .get.CIF, time=bh.k$time, cuminc=cif.k)
    RMT.k <-  sapply(time, .get.RMT, time = bh.k$time, cuminc=cif.k)

    res.a[[paste("Ev=", k, sep="")]] <- list(CumHaz=CumHaz, CIF=CIF.k, RMT=RMT.k)
  }
  res.a
} # end of .estimate.nonpar


.sequential.fit.nonpar <- function(df, T, E, A, C, wtype, cens, conf.level, bs, nbs.rep, seed)
{
  X <- df[[T]]
  E <- df[[E]]
  nobs <- length(X)
  trt <- df[[A]]
  time <- sort(unique(X[E!=cens]))
  E.set <- sort(unique(E))
  E.set <- E.set[E.set!=cens]

  res <- list(time=time) # this is the only place where I'll keep the time
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

  if (bs)
  {
    bs_seeds <- (1:nbs.rep) + seed

    # allocate memory for bs results:
    ntime <- length(res$time)
    bs.CumHaz <- bs.CIF <- bs.RMT <- list()
    for (k in E.set)
    {
      bs.CumHaz$trt.0[[paste("Ev=", k, sep="")]] <- bs.CumHaz$trt.1[[paste("Ev=", k, sep="")]] <-
        bs.CIF$trt.0[[paste("Ev=", k, sep="")]] <- bs.CIF$trt.1[[paste("Ev=", k, sep="")]] <-
        bs.RMT$trt.0[[paste("Ev=", k, sep="")]] <- bs.RMT$trt.1[[paste("Ev=", k, sep="")]] <-
        bs.CumHaz$logRatio[[paste("Ev=", k, sep="")]] <-
        bs.CIF$RD[[paste("Ev=", k, sep="")]] <- bs.CIF$RR[[paste("Ev=", k, sep="")]] <-
        bs.RMT$ATE[[paste("Ev=", k, sep="")]] <- matrix(nrow=nbs.rep, ncol=ntime)
    }

    pb <- txtProgressBar(min = 1, max = nbs.rep, style=3)
    # sequential bootstrap:
    for (i in 1:nbs.rep)
    {
      set.seed(bs_seeds[i])
      bs.w <- pmin(rexp(nobs,1), 5) # nobs = our sample size
      bs.w <- bs.w/mean(bs.w)

      if (wtype!="unadj")
      {
        bs.ps.fit <- get.weights(df, A, C, wtype, case.w = bs.w)
        bs.est.0 <- .estimate.nonpar(X=X[trt==0], E=E[trt==0],
                                     case.w=bs.w[trt==0]*bs.ps.fit$w[trt==0],
                                     cens=cens, time=time, E.set)
        bs.est.1 <- .estimate.nonpar(X=X[trt==1], E=E[trt==1],
                                     case.w=bs.w[trt==1]*bs.ps.fit$w[trt==1],
                                     cens=cens, time=time, E.set)
      }
      else
      {
        bs.est.0 <- .estimate.nonpar(X=X[trt==0], E=E[trt==0],
                                     case.w=bs.w[trt==0],
                                     cens=cens, time=time, E.set)
        bs.est.1 <- .estimate.nonpar(X=X[trt==1], E=E[trt==1],
                                     case.w=bs.w[trt==1],
                                     cens=cens, time=time, E.set)
      }

      # accumulate the results
      for (k in E.set)
      { # CumHaz
        bs.b <- log(bs.est.1[[paste("Ev=", k, sep="")]]$CumHaz) -
          log(bs.est.0[[paste("Ev=", k, sep="")]]$CumHaz) # calculate log(CumHaz1/CumHaz0)
        bs.CumHaz$trt.0[[paste("Ev=", k, sep="")]][i,] <- bs.est.0[[paste("Ev=", k, sep="")]]$CumHaz
        bs.CumHaz$trt.1[[paste("Ev=", k, sep="")]][i,] <- bs.est.1[[paste("Ev=", k, sep="")]]$CumHaz
        bs.CumHaz$logRatio[[paste("Ev=", k, sep="")]][i,] <- bs.b
        # CIF
        bs.CIF$trt.0[[paste("Ev=", k, sep="")]][i,] <- bs.est.0[[paste("Ev=", k, sep="")]]$CIF
        bs.CIF$trt.1[[paste("Ev=", k, sep="")]][i,] <- bs.est.1[[paste("Ev=", k, sep="")]]$CIF
        bs.CIF$RD[[paste("Ev=", k, sep="")]][i,] <- bs.est.1[[paste("Ev=", k, sep="")]]$CIF -
          bs.est.0[[paste("Ev=", k, sep="")]]$CIF
        bs.CIF$RR[[paste("Ev=", k, sep="")]][i,] <- bs.est.1[[paste("Ev=", k, sep="")]]$CIF/
          bs.est.0[[paste("Ev=", k, sep="")]]$CIF
        # RMT
        bs.RMT$trt.0[[paste("Ev=", k, sep="")]][i,] <- bs.est.0[[paste("Ev=", k, sep="")]]$RMT
        bs.RMT$trt.1[[paste("Ev=", k, sep="")]][i,] <- bs.est.1[[paste("Ev=", k, sep="")]]$RMT
        bs.RMT$ATE[[paste("Ev=", k, sep="")]][i,] <- bs.est.1[[paste("Ev=", k, sep="")]]$RMT -
          bs.est.0[[paste("Ev=", k, sep="")]]$RMT
      }
      #      if(i %% 25 == 0)
      #        cat("We are on replication ",i," of ",nbs.rep," replications.\n")
      setTxtProgressBar(pb, i)
    }
    close(pb)

    # summarize bs replications and save the results in 'res' object:
    # res is a list with 4 fields:
    # 1.time - a vector of times for which everything is estimated
    # 2.trt.0[[paste("Ev=", k, sep="")]]: CumHaz, CIF, RMT
    # 3.trt.1[[paste("Ev=", k, sep="")]]: CumHaz, CIF, RMT
    # 4. trt.eff[[paste("Ev=", k, sep="")]]: log.CumHazR, RD, RR, ATE.RMT
    alpha = 1-conf.level
    for (k in E.set)
    {
      # A. Cumulative Hazards::::::::::::::::::::::::::::::::
      # CI:
      res$trt.0[[paste("Ev=", k, sep="")]][["CumHaz.CI.L"]] <-
        apply(bs.CumHaz$trt.0[[paste("Ev=", k, sep="")]], 2, quantile, prob=alpha/2, na.rm=TRUE)
      res$trt.1[[paste("Ev=", k, sep="")]][["CumHaz.CI.L"]] <-
        apply(bs.CumHaz$trt.1[[paste("Ev=", k, sep="")]], 2, quantile, prob=alpha/2, na.rm=TRUE)
      res$trt.0[[paste("Ev=", k, sep="")]][["CumHaz.CI.U"]] <-
        apply(bs.CumHaz$trt.0[[paste("Ev=", k, sep="")]], 2, quantile, prob=1-alpha/2, na.rm=TRUE)
      res$trt.1[[paste("Ev=", k, sep="")]][["CumHaz.CI.U"]] <-
        apply(bs.CumHaz$trt.1[[paste("Ev=", k, sep="")]], 2, quantile, prob=1-alpha/2, na.rm=TRUE)
      res$trt.eff[[paste("Ev=", k, sep="")]][["log.CumHazR.CI.L"]] <-
        apply(bs.CumHaz$logRatio[[paste("Ev=", k, sep="")]], 2, quantile, prob=alpha/2, na.rm=TRUE)
      res$trt.eff[[paste("Ev=", k, sep="")]][["log.CumHazR.CI.U"]] <-
        apply(bs.CumHaz$logRatio[[paste("Ev=", k, sep="")]], 2, quantile, prob=1-alpha/2, na.rm=TRUE)
      # SE:
      res$trt.0[[paste("Ev=", k, sep="")]][["CumHaz.SE"]] <-
        apply(bs.CumHaz$trt.0[[paste("Ev=", k, sep="")]], 2, sd, na.rm=TRUE)
      res$trt.1[[paste("Ev=", k, sep="")]][["CumHaz.SE"]] <-
        apply(bs.CumHaz$trt.1[[paste("Ev=", k, sep="")]], 2, sd, na.rm=TRUE)
      res$trt.eff[[paste("Ev=", k, sep="")]][["log.CumHazR.SE"]] <-
        apply(bs.CumHaz$logRatio[[paste("Ev=", k, sep="")]], 2, sd, na.rm=TRUE)
      # bs.avg:
      res$trt.0[[paste("Ev=", k, sep="")]][["CumHaz.bs.avg"]] <-
        apply(bs.CumHaz$trt.0[[paste("Ev=", k, sep="")]], 2, mean, na.rm=TRUE)
      res$trt.1[[paste("Ev=", k, sep="")]][["CumHaz.bs.avg"]] <-
        apply(bs.CumHaz$trt.1[[paste("Ev=", k, sep="")]], 2, mean, na.rm=TRUE)
      res$trt.eff[[paste("Ev=", k, sep="")]][["log.CumHazR.bs.avg"]] <-
        apply(bs.CumHaz$logRatio[[paste("Ev=", k, sep="")]], 2, mean, na.rm=TRUE)

      # B. CIFs :::::::::::::::::::::::::::::::::::::::::::::::
      # CI:
      res$trt.0[[paste("Ev=", k, sep="")]][["CIF.CI.L"]] <-
        apply(bs.CIF$trt.0[[paste("Ev=", k, sep="")]], 2, quantile, prob=alpha/2, na.rm=TRUE)
      res$trt.1[[paste("Ev=", k, sep="")]][["CIF.CI.L"]] <-
        apply(bs.CIF$trt.1[[paste("Ev=", k, sep="")]], 2, quantile, prob=alpha/2, na.rm=TRUE)
      res$trt.0[[paste("Ev=", k, sep="")]][["CIF.CI.U"]] <-
        apply(bs.CIF$trt.0[[paste("Ev=", k, sep="")]], 2, quantile, prob=1-alpha/2, na.rm=TRUE)
      res$trt.1[[paste("Ev=", k, sep="")]][["CIF.CI.U"]] <-
        apply(bs.CIF$trt.1[[paste("Ev=", k, sep="")]], 2, quantile, prob=1-alpha/2, na.rm=TRUE)
      res$trt.eff[[paste("Ev=", k, sep="")]][["RD.CI.L"]] <-
        apply(bs.CIF$RD[[paste("Ev=", k, sep="")]], 2, quantile, prob=alpha/2, na.rm=TRUE)
      res$trt.eff[[paste("Ev=", k, sep="")]][["RD.CI.U"]] <-
        apply(bs.CIF$RD[[paste("Ev=", k, sep="")]], 2, quantile, prob=1-alpha/2, na.rm=TRUE)
      res$trt.eff[[paste("Ev=", k, sep="")]][["RR.CI.L"]] <-
        apply(bs.CIF$RR[[paste("Ev=", k, sep="")]], 2, quantile, prob=alpha/2, na.rm=TRUE)
      res$trt.eff[[paste("Ev=", k, sep="")]][["RR.CI.U"]] <-
        apply(bs.CIF$RR[[paste("Ev=", k, sep="")]], 2, quantile, prob=1-alpha/2, na.rm=TRUE)

      # SE:
      res$trt.0[[paste("Ev=", k, sep="")]][["CIF.SE"]] <-
        apply(bs.CIF$trt.0[[paste("Ev=", k, sep="")]], 2, sd, na.rm=TRUE)
      res$trt.1[[paste("Ev=", k, sep="")]][["CIF.SE"]] <-
        apply(bs.CIF$trt.1[[paste("Ev=", k, sep="")]], 2, sd, na.rm=TRUE)
      res$trt.eff[[paste("Ev=", k, sep="")]][["RD.SE"]] <-
        apply(bs.CIF$RD[[paste("Ev=", k, sep="")]], 2, sd, na.rm=TRUE)
      res$trt.eff[[paste("Ev=", k, sep="")]][["RR.SE"]] <-
        apply(bs.CIF$RR[[paste("Ev=", k, sep="")]], 2, sd, na.rm=TRUE)

      # bs.avg:
      res$trt.0[[paste("Ev=", k, sep="")]][["CIF.bs.avg"]] <-
        apply(bs.CIF$trt.0[[paste("Ev=", k, sep="")]], 2, mean, na.rm=TRUE)
      res$trt.1[[paste("Ev=", k, sep="")]][["CIF.bs.avg"]] <-
        apply(bs.CIF$trt.1[[paste("Ev=", k, sep="")]], 2, mean, na.rm=TRUE)
      res$trt.eff[[paste("Ev=", k, sep="")]][["RD.bs.avg"]] <-
        apply(bs.CIF$RD[[paste("Ev=", k, sep="")]], 2, mean, na.rm=TRUE)
      res$trt.eff[[paste("Ev=", k, sep="")]][["RR.bs.avg"]] <-
        apply(bs.CIF$RR[[paste("Ev=", k, sep="")]], 2, mean, na.rm=TRUE)

      # C. RMTs ::::::::::::::::::::::::::::::::::::::::::::::
      # CI:
      res$trt.0[[paste("Ev=", k, sep="")]][["RMT.CI.L"]] <-
        apply(bs.RMT$trt.0[[paste("Ev=", k, sep="")]], 2, quantile, prob=alpha/2, na.rm=TRUE)
      res$trt.1[[paste("Ev=", k, sep="")]][["RMT.CI.L"]] <-
        apply(bs.RMT$trt.1[[paste("Ev=", k, sep="")]], 2, quantile, prob=alpha/2, na.rm=TRUE)
      res$trt.0[[paste("Ev=", k, sep="")]][["RMT.CI.U"]] <-
        apply(bs.RMT$trt.0[[paste("Ev=", k, sep="")]], 2, quantile, prob=1-alpha/2, na.rm=TRUE)
      res$trt.1[[paste("Ev=", k, sep="")]][["RMT.CI.U"]] <-
        apply(bs.RMT$trt.1[[paste("Ev=", k, sep="")]], 2, quantile, prob=1-alpha/2, na.rm=TRUE)
      res$trt.eff[[paste("Ev=", k, sep="")]][["ATE.RMT.CI.L"]] <-
        apply(bs.RMT$ATE[[paste("Ev=", k, sep="")]], 2, quantile, prob=alpha/2, na.rm=TRUE)
      res$trt.eff[[paste("Ev=", k, sep="")]][["ATE.RMT.CI.U"]] <-
        apply(bs.RMT$ATE[[paste("Ev=", k, sep="")]], 2, quantile, prob=1-alpha/2, na.rm=TRUE)
      # SE:
      res$trt.0[[paste("Ev=", k, sep="")]][["RMT.SE"]] <-
        apply(bs.RMT$trt.0[[paste("Ev=", k, sep="")]], 2, sd, na.rm=TRUE)
      res$trt.1[[paste("Ev=", k, sep="")]][["RMT.SE"]] <-
        apply(bs.RMT$trt.1[[paste("Ev=", k, sep="")]], 2, sd, na.rm=TRUE)
      res$trt.eff[[paste("Ev=", k, sep="")]][["ATE.RMT.SE"]] <-
        apply(bs.RMT$ATE[[paste("Ev=", k, sep="")]], 2, sd, na.rm=TRUE)
      # bs.avg:
      res$trt.0[[paste("Ev=", k, sep="")]][["RMT.bs.avg"]] <-
        apply(bs.RMT$trt.0[[paste("Ev=", k, sep="")]], 2, mean, na.rm=TRUE)
      res$trt.1[[paste("Ev=", k, sep="")]][["RMT.bs.avg"]] <-
        apply(bs.RMT$trt.1[[paste("Ev=", k, sep="")]], 2, mean, na.rm=TRUE)
      res$trt.eff[[paste("Ev=", k, sep="")]][["ATE.RMT.bs.avg"]] <-
        apply(bs.RMT$ATE[[paste("Ev=", k, sep="")]], 2, mean, na.rm=TRUE)
    }
  }
  class(res) <- "cmprsk"
  return(res)
}

# assumes PH across 2 counterfactual worlds...
.estimate.cox <- function(X, E, trt, case.w, cens, time, E.set)
{
  res.a <- list() # a list of results; time is common for all the events
  # estimate and save base.Haz.k and HR.k:
  OS.1 <- OS.0 <- 0
  for (k in E.set)
  {
    Ek <- as.numeric(E==k)
    fit.cox.k <- coxph(Surv(time = X, event = Ek) ~ trt , weights = case.w)
    BaseCumHaz.k <- basehaz(fit.cox.k, centered = FALSE)
    res.a[[paste("Ev=", k, sep="")]] <- list(logHR=fit.cox.k$coef,
                                             bh=BaseCumHaz.k)
    OS.0 <- OS.0 - sapply(time, .get.CIF, BaseCumHaz.k$time, BaseCumHaz.k$haz)
    OS.1 <- OS.1 - sapply(time, .get.CIF, BaseCumHaz.k$time, BaseCumHaz.k$haz * exp(fit.cox.k$coef))
  }
  # OS=event-free survival in each of the trt arms at time=time
  res.a$OS <- list(time=time, OS.trt.0=exp(OS.0), OS.trt.1=exp(OS.1))
  res.a
} # end of .estimate.cox



.sequential.fit.cox <- function(df, T, E, A, C, wtype, cens, conf.level, bs, nbs.rep, seed)
{
  X <- df[[T]]
  E <- df[[E]]
  nobs <- length(X)
  trt <- df[[A]]
  time <- sort(unique(X[E!=cens]))
  E.set <- sort(unique(E))
  E.set <- E.set[E.set!=cens]

  res <- list(time=time)
  if (wtype!="unadj")
  {
    ps.fit <- get.weights(df, A, C, wtype)
    est <- .estimate.cox(X=X, E=E, trt=trt, case.w =ps.fit$w, cens=cens, time=time, E.set=E.set)
  }
  else
    est <- .estimate.cox(X=X, E=E, trt=trt, case.w=rep(1, nobs), cens=cens, time=time, E.set=E.set)

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

  if (bs)
  {
    bs_seeds <- (1:nbs.rep) + seed

    # allocate memory for bs results:
    ntime <- length(res$time)
    bs.CumHaz <- bs.CIF <- bs.RMT <- list()
    for (k in E.set)
    {
      bs.CumHaz$trt.0[[paste("Ev=", k, sep="")]] <- bs.CumHaz$trt.1[[paste("Ev=", k, sep="")]] <-
        bs.CIF$trt.0[[paste("Ev=", k, sep="")]] <- bs.CIF$trt.1[[paste("Ev=", k, sep="")]] <-
        bs.RMT$trt.0[[paste("Ev=", k, sep="")]] <- bs.RMT$trt.1[[paste("Ev=", k, sep="")]] <-
        bs.CIF$RD[[paste("Ev=", k, sep="")]] <- bs.CIF$RR[[paste("Ev=", k, sep="")]] <-
        bs.RMT$ATE[[paste("Ev=", k, sep="")]] <- matrix(nrow=nbs.rep, ncol=ntime)
      bs.CumHaz$logRatio[[paste("Ev=", k, sep="")]] <- vector("double", length=nbs.rep)
    }

    pb <- txtProgressBar(min = 1, max = nbs.rep, style=3)
    # sequential bootstrap:
    for (i in 1:nbs.rep)
    {
      set.seed(bs_seeds[i])
      bs.w <- pmin(rexp(nobs,1), 5) # nobs = our sample size
      bs.w <- bs.w/mean(bs.w)
      if (wtype!="unadj")
      {
        bs.ps.fit <- get.weights(df, A, C, wtype, case.w = bs.w)
        bs.est <- .estimate.cox(X=X, E=E, trt=trt, case.w =bs.w*bs.ps.fit$w, cens=cens, time=time, E.set=E.set)
      }
      else
        bs.est <- .estimate.cox(X=X, E=E, trt=trt, case.w =bs.w, cens=cens, time=time, E.set=E.set)

      # accumulate the results
      for (k in E.set)
      {
        # CumHaz:
        bs.time.k <- est[[paste("Ev=", k, sep="")]]$bh$time
        bs.CumHaz$trt.0[[paste("Ev=", k, sep="")]][i,] <-
          sapply(time, .get.CIF, bs.time.k,
                 bs.est[[paste("Ev=", k, sep="")]]$bh$haz)
        bs.CumHaz$trt.1[[paste("Ev=", k, sep="")]][i,] <-
          sapply(time, .get.CIF, bs.time.k,
                 bs.est[[paste("Ev=", k, sep="")]]$bh$haz * exp(bs.est[[paste("Ev=", k, sep="")]]$logHR))
        bs.CumHaz$logRatio[[paste("Ev=", k, sep="")]][i] <- bs.est[[paste("Ev=", k, sep="")]]$logHR

        # CIF
        bs.bh.k <- .base.haz.std(bs.est[[paste("Ev=", k, sep="")]]$bh)
        bs.cif.0.k <- cumsum(sapply(bs.bh.k$time, .get.S, bs.est$OS$time, bs.est$OS$OS.trt.0) * bs.bh.k$haz)
        bs.cif.1.k <- cumsum(sapply(bs.bh.k$time, .get.S, bs.est$OS$time, bs.est$OS$OS.trt.1) * bs.bh.k$haz *
                               exp(bs.est[[paste("Ev=", k, sep="")]]$logHR))
        bs.CIF$trt.0[[paste("Ev=", k, sep="")]][i,] <- sapply(time, .get.CIF, time=bs.bh.k$time,
                                                              cuminc=bs.cif.0.k)
        bs.CIF$trt.1[[paste("Ev=", k, sep="")]][i,] <- sapply(time, .get.CIF, time=bs.bh.k$time,
                                                              cuminc=bs.cif.1.k)
        bs.CIF$RD[[paste("Ev=", k, sep="")]][i,] <- bs.CIF$trt.1[[paste("Ev=", k, sep="")]][i,] -
          bs.CIF$trt.0[[paste("Ev=", k, sep="")]][i,]
        bs.CIF$RR[[paste("Ev=", k, sep="")]][i,] <- bs.CIF$trt.1[[paste("Ev=", k, sep="")]][i,] /
          bs.CIF$trt.0[[paste("Ev=", k, sep="")]][i,]

        # RMT
        bs.RMT$trt.0[[paste("Ev=", k, sep="")]][i,] <- sapply(time, .get.RMT, time=bs.bh.k$time,
                                                              cuminc=bs.cif.0.k)
        bs.RMT$trt.1[[paste("Ev=", k, sep="")]][i,] <- sapply(time, .get.RMT, time=bs.bh.k$time,
                                                              cuminc=bs.cif.1.k)
        bs.RMT$ATE[[paste("Ev=", k, sep="")]][i,] <- bs.RMT$trt.1[[paste("Ev=", k, sep="")]][i,] -
          bs.RMT$trt.0[[paste("Ev=", k, sep="")]][i,]
      }
      #     if (i %% 25 == 0)
      #       cat("We are on replication ",i," of ",nbs.rep," replications.\n")
      setTxtProgressBar(pb, i)
    }
    close(pb)

    # summarize bs replications and save the results in 'res' object:
    # res is a list with 4 fields:
    # 1.time - a vector of times for which everything is estimated
    # 2.trt.0[[paste("Ev=", k, sep="")]]: CumHaz, CIF, RMT
    # 3.trt.1[[paste("Ev=", k, sep="")]]: CumHaz, CIF, RMT
    # 4. trt.eff[[paste("Ev=", k, sep="")]]: log.CumHazR, RD, RR, ATE.RMT
    alpha = 1-conf.level
    for (k in E.set)
    {
      # A. Cumulative Hazards::::::::::::::::::::::::::::::::
      # CI:
      res$trt.0[[paste("Ev=", k, sep="")]][["CumHaz.CI.L"]] <-
        apply(bs.CumHaz$trt.0[[paste("Ev=", k, sep="")]], 2, quantile, prob=alpha/2, na.rm=TRUE)
      res$trt.1[[paste("Ev=", k, sep="")]][["CumHaz.CI.L"]] <-
        apply(bs.CumHaz$trt.1[[paste("Ev=", k, sep="")]], 2, quantile, prob=alpha/2, na.rm=TRUE)
      res$trt.0[[paste("Ev=", k, sep="")]][["CumHaz.CI.U"]] <-
        apply(bs.CumHaz$trt.0[[paste("Ev=", k, sep="")]], 2, quantile, prob=1-alpha/2, na.rm=TRUE)
      res$trt.1[[paste("Ev=", k, sep="")]][["CumHaz.CI.U"]] <-
        apply(bs.CumHaz$trt.1[[paste("Ev=", k, sep="")]], 2, quantile, prob=1-alpha/2, na.rm=TRUE)
      res$trt.eff[[paste("Ev=", k, sep="")]][["log.CumHazR.CI.L"]] <-
        quantile(bs.CumHaz$logRatio[[paste("Ev=", k, sep="")]], prob=alpha/2, na.rm=TRUE)
      res$trt.eff[[paste("Ev=", k, sep="")]][["log.CumHazR.CI.U"]] <-
        quantile(bs.CumHaz$logRatio[[paste("Ev=", k, sep="")]], prob=1-alpha/2, na.rm=TRUE)

      # SE:
      res$trt.0[[paste("Ev=", k, sep="")]][["CumHaz.SE"]] <-
        apply(bs.CumHaz$trt.0[[paste("Ev=", k, sep="")]], 2, sd, na.rm=TRUE)
      res$trt.1[[paste("Ev=", k, sep="")]][["CumHaz.SE"]] <-
        apply(bs.CumHaz$trt.1[[paste("Ev=", k, sep="")]], 2, sd, na.rm=TRUE)
      res$trt.eff[[paste("Ev=", k, sep="")]][["log.CumHazR.SE"]] <-
        sd(bs.CumHaz$logRatio[[paste("Ev=", k, sep="")]], na.rm=TRUE)

      # bs.avg:
      res$trt.0[[paste("Ev=", k, sep="")]][["CumHaz.bs.avg"]] <-
        apply(bs.CumHaz$trt.0[[paste("Ev=", k, sep="")]], 2, mean, na.rm=TRUE)
      res$trt.1[[paste("Ev=", k, sep="")]][["CumHaz.bs.avg"]] <-
        apply(bs.CumHaz$trt.1[[paste("Ev=", k, sep="")]], 2, mean, na.rm=TRUE)
      res$trt.eff[[paste("Ev=", k, sep="")]][["log.CumHazR.bs.avg"]] <-
        mean(bs.CumHaz$logRatio[[paste("Ev=", k, sep="")]], na.rm=TRUE)

      # B. CIFs :::::::::::::::::::::::::::::::::::::::::::::::
      # CI:
      res$trt.0[[paste("Ev=", k, sep="")]][["CIF.CI.L"]] <-
        apply(bs.CIF$trt.0[[paste("Ev=", k, sep="")]], 2, quantile, prob=alpha/2, na.rm=TRUE)
      res$trt.1[[paste("Ev=", k, sep="")]][["CIF.CI.L"]] <-
        apply(bs.CIF$trt.1[[paste("Ev=", k, sep="")]], 2, quantile, prob=alpha/2, na.rm=TRUE)
      res$trt.0[[paste("Ev=", k, sep="")]][["CIF.CI.U"]] <-
        apply(bs.CIF$trt.0[[paste("Ev=", k, sep="")]], 2, quantile, prob=1-alpha/2, na.rm=TRUE)
      res$trt.1[[paste("Ev=", k, sep="")]][["CIF.CI.U"]] <-
        apply(bs.CIF$trt.1[[paste("Ev=", k, sep="")]], 2, quantile, prob=1-alpha/2, na.rm=TRUE)
      res$trt.eff[[paste("Ev=", k, sep="")]][["RD.CI.L"]] <-
        apply(bs.CIF$RD[[paste("Ev=", k, sep="")]], 2, quantile, prob=alpha/2, na.rm=TRUE)
      res$trt.eff[[paste("Ev=", k, sep="")]][["RD.CI.U"]] <-
        apply(bs.CIF$RD[[paste("Ev=", k, sep="")]], 2, quantile, prob=1-alpha/2, na.rm=TRUE)
      res$trt.eff[[paste("Ev=", k, sep="")]][["RR.CI.L"]] <-
        apply(bs.CIF$RR[[paste("Ev=", k, sep="")]], 2, quantile, prob=alpha/2, na.rm=TRUE)
      res$trt.eff[[paste("Ev=", k, sep="")]][["RR.CI.U"]] <-
        apply(bs.CIF$RR[[paste("Ev=", k, sep="")]], 2, quantile, prob=1-alpha/2, na.rm=TRUE)

      # SE:
      res$trt.0[[paste("Ev=", k, sep="")]][["CIF.SE"]] <-
        apply(bs.CIF$trt.0[[paste("Ev=", k, sep="")]], 2, sd, na.rm=TRUE)
      res$trt.1[[paste("Ev=", k, sep="")]][["CIF.SE"]] <-
        apply(bs.CIF$trt.1[[paste("Ev=", k, sep="")]], 2, sd, na.rm=TRUE)
      res$trt.eff[[paste("Ev=", k, sep="")]][["RD.SE"]] <-
        apply(bs.CIF$RD[[paste("Ev=", k, sep="")]], 2, sd, na.rm=TRUE)
      res$trt.eff[[paste("Ev=", k, sep="")]][["RR.SE"]] <-
        apply(bs.CIF$RR[[paste("Ev=", k, sep="")]], 2, sd, na.rm=TRUE)

      # bs.avg:
      res$trt.0[[paste("Ev=", k, sep="")]][["CIF.bs.avg"]] <-
        apply(bs.CIF$trt.0[[paste("Ev=", k, sep="")]], 2, mean, na.rm=TRUE)
      res$trt.1[[paste("Ev=", k, sep="")]][["CIF.bs.avg"]] <-
        apply(bs.CIF$trt.1[[paste("Ev=", k, sep="")]], 2, mean, na.rm=TRUE)
      res$trt.eff[[paste("Ev=", k, sep="")]][["RD.bs.avg"]] <-
        apply(bs.CIF$RD[[paste("Ev=", k, sep="")]], 2, mean, na.rm=TRUE)
      res$trt.eff[[paste("Ev=", k, sep="")]][["RR.bs.avg"]] <-
        apply(bs.CIF$RR[[paste("Ev=", k, sep="")]], 2, mean, na.rm=TRUE)

      # C. RMTs ::::::::::::::::::::::::::::::::::::::::::::::
      # CI:
      res$trt.0[[paste("Ev=", k, sep="")]][["RMT.CI.L"]] <-
        apply(bs.RMT$trt.0[[paste("Ev=", k, sep="")]], 2, quantile, prob=alpha/2, na.rm=TRUE)
      res$trt.1[[paste("Ev=", k, sep="")]][["RMT.CI.L"]] <-
        apply(bs.RMT$trt.1[[paste("Ev=", k, sep="")]], 2, quantile, prob=alpha/2, na.rm=TRUE)
      res$trt.0[[paste("Ev=", k, sep="")]][["RMT.CI.U"]] <-
        apply(bs.RMT$trt.0[[paste("Ev=", k, sep="")]], 2, quantile, prob=1-alpha/2, na.rm=TRUE)
      res$trt.1[[paste("Ev=", k, sep="")]][["RMT.CI.U"]] <-
        apply(bs.RMT$trt.1[[paste("Ev=", k, sep="")]], 2, quantile, prob=1-alpha/2, na.rm=TRUE)
      res$trt.eff[[paste("Ev=", k, sep="")]][["ATE.RMT.CI.L"]] <-
        apply(bs.RMT$ATE[[paste("Ev=", k, sep="")]], 2, quantile, prob=alpha/2, na.rm=TRUE)
      res$trt.eff[[paste("Ev=", k, sep="")]][["ATE.RMT.CI.U"]] <-
        apply(bs.RMT$ATE[[paste("Ev=", k, sep="")]], 2, quantile, prob=1-alpha/2, na.rm=TRUE)
      # SE:
      res$trt.0[[paste("Ev=", k, sep="")]][["RMT.SE"]] <-
        apply(bs.RMT$trt.0[[paste("Ev=", k, sep="")]], 2, sd, na.rm=TRUE)
      res$trt.1[[paste("Ev=", k, sep="")]][["RMT.SE"]] <-
        apply(bs.RMT$trt.1[[paste("Ev=", k, sep="")]], 2, sd, na.rm=TRUE)
      res$trt.eff[[paste("Ev=", k, sep="")]][["ATE.RMT.SE"]] <-
        apply(bs.RMT$ATE[[paste("Ev=", k, sep="")]], 2, sd, na.rm=TRUE)
      # bs.avg:
      res$trt.0[[paste("Ev=", k, sep="")]][["RMT.bs.avg"]] <-
        apply(bs.RMT$trt.0[[paste("Ev=", k, sep="")]], 2, mean, na.rm=TRUE)
      res$trt.1[[paste("Ev=", k, sep="")]][["RMT.bs.avg"]] <-
        apply(bs.RMT$trt.1[[paste("Ev=", k, sep="")]], 2, mean, na.rm=TRUE)
      res$trt.eff[[paste("Ev=", k, sep="")]][["ATE.RMT.bs.avg"]] <-
        apply(bs.RMT$ATE[[paste("Ev=", k, sep="")]], 2, mean, na.rm=TRUE)
    }
  }
  class(res) <- "cmprsk"
  return(res)
}

#' Returns point estimates corresponding to a specific time point
#'
#' @param cmprsk.obj a "cmprsk" object returned by one of the function \code{fit.cox} or \code{fit.nonpar}
#' @param timepoint a scalar value of the time point of interest
#'
#' @return  A list with the following fields:
#' \tabular{ll}{
#' \code{time}  \tab \cr a scalar timepoint passed into the function \cr
#' \code{trt.0} \tab \cr a list of estimates of the absolute counterfactual parameters
#' corresponding to \code{A}=0 and the type of event \code{E}. \code{trt.0} has the number of
#'  fields as the number of different types of events in the data set.
#' For each type of event there is a list of estimates:  \cr }
#' \itemize{
#' \item{\code{CumHaz}} point cumulative hazard estimates
#' \item{\code{CIF}} point cumulative incidence function estimated
#' \item{\code{RMT}} point estimates of the restricted mean time}
#' \tabular{ll}{
#' \code{trt.1} \tab \cr a list of estimates of the absolute counterfactual parameters
#' corresponding to \code{A}=1 and the type of event \code{E}. \code{trt.1} has the number of
#'  fields as the number of different types of events  in the data set.
#' For each type of event there is a list of estimates:  \cr }
#' \itemize{
#' \item{\code{CumHaz}} point cumulative hazard estimates
#' \item{\code{CIF}} point cumulative incidence function estimated
#' \item{\code{RMT}} point estimates of the restricted mean time}
#' \tabular{ll}{
#' \code{trt.eff} \tab \cr a list of estimates of the treatment effect measures
#' corresponding to the type of event \code{E}. \code{trt.eff} has the number of
#'  fields as the number of different types of events  in the data set.
#' For each type of event there is a list of estimates: \cr}
#' \itemize{
#' \item{\code{log.CumHazR}} a point estimate of the log of the ratio of hazards in two treatment arms
#' \item{\code{RD}} a point estimate of the risk difference between two treatment arms
#' \item{\code{RR}} a point estimate of the risk ratio between two treatment arms
#' \item{\code{ATE.RMT}} a point estimate of the restricted mean time difference
#'  between two treatment arms}
#'
#' @seealso \code{\link{fit.cox}}, \code{\link{fit.nonpar}}, \code{\link{causalCmprsk}}
#' @examples
#' # please see our package vignette for practical examples
#'
#' @export
get.pointEst <- function(cmprsk.obj, timepoint) # assumes timepoint is a scalar
{
  if (class(cmprsk.obj)!="cmprsk")
  {
    cat("Error. 'cmprsk.obj' is not of 'cmprsk' type.")
    return(NULL)
  }
  if (length(timepoint)>1)
  {
    cat("Error. 'timepoint' is supposed to be a scalar.")
    return(NULL)
  }
  point.res <- list()
  point.res$time <- timepoint
  K <-  length(cmprsk.obj$trt.0) # number of events
  for (k in 1:K)
  {
    point.res$trt.0[[paste("Ev=", k, sep="")]] <- list()
    point.res$trt.1[[paste("Ev=", k, sep="")]] <- list()
    point.res$trt.eff[[paste("Ev=", k, sep="")]] <- list()
    #trt.0
    point.res$trt.0[[paste("Ev=", k, sep="")]]$CumHaz <- .get.CIF(timepoint, cmprsk.obj$time, cmprsk.obj$trt.0[[paste("Ev=", k, sep="")]]$CumHaz)
    point.res$trt.0[[paste("Ev=", k, sep="")]]$CumHaz.CI.L <- .get.CIF(timepoint, cmprsk.obj$time, cmprsk.obj$trt.0[[paste("Ev=", k, sep="")]]$CumHaz.CI.L)
    point.res$trt.0[[paste("Ev=", k, sep="")]]$CumHaz.CI.U <- .get.CIF(timepoint, cmprsk.obj$time, cmprsk.obj$trt.0[[paste("Ev=", k, sep="")]]$CumHaz.CI.U)

    point.res$trt.0[[paste("Ev=", k, sep="")]]$CIF <- .get.CIF(timepoint, cmprsk.obj$time, cmprsk.obj$trt.0[[paste("Ev=", k, sep="")]]$CIF)
    point.res$trt.0[[paste("Ev=", k, sep="")]]$CIF.CI.L <- .get.CIF(timepoint, cmprsk.obj$time, cmprsk.obj$trt.0[[paste("Ev=", k, sep="")]]$CIF.CI.L)
    point.res$trt.0[[paste("Ev=", k, sep="")]]$CIF.CI.U <- .get.CIF(timepoint, cmprsk.obj$time, cmprsk.obj$trt.0[[paste("Ev=", k, sep="")]]$CIF.CI.U)

    point.res$trt.0[[paste("Ev=", k, sep="")]]$RMT <- .get.CIF(timepoint, cmprsk.obj$time, cmprsk.obj$trt.0[[paste("Ev=", k, sep="")]]$RMT)
    point.res$trt.0[[paste("Ev=", k, sep="")]]$RMT.CI.L <- .get.CIF(timepoint, cmprsk.obj$time, cmprsk.obj$trt.0[[paste("Ev=", k, sep="")]]$RMT.CI.L)
    point.res$trt.0[[paste("Ev=", k, sep="")]]$RMT.CI.U <- .get.CIF(timepoint, cmprsk.obj$time, cmprsk.obj$trt.0[[paste("Ev=", k, sep="")]]$RMT.CI.U)

    # trt.1
    point.res$trt.1[[paste("Ev=", k, sep="")]]$CumHaz <- .get.CIF(timepoint, cmprsk.obj$time, cmprsk.obj$trt.1[[paste("Ev=", k, sep="")]]$CumHaz)
    point.res$trt.1[[paste("Ev=", k, sep="")]]$CumHaz.CI.L <- .get.CIF(timepoint, cmprsk.obj$time, cmprsk.obj$trt.1[[paste("Ev=", k, sep="")]]$CumHaz.CI.L)
    point.res$trt.1[[paste("Ev=", k, sep="")]]$CumHaz.CI.U <- .get.CIF(timepoint, cmprsk.obj$time, cmprsk.obj$trt.1[[paste("Ev=", k, sep="")]]$CumHaz.CI.U)

    point.res$trt.1[[paste("Ev=", k, sep="")]]$CIF <- .get.CIF(timepoint, cmprsk.obj$time, cmprsk.obj$trt.1[[paste("Ev=", k, sep="")]]$CIF)
    point.res$trt.1[[paste("Ev=", k, sep="")]]$CIF.CI.L <- .get.CIF(timepoint, cmprsk.obj$time, cmprsk.obj$trt.1[[paste("Ev=", k, sep="")]]$CIF.CI.L)
    point.res$trt.1[[paste("Ev=", k, sep="")]]$CIF.CI.U <- .get.CIF(timepoint, cmprsk.obj$time, cmprsk.obj$trt.1[[paste("Ev=", k, sep="")]]$CIF.CI.U)

    point.res$trt.1[[paste("Ev=", k, sep="")]]$RMT <- .get.CIF(timepoint, cmprsk.obj$time, cmprsk.obj$trt.1[[paste("Ev=", k, sep="")]]$RMT)
    point.res$trt.1[[paste("Ev=", k, sep="")]]$RMT.CI.L <- .get.CIF(timepoint, cmprsk.obj$time, cmprsk.obj$trt.1[[paste("Ev=", k, sep="")]]$RMT.CI.L)
    point.res$trt.1[[paste("Ev=", k, sep="")]]$RMT.CI.U <- .get.CIF(timepoint, cmprsk.obj$time, cmprsk.obj$trt.1[[paste("Ev=", k, sep="")]]$RMT.CI.U)

    # trt.eff
    if (length(cmprsk.obj$trt.eff[[paste("Ev=", k, sep="")]]$log.CumHazR)!=1)
    {
      point.res$trt.eff[[paste("Ev=", k, sep="")]]$log.CumHazR <- .get.CIF(timepoint, cmprsk.obj$time, cmprsk.obj$trt.eff[[paste("Ev=", k, sep="")]]$log.CumHazR)
      point.res$trt.eff[[paste("Ev=", k, sep="")]]$log.CumHazR.CI.L <- .get.CIF(timepoint, cmprsk.obj$time, cmprsk.obj$trt.eff[[paste("Ev=", k, sep="")]]$log.CumHazR.CI.L)
      point.res$trt.eff[[paste("Ev=", k, sep="")]]$log.CumHazR.CI.U <- .get.CIF(timepoint, cmprsk.obj$time, cmprsk.obj$trt.eff[[paste("Ev=", k, sep="")]]$log.CumHazR.CI.U)
    }
    else
    {
      point.res$trt.eff[[paste("Ev=", k, sep="")]]$CumHaz <- cmprsk.obj$trt.eff[[paste("Ev=", k, sep="")]]$log.CumHazR
      point.res$trt.eff[[paste("Ev=", k, sep="")]]$CumHaz.CI.L <- cmprsk.obj$trt.eff[[paste("Ev=", k, sep="")]]$log.CumHazR.CI.L
      point.res$trt.eff[[paste("Ev=", k, sep="")]]$CumHaz.CI.U <- cmprsk.obj$trt.eff[[paste("Ev=", k, sep="")]]$log.CumHazR.CI.U
    }

    point.res$trt.eff[[paste("Ev=", k, sep="")]]$RD <- .get.CIF(timepoint, cmprsk.obj$time, cmprsk.obj$trt.eff[[paste("Ev=", k, sep="")]]$RD)
    point.res$trt.eff[[paste("Ev=", k, sep="")]]$RD.CI.L <- .get.CIF(timepoint, cmprsk.obj$time, cmprsk.obj$trt.eff[[paste("Ev=", k, sep="")]]$RD.CI.L)
    point.res$trt.eff[[paste("Ev=", k, sep="")]]$RD.CI.U <- .get.CIF(timepoint, cmprsk.obj$time, cmprsk.obj$trt.eff[[paste("Ev=", k, sep="")]]$RD.CI.U)

    point.res$trt.eff[[paste("Ev=", k, sep="")]]$RR <- .get.CIF(timepoint, cmprsk.obj$time, cmprsk.obj$trt.eff[[paste("Ev=", k, sep="")]]$RR)
    point.res$trt.eff[[paste("Ev=", k, sep="")]]$RR.CI.L <- .get.CIF(timepoint, cmprsk.obj$time, cmprsk.obj$trt.eff[[paste("Ev=", k, sep="")]]$RR.CI.L)
    point.res$trt.eff[[paste("Ev=", k, sep="")]]$RR.CI.U <- .get.CIF(timepoint, cmprsk.obj$time, cmprsk.obj$trt.eff[[paste("Ev=", k, sep="")]]$RR.CI.U)

    point.res$trt.eff[[paste("Ev=", k, sep="")]]$ATE.RMT <- .get.CIF(timepoint, cmprsk.obj$time, cmprsk.obj$trt.eff[[paste("Ev=", k, sep="")]]$ATE.RMT)
    point.res$trt.eff[[paste("Ev=", k, sep="")]]$ATE.RMT.CI.L <- .get.CIF(timepoint, cmprsk.obj$time, cmprsk.obj$trt.eff[[paste("Ev=", k, sep="")]]$ATE.RMT.CI.L)
    point.res$trt.eff[[paste("Ev=", k, sep="")]]$ATE.RMT.CI.U <- .get.CIF(timepoint, cmprsk.obj$time, cmprsk.obj$trt.eff[[paste("Ev=", k, sep="")]]$ATE.RMT.CI.U)

  }
  return(point.res)
}

.nonpar.run <- function(df, T, E, A, C, wtype, cens, E.set,time,trt,nobs,X,case.w)
{
  if (wtype!="unadj")
  {
    ps.fit <- get.weights(df, A, C, wtype, case.w = case.w)
    est.0 <- .estimate.nonpar(X=X[trt==0], E=E[trt==0],
                              case.w = case.w[trt==0]*ps.fit$w[trt==0],
                              cens=cens, time=time, E.set)
    est.1 <- .estimate.nonpar(X=X[trt==1], E=E[trt==1],
                              case.w = case.w[trt==1]*ps.fit$w[trt==1],
                              cens=cens, time=time, E.set)
  }
  else
  {
    est.0 <- .estimate.nonpar(X=X[trt==0], E=E[trt==0],
                              case.w = case.w[trt==0], cens=cens, time=time, E.set)
    est.1 <- .estimate.nonpar(X=X[trt==1], E=E[trt==1],
                              case.w = case.w[trt==1], cens=cens, time=time, E.set)
  }
  res <- list(time=time)
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

  res <- list(time=time)
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

# Parallel Cox Fit
.parallel.fit.cox <- function(df, T, E, A, C, wtype, cens, conf.level, bs, nbs.rep, seed)
{
  X <- df[[T]]
  E <- df[[E]]
  nobs <- length(X)
  trt <- df[[A]]
  time <- sort(unique(X[E!=cens]))
  E.set <- sort(unique(E))
  E.set <- E.set[E.set!=cens]

  #  res <- list(time=time)
  res <- .cox.run(df, T, E, A, C, wtype, cens, E.set,time,trt,nobs,X, case.w = rep(1,nobs))
  # res is a list with 4 fields: time, trt.0, trt.0, trt.eff
  if (bs)
  {
    bs_seeds <- seq(1,nbs.rep,1) + seed
    # allocate memory for bs results:
    ntime <- length(res$time)
    bs.CumHaz <- bs.CIF <- bs.RMT <- list()
    for (k in E.set)
    {
      bs.CumHaz$trt.0[[paste("Ev=", k, sep="")]] <- bs.CumHaz$trt.1[[paste("Ev=", k, sep="")]] <-
        bs.CIF$trt.0[[paste("Ev=", k, sep="")]] <- bs.CIF$trt.1[[paste("Ev=", k, sep="")]] <-
        bs.RMT$trt.0[[paste("Ev=", k, sep="")]] <- bs.RMT$trt.1[[paste("Ev=", k, sep="")]] <-
        bs.CIF$RD[[paste("Ev=", k, sep="")]] <- bs.CIF$RR[[paste("Ev=", k, sep="")]] <-
        bs.RMT$ATE[[paste("Ev=", k, sep="")]] <- matrix(nrow=nbs.rep, ncol=ntime)
      bs.CumHaz$logRatio[[paste("Ev=", k, sep="")]] <- vector("double", length=nbs.rep)
    }
    i=0 # to remove the package NOTE
    pb <- txtProgressBar(min = 1, max = nbs.rep, style=3)
    bs_aggregates <- foreach(i = 1:nbs.rep, .export=c(".cox.run", "get.weights", ".estimate.cox", ".base.haz.std",
                                                      ".get.CIF",".get.S", ".get.RMT"),
                             .packages=c("survival")) %dopar%
      {
        set.seed(bs_seeds[i])
        bs.w <- pmin(rexp(nobs,1), 5) # nobs = our sample size
        bs.w <- bs.w/mean(bs.w)
        bs_aggregates <- .cox.run(df, T, E, A, C, wtype, cens, E.set,time,trt,nobs,X,case.w = bs.w)
        #        if(i %% 25 == 0)
        #         cat(paste0('We are on replication ',i,' of ',nbs.rep,' replications.'),file= stdout())
        setTxtProgressBar(pb, i)
        bs_aggregates
      }
    close(pb)
    # summarize bs replications and save the results in 'res' object:
    # res is a list with 4 fields:
    # 1.time - a vector of times for which everything is estimated
    # 2.trt.0[[paste("Ev=", k, sep="")]]: CumHaz, CIF, RMT
    # 3.trt.1[[paste("Ev=", k, sep="")]]: CumHaz, CIF, RMT
    # 4. trt.eff[[paste("Ev=", k, sep="")]]: log.CumHazR, RD, RR, ATE.RMT

    # aggregate all bootstrap replications
    for (k in E.set){
      bs.CumHaz$trt.0[[paste("Ev=", k, sep="")]] <- t(rbindlist(list(purrr::map(bs_aggregates,~.x[['trt.0']][[paste("Ev=", k, sep="")]][['CumHaz']]))))
      bs.CIF$trt.0[[paste("Ev=", k, sep="")]] <- t(rbindlist(list(purrr::map(bs_aggregates,~.x[['trt.0']][[paste("Ev=", k, sep="")]][['CIF']]))))
      bs.RMT$trt.0[[paste("Ev=", k, sep="")]] <- t(rbindlist(list(purrr::map(bs_aggregates,~.x[['trt.0']][[paste("Ev=", k, sep="")]][['RMT']]))))
      bs.CumHaz$trt.1[[paste("Ev=", k, sep="")]] <- t(rbindlist(list(purrr::map(bs_aggregates,~.x[['trt.1']][[paste("Ev=", k, sep="")]][['CumHaz']]))))
      bs.CIF$trt.1[[paste("Ev=", k, sep="")]] <- t(rbindlist(list(purrr::map(bs_aggregates,~.x[['trt.1']][[paste("Ev=", k, sep="")]][['CIF']]))))
      bs.RMT$trt.1[[paste("Ev=", k, sep="")]] <- t(rbindlist(list(purrr::map(bs_aggregates,~.x[['trt.1']][[paste("Ev=", k, sep="")]][['RMT']]))))
      bs.CumHaz$logRatio[[paste("Ev=", k, sep="")]] <- t(rbindlist(list((purrr::map(bs_aggregates,~.x[['trt.eff']][[paste("Ev=", k, sep="")]][['log.CumHazR']])))))
      bs.CIF$RD[[paste("Ev=", k, sep="")]] <- t(rbindlist(list(purrr::map(bs_aggregates,~.x[['trt.eff']][[paste("Ev=", k, sep="")]][['RD']]))))
      bs.CIF$RR[[paste("Ev=", k, sep="")]] <- t(rbindlist(list(purrr::map(bs_aggregates,~.x[['trt.eff']][[paste("Ev=", k, sep="")]][['RR']]))))
      bs.RMT$ATE[[paste("Ev=", k, sep="")]] <- t(rbindlist(list(purrr::map(bs_aggregates,~.x[['trt.eff']][[paste("Ev=", k, sep="")]][['ATE.RMT']]))))
    }

    alpha = 1-conf.level
    for (k in E.set)
    {
      # A. Cumulative Hazards::::::::::::::::::::::::::::::::
      # CI:
      res$trt.0[[paste("Ev=", k, sep="")]][["CumHaz.CI.L"]] <-
        apply(bs.CumHaz$trt.0[[paste("Ev=", k, sep="")]], 2, quantile, prob=alpha/2, na.rm=TRUE)
      res$trt.1[[paste("Ev=", k, sep="")]][["CumHaz.CI.L"]] <-
        apply(bs.CumHaz$trt.1[[paste("Ev=", k, sep="")]], 2, quantile, prob=alpha/2, na.rm=TRUE)
      res$trt.0[[paste("Ev=", k, sep="")]][["CumHaz.CI.U"]] <-
        apply(bs.CumHaz$trt.0[[paste("Ev=", k, sep="")]], 2, quantile, prob=1-alpha/2, na.rm=TRUE)
      res$trt.1[[paste("Ev=", k, sep="")]][["CumHaz.CI.U"]] <-
        apply(bs.CumHaz$trt.1[[paste("Ev=", k, sep="")]], 2, quantile, prob=1-alpha/2, na.rm=TRUE)
      res$trt.eff[[paste("Ev=", k, sep="")]][["log.CumHazR.CI.L"]] <-
        quantile(bs.CumHaz$logRatio[[paste("Ev=", k, sep="")]], prob=alpha/2, na.rm=TRUE)
      res$trt.eff[[paste("Ev=", k, sep="")]][["log.CumHazR.CI.U"]] <-
        quantile(bs.CumHaz$logRatio[[paste("Ev=", k, sep="")]], prob=1-alpha/2, na.rm=TRUE)

      # SE:
      res$trt.0[[paste("Ev=", k, sep="")]][["CumHaz.SE"]] <-
        apply(bs.CumHaz$trt.0[[paste("Ev=", k, sep="")]], 2, sd, na.rm=TRUE)
      res$trt.1[[paste("Ev=", k, sep="")]][["CumHaz.SE"]] <-
        apply(bs.CumHaz$trt.1[[paste("Ev=", k, sep="")]], 2, sd, na.rm=TRUE)
      res$trt.eff[[paste("Ev=", k, sep="")]][["log.CumHazR.SE"]] <-
        sd(bs.CumHaz$logRatio[[paste("Ev=", k, sep="")]], na.rm=TRUE)

      # bs.avg:
      res$trt.0[[paste("Ev=", k, sep="")]][["CumHaz.bs.avg"]] <-
        apply(bs.CumHaz$trt.0[[paste("Ev=", k, sep="")]], 2, mean, na.rm=TRUE)
      res$trt.1[[paste("Ev=", k, sep="")]][["CumHaz.bs.avg"]] <-
        apply(bs.CumHaz$trt.1[[paste("Ev=", k, sep="")]], 2, mean, na.rm=TRUE)
      res$trt.eff[[paste("Ev=", k, sep="")]][["log.CumHazR.bs.avg"]] <-
        mean(bs.CumHaz$logRatio[[paste("Ev=", k, sep="")]], na.rm=TRUE)

      # B. CIFs :::::::::::::::::::::::::::::::::::::::::::::::
      # CI:
      res$trt.0[[paste("Ev=", k, sep="")]][["CIF.CI.L"]] <-
        apply(bs.CIF$trt.0[[paste("Ev=", k, sep="")]], 2, quantile, prob=alpha/2, na.rm=TRUE)
      res$trt.1[[paste("Ev=", k, sep="")]][["CIF.CI.L"]] <-
        apply(bs.CIF$trt.1[[paste("Ev=", k, sep="")]], 2, quantile, prob=alpha/2, na.rm=TRUE)
      res$trt.0[[paste("Ev=", k, sep="")]][["CIF.CI.U"]] <-
        apply(bs.CIF$trt.0[[paste("Ev=", k, sep="")]], 2, quantile, prob=1-alpha/2, na.rm=TRUE)
      res$trt.1[[paste("Ev=", k, sep="")]][["CIF.CI.U"]] <-
        apply(bs.CIF$trt.1[[paste("Ev=", k, sep="")]], 2, quantile, prob=1-alpha/2, na.rm=TRUE)
      res$trt.eff[[paste("Ev=", k, sep="")]][["RD.CI.L"]] <-
        apply(bs.CIF$RD[[paste("Ev=", k, sep="")]], 2, quantile, prob=alpha/2, na.rm=TRUE)
      res$trt.eff[[paste("Ev=", k, sep="")]][["RD.CI.U"]] <-
        apply(bs.CIF$RD[[paste("Ev=", k, sep="")]], 2, quantile, prob=1-alpha/2, na.rm=TRUE)
      res$trt.eff[[paste("Ev=", k, sep="")]][["RR.CI.L"]] <-
        apply(bs.CIF$RR[[paste("Ev=", k, sep="")]], 2, quantile, prob=alpha/2, na.rm=TRUE)
      res$trt.eff[[paste("Ev=", k, sep="")]][["RR.CI.U"]] <-
        apply(bs.CIF$RR[[paste("Ev=", k, sep="")]], 2, quantile, prob=1-alpha/2, na.rm=TRUE)

      # SE:
      res$trt.0[[paste("Ev=", k, sep="")]][["CIF.SE"]] <-
        apply(bs.CIF$trt.0[[paste("Ev=", k, sep="")]], 2, sd, na.rm=TRUE)
      res$trt.1[[paste("Ev=", k, sep="")]][["CIF.SE"]] <-
        apply(bs.CIF$trt.1[[paste("Ev=", k, sep="")]], 2, sd, na.rm=TRUE)
      res$trt.eff[[paste("Ev=", k, sep="")]][["RD.SE"]] <-
        apply(bs.CIF$RD[[paste("Ev=", k, sep="")]], 2, sd, na.rm=TRUE)
      res$trt.eff[[paste("Ev=", k, sep="")]][["RR.SE"]] <-
        apply(bs.CIF$RR[[paste("Ev=", k, sep="")]], 2, sd, na.rm=TRUE)

      # bs.avg:
      res$trt.0[[paste("Ev=", k, sep="")]][["CIF.bs.avg"]] <-
        apply(bs.CIF$trt.0[[paste("Ev=", k, sep="")]], 2, mean, na.rm=TRUE)
      res$trt.1[[paste("Ev=", k, sep="")]][["CIF.bs.avg"]] <-
        apply(bs.CIF$trt.1[[paste("Ev=", k, sep="")]], 2, mean, na.rm=TRUE)
      res$trt.eff[[paste("Ev=", k, sep="")]][["RD.bs.avg"]] <-
        apply(bs.CIF$RD[[paste("Ev=", k, sep="")]], 2, mean, na.rm=TRUE)
      res$trt.eff[[paste("Ev=", k, sep="")]][["RR.bs.avg"]] <-
        apply(bs.CIF$RR[[paste("Ev=", k, sep="")]], 2, mean, na.rm=TRUE)

      # C. RMTs ::::::::::::::::::::::::::::::::::::::::::::::
      # CI:
      res$trt.0[[paste("Ev=", k, sep="")]][["RMT.CI.L"]] <-
        apply(bs.RMT$trt.0[[paste("Ev=", k, sep="")]], 2, quantile, prob=alpha/2, na.rm=TRUE)
      res$trt.1[[paste("Ev=", k, sep="")]][["RMT.CI.L"]] <-
        apply(bs.RMT$trt.1[[paste("Ev=", k, sep="")]], 2, quantile, prob=alpha/2, na.rm=TRUE)
      res$trt.0[[paste("Ev=", k, sep="")]][["RMT.CI.U"]] <-
        apply(bs.RMT$trt.0[[paste("Ev=", k, sep="")]], 2, quantile, prob=1-alpha/2, na.rm=TRUE)
      res$trt.1[[paste("Ev=", k, sep="")]][["RMT.CI.U"]] <-
        apply(bs.RMT$trt.1[[paste("Ev=", k, sep="")]], 2, quantile, prob=1-alpha/2, na.rm=TRUE)
      res$trt.eff[[paste("Ev=", k, sep="")]][["ATE.RMT.CI.L"]] <-
        apply(bs.RMT$ATE[[paste("Ev=", k, sep="")]], 2, quantile, prob=alpha/2, na.rm=TRUE)
      res$trt.eff[[paste("Ev=", k, sep="")]][["ATE.RMT.CI.U"]] <-
        apply(bs.RMT$ATE[[paste("Ev=", k, sep="")]], 2, quantile, prob=1-alpha/2, na.rm=TRUE)
      # SE:
      res$trt.0[[paste("Ev=", k, sep="")]][["RMT.SE"]] <-
        apply(bs.RMT$trt.0[[paste("Ev=", k, sep="")]], 2, sd, na.rm=TRUE)
      res$trt.1[[paste("Ev=", k, sep="")]][["RMT.SE"]] <-
        apply(bs.RMT$trt.1[[paste("Ev=", k, sep="")]], 2, sd, na.rm=TRUE)
      res$trt.eff[[paste("Ev=", k, sep="")]][["ATE.RMT.SE"]] <-
        apply(bs.RMT$ATE[[paste("Ev=", k, sep="")]], 2, sd, na.rm=TRUE)
      # bs.avg:
      res$trt.0[[paste("Ev=", k, sep="")]][["RMT.bs.avg"]] <-
        apply(bs.RMT$trt.0[[paste("Ev=", k, sep="")]], 2, mean, na.rm=TRUE)
      res$trt.1[[paste("Ev=", k, sep="")]][["RMT.bs.avg"]] <-
        apply(bs.RMT$trt.1[[paste("Ev=", k, sep="")]], 2, mean, na.rm=TRUE)
      res$trt.eff[[paste("Ev=", k, sep="")]][["ATE.RMT.bs.avg"]] <-
        apply(bs.RMT$ATE[[paste("Ev=", k, sep="")]], 2, mean, na.rm=TRUE)
    }
    class(res) <- "cmprsk"
    return(res)
  }
}

.parallel.fit.nonpar <- function(df, T, E, A, C, wtype, cens, conf.level, bs, nbs.rep, seed)
{
  X <- df[[T]]
  E <- df[[E]]
  nobs <- length(X)
  trt <- df[[A]]
  time <- sort(unique(X[E!=cens]))
  E.set <- sort(unique(E))
  E.set <- E.set[E.set!=cens]

  #  res <- list(time=time) # this is the only place where I'll keep the time
  res <- .nonpar.run(df, T, E, A, C, wtype, cens, E.set,time,trt,nobs,X, case.w = rep(1,nobs))

  # res is a list with 4 fields: time, trt.0, trt.0, trt.eff
  if (bs)
  {
    bs_seeds <- (1:nbs.rep) + seed
    # allocate memory for bs results:
    ntime <- length(res$time)
    bs.CumHaz <- bs.CIF <- bs.RMT <- list()
    for (k in E.set)
    {
      bs.CumHaz$trt.0[[paste("Ev=", k, sep="")]] <- bs.CumHaz$trt.1[[paste("Ev=", k, sep="")]] <-
        bs.CIF$trt.0[[paste("Ev=", k, sep="")]] <- bs.CIF$trt.1[[paste("Ev=", k, sep="")]] <-
        bs.RMT$trt.0[[paste("Ev=", k, sep="")]] <- bs.RMT$trt.1[[paste("Ev=", k, sep="")]] <-
        bs.CumHaz$logRatio[[paste("Ev=", k, sep="")]] <-
        bs.CIF$RD[[paste("Ev=", k, sep="")]] <- bs.CIF$RR[[paste("Ev=", k, sep="")]] <-
        bs.RMT$ATE[[paste("Ev=", k, sep="")]] <- matrix(nrow=nbs.rep, ncol=ntime)
    }
    i=0 # to remove the package NOTE
    # Create the progress bar.
    pb <- txtProgressBar(min = 1, max = nbs.rep, style=3, file=stdout())
    bs_aggregates <- foreach(i = 1:nbs.rep, .export=c(".nonpar.run", "get.weights", ".estimate.nonpar", ".base.haz.std",
                                                      ".get.CIF",".get.S", ".get.RMT"),
                             .packages=c("survival")) %dopar%
      {
        set.seed(bs_seeds[i])
        bs.w <- pmin(rexp(nobs,1), 5) # nobs = our sample size
        bs.w <- bs.w/mean(bs.w)
        bs_aggregates <- .nonpar.run(df, T, E, A, C, wtype, cens, E.set,time,trt,nobs,X,case.w = bs.w)
        setTxtProgressBar(pb, i)
        bs_aggregates
      }          # df, T, E, A, C, wtype, cens, E.set,time,trt,nobs,X,case.w
    #        if(i %% 25 == 0)
    #          cat(paste0('We are on replication ',i,' of ',nbs.rep,' replications.'),file= stdout())
    close(pb)

    for (k in E.set){
      bs.CumHaz$trt.0[[paste("Ev=", k, sep="")]] <- t(rbindlist(list(purrr::map(bs_aggregates,~.x[['trt.0']][[paste("Ev=", k, sep="")]][['CumHaz']]))))
      bs.CIF$trt.0[[paste("Ev=", k, sep="")]] <- t(rbindlist(list(purrr::map(bs_aggregates,~.x[['trt.0']][[paste("Ev=", k, sep="")]][['CIF']]))))
      bs.RMT$trt.0[[paste("Ev=", k, sep="")]] <- t(rbindlist(list(purrr::map(bs_aggregates,~.x[['trt.0']][[paste("Ev=", k, sep="")]][['RMT']]))))
      bs.CumHaz$trt.1[[paste("Ev=", k, sep="")]] <- t(rbindlist(list(purrr::map(bs_aggregates,~.x[['trt.1']][[paste("Ev=", k, sep="")]][['CumHaz']]))))
      bs.CIF$trt.1[[paste("Ev=", k, sep="")]] <- t(rbindlist(list(purrr::map(bs_aggregates,~.x[['trt.1']][[paste("Ev=", k, sep="")]][['CIF']]))))
      bs.RMT$trt.1[[paste("Ev=", k, sep="")]] <- t(rbindlist(list(purrr::map(bs_aggregates,~.x[['trt.1']][[paste("Ev=", k, sep="")]][['RMT']]))))
      bs.CumHaz$logRatio[[paste("Ev=", k, sep="")]] <- t(rbindlist(list(purrr::map(bs_aggregates,~.x[['trt.eff']][[paste("Ev=", k, sep="")]][['log.CumHazR']]))))
      bs.CIF$RD[[paste("Ev=", k, sep="")]] <- t(rbindlist(list(purrr::map(bs_aggregates,~.x[['trt.eff']][[paste("Ev=", k, sep="")]][['RD']]))))
      bs.CIF$RR[[paste("Ev=", k, sep="")]] <- t(rbindlist(list(purrr::map(bs_aggregates,~.x[['trt.eff']][[paste("Ev=", k, sep="")]][['RR']]))))
      bs.RMT$ATE[[paste("Ev=", k, sep="")]] <- t(rbindlist(list(purrr::map(bs_aggregates,~.x[['trt.eff']][[paste("Ev=", k, sep="")]][['ATE.RMT']]))))
    }


    # summarize bs replications and save the results in 'res' object:
    # res is a list with 4 fields:
    # 1.time - a vector of times for which everything is estimated
    # 2.trt.0[[paste("Ev=", k, sep="")]]: CumHaz, CIF, RMT
    # 3.trt.1[[paste("Ev=", k, sep="")]]: CumHaz, CIF, RMT
    # 4. trt.eff[[paste("Ev=", k, sep="")]]: log.CumHazR, RD, RR, ATE.RMT
    alpha = 1-conf.level
    for (k in E.set)
    {
      # A. Cumulative Hazards::::::::::::::::::::::::::::::::
      # CI:
      res$trt.0[[paste("Ev=", k, sep="")]][["CumHaz.CI.L"]] <-
        apply(bs.CumHaz$trt.0[[paste("Ev=", k, sep="")]], 2, quantile, prob=alpha/2, na.rm=TRUE)
      res$trt.1[[paste("Ev=", k, sep="")]][["CumHaz.CI.L"]] <-
        apply(bs.CumHaz$trt.1[[paste("Ev=", k, sep="")]], 2, quantile, prob=alpha/2, na.rm=TRUE)
      res$trt.0[[paste("Ev=", k, sep="")]][["CumHaz.CI.U"]] <-
        apply(bs.CumHaz$trt.0[[paste("Ev=", k, sep="")]], 2, quantile, prob=1-alpha/2, na.rm=TRUE)
      res$trt.1[[paste("Ev=", k, sep="")]][["CumHaz.CI.U"]] <-
        apply(bs.CumHaz$trt.1[[paste("Ev=", k, sep="")]], 2, quantile, prob=1-alpha/2, na.rm=TRUE)
      res$trt.eff[[paste("Ev=", k, sep="")]][["log.CumHazR.CI.L"]] <-
        apply(bs.CumHaz$logRatio[[paste("Ev=", k, sep="")]], 2, quantile, prob=alpha/2, na.rm=TRUE)
      res$trt.eff[[paste("Ev=", k, sep="")]][["log.CumHazR.CI.U"]] <-
        apply(bs.CumHaz$logRatio[[paste("Ev=", k, sep="")]], 2, quantile, prob=1-alpha/2, na.rm=TRUE)
      # SE:
      res$trt.0[[paste("Ev=", k, sep="")]][["CumHaz.SE"]] <-
        apply(bs.CumHaz$trt.0[[paste("Ev=", k, sep="")]], 2, sd, na.rm=TRUE)
      res$trt.1[[paste("Ev=", k, sep="")]][["CumHaz.SE"]] <-
        apply(bs.CumHaz$trt.1[[paste("Ev=", k, sep="")]], 2, sd, na.rm=TRUE)
      res$trt.eff[[paste("Ev=", k, sep="")]][["log.CumHazR.SE"]] <-
        apply(bs.CumHaz$logRatio[[paste("Ev=", k, sep="")]], 2, sd, na.rm=TRUE)
      # bs.avg:
      res$trt.0[[paste("Ev=", k, sep="")]][["CumHaz.bs.avg"]] <-
        apply(bs.CumHaz$trt.0[[paste("Ev=", k, sep="")]], 2, mean, na.rm=TRUE)
      res$trt.1[[paste("Ev=", k, sep="")]][["CumHaz.bs.avg"]] <-
        apply(bs.CumHaz$trt.1[[paste("Ev=", k, sep="")]], 2, mean, na.rm=TRUE)
      res$trt.eff[[paste("Ev=", k, sep="")]][["log.CumHazR.bs.avg"]] <-
        apply(bs.CumHaz$logRatio[[paste("Ev=", k, sep="")]], 2, mean, na.rm=TRUE)

      # B. CIFs :::::::::::::::::::::::::::::::::::::::::::::::
      # CI:
      res$trt.0[[paste("Ev=", k, sep="")]][["CIF.CI.L"]] <-
        apply(bs.CIF$trt.0[[paste("Ev=", k, sep="")]], 2, quantile, prob=alpha/2, na.rm=TRUE)
      res$trt.1[[paste("Ev=", k, sep="")]][["CIF.CI.L"]] <-
        apply(bs.CIF$trt.1[[paste("Ev=", k, sep="")]], 2, quantile, prob=alpha/2, na.rm=TRUE)
      res$trt.0[[paste("Ev=", k, sep="")]][["CIF.CI.U"]] <-
        apply(bs.CIF$trt.0[[paste("Ev=", k, sep="")]], 2, quantile, prob=1-alpha/2, na.rm=TRUE)
      res$trt.1[[paste("Ev=", k, sep="")]][["CIF.CI.U"]] <-
        apply(bs.CIF$trt.1[[paste("Ev=", k, sep="")]], 2, quantile, prob=1-alpha/2, na.rm=TRUE)
      res$trt.eff[[paste("Ev=", k, sep="")]][["RD.CI.L"]] <-
        apply(bs.CIF$RD[[paste("Ev=", k, sep="")]], 2, quantile, prob=alpha/2, na.rm=TRUE)
      res$trt.eff[[paste("Ev=", k, sep="")]][["RD.CI.U"]] <-
        apply(bs.CIF$RD[[paste("Ev=", k, sep="")]], 2, quantile, prob=1-alpha/2, na.rm=TRUE)
      res$trt.eff[[paste("Ev=", k, sep="")]][["RR.CI.L"]] <-
        apply(bs.CIF$RR[[paste("Ev=", k, sep="")]], 2, quantile, prob=alpha/2, na.rm=TRUE)
      res$trt.eff[[paste("Ev=", k, sep="")]][["RR.CI.U"]] <-
        apply(bs.CIF$RR[[paste("Ev=", k, sep="")]], 2, quantile, prob=1-alpha/2, na.rm=TRUE)

      # SE:
      res$trt.0[[paste("Ev=", k, sep="")]][["CIF.SE"]] <-
        apply(bs.CIF$trt.0[[paste("Ev=", k, sep="")]], 2, sd, na.rm=TRUE)
      res$trt.1[[paste("Ev=", k, sep="")]][["CIF.SE"]] <-
        apply(bs.CIF$trt.1[[paste("Ev=", k, sep="")]], 2, sd, na.rm=TRUE)
      res$trt.eff[[paste("Ev=", k, sep="")]][["RD.SE"]] <-
        apply(bs.CIF$RD[[paste("Ev=", k, sep="")]], 2, sd, na.rm=TRUE)
      res$trt.eff[[paste("Ev=", k, sep="")]][["RR.SE"]] <-
        apply(bs.CIF$RR[[paste("Ev=", k, sep="")]], 2, sd, na.rm=TRUE)

      # bs.avg:
      res$trt.0[[paste("Ev=", k, sep="")]][["CIF.bs.avg"]] <-
        apply(bs.CIF$trt.0[[paste("Ev=", k, sep="")]], 2, mean, na.rm=TRUE)
      res$trt.1[[paste("Ev=", k, sep="")]][["CIF.bs.avg"]] <-
        apply(bs.CIF$trt.1[[paste("Ev=", k, sep="")]], 2, mean, na.rm=TRUE)
      res$trt.eff[[paste("Ev=", k, sep="")]][["RD.bs.avg"]] <-
        apply(bs.CIF$RD[[paste("Ev=", k, sep="")]], 2, mean, na.rm=TRUE)
      res$trt.eff[[paste("Ev=", k, sep="")]][["RR.bs.avg"]] <-
        apply(bs.CIF$RR[[paste("Ev=", k, sep="")]], 2, mean, na.rm=TRUE)

      # C. RMTs ::::::::::::::::::::::::::::::::::::::::::::::
      # CI:
      res$trt.0[[paste("Ev=", k, sep="")]][["RMT.CI.L"]] <-
        apply(bs.RMT$trt.0[[paste("Ev=", k, sep="")]], 2, quantile, prob=alpha/2, na.rm=TRUE)
      res$trt.1[[paste("Ev=", k, sep="")]][["RMT.CI.L"]] <-
        apply(bs.RMT$trt.1[[paste("Ev=", k, sep="")]], 2, quantile, prob=alpha/2, na.rm=TRUE)
      res$trt.0[[paste("Ev=", k, sep="")]][["RMT.CI.U"]] <-
        apply(bs.RMT$trt.0[[paste("Ev=", k, sep="")]], 2, quantile, prob=1-alpha/2, na.rm=TRUE)
      res$trt.1[[paste("Ev=", k, sep="")]][["RMT.CI.U"]] <-
        apply(bs.RMT$trt.1[[paste("Ev=", k, sep="")]], 2, quantile, prob=1-alpha/2, na.rm=TRUE)
      res$trt.eff[[paste("Ev=", k, sep="")]][["ATE.RMT.CI.L"]] <-
        apply(bs.RMT$ATE[[paste("Ev=", k, sep="")]], 2, quantile, prob=alpha/2, na.rm=TRUE)
      res$trt.eff[[paste("Ev=", k, sep="")]][["ATE.RMT.CI.U"]] <-
        apply(bs.RMT$ATE[[paste("Ev=", k, sep="")]], 2, quantile, prob=1-alpha/2, na.rm=TRUE)
      # SE:
      res$trt.0[[paste("Ev=", k, sep="")]][["RMT.SE"]] <-
        apply(bs.RMT$trt.0[[paste("Ev=", k, sep="")]], 2, sd, na.rm=TRUE)
      res$trt.1[[paste("Ev=", k, sep="")]][["RMT.SE"]] <-
        apply(bs.RMT$trt.1[[paste("Ev=", k, sep="")]], 2, sd, na.rm=TRUE)
      res$trt.eff[[paste("Ev=", k, sep="")]][["ATE.RMT.SE"]] <-
        apply(bs.RMT$ATE[[paste("Ev=", k, sep="")]], 2, sd, na.rm=TRUE)
      # bs.avg:
      res$trt.0[[paste("Ev=", k, sep="")]][["RMT.bs.avg"]] <-
        apply(bs.RMT$trt.0[[paste("Ev=", k, sep="")]], 2, mean, na.rm=TRUE)
      res$trt.1[[paste("Ev=", k, sep="")]][["RMT.bs.avg"]] <-
        apply(bs.RMT$trt.1[[paste("Ev=", k, sep="")]], 2, mean, na.rm=TRUE)
      res$trt.eff[[paste("Ev=", k, sep="")]][["ATE.RMT.bs.avg"]] <-
        apply(bs.RMT$ATE[[paste("Ev=", k, sep="")]], 2, mean, na.rm=TRUE)
    }
  }
  class(res) <- "cmprsk"
  return(res)
}

#' Cox-based estimation of ATE corresponding to the target population
#'
#' Implements Cox-based estimation  of ATE assuming a structural proportional hazards model for two potential outcomes.
#' It provides three measures of treatment effects on time-to-event outcomes:
#' (1) cause-specific hazard ratios which are time-dependent measures under a nonparametric model,
#' (2) risk-based measures such as cause-specific risk differences and cause-specific risk ratios, and
#' (3) restricted-mean-time differences which quantify how much time on average was lost (or gained)
#' due to treatment by some specified time point.
#' Please see our package vignette for more details.
#'
#'
#' @param df a data frame that includes time-to-event \code{T}, type of event \code{E},
#' a treatment indicator \code{A}  and covariates \code{C}.
#' @param T a character string specifying the name of the time-to-event variable in \code{df}.
#' @param E a character string specifying the name of the "event type" variable in \code{df}.
#' @param A a character specifying the name of the treatment/exposure variable.
#' It is assumed that \code{A} is a numeric binary indicator with 0/1 values, where \code{A}=1
#' is assumed a treatment group, and \code{A}=0 a control group.
#' @param C a vector of character strings with variable names (potential confounders)
#' in the logistic regression model for Propensity Scores, i.e. P(A=1|C=c).
#' The default value of \code{C} is NULL corresponding to \code{wtype}="unadj"
#' that will estimate treatment effects in the raw (observed) data.
#' @param wtype a character string variable indicating the type of weights that will define the target
#' population for which the ATE will be estimated.
#' The default is "unadj" - this will not adjust for possible
#' treatment selection bias and will not use propensity scores weighting. It can be used, for example,
#' in data from a randized controlled trial (RCT) where there is no need for emulation of baseline randomization.
#' Other possible values are "stab.ATE", "ATE", "ATT", "ATC" and "overlap".
#' See Table 1 from Li, Morgan, and Zaslavsky (2018).
#' "stab.ATE" is defined as P(A=a)/P(A=a|C=c) - see Hernan et al. (2000).
#' @param cens an integer value in \code{E} that corresponds to censoring times recorded in \code{T}.
#' By default \code{fit.nonpar} assumes \code{cens}=0
#' @param conf.level the confidence level that will be used in the bootstrap confidence intervals.
#' The default is 0.95
#' @param bs a logical flag indicating whether to perform bootstrap in order to obtain confidence intervals. There are no
#' analytical confidence intervals in \code{fit.nonpar}
#' @param nbs.rep number of bootstrap replications
#' @param seed the random seed for the bootstrap, in order to make the results reproducible
#' @param parallel a logical flag indicating whether to perform bootstrap sequentially or in parallel
#' using 2 cores simultaneously. The default value is FALSE.
#'
#' @return  A list with the following fields:
#' \tabular{ll}{
#' \code{time}  \tab \cr a vector of time points for which all the parameters are estimated \cr
#' \code{trt.0} \tab \cr a list of estimates of the absolute counterfactual parameters
#' corresponding to \code{A}=0 and the type of event \code{E}. \code{trt.0} has the number of
#'  fields as the number of different types of events in the data set.
#' For each type of event there is a list of estimates:  \cr }
#' \itemize{
#' \item{\code{CumHaz}} a vector of cumulative hazard estimates
#' \item{\code{CIF}} a vector of cumulative incidence functions
#' \item{\code{RMT}} a vector of restricted mean time estimates}
#' \tabular{ll}{
#' \code{trt.1} \tab \cr a list of estimates of the absolute counterfactual parameters
#' corresponding to \code{A}=1 and the type of event \code{E}. \code{trt.1} has the number of
#'  fields as the number of different types of events  in the data set.
#' For each type of event there is a list of estimates:  \cr }
#' \itemize{
#' \item{\code{CumHaz}} a vector of cumulative hazard estimates
#' \item{\code{CIF}} a vector of cumulative incidence functions
#' \item{\code{RMT}} a vector of restricted mean time estimates}
#' \tabular{ll}{
#' \code{trt.eff} \tab \cr a list of estimates of the treatment effect measures
#' corresponding to the type of event \code{E}. \code{trt.eff} has the number of
#'  fields as the number of different types of events  in the data set.
#' For each type of event there is a list of estimates: \cr}
#' \itemize{
#' \item{\code{log.CumHazR}} an estimate of the log of the hazard ratio. It is a scalar since the Cox model is assumed.
#' \item{\code{RD}} a vector of time-varying risk difference between two treatment arms
#' \item{\code{RR}} a vector of time-varying risk ratio between two treatment arms
#' \item{\code{ATE.RMT}} a vector of the time-varying restricted mean time difference
#'  between two treatment arms}
#'
#' @seealso \code{\link{fit.nonpar}}, \code{\link{get.pointEst}}, \code{\link{causalCmprsk}}
#' @examples
#' # please see our package vignette for practical examples
#'
#' @references F. Li, K.L. Morgan, and A.M. Zaslavsky. 2018. Balancing Covariates via Propensity Score Weighting. Journal of the American Statistical Association 113 (521): 390–400.
#' @references M. Hernan, K.L. Morgan, and A.M. Zaslavsky. 2000. Balancing Covariates via Propensity Score Weighting. Journal of the American Statistical Association 113 (521): 390–400.
#'
#'
#' @export
fit.cox <- function(df, T, E, A, C=NULL, wtype="unadj", cens=0, conf.level=0.95, bs=FALSE, nbs.rep=400, seed=17, parallel = FALSE){

  get_os <- function(){
    platform <- .Platform$OS.type
    return(platform)
  }
  if(bs == TRUE && parallel == TRUE){

    # if windows (snow)
    if(grepl('win',get_os(),ignore.case = T)) {
      #      num_cores <- round(parallel::detectCores()/2,digits = 0)
      num_cores <- 2
      cluster <- parallel::makeCluster(num_cores)
      doParallel::registerDoParallel(cl = cluster)
    }

    # if unix (multicore)
    if(grepl('unix|linux',get_os(),ignore.case = T)){
      num_cores <- 2 # https://stackoverflow.com/questions/41307178/error-processing-vignette-failed-with-diagnostics-4-simultaneous-processes-spa
      cluster <- parallel::makeCluster(num_cores)
      doParallel::registerDoParallel(cl=cluster)
    }

    .parallel.fit.cox(df = df, T = T, E = E, A = A, C = C, wtype = wtype, cens = cens, conf.level = conf.level, bs = bs, nbs.rep = nbs.rep, seed = seed)
  }
  else
    .sequential.fit.cox(df = df, T = T, E = E, A = A, C = C, wtype = wtype, cens = cens, conf.level = conf.level, bs = bs, nbs.rep = nbs.rep, seed = seed)

}

#' Nonparametric estimation of ATE corresponding to the target population
#'
#'
#' Implements nonparametric estimation (based on the weighted Aalen-Johansen estimator) of ATE meaning that
#' it does not assume any model for potential outcomes.
#' It provides three measures of treatment effects on time-to-event outcomes:
#' (1) cause-specific hazard ratios which are time-dependent measures under a nonparametric model,
#' (2) risk-based measures such as cause-specific risk differences and cause-specific risk ratios, and
#' (3) restricted-mean-time differences which quantify how much time on average was lost (or gained)
#' due to treatment by some specified time point.
#' Please see our package vignette for more details.
#'
#'
#' @param df a data frame that includes time-to-event \code{T}, type of event \code{E},
#' a treatment indicator \code{A}  and covariates \code{C}.
#' @param T a character string specifying the name of the time-to-event variable in \code{df}.
#' @param E a character string specifying the name of the "event type" variable in \code{df}.
#' @param A a character specifying the name of the treatment/exposure variable.
#' It is assumed that \code{A} is a numeric binary indicator with 0/1 values, where \code{A}=1
#' is assumed a treatment group, and \code{A}=0 a control group.
#' @param C a vector of character strings with variable names (potential confounders)
#' in the logistic regression model for Propensity Scores, i.e. P(A=1|C=c).
#' The default value of \code{C} is NULL corresponding to \code{wtype}="unadj"
#' that will estimate treatment effects in the raw (observed) data.
#' @param wtype a character string variable indicating the type of weights that will define the target
#' population for which the ATE will be estimated.
#' The default is "unadj" - this will not adjust for possible
#' treatment selection bias and will not use propensity scores weighting. It can be used, for example,
#' in data from a randized controlled trial (RCT) where there is no need for emulation of baseline randomization.
#' Other possible values are "stab.ATE", "ATE", "ATT", "ATC" and "overlap".
#' See Table 1 from Li, Morgan, and Zaslavsky (2018).
#' "stab.ATE" is defined as P(A=a)/P(A=a|C=c) - see Hernan et al. (2000).
#' @param cens an integer value in \code{E} that corresponds to censoring times recorded in \code{T}.
#' By default \code{fit.nonpar} assumes \code{cens}=0
#' @param conf.level the confidence level that will be used in the bootstrap confidence intervals.
#' The default is 0.95
#' @param bs a logical flag indicating whether to perform bootstrap in order to obtain confidence intervals. There are no
#' analytical confidence intervals in \code{fit.nonpar}
#' @param nbs.rep number of bootstrap replications
#' @param seed the random seed for the bootstrap, in order to make the results reproducible
#' @param parallel a logical flag indicating whether to perform bootstrap sequentially or in parallel
#' using 2 cores simultaneously. The default value is FALSE.
#'
#' @return  A list with the following fields:
#' \tabular{ll}{
#' \code{time}  \tab \cr a vector of time points for which all the parameters are estimated \cr
#' \code{trt.0} \tab \cr a list of estimates of the absolute counterfactual parameters
#' corresponding to \code{A}=0 and the type of event \code{E}. \code{trt.0} has the number of
#'  fields as the number of different types of events in the data set.
#' For each type of event there is a list of estimates:  \cr }
#' \itemize{
#' \item{\code{CumHaz}} a vector of cumulative hazard estimates
#' \item{\code{CIF}} a vector of cumulative incidence functions
#' \item{\code{RMT}} a vector of restricted mean time estimates}
#' \tabular{ll}{
#' \code{trt.1} \tab \cr a list of estimates of the absolute counterfactual parameters
#' corresponding to \code{A}=1 and the type of event \code{E}. \code{trt.1} has the number of
#'  fields as the number of different types of events  in the data set.
#' For each type of event there is a list of estimates:  \cr }
#' \itemize{
#' \item{\code{CumHaz}} a vector of cumulative hazard estimates
#' \item{\code{CIF}} a vector of cumulative incidence functions
#' \item{\code{RMT}} a vector of restricted mean time estimates}
#' \tabular{ll}{
#' \code{trt.eff} \tab \cr a list of estimates of the treatment effect measures
#' corresponding to the type of event \code{E}. \code{trt.eff} has the number of
#'  fields as the number of different types of events  in the data set.
#' For each type of event there is a list of estimates: \cr}
#' \itemize{
#' \item{\code{log.CumHazR}} a vector of the log of the time-varying ratio of hazards in two treatment arms
#' \item{\code{RD}} a vector of time-varying risk difference between two treatment arms
#' \item{\code{RR}} a vector of time-varying risk ratio between two treatment arms
#' \item{\code{ATE.RMT}} a vector of the time-varying restricted mean time difference
#'  between two treatment arms}
#'
#' @seealso \code{\link{fit.cox}}, \code{\link{get.pointEst}}, \code{\link{causalCmprsk}}
#' @examples
#' # please see our package vignette for practical examples
#'
#' @references F. Li, K.L. Morgan, and A.M. Zaslavsky. 2018. Balancing Covariates via Propensity Score Weighting. Journal of the American Statistical Association 113 (521): 390–400.
#' @references M. Hernan, K.L. Morgan, and A.M. Zaslavsky. 2000. Balancing Covariates via Propensity Score Weighting. Journal of the American Statistical Association 113 (521): 390–400.
#'
#' @export
fit.nonpar <- function(df, T, E, A, C=NULL, wtype="unadj", cens=0, conf.level=0.95, bs=FALSE, nbs.rep=400, seed=17, parallel = FALSE){

  get_os <- function(){
    platform <- .Platform$OS.type
    return(platform)
  }
  if(bs == TRUE && parallel == TRUE){

    # if windows (snow)
    if(grepl('win',get_os(),ignore.case = T)) {
#      num_cores <- round(parallel::detectCores()/2,digits = 0)
      num_cores <- 2
      cluster <- parallel::makeCluster(num_cores)
      doParallel::registerDoParallel(cl = cluster)
    }

    # if unix (multicore)
    if(grepl('unix|linux',get_os(),ignore.case = T)){
      num_cores <- 2 # https://stackoverflow.com/questions/41307178/error-processing-vignette-failed-with-diagnostics-4-simultaneous-processes-spa
      cluster <- parallel::makeCluster(num_cores)
      doParallel::registerDoParallel(cl=cluster)
    }

    .parallel.fit.nonpar(df = df, T = T, E = E, A = A, C = C, wtype = wtype, cens = cens, conf.level = conf.level, bs = bs, nbs.rep = nbs.rep, seed = seed)
  }
  else
    .sequential.fit.nonpar(df = df, T = T, E = E, A = A, C = C, wtype = wtype, cens = cens, conf.level = conf.level, bs = bs, nbs.rep = nbs.rep, seed = seed)

}

