

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

# get.weights fits PS model, returns both PSs and requested weights. PSs estimates can be used for further diagnostics - GOF, positivity and covariates balance - see below.
# df - a data frame containing the variables in the propensity score model, i.e.containing the treatment variable and confounders.
# A - a character specifying the name of the outcome variable in df. treatment/exposure variable. It is assumed that A is a numeric binary indicator with 0/1 values, where A=1 is assumed a treatment group, and A=0 a control group.
# C - a vector of character strings with variable names (potential confounders) in the logistic regression model for Propensity Scores, i.e. P(A=1|C=c).

# wtype - a character variable indicating the type of weight.
#The default is "stab.ATE" defined as w_stab=P(A=a)/P(A=a|C=c) - see Hernan et al. (2000).
#Other possible values are "ATE", "ATT", "ATC" and "overlap". See Table 1 from Li, Morgan, and Zaslavsky (2018).

# case.w - should it be a part of a data frame? or an external vector?


#' Fitting propensity scores model and estimating weights
#'
#' Description...
#'
#' @param case.w a vector of sampling weights
#'
#'
#' @references F. Li, K.L. Morgan, and A.M. Zaslavsky. 2018. Balancing Covariates via Propensity Score Weighting. Journal of the American Statistical Association 113 (521): 390–400.
#' @references M. Hernan, K.L. Morgan, and A.M. Zaslavsky. 2000. Balancing Covariates via Propensity Score Weighting. Journal of the American Statistical Association 113 (521): 390–400.
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

# wtype - a character variable indicating the type of weight.
# The default is "stab.ATE" defined as w_stab=P(A=a)/P(A=a|C=c) - see Hernan et al. (2000).
# Other possible values are "ATE", "ATT", "ATC" and "overlap". See Table 1 from Li, Morgan, and Zaslavsky (2018).

#' Obtaining number-at-risk statistic in the raw and weighted data
#'
#'
#' Description ...
#'
#' @param df a data frame with ...
#'
#' @examples
#' # please see our package vignette for practical examples
#'
#' @export
get.numAtRisk <- function(df, T, E, A, C, wtype="stab.ATE", cens=0)
{
  X <- df[[T]]
  E <- df[[E]]
  nobs <- length(X)
  trt <- df[[A]]
  ps.fit <- get.weights(df=df, A=A, C=C, wtype=wtype)

  res <- list()
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


# T - a character string specifying the name of the time-to-event variable in df
# E - a character string specifying the name of the "event type" variable in df
# A - a character specifying the name of the outcome variable in df. treatment/exposure variable. It is assumed that A is a numeric binary indicator with 0/1 values, where A=1 is assumed a treatment group, and A=0 a control group.
# C - a vector of character strings with variable names (potential confounders) in the logistic regression model for Propensity Scores, i.e. P(A=1|C=c).
# seed - for bs
# cens - censoring value
# wtype - a character variable indicating the type of weight.
# The default is "stab.ATE" defined as w_stab=P(A=a)/P(A=a|C=c) - see Hernan et al. (2000).
# Other possible values are "ATE", "ATT", "ATC" and "overlap". See Table 1 from Li, Morgan, and Zaslavsky (2018).
# There is also an option of "unajusted" estimation, set with wtype="unadj", that does not adjust for possible
# treatment selection bias and does not use propensity scores weighting. It can be used, for example,
# in data from an RCT where there is no need for emulation of baseline randomization.


#' Nonparametric estimation of ATE corresponding to the target population
#'
#'
#' Description ...
#'
#'
#' @param df a data frame with time to event, type of event, treatment indicator and covariates.
#' @param T a character string specifying the name of the time-to-event variable in \code{df}.
#' @param E a character string specifying the name of the "event type" variable in \code{df}.
#' @param A a character specifying the name of the treatment/exposure variable.
#' It is assumed that \code{A} is a numeric binary indicator with 0/1 values, where \code{A}=1
#' is assumed a treatment group, and \code{A}=0 a control group.
#' @param C a vector of character strings with variable names (potential confounders)
#' in the logistic regression model for Propensity Scores, i.e. P(A=1|C=c).
#' @param seed for the bootstrap
#' @param cens integer value in \code{E} that corresponds to censorings resorded in \code{T}.
#' @param wtype a character variable indicating the type of weight.
#' The default is "stab.ATE" defined as w_stab=P(A=a)/P(A=a|C=c) - see Hernan et al. (2000).
#' Other possible values are "ATE", "ATT", "ATC" and "overlap". See Table 1 from Li, Morgan, and Zaslavsky (2018).
#' There is also an option of "unajusted" estimation, by settinh wtype="unadj" - this will not adjust for possible
#' treatment selection bias and will not use propensity scores weighting. It can be used, for example,
#' in data from an RCT where there is no need for emulation of baseline randomization.
#'
#' @return  A list with components:
#'
#' @examples
#' # please see our package vignette for practical examples
#'
#' @export
fit.nonpar <- function(df, T, E, A, C, wtype="stab.ATE", bs=FALSE, nbs.rep=400, seed=17, cens=0, conf.level=0.95)
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

#  set.seed(seed)

  if (bs)
  {
    # allocate memory for bs results:
    bs_seeds <- (1:nbs.rep) + seed
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


# T - a character string specifying the name of the time-to-event variable in df
# E - a character string specifying the name of the "event type" variable in df
# A - a character specifying the name of the outcome variable in df. treatment/exposure variable. It is assumed that A is a numeric binary indicator with 0/1 values, where A=1 is assumed a treatment group, and A=0 a control group.
# C - a vector of character strings with variable names (potential confounders) in the logistic regression model for Propensity Scores, i.e. P(A=1|C=c).
# seed - for bs
# cens - censoring value


#' Cox-based estimation of ATE corresponding to the target population
#'
#'
#' Description ...
#'
#'
#' @param df a data frame with time to event, type of event, treatment indicator and covariates.
#' @param T a character string specifying the name of the time-to-event variable in \code{df}.
#' @param E a character string specifying the name of the "event type" variable in \code{df}.
#' @param A a character specifying the name of the treatment/exposure variable.
#' It is assumed that \code{A} is a numeric binary indicator with 0/1 values, where \code{A}=1
#' is assumed a treatment group, and \code{A}=0 a control group.
#' @param C a vector of character strings with variable names (potential confounders)
#' in the logistic regression model for Propensity Scores, i.e. P(A=1|C=c).
#' @param seed for the bootstrap
#' @param cens integer value in \code{E} that corresponds to censorings resorded in \code{T}.
#' @param wtype a character variable indicating the type of weight.
#' The default is "stab.ATE" defined as w_stab=P(A=a)/P(A=a|C=c) - see Hernan et al. (2000).
#' Other possible values are "ATE", "ATT", "ATC" and "overlap". See Table 1 from Li, Morgan, and Zaslavsky (2018).
#' There is also an option of "unajusted" estimation, by settinh wtype="unadj" - this will not adjust for possible
#' treatment selection bias and will not use propensity scores weighting. It can be used, for example,
#' in data from an RCT where there is no need for emulation of baseline randomization.
#'
#' @return  A list with components:
#'
#' @examples
#' # please see our package vignette for practical examples
#'
#' @export
fit.cox <- function(df, T, E, A, C, wtype="stab.ATE", bs=FALSE, nbs.rep=400, seed=17, cens=0, conf.level=0.95)
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

  set.seed(seed)
  if (bs)
  {
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
      bs.CumHaz$logRatio[[paste("Ev=", k, sep="")]] <- vector("double", len=nbs.rep)
    }
    # sequential bootstrap:
    for (i in 1:nbs.rep)
    {
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

#' Returns point estimates corresponding to a time point
#'
#'
#' @param cmprsk.obj a "cmprsk" object returned by one of the function \code{fit.cox} or \code{fit.nonpar}
#' @param timepoint a scalar value of the time point of interest
#'
#' @return  A list with components: ...
#'
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
    point.res$trt.0[[k]] <- list()
    point.res$trt.1[[k]] <- list()
    point.res$trt.eff[[k]] <- list()
    #trt.0
    point.res$trt.0[[k]]$CumHaz <- .get.CIF(timepoint, cmprsk.obj$time, cmprsk.obj$trt.0[[k]]$CumHaz)
    point.res$trt.0[[k]]$CumHaz.CI.L <- .get.CIF(timepoint, cmprsk.obj$time, cmprsk.obj$trt.0[[k]]$CumHaz.CI.L)
    point.res$trt.0[[k]]$CumHaz.CI.U <- .get.CIF(timepoint, cmprsk.obj$time, cmprsk.obj$trt.0[[k]]$CumHaz.CI.U)

    point.res$trt.0[[k]]$CIF <- .get.CIF(timepoint, cmprsk.obj$time, cmprsk.obj$trt.0[[k]]$CIF)
    point.res$trt.0[[k]]$CIF.CI.L <- .get.CIF(timepoint, cmprsk.obj$time, cmprsk.obj$trt.0[[k]]$CIF.CI.L)
    point.res$trt.0[[k]]$CIF.CI.U <- .get.CIF(timepoint, cmprsk.obj$time, cmprsk.obj$trt.0[[k]]$CIF.CI.U)

    point.res$trt.0[[k]]$RMT <- .get.CIF(timepoint, cmprsk.obj$time, cmprsk.obj$trt.0[[k]]$RMT)
    point.res$trt.0[[k]]$RMT.CI.L <- .get.CIF(timepoint, cmprsk.obj$time, cmprsk.obj$trt.0[[k]]$RMT.CI.L)
    point.res$trt.0[[k]]$RMT.CI.U <- .get.CIF(timepoint, cmprsk.obj$time, cmprsk.obj$trt.0[[k]]$RMT.CI.U)

    # trt.1
    point.res$trt.1[[k]]$CumHaz <- .get.CIF(timepoint, cmprsk.obj$time, cmprsk.obj$trt.1[[k]]$CumHaz)
    point.res$trt.1[[k]]$CumHaz.CI.L <- .get.CIF(timepoint, cmprsk.obj$time, cmprsk.obj$trt.1[[k]]$CumHaz.CI.L)
    point.res$trt.1[[k]]$CumHaz.CI.U <- .get.CIF(timepoint, cmprsk.obj$time, cmprsk.obj$trt.1[[k]]$CumHaz.CI.U)

    point.res$trt.1[[k]]$CIF <- .get.CIF(timepoint, cmprsk.obj$time, cmprsk.obj$trt.1[[k]]$CIF)
    point.res$trt.1[[k]]$CIF.CI.L <- .get.CIF(timepoint, cmprsk.obj$time, cmprsk.obj$trt.1[[k]]$CIF.CI.L)
    point.res$trt.1[[k]]$CIF.CI.U <- .get.CIF(timepoint, cmprsk.obj$time, cmprsk.obj$trt.1[[k]]$CIF.CI.U)

    point.res$trt.1[[k]]$RMT <- .get.CIF(timepoint, cmprsk.obj$time, cmprsk.obj$trt.1[[k]]$RMT)
    point.res$trt.1[[k]]$RMT.CI.L <- .get.CIF(timepoint, cmprsk.obj$time, cmprsk.obj$trt.1[[k]]$RMT.CI.L)
    point.res$trt.1[[k]]$RMT.CI.U <- .get.CIF(timepoint, cmprsk.obj$time, cmprsk.obj$trt.1[[k]]$RMT.CI.U)

    # trt.eff
    if (length(point.res$trt.eff[[k]]$log.CumHazR)!=1)
    {
      point.res$trt.eff[[k]]$log.CumHazR <- .get.CIF(timepoint, cmprsk.obj$time, cmprsk.obj$trt.eff[[k]]$log.CumHazR)
      point.res$trt.eff[[k]]$log.CumHazR.CI.L <- .get.CIF(timepoint, cmprsk.obj$time, cmprsk.obj$trt.eff[[k]]$log.CumHazR.CI.L)
      point.res$trt.eff[[k]]$log.CumHazR.CI.U <- .get.CIF(timepoint, cmprsk.obj$time, cmprsk.obj$trt.eff[[k]]$log.CumHazR.CI.U)
    }
    else
    {
      point.res$trt.eff[[k]]$CumHaz <- cmprsk.obj$trt.eff[[k]]$log.CumHazR
      point.res$trt.eff[[k]]$CumHaz.CI.L <- cmprsk.obj$trt.eff[[k]]$log.CumHazR.CI.L
      point.res$trt.eff[[k]]$CumHaz.CI.U <- cmprsk.obj$trt.eff[[k]]$log.CumHazR.CI.U
    }

    point.res$trt.eff[[k]]$RD <- .get.CIF(timepoint, cmprsk.obj$time, cmprsk.obj$trt.eff[[k]]$RD)
    point.res$trt.eff[[k]]$RD.CI.L <- .get.CIF(timepoint, cmprsk.obj$time, cmprsk.obj$trt.eff[[k]]$RD.CI.L)
    point.res$trt.eff[[k]]$RD.CI.U <- .get.CIF(timepoint, cmprsk.obj$time, cmprsk.obj$trt.eff[[k]]$RD.CI.U)

    point.res$trt.eff[[k]]$RR <- .get.CIF(timepoint, cmprsk.obj$time, cmprsk.obj$trt.eff[[k]]$RR)
    point.res$trt.eff[[k]]$RR.CI.L <- .get.CIF(timepoint, cmprsk.obj$time, cmprsk.obj$trt.eff[[k]]$RR.CI.L)
    point.res$trt.eff[[k]]$RR.CI.U <- .get.CIF(timepoint, cmprsk.obj$time, cmprsk.obj$trt.eff[[k]]$RR.CI.U)

    point.res$trt.eff[[k]]$ATE.RMT <- .get.CIF(timepoint, cmprsk.obj$time, cmprsk.obj$trt.eff[[k]]$ATE.RMT)
    point.res$trt.eff[[k]]$ATE.RMT.CI.L <- .get.CIF(timepoint, cmprsk.obj$time, cmprsk.obj$trt.eff[[k]]$ATE.RMT.CI.L)
    point.res$trt.eff[[k]]$ATE.RMT.CI.U <- .get.CIF(timepoint, cmprsk.obj$time, cmprsk.obj$trt.eff[[k]]$ATE.RMT.CI.U)

  }

  return(point.res)
}




# Parallel Cox Fit

boot_comb <- function(x, ...) {
  # takes each bootstrap iteration and combines it with some logic
}

get_os <- function(){
platform <- .Platform$OS.type
return(platform)
}

# if windows (snow)
if(grepl('win',get_os(),ignore.case = T)) {
num_cores <- round(detectCores()/2,digits = 0)
cluster <- makeCluster(num_cores)
registerDoParallel(cl = cluster)
}

# if unix (multicore)
if(grepl('unix|linux',get_os(),ignore.case = T)){
registerDoParallel()
}

df <- data.frame(trt = rbinom(100,1,.5),X = rexp(100,1),E = {rbinom(100,1,.3) + 1},C = rbinom(100,1,.7))

# Cox-based estimation:
#res.cox.stab.ATE <- fit.cox(df=rhc, T="T", E="E", A="RHC", C=covs.names, wtype="stab.ATE", bs=TRUE, nbs.rep=200, seed=17, cens=0, conf.level=0.95)

demo_parallel_boot <- function(i,df, T, E, A, C, wtype,E.set,cens,time) {
  X <- df[[T]]
  E <- df[[E]]
  nobs <- length(X)
  trt <- df[[A]]
  bs.w <- pmin(rexp(nobs,1), 5) # nobs = our sample size
  bs.w <- bs.w/mean(bs.w)
    if (wtype!="unadj"){
        bs.ps.fit <- get.weights(df,A,C,wtype, case.w = bs.w)
        bs.est <- .estimate.cox(X=X, E=E, trt=trt, case.w =bs.w*bs.ps.fit$w, cens=cens, time=time, E.set=E.set)
      }
    else
      bs.est <- .estimate.cox(X=X, E=E, trt=df$trt, case.w =bs.w, cens=cens, time=time, E.set=E.set)
}


my_vector <- foreach(i = 1:10) %dopar% demo_parallel_boot(i,df,'X','E','trt',c('C'),'stab.ATE',c(1,2),cens = 0,time = sort(df$X))

# accumulate results

### inside event 1, we have list of 2: logHR is scalar, bh is a dataframe with 2 cols, hazard and time
### combine logHR into  a vector
### combine bh, ignore time, combine hazards row for every iteration (rbind)

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
      bs.CumHaz$logRatio[[paste("Ev=", k, sep="")]] <- vector("double", len=nbs.rep)
    }

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



parallel.fit.cox <- function(df, T, E, A, C, wtype="stab.ATE", bs=FALSE, nbs.rep=400, seed=17, cens=0, conf.level=0.95)
{
  X <- df[[T]]
  E <- df[[E]]
  nobs <- length(X)
  trt <- df[[A]]
  time <- sort(unique(X[E!=cens]))
  E.set <- sort(unique(E))
  E.set <- E.set[E.set!=cens]

  res <- list(time=time)
  res <- .cox.run(df, T, E, A, C, wtype, cens, E.set,time,trt,nobs,X,case.w = rep(1,nobs))
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
      bs.CumHaz$logRatio[[paste("Ev=", k, sep="")]] <- vector("double", len=nbs.rep)
    }
    # sequential bootstrap:
    foreach(i = 1:nbs.rep) %dopar%
    {
      set.seed(bs_seeds[i])
      bs.w <- pmin(rexp(nobs,1), 5) # nobs = our sample size
      bs.w <- bs.w/mean(bs.w)
      bs_aggregates <- .cox.run(df, T, E, A, C, wtype, cens, E.set,time,trt,nobs,X,case.w = bs.w)
    }
    # summarize bs replications and save the results in 'res' object:
    # res is a list with 4 fields:
    # 1.time - a vector of times for which everything is estimated
    # 2.trt.0[[paste("Ev=", k, sep="")]]: CumHaz, CIF, RMT
    # 3.trt.1[[paste("Ev=", k, sep="")]]: CumHaz, CIF, RMT
    # 4. trt.eff[[paste("Ev=", k, sep="")]]: log.CumHazR, RD, RR, ATE.RMT

    # aggregate all bootstrap replications
    for (k in E.set){
    bs.CumHaz$trt.0[[paste("Ev=", k, sep="")]] <- t(rbindlist(list(map(bs_aggregates,~.x[['trt.0']][[paste("Ev=", k, sep="")]][['CumHaz']]))))
    bs.CIF$trt.0[[paste("Ev=", k, sep="")]] <- t(rbindlist(list(map(bs_aggregates,~.x[['trt.0']][[paste("Ev=", k, sep="")]][['CIF']]))))
    bs.RMT$trt.0[[paste("Ev=", k, sep="")]] <- t(rbindlist(list(map(bs_aggregates,~.x[['trt.0']][[paste("Ev=", k, sep="")]][['RMT']]))))
    bs.CumHaz$trt.1[[paste("Ev=", k, sep="")]] <- t(rbindlist(list(map(bs_aggregates,~.x[['trt.1']][[paste("Ev=", k, sep="")]][['CumHaz']]))))
    bs.CIF$trt.1[[paste("Ev=", k, sep="")]] <- t(rbindlist(list(map(bs_aggregates,~.x[['trt.1']][[paste("Ev=", k, sep="")]][['CIF']]))))
    bs.RMT$trt.1[[paste("Ev=", k, sep="")]] <- t(rbindlist(list(map(bs_aggregates,~.x[['trt.1']][[paste("Ev=", k, sep="")]][['RMT']]))))
    bs.CumHaz$logRatio[[paste("Ev=", k, sep="")]] <- t(rbindlist(map(bs_aggregates,~.x[['trt.eff']][[paste("Ev=", k, sep="")]][['log.CumHazR']])))
    bs.CIF$RD[[paste("Ev=", k, sep="")]] <- t(rbindlist(list(map(bs_aggregates,~.x[['trt.eff']][[paste("Ev=", k, sep="")]][['RD']]))))
    bs.CIF$RR[[paste("Ev=", k, sep="")]] <- t(rbindlist(list(map(bs_aggregates,~.x[['trt.eff']][[paste("Ev=", k, sep="")]][['RR']]))))
    bs.RMT$ATE[[paste("Ev=", k, sep="")]] <- t(rbindlist(list(map(bs_aggregates,~.x[['trt.eff']][[paste("Ev=", k, sep="")]][['ATE.RMT']]))))
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

parallel.fit.nonpar <- function(df, T, E, A, C, wtype="stab.ATE", bs=FALSE, nbs.rep=400, seed=17, cens=0, conf.level=0.95)
{
  X <- df[[T]]
  E <- df[[E]]
  nobs <- length(X)
  trt <- df[[A]]
  time <- sort(unique(X[E!=cens]))
  E.set <- sort(unique(E))
  E.set <- E.set[E.set!=cens]

  res <- list(time=time) # this is the only place where I'll keep the time
  res <- .nonpar.run(df, T, E, A, C, wtype, cens, E.set,time,trt,nobs,X,case.w)

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
    bs_seeds <- seq(1,nbs.rep,1) + seed
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
    # sequential bootstrap:
    foreach(i = 1:nbs.rep) %dopar%
    {
      set.seed(bs_seeds[i])
      bs.w <- pmin(rexp(nobs,1), 5) # nobs = our sample size
      bs.w <- bs.w/mean(bs.w)
      bs_aggregates <- .nonpar.run(df, T, E, A, C, wtype, cens, E.set,time,trt,nobs,X,case.w = bs.w)
    }

    for (k in E.set){
    bs.CumHaz$trt.0[[paste("Ev=", k, sep="")]] <- t(rbindlist(list(map(bs_aggregates,~.x[['trt.0']][[paste("Ev=", k, sep="")]][['CumHaz']]))))
    bs.CIF$trt.0[[paste("Ev=", k, sep="")]] <- t(rbindlist(list(map(bs_aggregates,~.x[['trt.0']][[paste("Ev=", k, sep="")]][['CIF']]))))
    bs.RMT$trt.0[[paste("Ev=", k, sep="")]] <- t(rbindlist(list(map(bs_aggregates,~.x[['trt.0']][[paste("Ev=", k, sep="")]][['RMT']]))))
    bs.CumHaz$trt.1[[paste("Ev=", k, sep="")]] <- t(rbindlist(list(map(bs_aggregates,~.x[['trt.1']][[paste("Ev=", k, sep="")]][['CumHaz']]))))
    bs.CIF$trt.0[[paste("Ev=", k, sep="")]] <- t(rbindlist(list(map(bs_aggregates,~.x[['trt.1']][[paste("Ev=", k, sep="")]][['CIF']]))))
    bs.RMT$trt.0[[paste("Ev=", k, sep="")]] <- t(rbindlist(list(map(bs_aggregates,~.x[['trt.1']][[paste("Ev=", k, sep="")]][['RMT']]))))
    bs.CumHaz$logRatio[[paste("Ev=", k, sep="")]] <- t(rbindlist(list(map(bs_aggregates,~.x[['trt.eff']][[paste("Ev=", k, sep="")]][['log.CumHazR']]))))
    bs.CIF$RD[[paste("Ev=", k, sep="")]] <- t(rbindlist(list(map(bs_aggregates,~.x[['trt.eff']][[paste("Ev=", k, sep="")]][['RD']]))))
    bs.CIF$RR[[paste("Ev=", k, sep="")]] <- t(rbindlist(list(map(bs_aggregates,~.x[['trt.eff']][[paste("Ev=", k, sep="")]][['RR']]))))
    bs.RMT$ATE[[paste("Ev=", k, sep="")]] <- t(rbindlist(list(map(bs_aggregates,~.x[['trt.eff']][[paste("Ev=", k, sep="")]][['ATE.RMT']]))))
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








