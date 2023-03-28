.get.S <- function(t, time, surv)	# assumes that t is a scalar
{
  if (length(t)!=1) stop("Error in .get.S: t is not a scalar")
  n <- length(time)
  if (t<time[1] && abs(t-time[1])>100*.Machine$double.eps) return(1)
  if (t>time[n] && abs(t-time[n])>100*.Machine$double.eps) return(surv[n])
  return( surv[sum(time<=t | abs(time-t)<100*.Machine$double.eps)] )
}

# obtain CIF or CumHaz for arbitrary t:
.get.CIF <- function(t, time, cuminc)	# assumes that t is a scalar
{
  if (length(t)!=1) stop("Error in .get.S: t is not a scalar")
  n <- length(time)
  if (t<time[1] && abs(t-time[1])>100*.Machine$double.eps) return(0)
  if (t>time[n] && abs(t-time[n])>100*.Machine$double.eps) return(cuminc[n])
  return( cuminc[sum(time<=t | abs(time-t)<100*.Machine$double.eps)] )
}

# RMT is an integral of CIF:
.get.RMT <- function(t, time, cuminc)	# assumes that t is a scalar
{
  if (length(t)!=1) stop("Error in .get.S: t is not a scalar")
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
#' The default is "stab.ATE" defined as P(A=a)/P(A=a|C=c) - see Hernán et al. (2000).
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
#' @references M.A. Hernán, B. Brumback, and J.M. Robins. 2000. Marginal structural models and to estimate the causal effect of zidovudine on the survival of HIV-positive men. Epidemiology, 11 (5): 561-570.
#'
#' @seealso \code{\link{fit.nonpar}}, \code{\link{fit.cox}}, \code{\link{causalCmprsk}}
#' @examples
#' # create a data set
#' n <- 1000
#' set.seed(7)
#' c1 <- runif(n)
#' c2 <- as.numeric(runif(n)< 0.2)
#' set.seed(77)
#' cf.m.T1 <- rweibull(n, shape=1, scale=exp(-(-1 + 2*c1)))
#' cf.m.T2 <-  rweibull(n, shape=1, scale=exp(-(1 + 1*c2)))
#' cf.m.T <- pmin( cf.m.T1, cf.m.T2)
#' cf.m.E <- rep(0, n)
#' cf.m.E[cf.m.T1<=cf.m.T2] <- 1
#' cf.m.E[cf.m.T2<cf.m.T1] <- 2
#' set.seed(77)
#' cf.s.T1 <- rweibull(n, shape=1, scale=exp(-1*c1 ))
#' cf.s.T2 <-  rweibull(n, shape=1, scale=exp(-2*c2))
#' cf.s.T <- pmin( cf.s.T1, cf.s.T2)
#' cf.s.E <- rep(0, n)
#' cf.s.E[cf.s.T1<=cf.s.T2] <- 1
#' cf.s.E[cf.s.T2<cf.s.T1] <- 2
#' exp.z <- exp(0.5 + 1*c1 - 1*c2)
#' pr <- exp.z/(1+exp.z)
#' TRT <- ifelse(runif(n)< pr, 1, 0)
#' X <- ifelse(TRT==1, cf.m.T, cf.s.T)
#' E <- ifelse(TRT==1, cf.m.E, cf.s.E)
#' covs.names <- c("c1", "c2")
#' data <- data.frame(X=X, E=E, TRT=TRT, c1=c1, c2=c2)
#' wei <- get.weights(data, "TRT", covs.names, wtype = "overlap")
#' hist(wei$ps[data$TRT==1], col="red", breaks = seq(0,1,0.05))
#' par(new=TRUE)
#' hist(wei$ps[data$TRT==0], col="blue", breaks = seq(0,1,0.05))
#'
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
#' @param df a data frame that includes time-to-event \code{X}, type of event \code{E},
#' a treatment indicator \code{A}  and covariates \code{C}.
#' @param X a character string specifying the name of the time-to-event variable in \code{df}.
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
#' "stab.ATE" is defined as P(A=a)/P(A=a|C=c) - see Hernán et al. (2000).
#' @param cens an integer value in \code{E} that corresponds to censoring times recorded in \code{X}.
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
#' @examples
#' # create a data set
#' n <- 1000
#' set.seed(7)
#' c1 <- runif(n)
#' c2 <- as.numeric(runif(n)< 0.2)
#' set.seed(77)
#' cf.m.T1 <- rweibull(n, shape=1, scale=exp(-(-1 + 2*c1)))
#' cf.m.T2 <-  rweibull(n, shape=1, scale=exp(-(1 + 1*c2)))
#' cf.m.T <- pmin( cf.m.T1, cf.m.T2)
#' cf.m.E <- rep(0, n)
#' cf.m.E[cf.m.T1<=cf.m.T2] <- 1
#' cf.m.E[cf.m.T2<cf.m.T1] <- 2
#' set.seed(77)
#' cf.s.T1 <- rweibull(n, shape=1, scale=exp(-1*c1 ))
#' cf.s.T2 <-  rweibull(n, shape=1, scale=exp(-2*c2))
#' cf.s.T <- pmin( cf.s.T1, cf.s.T2)
#' cf.s.E <- rep(0, n)
#' cf.s.E[cf.s.T1<=cf.s.T2] <- 1
#' cf.s.E[cf.s.T2<cf.s.T1] <- 2
#' exp.z <- exp(0.5 + 1*c1 - 1*c2)
#' pr <- exp.z/(1+exp.z)
#' TRT <- ifelse(runif(n)< pr, 1, 0)
#' X <- ifelse(TRT==1, cf.m.T, cf.s.T)
#' E <- ifelse(TRT==1, cf.m.E, cf.s.E)
#' covs.names <- c("c1", "c2")
#' data <- data.frame(X=X, E=E, TRT=TRT, c1=c1, c2=c2)
#'
#' num.atrisk <- get.numAtRisk(data, "X", "E", "TRT", C=covs.names, wtype="overlap", cens=0)
#' plot(num.atrisk$trt.1$time[num.atrisk$trt.1$sample=="Weighted"],
#'      num.atrisk$trt.1$num[num.atrisk$trt.1$sample=="Weighted"], col="red", type="s",
#'      xlab="time", ylab="number at risk",
#'      main="Number at risk in TRT=1", ylim=c(0, max(num.atrisk$trt.1$num)))
#' lines(num.atrisk$trt.1$time[num.atrisk$trt.1$sample=="Unadjusted"],
#'       num.atrisk$trt.1$num[num.atrisk$trt.1$sample=="Unadjusted"], col="blue", type="s")
#' legend("topright", legend=c("Weighted", "Unadjusted"), lty=1:1,  col=c("red", "blue"))
#' plot(num.atrisk$trt.0$time[num.atrisk$trt.0$sample=="Weighted"],
#'      num.atrisk$trt.0$num[num.atrisk$trt.0$sample=="Weighted"], col="red", type="s",
#'      xlab="time", ylab="number at risk",
#'      main="Number at risk in TRT=0", ylim=c(0, max(num.atrisk$trt.0$num)))
#' lines(num.atrisk$trt.0$time[num.atrisk$trt.0$sample=="Unadjusted"],
#'       num.atrisk$trt.0$num[num.atrisk$trt.0$sample=="Unadjusted"], col="blue", type="s")
#' legend("topright", legend=c("Weighted", "Unadjusted"), lty=1:1,  col=c("red", "blue"))
#'
#' @export
get.numAtRisk <- function(df, X, E, A, C=NULL, wtype="unadj", cens=0)
{
  X <- df[[X]]
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


.sequential.fit.nonpar <- function(df, X, E, A, C, wtype, cens, conf.level, bs, nbs.rep, seed, verbose)
{
  X <- df[[X]]
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

    if (verbose) pb <- txtProgressBar(min = 1, max = nbs.rep, style=3)
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
      if (verbose) setTxtProgressBar(pb, i)
    }
    if (verbose) close(pb)

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

      res$trt.0[[paste("Ev=", k, sep="")]][["bs.CumHaz"]] <- bs.CumHaz$trt.0[[paste("Ev=", k, sep="")]]
      res$trt.1[[paste("Ev=", k, sep="")]][["bs.CumHaz"]] <- bs.CumHaz$trt.1[[paste("Ev=", k, sep="")]]
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


.sequential.fit.cox <- function(df, X, E, A, C, wtype, cens, conf.level, bs, nbs.rep, seed, verbose)
{
  X <- df[[X]]
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
    log.CumHazR <- est[[paste("Ev=", k, sep="")]]$logHR
    names(log.CumHazR) <- NULL
    res$trt.eff[[paste("Ev=", k, sep="")]] <- list(log.CumHazR=log.CumHazR,
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

    if (verbose) pb <- txtProgressBar(min = 1, max = nbs.rep, style=3)
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
      if (verbose) setTxtProgressBar(pb, i)
    }
    if (verbose) close(pb)

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
      log.CumHazR.CI.L <- quantile(bs.CumHaz$logRatio[[paste("Ev=", k, sep="")]], prob=alpha/2, na.rm=TRUE)
      names(log.CumHazR.CI.L) <- NULL
      res$trt.eff[[paste("Ev=", k, sep="")]][["log.CumHazR.CI.L"]] <- log.CumHazR.CI.L
      log.CumHazR.CI.U <- quantile(bs.CumHaz$logRatio[[paste("Ev=", k, sep="")]], prob=1-alpha/2, na.rm=TRUE)
      names(log.CumHazR.CI.U) <- NULL
      res$trt.eff[[paste("Ev=", k, sep="")]][["log.CumHazR.CI.U"]] <- log.CumHazR.CI.U

      # SE:
      res$trt.0[[paste("Ev=", k, sep="")]][["CumHaz.SE"]] <-
        apply(bs.CumHaz$trt.0[[paste("Ev=", k, sep="")]], 2, sd, na.rm=TRUE)
      res$trt.1[[paste("Ev=", k, sep="")]][["CumHaz.SE"]] <-
        apply(bs.CumHaz$trt.1[[paste("Ev=", k, sep="")]], 2, sd, na.rm=TRUE)
      SD <- sd(bs.CumHaz$logRatio[[paste("Ev=", k, sep="")]], na.rm=TRUE)
      res$trt.eff[[paste("Ev=", k, sep="")]][["log.CumHazR.SE"]] <- SD

      # pvalue for log.CumHazR
      pvalue <- 2*pnorm(-abs(res$trt.eff[[paste("Ev=", k, sep="")]][["log.CumHazR"]])/SD)
      res$trt.eff[[paste("Ev=", k, sep="")]][["log.CumHazR.pvalue"]] <- pvalue

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

      res$trt.0[[paste("Ev=", k, sep="")]][["bs.CumHaz"]] <- bs.CumHaz$trt.0[[paste("Ev=", k, sep="")]]
      res$trt.1[[paste("Ev=", k, sep="")]][["bs.CumHaz"]] <- bs.CumHaz$trt.1[[paste("Ev=", k, sep="")]]
    }
  }
  class(res) <- "cmprsk"
  return(res)
}

#' Returns point estimates and \code{conf.level}\% confidence intervals corresponding to a specific time point
#'
#' The confidence interval returned by this function corresponds to the value \code{conf.level} passed to the function
#' \code{fit.cox} or \code{fit.nonpar}. The first input argument \code{cmprsk.obj} is a result corresponding to \code{conf.level}.
#'
#' @param cmprsk.obj a \code{cmprsk} object returned by one of the functions \code{fit.cox} or \code{fit.nonpar}
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
#' \item{\code{CumHaz}} a point estimate of the cumulative hazard
#' \item{\code{CumHaz.CI.L}} a bootstrap-based quantile estimate of a lower bound of a \code{conf.level}\% confidence
#' interval for the cumulative hazard
#' \item{\code{CumHaz.CI.U}} a bootstrap-based quantile estimate of an upper bound of a \code{conf.level}\% confidence
#' interval for the cumulative hazard
#' \item{\code{CIF}} a point estimate of the cumulative incidence function
#' \item{\code{CIF.CI.L}} a bootstrap-based quantile estimate of a lower bound of a \code{conf.level}\% confidence
#' interval for the  cumulative incidence function
#' \item{\code{CIF.CI.U}} a bootstrap-based quantile estimate of an upper bound of a \code{conf.level}\% confidence
#' interval for the cumulative incidence function
#' \item{\code{RMT}} a point estimate of the restricted mean time
#' \item{\code{RMT.CI.L}} a bootstrap-based quantile estimate of a lower bound of a \code{conf.level}\% confidence
#' interval for the restricted mean time
#' \item{\code{RMT.CI.U}} a bootstrap-based quantile estimate of an upper bound of a \code{conf.level}\% confidence
#' interval for the restricted mean time}
#' \tabular{ll}{
#' \code{trt.1} \tab \cr a list of estimates of the absolute counterfactual parameters
#' corresponding to \code{A}=1 and the type of event \code{E}. \code{trt.1} has the number of
#'  fields as the number of different types of events  in the data set.
#' For each type of event there is a list of estimates:  \cr }
#' \itemize{
#' \item{\code{CumHaz}} a point estimate of the cumulative hazard
#' \item{\code{CumHaz.CI.L}} a bootstrap-based quantile estimate of a lower bound of a \code{conf.level}\% confidence
#' interval for the cumulative hazard
#' \item{\code{CumHaz.CI.U}} a bootstrap-based quantile estimate of an upper bound of a \code{conf.level}\% confidence
#' interval for the cumulative hazard
#' \item{\code{CIF}} a point estimate of the cumulative incidence function
#' \item{\code{CIF.CI.L}} a bootstrap-based quantile estimate of a lower bound of a \code{conf.level}\% confidence
#' interval for the  cumulative incidence function
#' \item{\code{CIF.CI.U}} a bootstrap-based quantile estimate of an upper bound of a \code{conf.level}\% confidence
#' interval for the cumulative incidence function
#' \item{\code{RMT}} a point estimate of the restricted mean time
#' \item{\code{RMT.CI.L}} a bootstrap-based quantile estimate of a lower bound of a \code{conf.level}\% confidence
#' interval for the restricted mean time
#' \item{\code{RMT.CI.U}} a bootstrap-based quantile estimate of an upper bound of a \code{conf.level}\% confidence
#' interval for the restricted mean time}
#' \tabular{ll}{
#' \code{trt.eff} \tab \cr a list of estimates of the treatment effect measures
#' corresponding to the type of event \code{E}. \code{trt.eff} has the number of
#'  fields as the number of different types of events  in the data set.
#' For each type of event there is a list of estimates: \cr}
#' \itemize{
#' \item{\code{log.CumHazR}} a point estimate of the log of the ratio of hazards between two treatment arms
#' \item{\code{log.CumHazR.CI.L}} a bootstrap-based quantile estimate of a lower bound of a \code{conf.level}\% confidence
#' interval for the log of the ratio of hazards between two treatment arms
#' \item{\code{log.CumHazR.CI.U}} a bootstrap-based quantile estimate of an upper bound of a \code{conf.level}\% confidence
#' interval for the log of the ratio of hazards between two treatment arms
#' \item{\code{RD}} a point estimate of the risk difference between two treatment arms
#' \item{\code{RD.CI.L}} a bootstrap-based quantile estimate of a lower bound of a \code{conf.level}\% confidence
#' interval for the risk difference between two treatment arms
#' \item{\code{RD.CI.U}} a bootstrap-based quantile estimate of an upper bound of a \code{conf.level}\% confidence
#' interval for the risk difference between two treatment arms
#' \item{\code{RR}} a point estimate of the risk ratio between two treatment arms
#' \item{\code{RR.CI.L}} a bootstrap-based quantile estimate of a lower bound of a \code{conf.level}\% confidence
#' interval for the risk ratio between two treatment arms
#' \item{\code{RR.CI.U}} a bootstrap-based quantile estimate of an upper bound of a \code{conf.level}\% confidence
#' interval for the risk ratio between two treatment arms
#' \item{\code{ATE.RMT}} a point estimate of the restricted mean time difference
#'  between two treatment arms
#' \item{\code{ATE.RMT.CI.L}} a bootstrap-based quantile estimate of a lower bound of a \code{conf.level}\% confidence
#' interval for the restricted mean time difference between two treatment arms
#' \item{\code{ATE.RMT.CI.U}} a bootstrap-based quantile estimate of an upper bound of a \code{conf.level}\% confidence
#' interval for the restricted mean time difference between two treatment arms}
#'
#' @seealso \code{\link{fit.cox}}, \code{\link{fit.nonpar}}, \code{\link{causalCmprsk}}
#' @examples
#' # create a data set
#' n <- 1000
#' set.seed(7)
#' c1 <-  runif(n)
#' c2 <- as.numeric(runif(n)< 0.2)
#' set.seed(77)
#' cf.m.T1 <- rweibull(n, shape=1, scale=exp(-(-1 + 2*c1)))
#' cf.m.T2 <-  rweibull(n, shape=1, scale=exp(-(1 + 1*c2)))
#' cf.m.T <- pmin( cf.m.T1, cf.m.T2)
#' cf.m.E <- rep(0, n)
#' cf.m.E[cf.m.T1<=cf.m.T2] <- 1
#' cf.m.E[cf.m.T2<cf.m.T1] <- 2
#' set.seed(77)
#' cf.s.T1 <- rweibull(n, shape=1, scale=exp(-1*c1 ))
#' cf.s.T2 <-  rweibull(n, shape=1, scale=exp(-2*c2))
#' cf.s.T <- pmin( cf.s.T1, cf.s.T2)
#' cf.s.E <- rep(0, n)
#' cf.s.E[cf.s.T1<=cf.s.T2] <- 1
#' cf.s.E[cf.s.T2<cf.s.T1] <- 2
#' exp.z <- exp(0.5 + 1*c1 - 1*c2)
#' pr <- exp.z/(1+exp.z)
#' TRT <- ifelse(runif(n)< pr, 1, 0)
#' X <- ifelse(TRT==1, cf.m.T, cf.s.T)
#' E <- ifelse(TRT==1, cf.m.E, cf.s.E)
#' covs.names <- c("c1", "c2")
#' data <- data.frame(X=X, E=E, TRT=TRT, c1=c1, c2=c2)
#' wei <- get.weights(data, "TRT", covs.names, wtype = "overlap")
#' hist(wei$ps[data$TRT==1], col="red", breaks = seq(0,1,0.05))
#' par(new=TRUE)
#' hist(wei$ps[data$TRT==0], col="blue", breaks = seq(0,1,0.05))
#' # Nonparametric estimation:
#' res.ATE <- fit.nonpar(df=data, X="X", E="E", A="TRT", C=covs.names, wtype="stab.ATE")
#' nonpar.pe <- get.pointEst(res.ATE, 0.5)
#' nonpar.pe$trt.eff[[1]]$RD
#' # Cox-based estimation:
#' res.cox.ATE <- fit.cox(df=data, X="X", E="E", A="TRT", C=covs.names, wtype="stab.ATE")
#' cox.pe <- get.pointEst(res.cox.ATE, 0.5)
#' cox.pe$trt.eff[[1]]$RD
#'
#' # please see our package vignette for practical examples
#'
#' @export
get.pointEst <- function(cmprsk.obj, timepoint) # assumes timepoint is a scalar
{
  if (! methods::is(cmprsk.obj, "cmprsk"))
  {
    stop("Error. 'cmprsk.obj' is not of 'cmprsk' type.\n")
  }
  if (length(timepoint)>1)
  {
    stop("Error. 'timepoint' is supposed to be a scalar.\n")
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

.nonpar.run <- function(df, X, E, A, C, wtype, cens, E.set,time,trt,nobs, case.w)
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

.cox.run <- function(df, X, E, A, C, wtype, cens, E.set,time,trt,nobs,case.w)
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

    log.CumHazR=est[[paste("Ev=", k, sep="")]]$logHR
    names(log.CumHazR) <- NULL
    res$trt.eff[[paste("Ev=", k, sep="")]] <- list(log.CumHazR=log.CumHazR,
                                                   RD=RD, RR=RR,
                                                   ATE.RMT=ATE.RMT)
  } # res is a list with 4 fields: time, trt.0, trt.0, trt.eff
  res
}

# Parallel Cox Fit
.parallel.fit.cox <- function(df, X, E, A, C, wtype, cens, conf.level, bs, nbs.rep, seed, verbose)
{
  X <- df[[X]]
  E <- df[[E]]
  nobs <- length(X)
  trt <- df[[A]]
  time <- sort(unique(X[E!=cens]))
  E.set <- sort(unique(E))
  E.set <- E.set[E.set!=cens]

  #  res <- list(time=time)
  res <- .cox.run(df, X, E, A, C, wtype, cens, E.set,time, trt,nobs, case.w = rep(1,nobs))
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
    if (verbose) pb <- txtProgressBar(min = 1, max = nbs.rep, style=3)
    bs_aggregates <- foreach(i = 1:nbs.rep, .export=c(".cox.run", "get.weights", ".estimate.cox",
                                                      ".base.haz.std",
                                                      ".get.CIF",".get.S", ".get.RMT"),
                             .packages=c("survival")) %dopar%
      {
        set.seed(bs_seeds[i])
        bs.w <- pmin(rexp(nobs,1), 5) # nobs = our sample size
        bs.w <- bs.w/mean(bs.w)
        bs_aggregates <- .cox.run(df, X, E, A, C, wtype, cens, E.set,time,trt,nobs, case.w = bs.w)
        if (verbose) setTxtProgressBar(pb, i)
        bs_aggregates # by default the results are combined in a list
      }
    if (verbose) close(pb)
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
      log.CumHazR.CI.L <- quantile(bs.CumHaz$logRatio[[paste("Ev=", k, sep="")]], prob=alpha/2, na.rm=TRUE)
      names(log.CumHazR.CI.L) <- NULL
      res$trt.eff[[paste("Ev=", k, sep="")]][["log.CumHazR.CI.L"]] <- log.CumHazR.CI.L
      log.CumHazR.CI.U <- quantile(bs.CumHaz$logRatio[[paste("Ev=", k, sep="")]], prob=1-alpha/2, na.rm=TRUE)
      names(log.CumHazR.CI.U) <- NULL
      res$trt.eff[[paste("Ev=", k, sep="")]][["log.CumHazR.CI.U"]] <- log.CumHazR.CI.U

      # SE:
      res$trt.0[[paste("Ev=", k, sep="")]][["CumHaz.SE"]] <-
        apply(bs.CumHaz$trt.0[[paste("Ev=", k, sep="")]], 2, sd, na.rm=TRUE)
      res$trt.1[[paste("Ev=", k, sep="")]][["CumHaz.SE"]] <-
        apply(bs.CumHaz$trt.1[[paste("Ev=", k, sep="")]], 2, sd, na.rm=TRUE)
      SD <- sd(bs.CumHaz$logRatio[[paste("Ev=", k, sep="")]], na.rm=TRUE)
      res$trt.eff[[paste("Ev=", k, sep="")]][["log.CumHazR.SE"]] <- SD

      # pvalue
      pvalue <- 2*pnorm(-abs(res$trt.eff[[paste("Ev=", k, sep="")]][["log.CumHazR"]])/SD)
      res$trt.eff[[paste("Ev=", k, sep="")]][["log.CumHazR.pvalue"]] <- pvalue

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

      res$trt.0[[paste("Ev=", k, sep="")]][["bs.CumHaz"]] <- bs.CumHaz$trt.0[[paste("Ev=", k, sep="")]]
      res$trt.1[[paste("Ev=", k, sep="")]][["bs.CumHaz"]] <- bs.CumHaz$trt.1[[paste("Ev=", k, sep="")]]
    }
    class(res) <- "cmprsk"
    return(res)
  }
}

.parallel.fit.nonpar <- function(df, X, E, A, C, wtype, cens, conf.level, bs, nbs.rep, seed, verbose)
{
  X <- df[[X]]
  E <- df[[E]]
  nobs <- length(X)
  trt <- df[[A]]
  time <- sort(unique(X[E!=cens]))
  E.set <- sort(unique(E))
  E.set <- E.set[E.set!=cens]

  #  res <- list(time=time) # this is the only place where I'll keep the time
  res <- .nonpar.run(df, X, E, A, C, wtype, cens, E.set,time,trt,nobs, case.w = rep(1,nobs))

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
    if (verbose) pb <- txtProgressBar(min = 1, max = nbs.rep, style=3, file=stdout())
    bs_aggregates <- foreach(i = 1:nbs.rep, .export=c(".nonpar.run", "get.weights",
                                                      ".estimate.nonpar", ".base.haz.std",
                                                      ".get.CIF",".get.S", ".get.RMT"),
                             .packages=c("survival")) %dopar%
      {
        set.seed(bs_seeds[i])
        bs.w <- pmin(rexp(nobs,1), 5) # nobs = our sample size
        bs.w <- bs.w/mean(bs.w)
        bs_aggregates <- .nonpar.run(df, X, E, A, C, wtype, cens, E.set,time,trt,nobs, case.w = bs.w)
        if (verbose) setTxtProgressBar(pb, i)
        bs_aggregates
      }          # df, T, E, A, C, wtype, cens, E.set,time,trt,nobs,X,case.w
    if (verbose) close(pb)

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

      res$trt.0[[paste("Ev=", k, sep="")]][["bs.CumHaz"]] <- bs.CumHaz$trt.0[[paste("Ev=", k, sep="")]]
      res$trt.1[[paste("Ev=", k, sep="")]][["bs.CumHaz"]] <- bs.CumHaz$trt.1[[paste("Ev=", k, sep="")]]
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
#' @param df a data frame that includes time-to-event \code{X}, type of event \code{E},
#' a treatment indicator \code{A}  and covariates \code{C}.
#' @param X a character string specifying the name of the time-to-event variable in \code{df}.
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
#' in data from a randomized controlled trial (RCT) where there is no need for emulation of baseline randomization.
#' Other possible values are "stab.ATE", "ATE", "ATT", "ATC" and "overlap".
#' See Table 1 from Li, Morgan, and Zaslavsky (2018).
#' "stab.ATE" is defined as P(A=a)/P(A=a|C=c) - see Hernán et al. (2000).
#' @param cens an integer value in \code{E} that corresponds to censoring times recorded in \code{X}.
#' By default \code{fit.nonpar} assumes \code{cens}=0
#' @param conf.level the confidence level that will be used in the bootstrap confidence intervals.
#' The default is 0.95
#' @param bs a logical flag indicating whether to perform bootstrap in order to obtain confidence intervals. There are no
#' analytical confidence intervals in \code{fit.nonpar}
#' @param nbs.rep number of bootstrap replications
#' @param seed the random seed for the bootstrap, in order to make the results reproducible
#' @param parallel a logical flag indicating whether to perform bootstrap sequentially or in parallel,
#' using several cores simultaneously. The default value is FALSE. In parallel execution, the number
#' of available cores is detected, and the parallel jobs are assigned to the number of
#' detected available cores minus one.
#' @param verbose a logical flag indicating whether to show a progress of bootstrap.
#' The progress bar is shown only for sequential bootstrap computation.
#' The default value is FALSE.
#'
#' @return  A list of class \code{cmprsk} with the following fields:
#' \tabular{ll}{
#' \code{time}  \tab \cr a vector of time points for which all the parameters are estimated \cr
#' \code{trt.0} \tab \cr a list of estimates of the counterfactual parameters
#' corresponding to \code{A}=0 and the type of event \code{E}. \code{trt.0}
#' has K
#'  fields as the number of competing events in the data set.
#' For each competing risk there is a list of point estimates, their standard errors and
#' \code{conf.level}\% confidence intervals:  \cr }
#' \itemize{
#' \item{\code{CumHaz}} a vector of cumulative hazard estimates
#' \item{\code{CIF}} a vector of cumulative incidence functions (CIF)
#' \item{\code{RMT}} a vector of restricted mean time (RMT) estimates
#' \item{\code{CumHaz.CI.L}} a vector of bootstrap-based quantile estimate of
#'  lower confidence limits for cumulative hazard estimates
#' \item{\code{CumHaz.CI.U}} a vector of bootstrap-based quantile estimate of
#'  upper confidence limits for cumulative hazard estimates
#' \item{\code{CumHaz.SE}} a vector of the bootstrap-based estimated standard errors
#' of cumulative hazard estimates
#' \item{\code{CIF.CI.L}} a vector of bootstrap-based quantile estimate of
#'  lower confidence limits for CIF estimates
#' \item{\code{CIF.CI.U}} a vector of bootstrap-based quantile estimate of
#'  upper confidence limits for CIF estimates
#' \item{\code{CIF.SE}} a vector of bootstrap-based estimated standard error
#' of CIF estimates
#' \item{\code{RMT.CI.L}} a vector of bootstrap-based quantile estimate of
#'  lower confidence limits for RMT estimates
#' \item{\code{RMT.CI.U}} a vector of bootstrap-based quantile estimate of
#'  upper confidence limits for RMT estimates
#' \item{\code{RMT.SE}} a vector of the bootstrap-based estimated standard errors
#' of RMT estimates
#' \item{\code{bs.CumHaz}} a matrix of dimension \code{nbs.rep} by the length of \code{time} vector,
#' with cumulative hazard estimates for \code{nbs.rep} bootstrap samples
#' }
#' \tabular{ll}{
#' \code{trt.1} \tab \cr a list of estimates of the counterfactual parameters
#' corresponding to \code{A}=1 and the type of event \code{E}. \code{trt.1} has K
#'  fields as the number of competing events (risks) in the data set.
#' For each competing risk there is a list of point estimates:  \cr }
#' \itemize{
#' \item{\code{CumHaz}} a vector of cumulative hazard estimates
#' \item{\code{CIF}} a vector of cumulative incidence functions
#' \item{\code{RMT}} a vector of restricted mean time estimates
#' \item{\code{CumHaz.CI.L}} a vector of bootstrap-based quantile estimate of
#'  lower confidence limits for cumulative hazard estimates
#' \item{\code{CumHaz.CI.U}} a vector of bootstrap-based quantile estimate of
#'  upper confidence limits for cumulative hazard estimates
#' \item{\code{CumHaz.SE}} a vector of the bootstrap-based estimated standard errors
#' of cumulative hazard estimates
#' \item{\code{CIF.CI.L}} a vector of bootstrap-based quantile estimate of
#'  lower confidence limits for CIF estimates
#' \item{\code{CIF.CI.U}} a vector of bootstrap-based quantile estimate of
#'  upper confidence limits for CIF estimates
#' \item{\code{CIF.SE}} a vector of bootstrap-based estimated standard error
#' for CIF estimates
#' \item{\code{RMT.CI.L}} a vector of bootstrap-based quantile estimate of
#'  lower confidence limits for RMT estimates
#' \item{\code{RMT.CI.U}} a vector of bootstrap-based quantile estimate of
#'  upper confidence limits for RMT estimates
#' \item{\code{RMT.SE}} a vector of the bootstrap-based estimated standard errors
#' of the RMT estimates
#' \item{\code{bs.CumHaz}} a matrix of dimension \code{nbs.rep} by the length of \code{time} vector,
#' with cumulative hazard estimates for \code{nbs.rep} bootstrap samples
#' }
#' \tabular{ll}{
#' \code{trt.eff} \tab \cr a list of estimates of the treatment effect measures
#' corresponding to the type of event \code{E}. \code{trt.eff} has the number of
#'  fields as the number of different types of events (risks)  in the data set.
#' For each competing risk there is a list of estimates: \cr}
#' \itemize{
#' \item{\code{log.CumHazR}} an estimate of the log of the hazard ratio.
#' It is a scalar since the Cox model is assumed.
#' \item{\code{RD}} a vector of time-varying Risk Difference between two treatment arms
#' \item{\code{RR}} a vector of time-varying Risk Ratio between two treatment arms
#' \item{\code{ATE.RMT}} a vector of the time-varying Restricted Mean Time Difference
#'  between two treatment arms
#' \item{\code{log.CumHazR.CI.L}} a bootstrap-based quantile estimate of the
#'  lower confidence limit of \code{log.CumHazR}
#' \item{\code{log.CumHazR.CI.U}} a bootstrap-based quantile estimate of the
#'  upper confidence limit of \code{log.CumHazR}
#' \item{\code{log.CumHazR.SE}} a bootstrap-based estimated standard error
#' of \code{log.CumHazR}
#' \item{\code{log.CumHazR.pvalue}} p-value from a Wald test of a two-sided hypothesis H0: HR(A=1)/HR(A=0)=1
#' \item{\code{RD.CI.L}} a vector of bootstrap-based quantile estimates of the
#'  lower confidence limits of the Risk Difference estimates
#' \item{\code{RD.CI.U}} a vector of bootstrap-based quantile estimate of the
#'  upper confidence limits of the Risk Difference estimates
#' \item{\code{RD.SE}} a vector of the bootstrap-based estimated standard errors
#' of the Risk Difference
#' \item{\code{RR.CI.L}} a vector of bootstrap-based quantile estimates of the
#'  lower confidence limits of the Risk Ratio estimates
#' \item{\code{RR.CI.U}} a vector of bootstrap-based quantile estimate of the
#'  upper confidence limits of the Risk Ratio estimates
#' \item{\code{RR.SE}} a vector of the bootstrap-based estimated standard errors
#' of the Risk Ratio
#' \item{\code{ATE.RMT.CI.L}} a vector of bootstrap-based quantile estimate of
#'  lower confidence limits for the RMT difference estimates
#' \item{\code{ATE.RMT.CI.U}} a vector of bootstrap-based quantile estimate of
#'  upper confidence limits for the RMT difference estimates
#' \item{\code{ATE.RMT.SE}} a vector of bootstrap-based estimated standard errors
#' of the RMT difference estimates
#'  }
#'
#' @seealso \code{\link{fit.nonpar}}, \code{\link{get.pointEst}}, \code{\link{causalCmprsk}}
#'
#' @examples
#' # create a data set
#' n <- 1000
#' set.seed(7)
#' c1 <- runif(n)
#' c2 <- as.numeric(runif(n)< 0.2)
#' set.seed(77)
#' cf.m.T1 <- rweibull(n, shape=1, scale=exp(-(-1 + 2*c1)))
#' cf.m.T2 <-  rweibull(n, shape=1, scale=exp(-(1 + 1*c2)))
#' cf.m.T <- pmin( cf.m.T1, cf.m.T2)
#' cf.m.E <- rep(0, n)
#' cf.m.E[cf.m.T1<=cf.m.T2] <- 1
#' cf.m.E[cf.m.T2<cf.m.T1] <- 2
#' set.seed(77)
#' cf.s.T1 <- rweibull(n, shape=1, scale=exp(-1*c1 ))
#' cf.s.T2 <-  rweibull(n, shape=1, scale=exp(-2*c2))
#' cf.s.T <- pmin( cf.s.T1, cf.s.T2)
#' cf.s.E <- rep(0, n)
#' cf.s.E[cf.s.T1<=cf.s.T2] <- 1
#' cf.s.E[cf.s.T2<cf.s.T1] <- 2
#' exp.z <- exp(0.5 + 1*c1 - 1*c2)
#' pr <- exp.z/(1+exp.z)
#' TRT <- ifelse(runif(n)< pr, 1, 0)
#' X <- ifelse(TRT==1, cf.m.T, cf.s.T)
#' E <- ifelse(TRT==1, cf.m.E, cf.s.E)
#' covs.names <- c("c1", "c2")
#' data <- data.frame(X=X, E=E, TRT=TRT, c1=c1, c2=c2)
#' wei <- get.weights(data, "TRT", covs.names, wtype = "overlap")
#' hist(wei$ps[data$TRT==1], col="red", breaks = seq(0,1,0.05))
#' hist(wei$ps[data$TRT==0], col="blue", breaks = seq(0,1,0.05))
#' # Cox-based estimation:
#' res.cox.ATE <- fit.cox(df=data, X="X", E="E", A="TRT", C=covs.names, wtype="stab.ATE")
#' cox.pe <- get.pointEst(res.cox.ATE, 0.5)
#' cox.pe$trt.eff[[1]]$RD
#'
#' # please see our package vignette for practical examples
#'
#' @references F. Li, K.L. Morgan, and A.M. Zaslavsky. 2018. Balancing Covariates via Propensity Score Weighting. Journal of the American Statistical Association, 113 (521): 390–400.
#' @references M.A. Hernán, B. Brumback, and J.M. Robins. 2000. Marginal structural models and to estimate the causal effect of zidovudine on the survival of HIV-positive men. Epidemiology, 11 (5): 561-570.
#'
#'
#' @export
fit.cox <- function(df, X, E, A, C=NULL, wtype="unadj", cens=0, conf.level=0.95, bs=FALSE, nbs.rep=400, seed=17, parallel = FALSE, verbose=FALSE){

  get_os <- function(){
    platform <- .Platform$OS.type
    return(platform)
  }
  if(bs == TRUE && parallel == TRUE){

    # if windows (snow)
    if(grepl('win',get_os(),ignore.case = TRUE)) {
      #      num_cores <- round(parallel::detectCores()/2,digits = 0)
      chk <- Sys.getenv("_R_CHECK_LIMIT_CORES_", "")
      if (nzchar(chk) && chk == "TRUE")
        { # use 2 cores in CRAN/Travis/AppVeyor
        num_cores <- 2L
        }
      else
        { # use all cores in devtools::test()
          num_cores <- parallel::detectCores() - 1
      }
      cluster <- parallel::makeCluster(spec=num_cores)
      doParallel::registerDoParallel(cl = cluster)
    }

    # if unix (multicore)
    if(grepl('unix|linux',get_os(),ignore.case = TRUE)){
      num_cores <- 2 # https://stackoverflow.com/questions/41307178/error-processing-vignette-failed-with-diagnostics-4-simultaneous-processes-spa
      cluster <- parallel::makeCluster(spec=num_cores)
      doParallel::registerDoParallel(cl=cluster)
    }

    res <- .parallel.fit.cox(df = df, X = X, E = E, A = A, C = C, wtype = wtype, cens = cens, conf.level = conf.level, bs = bs, nbs.rep = nbs.rep, seed = seed, verbose=verbose)
    parallel::stopCluster(cluster)
  }
  else
  res <- .sequential.fit.cox(df = df, X = X, E = E, A = A, C = C, wtype = wtype, cens = cens, conf.level = conf.level, bs = bs, nbs.rep = nbs.rep, seed = seed, verbose=verbose)
  res
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
#' @param df a data frame that includes time-to-event \code{X}, type of event \code{E},
#' a treatment indicator \code{A}  and covariates \code{C}.
#' @param X a character string specifying the name of the time-to-event variable in \code{df}.
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
#' "stab.ATE" is defined as P(A=a)/P(A=a|C=c) - see Hernán et al. (2000).
#' @param cens an integer value in \code{E} that corresponds to censoring times recorded in \code{X}.
#' By default \code{fit.nonpar} assumes \code{cens}=0
#' @param conf.level the confidence level that will be used in the bootstrap confidence intervals.
#' The default is 0.95
#' @param bs a logical flag indicating whether to perform bootstrap in order to obtain confidence intervals. There are no
#' analytical confidence intervals in \code{fit.nonpar}
#' @param nbs.rep number of bootstrap replications
#' @param seed the random seed for the bootstrap, in order to make the results reproducible
#' @param parallel a logical flag indicating whether to perform bootstrap sequentially or in parallel,
#' using several cores simultaneously. The default value is FALSE. In parallel execution, the number
#' of available cores is detected, and the parallel jobs are assigned to the number of
#' detected available cores minus one.
#' @param verbose a logical flag indicating whether to show a progress of bootstrap.
#' The progress bar is shown only for sequential bootstrap computation.
#' The default value is FALSE.
#'
#' @return  A list of class \code{cmprsk} with the following fields:
#' \tabular{ll}{
#' \code{time}  \tab \cr a vector of time points for which all the parameters are estimated \cr
#' \code{trt.0} \tab \cr a list of estimates of the absolute counterfactual parameters
#' corresponding to \code{A}=0 and the type of event \code{E}. \code{trt.0} has the number of
#'  fields as the number of different types of events in the data set.
#' For each type of event there is a list of estimates:  \cr }
#' \itemize{
#' \item{\code{CumHaz}} a vector of cumulative hazard estimates
#' \item{\code{CIF}} a vector of cumulative incidence functions (CIF)
#' \item{\code{RMT}} a vector of restricted mean time (RMT) estimates
#' \item{\code{CumHaz.CI.L}} a vector of bootstrap-based quantile estimate of
#'  lower confidence limits for cumulative hazard estimates
#' \item{\code{CumHaz.CI.U}} a vector of bootstrap-based quantile estimate of
#'  upper confidence limits for cumulative hazard estimates
#' \item{\code{CumHaz.SE}} a vector of the bootstrap-based estimated standard errors
#' of cumulative hazard estimates
#' \item{\code{CIF.CI.L}} a vector of bootstrap-based quantile estimate of
#'  lower confidence limits for CIF estimates
#' \item{\code{CIF.CI.U}} a vector of bootstrap-based quantile estimate of
#'  upper confidence limits for CIF estimates
#' \item{\code{CIF.SE}} a vector of bootstrap-based estimated standard error
#' of CIF estimates
#' \item{\code{RMT.CI.L}} a vector of bootstrap-based quantile estimate of
#'  lower confidence limits for RMT estimates
#' \item{\code{RMT.CI.U}} a vector of bootstrap-based quantile estimate of
#'  upper confidence limits for RMT estimates
#' \item{\code{RMT.SE}} a vector of the bootstrap-based estimated standard errors
#' of RMT estimates
#' \item{\code{bs.CumHaz}} a matrix of dimension \code{nbs.rep} by the length of \code{time} vector,
#' with cumulative hazard estimates for \code{nbs.rep} bootstrap samples
#' }
#' \tabular{ll}{
#' \code{trt.1} \tab \cr a list of estimates of the absolute counterfactual parameters
#' corresponding to \code{A}=1 and the type of event \code{E}. \code{trt.1} has the number of
#'  fields as the number of different types of events  in the data set.
#' For each type of event there is a list of estimates:  \cr }
#' \itemize{
#' \item{\code{CumHaz}} a vector of cumulative hazard estimates
#' \item{\code{CIF}} a vector of cumulative incidence functions
#' \item{\code{RMT}} a vector of restricted mean time estimates
#' \item{\code{CumHaz.CI.L}} a vector of bootstrap-based quantile estimate of
#'  lower confidence limits for cumulative hazard estimates
#' \item{\code{CumHaz.CI.U}} a vector of bootstrap-based quantile estimate of
#'  upper confidence limits for cumulative hazard estimates
#' \item{\code{CumHaz.SE}} a vector of the bootstrap-based estimated standard errors
#' of cumulative hazard estimates
#' \item{\code{CIF.CI.L}} a vector of bootstrap-based quantile estimate of
#'  lower confidence limits for CIF estimates
#' \item{\code{CIF.CI.U}} a vector of bootstrap-based quantile estimate of
#'  upper confidence limits for CIF estimates
#' \item{\code{CIF.SE}} a vector of bootstrap-based estimated standard error
#' for CIF estimates
#' \item{\code{RMT.CI.L}} a vector of bootstrap-based quantile estimate of
#'  lower confidence limits for RMT estimates
#' \item{\code{RMT.CI.U}} a vector of bootstrap-based quantile estimate of
#'  upper confidence limits for RMT estimates
#' \item{\code{RMT.SE}} a vector of the bootstrap-based estimated standard errors
#' of the RMT estimates
#' \item{\code{bs.CumHaz}} a matrix of dimension \code{nbs.rep} by the length of \code{time} vector,
#' with cumulative hazard estimates for \code{nbs.rep} bootstrap samples
#' }
#' \tabular{ll}{
#' \code{trt.eff} \tab \cr a list of estimates of the treatment effect measures
#' corresponding to the type of event \code{E}. \code{trt.eff} has the number of
#'  fields as the number of different types of events  in the data set.
#' For each type of event there is a list of estimates: \cr}
#' \itemize{
#' \item{\code{log.CumHazR}} a vector of the log of the time-varying ratio of hazards in two treatment arms
#' \item{\code{RD}} a vector of time-varying Risk Difference between two treatment arms
#' \item{\code{RR}} a vector of time-varying Risk Ratio between two treatment arms
#' \item{\code{ATE.RMT}} a vector of the time-varying Restricted Mean Time Difference
#'  between two treatment arms
#' \item{\code{log.CumHazR.CI.L}} a vector of bootstrap-based quantile estimates of the
#'  lower confidence limits of \code{log.CumHazR}
#' \item{\code{log.CumHazR.CI.U}} a vector of bootstrap-based quantile estimates of the
#'  upper confidence limits of \code{log.CumHazR}
#' \item{\code{log.CumHazR.SE}} a vector of bootstrap-based estimated standard errors
#' of \code{log.CumHazR}
#' \item{\code{RD.CI.L}} a vector of bootstrap-based quantile estimates of the
#'  lower confidence limits of the Risk Difference estimates
#' \item{\code{RD.CI.U}} a vector of bootstrap-based quantile estimate of the
#'  upper confidence limits of the Risk Difference estimates
#' \item{\code{RD.SE}} a vector of the bootstrap-based estimated standard errors
#' of the Risk Difference
#' \item{\code{RR.CI.L}} a vector of bootstrap-based quantile estimates of the
#'  lower confidence limits of the Risk Ratio estimates
#' \item{\code{RR.CI.U}} a vector of bootstrap-based quantile estimate of the
#'  upper confidence limits of the Risk Ratio estimates
#' \item{\code{RR.SE}} a vector of the bootstrap-based estimated standard errors
#' of the Risk Ratio
#' \item{\code{ATE.RMT.CI.L}} a vector of bootstrap-based quantile estimate of
#'  lower confidence limits for the RMT difference estimates
#' \item{\code{ATE.RMT.CI.U}} a vector of bootstrap-based quantile estimate of
#'  upper confidence limits for the RMT difference estimates
#' \item{\code{ATE.RMT.SE}} a vector of bootstrap-based estimated standard errors
#' of the RMT difference estimates
#'  }
#'
#' @seealso \code{\link{fit.cox}}, \code{\link{get.pointEst}}, \code{\link{causalCmprsk}}
#'
#' @examples
#' # create a data set
#' n <- 1000
#' set.seed(7)
#' c1 <-  runif(n)
#' c2 <- as.numeric(runif(n)< 0.2)
#' set.seed(77)
#' cf.m.T1 <- rweibull(n, shape=1, scale=exp(-(-1 + 2*c1)))
#' cf.m.T2 <-  rweibull(n, shape=1, scale=exp(-(1 + 1*c2)))
#' cf.m.T <- pmin( cf.m.T1, cf.m.T2)
#' cf.m.E <- rep(0, n)
#' cf.m.E[cf.m.T1<=cf.m.T2] <- 1
#' cf.m.E[cf.m.T2<cf.m.T1] <- 2
#' set.seed(77)
#' cf.s.T1 <- rweibull(n, shape=1, scale=exp(-1*c1 ))
#' cf.s.T2 <-  rweibull(n, shape=1, scale=exp(-2*c2))
#' cf.s.T <- pmin( cf.s.T1, cf.s.T2)
#' cf.s.E <- rep(0, n)
#' cf.s.E[cf.s.T1<=cf.s.T2] <- 1
#' cf.s.E[cf.s.T2<cf.s.T1] <- 2
#' exp.z <- exp(0.5 + 1*c1 - 1*c2)
#' pr <- exp.z/(1+exp.z)
#' TRT <- ifelse(runif(n)< pr, 1, 0)
#' X <- ifelse(TRT==1, cf.m.T, cf.s.T)
#' E <- ifelse(TRT==1, cf.m.E, cf.s.E)
#' covs.names <- c("c1", "c2")
#' data <- data.frame(X=X, E=E, TRT=TRT, c1=c1, c2=c2)
#' wei <- get.weights(data, "TRT", covs.names, wtype = "overlap")
#' hist(wei$ps[data$TRT==1], col="red", breaks = seq(0,1,0.05))
#' hist(wei$ps[data$TRT==0], col="blue", breaks = seq(0,1,0.05))
#' # Nonparametric estimation:
#' res.ATE <- fit.nonpar(df=data, X="X", E="E", A="TRT", C=covs.names, wtype="stab.ATE")
#' nonpar.pe <- get.pointEst(res.ATE, 0.5)
#' nonpar.pe$trt.eff[[1]]$RD
#'
#' # please see our package vignette for practical examples
#'
#' @references F. Li, K.L. Morgan, and A.M. Zaslavsky. 2018. Balancing Covariates via Propensity Score Weighting. Journal of the American Statistical Association 113 (521): 390–400.
#' @references M.A. Hernán, B. Brumback, and J.M. Robins. 2000. Marginal structural models and to estimate the causal effect of zidovudine on the survival of HIV-positive men. Epidemiology, 11 (5): 561-570.
#'
#' @export
fit.nonpar <- function(df, X, E, A, C=NULL, wtype="unadj", cens=0, conf.level=0.95, bs=FALSE, nbs.rep=400, seed=17, parallel = FALSE, verbose=FALSE){

  get_os <- function(){
    platform <- .Platform$OS.type
    return(platform)
  }
  if(bs == TRUE && parallel == TRUE){

    # if windows (snow)
    if(grepl('win',get_os(),ignore.case = TRUE)) {
#      num_cores <- round(parallel::detectCores()/2,digits = 0)
      chk <- Sys.getenv("_R_CHECK_LIMIT_CORES_", "")
      if (nzchar(chk) && chk == "TRUE")
      { # use 2 cores in CRAN/Travis/AppVeyor
        num_cores <- 2L
      }
      else
      { # use all cores in devtools::test()
        num_cores <- parallel::detectCores() - 1
      }
      cluster <- parallel::makeCluster(spec=num_cores)
      doParallel::registerDoParallel(cl=cluster)
    }

    # if unix (multicore)
    if(grepl('unix|linux',get_os(),ignore.case = TRUE)){
      num_cores <- 2 # https://stackoverflow.com/questions/41307178/error-processing-vignette-failed-with-diagnostics-4-simultaneous-processes-spa
      cluster <- parallel::makeCluster(spec=num_cores)
      doParallel::registerDoParallel(cl=cluster)
    }
    res <- .parallel.fit.nonpar(df = df, X = X, E = E, A = A, C = C, wtype = wtype, cens = cens, conf.level = conf.level, bs = bs, nbs.rep = nbs.rep, seed = seed, verbose=verbose)
    parallel::stopCluster(cl=cluster)
  }
  else
    res <- .sequential.fit.nonpar(df = df, X = X, E = E, A = A, C = C, wtype = wtype, cens = cens, conf.level = conf.level, bs = bs, nbs.rep = nbs.rep, seed = seed, verbose=verbose)
  res
}

#' Summary of Event-specific Cumulative Hazards, Cumulative Incidence Functions and Various Treatment Effects
#'
#' Returns an object of class \code{data.frame} containing the summary extracted from the \code{cmprsk} object.
#'
#' @param object an object of class \code{cmprsk} (output from \code{fit.nonpar} or \code{fit.cox} functions)
#' @param event an integer number (a code) of an event of interest
#' @param estimand a character string naming the type of estimand to extract from \code{object}. \code{estimand} can be one of the following: "CumHaz" (Cumulative Hazard function),
#' "CIF" (Cumulative Incidence Function), "RMT" (Restricted Mean Time),
#' "logHR" (logarithm of the ratio of Cumulative Hazards in two treatment arms),
#' "RD" (Risk Difference, or the difference between the CIFs in two treatment arms),
#' "RR" (Risk Ratio, or the ratio of CIFs in two treatment arms),
#' "ATE.RMT" (Restricted mean time gained/lost due to treatment, or the difference between RMTs in two treatment arms).
#' The default value is "CIF".
#' @param ... This is not currently used, included for future methods.
#'
#' @return  \code{summary.cmprsk} returns a \code{data.frame} object with 7 or 6 columns:
#' the time vector, an indicator of the treatment arm
#' (if the requested \code{estimand} is one of c("logHR", "RD", "RR", "ATE.RMT"), this column is omitted),
#' an indicator of the type of event,
#' the point estimate for the requested \code{estimand}, the lower and upper bounds of the
#' confidence interval (for \code{conf.level} \% of the confidence level),
#' and the standard error of the point estimate. For example, if \code{estimand="CIF"},
#' the returned \code{data.frame} will include the following columns:
#' \code{time}, \code{TRT}, \code{Event}, \code{CIF}, \code{CIL.CIF}, \code{CIU.CIF}, \code{SE.CIF}.

#' @seealso \code{\link{fit.cox}}, \code{\link{fit.nonpar}}, \code{\link{causalCmprsk}}
#'
#' @examples
#' # create a data set
#' n <- 1000
#' set.seed(7)
#' c1 <-  runif(n)
#' c2 <- as.numeric(runif(n)< 0.2)
#' set.seed(77)
#' cf.m.T1 <- rweibull(n, shape=1, scale=exp(-(-1 + 2*c1)))
#' cf.m.T2 <-  rweibull(n, shape=1, scale=exp(-(1 + 1*c2)))
#' cf.m.T <- pmin( cf.m.T1, cf.m.T2)
#' cf.m.E <- rep(0, n)
#' cf.m.E[cf.m.T1<=cf.m.T2] <- 1
#' cf.m.E[cf.m.T2<cf.m.T1] <- 2
#' set.seed(77)
#' cf.s.T1 <- rweibull(n, shape=1, scale=exp(-1*c1 ))
#' cf.s.T2 <-  rweibull(n, shape=1, scale=exp(-2*c2))
#' cf.s.T <- pmin( cf.s.T1, cf.s.T2)
#' cf.s.E <- rep(0, n)
#' cf.s.E[cf.s.T1<=cf.s.T2] <- 1
#' cf.s.E[cf.s.T2<cf.s.T1] <- 2
#' exp.z <- exp(0.5 + c1 - c2)
#' pr <- exp.z/(1+exp.z)
#' TRT <- ifelse( runif(n)< pr, 1, 0)
#' X <- ifelse(TRT==1, cf.m.T, cf.s.T)
#' E <- ifelse(TRT==1, cf.m.E, cf.s.E)
#' covs.names <- c("c1", "c2")
#' data <- data.frame(X=X, E=E, TRT=TRT, c1=c1, c2=c2)
#' # Nonparametric estimation:
#' res.ATE <- fit.nonpar(df=data, X="X", E="E", A="TRT", C=covs.names, wtype="stab.ATE")
#' # summarizing results on the Risk Difference for event=2
#' fit.summary <- summary(object=res.ATE, event = 2, estimand="RD")
#' head(fit.summary)
#' # summarizing results on the CIFs for event=1
#' fit.summary <- summary(object=res.ATE, event = 1, estimand="CIF")
#' head(fit.summary)
#'
#' @references M.-L. Charpignon, B. Vakulenko-Lagun, B. Zheng, C. Magdamo, B. Su, K.E. Evans, S. Rodriguez, et al. 2022. Causal inference in medical records and complementary systems pharmacology for metformin drug repurposing towards dementia. Nature Communications 13:7652.
#'
#' @export
summary.cmprsk <- function(object, event, estimand="CIF", ...)
{
  # S3 method for class 'cmprsk'
  if (length(object$trt.0[[paste("Ev=", 1, sep="")]]) == 3)
  {
    if (estimand %in% c("CumHaz", "CIF", "RMT"))
    {
      df <- rbind( data.frame(time=object$time, TRT=1, Event=event,
                              CIF=object$trt.1[[paste("Ev=", event, sep="")]]$CIF,
                              RMT=object$trt.1[[paste("Ev=", event, sep="")]]$RMT,
                              CumHaz=object$trt.1[[paste("Ev=", event, sep="")]]$CumHaz
      ),
      data.frame(time=object$time, TRT=0, Event=event,
                 CIF=object$trt.0[[paste("Ev=", event, sep="")]]$CIF,
                 RMT=object$trt.0[[paste("Ev=", event, sep="")]]$RMT,
                 CumHaz=object$trt.0[[paste("Ev=", event, sep="")]]$CumHaz)
      )
      if (estimand=="CumHaz")
        df <- df[, c("time", "TRT", "Event", "CumHaz")]
      if (estimand=="CIF")
        df <- df[, c("time", "TRT", "Event", "CIF")]
      if (estimand=="RMT")
        df <- df[, c("time", "TRT", "Event", "RMT")]
    }
    else if (estimand %in% c("logHR", "RD", "RR", "ATE.RMT"))
    {

      df <- data.frame(time=object$time, Event=event,
                       logCumHazR=object$trt.eff[[paste("Ev=", event, sep="")]]$log.CumHazR,
                       RD=object$trt.eff[[paste("Ev=", event, sep="")]]$RD,
                       RR=object$trt.eff[[paste("Ev=", event, sep="")]]$RR,
                       ATE.RMT=object$trt.eff[[paste("Ev=", event, sep="")]]$ATE.RMT
      )
      if (estimand=="logHR")
        df <- df[, c("time", "Event", "logCumHazR")]
      if (estimand=="RD")
        df <- df[, c("time", "Event", "RD")]
      if (estimand=="RR")
        df <- df[, c("time", "Event", "RR")]
      if (estimand=="ATE.RMT")
        df <- df[, c("time", "Event", "ATE.RMT")]
    }
    else # error message
    {
      stop("'estimand' is not one of the following: c('CumHaz', 'CIF', 'RMT', 'logHR', 'RD', 'RR', 'ATE.RMT'). \n")
    }
  } else
  {
  if (estimand %in% c("CumHaz", "CIF", "RMT"))
  {
    df <- rbind( data.frame(time=object$time, TRT=1, Event=event,
                            CIF=object$trt.1[[paste("Ev=", event, sep="")]]$CIF,
                            RMT=object$trt.1[[paste("Ev=", event, sep="")]]$RMT,
                            CumHaz=object$trt.1[[paste("Ev=", event, sep="")]]$CumHaz,
                            CIL.CIF=object$trt.1[[paste("Ev=", event, sep="")]]$CIF.CI.L,
                            CIU.CIF=object$trt.1[[paste("Ev=", event, sep="")]]$CIF.CI.U,
                            SE.CIF=object$trt.1[[paste("Ev=", event, sep="")]]$CIF.SE,
                            CIL.RMT=object$trt.1[[paste("Ev=", event, sep="")]]$RMT.CI.L,
                            CIU.RMT=object$trt.1[[paste("Ev=", event, sep="")]]$RMT.CI.U,
                            SE.RMT=object$trt.1[[paste("Ev=", event, sep="")]]$RMT.SE,
                            CIL.CumHaz=object$trt.1[[paste("Ev=", event, sep="")]]$CumHaz.CI.L,
                            CIU.CumHaz=object$trt.1[[paste("Ev=", event, sep="")]]$CumHaz.CI.U,
                            SE.CumHaz=object$trt.1[[paste("Ev=", event, sep="")]]$CumHaz.SE
    ),
    data.frame(time=object$time, TRT=0, Event=event,
               CIF=object$trt.0[[paste("Ev=", event, sep="")]]$CIF,
               RMT=object$trt.0[[paste("Ev=", event, sep="")]]$RMT,
               CumHaz=object$trt.0[[paste("Ev=", event, sep="")]]$CumHaz,
               CIL.CIF=object$trt.0[[paste("Ev=", event, sep="")]]$CIF.CI.L,
               CIU.CIF=object$trt.0[[paste("Ev=", event, sep="")]]$CIF.CI.U,
               SE.CIF=object$trt.0[[paste("Ev=", event, sep="")]]$CIF.SE,
               CIL.RMT=object$trt.0[[paste("Ev=", event, sep="")]]$RMT.CI.L,
               CIU.RMT=object$trt.0[[paste("Ev=", event, sep="")]]$RMT.CI.U,
               SE.RMT=object$trt.0[[paste("Ev=", event, sep="")]]$RMT.SE,
               CIL.CumHaz=object$trt.0[[paste("Ev=", event, sep="")]]$CumHaz.CI.L,
               CIU.CumHaz=object$trt.0[[paste("Ev=", event, sep="")]]$CumHaz.CI.U,
               SE.CumHaz=object$trt.0[[paste("Ev=", event, sep="")]]$CumHaz.SE)
    )
    if (estimand=="CumHaz")
      df <- df[, c("time", "TRT", "Event", "CumHaz", "CIL.CumHaz", "CIU.CumHaz", "SE.CumHaz")]
    if (estimand=="CIF")
      df <- df[, c("time", "TRT", "Event", "CIF", "CIL.CIF", "CIU.CIF", "SE.CIF")]
    if (estimand=="RMT")
      df <- df[, c("time", "TRT", "Event", "RMT", "CIL.RMT", "CIU.RMT", "SE.RMT")]
  }
  else
    if (estimand %in% c("logHR", "RD", "RR", "ATE.RMT"))
  {
    df <- data.frame(time=object$time, Event=event,
                     logCumHazR=object$trt.eff[[paste("Ev=", event, sep="")]]$log.CumHazR,
                     CIL.logCumHazR=object$trt.eff[[paste("Ev=", event, sep="")]]$log.CumHazR.CI.L,
                     CIU.logCumHazR=object$trt.eff[[paste("Ev=", event, sep="")]]$log.CumHazR.CI.U,
                     SE.logCumHazR=object$trt.eff[[paste("Ev=", event, sep="")]]$log.CumHazR.SE,
                     RD=object$trt.eff[[paste("Ev=", event, sep="")]]$RD,
                     CIL.RD=object$trt.eff[[paste("Ev=", event, sep="")]]$RD.CI.L,
                     CIU.RD=object$trt.eff[[paste("Ev=", event, sep="")]]$RD.CI.U,
                     SE.RD=object$trt.eff[[paste("Ev=", event, sep="")]]$RD.SE,
                     RR=object$trt.eff[[paste("Ev=", event, sep="")]]$RR,
                     CIL.RR=object$trt.eff[[paste("Ev=", event, sep="")]]$RR.CI.L,
                     CIU.RR=object$trt.eff[[paste("Ev=", event, sep="")]]$RR.CI.U,
                     SE.RR=object$trt.eff[[paste("Ev=", event, sep="")]]$RR.SE,
                     ATE.RMT=object$trt.eff[[paste("Ev=", event, sep="")]]$ATE.RMT,
                     CIL.ATE.RMT=object$trt.eff[[paste("Ev=", event, sep="")]]$ATE.RMT.CI.L,
                     CIU.ATE.RMT=object$trt.eff[[paste("Ev=", event, sep="")]]$ATE.RMT.CI.U,
                     SE.ATE.RMT=object$trt.eff[[paste("Ev=", event, sep="")]]$ATE.RMT.SE)
    if (estimand=="logHR")
      df <- df[, c("time", "Event", "logCumHazR", "CIL.logCumHazR", "CIU.logCumHazR", "SE.logCumHazR")]
    if (estimand=="RD")
      df <- df[, c("time", "Event", "RD", "CIL.RD", "CIU.RD", "SE.RD")]
    if (estimand=="RR")
      df <- df[, c("time", "Event", "RR", "CIL.RR", "CIU.RR", "SE.RR")]
    if (estimand=="ATE.RMT")
      df <- df[, c("time", "Event", "ATE.RMT", "CIL.ATE.RMT", "CIU.ATE.RMT", "SE.ATE.RMT")]
  }
  else # error message
  {
    stop("'estimand' is not one of the following: c('CumHaz', 'CIF', 'RMT', 'logHR', 'RD', 'RR', 'ATE.RMT'). \n")
  }
  }
  return(df)
}


