% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/causalCmprsk-functions.R
\name{get.pointEst}
\alias{get.pointEst}
\title{Returns point estimates corresponding to a specific time point}
\usage{
get.pointEst(cmprsk.obj, timepoint)
}
\arguments{
\item{cmprsk.obj}{a "cmprsk" object returned by one of the function \code{fit.cox} or \code{fit.nonpar}}

\item{timepoint}{a scalar value of the time point of interest}
}
\value{
A list with the following fields:
\tabular{ll}{
\code{time}  \tab \cr a scalar timepoint passed into the function \cr
\code{trt.0} \tab \cr a list of estimates of the absolute counterfactual parameters
corresponding to \code{A}=0 and the type of event \code{E}. \code{trt.0} has the number of
 fields as the number of different types of events in the data set.
For each type of event there is a list of estimates:  \cr }
\itemize{
\item{\code{CumHaz}} point cumulative hazard estimates
\item{\code{CIF}} point cumulative incidence function estimated
\item{\code{RMT}} point estimates of the restricted mean time}
\tabular{ll}{
\code{trt.1} \tab \cr a list of estimates of the absolute counterfactual parameters
corresponding to \code{A}=1 and the type of event \code{E}. \code{trt.1} has the number of
 fields as the number of different types of events  in the data set.
For each type of event there is a list of estimates:  \cr }
\itemize{
\item{\code{CumHaz}} point cumulative hazard estimates
\item{\code{CIF}} point cumulative incidence function estimated
\item{\code{RMT}} point estimates of the restricted mean time}
\tabular{ll}{
\code{trt.eff} \tab \cr a list of estimates of the treatment effect measures
corresponding to the type of event \code{E}. \code{trt.eff} has the number of
 fields as the number of different types of events  in the data set.
For each type of event there is a list of estimates: \cr}
\itemize{
\item{\code{log.CumHazR}} a point estimate of the log of the ratio of hazards in two treatment arms
\item{\code{RD}} a point estimate of the risk difference between two treatment arms
\item{\code{RR}} a point estimate of the risk ratio between two treatment arms
\item{\code{ATE.RMT}} a point estimate of the restricted mean time difference
 between two treatment arms}
}
\description{
Returns point estimates corresponding to a specific time point
}
\examples{
# create a data set
n <- 1000
set.seed(7)
c1 <-  runif(n)
c2 <- as.numeric(runif(n)< 0.2)
set.seed(77)
cf.m.T1 <- rweibull(n, shape=1, scale=exp(-(-1 + 2*c1)))
cf.m.T2 <-  rweibull(n, shape=1, scale=exp(-(1 + 1*c2)))
cf.m.T <- pmin( cf.m.T1, cf.m.T2)
cf.m.E <- rep(0, n)
cf.m.E[cf.m.T1<=cf.m.T2] <- 1
cf.m.E[cf.m.T2<cf.m.T1] <- 2
set.seed(77)
cf.s.T1 <- rweibull(n, shape=1, scale=exp(-1*c1 ))
cf.s.T2 <-  rweibull(n, shape=1, scale=exp(-2*c2))
cf.s.T <- pmin( cf.s.T1, cf.s.T2)
cf.s.E <- rep(0, n)
cf.s.E[cf.s.T1<=cf.s.T2] <- 1
cf.s.E[cf.s.T2<cf.s.T1] <- 2
exp.z <- exp(0.5 + 1*c1 - 1*c2)
pr <- exp.z/(1+exp.z)
TRT <- ifelse(runif(n)< pr, 1, 0)
X <- ifelse(TRT==1, cf.m.T, cf.s.T)
E <- ifelse(TRT==1, cf.m.E, cf.s.E)
covs.names <- c("c1", "c2")
data <- data.frame(X=X, E=E, TRT=TRT, c1=c1, c2=c2)
wei <- get.weights(data, "TRT", covs.names, wtype = "overlap")
hist(wei$ps[data$TRT==1], col="red", breaks = seq(0,1,0.05))
par(new=TRUE)
hist(wei$ps[data$TRT==0], col="blue", breaks = seq(0,1,0.05))
# Nonparametric estimation:
res.ATE <- fit.nonpar(df=data, X="X", E="E", A="TRT", C=covs.names, wtype="stab.ATE")
nonpar.pe <- get.pointEst(res.ATE, 0.5)
nonpar.pe$trt.eff[[1]]$RD
# Cox-based estimation:
res.cox.ATE <- fit.cox(df=data, X="X", E="E", A="TRT", C=covs.names, wtype="stab.ATE")
cox.pe <- get.pointEst(res.cox.ATE, 0.5)
cox.pe$trt.eff[[1]]$RD

# please see our package vignette for practical examples

}
\seealso{
\code{\link{fit.cox}}, \code{\link{fit.nonpar}}, \code{\link{causalCmprsk}}
}
