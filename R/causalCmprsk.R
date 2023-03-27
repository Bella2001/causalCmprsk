#' Estimation of Average Treatment Effects (ATE) of Point Intervention on Time-to-Event Outcomes with Competing Risks
#'
#' The package accompanies the paper of Charpignon et al. (2022).
#' It can be applied to data with any number of competing events, including the case of only one type of event.
#' The method uses propensity scores weighting for emulation of baseline randomization.
#' The package implements different types of weights: ATE, stabilized ATE,
#' ATT, ATC and overlap weights, as described in Li et al. (2018),
#' and different treatment effect measures (hazard ratios, risk differences, risk ratios,
#'  and restricted mean time differences).
#'
#' The \pkg{causalCmprsk} package provides two main functions:
#' \code{\link{fit.cox}} that assumes Cox proportional hazards structural models for cause-specific hazards,
#'  and \code{\link{fit.nonpar}} that does not assume any model for potential outcomes.
#'  The function \code{\link{get.weights}} returns estimated weights that are aimed for
#'  emulation of a baseline randomization in observational data where the treatment was not assigned randomly, and where conditional exchangeability is assumed.
#'  The function \code{\link{get.pointEst}} extracts a point estimate corresponding to a specific time point
#'  from the time-varying functionals returned by \code{\link{fit.cox}} and \code{\link{fit.nonpar}}.
#'  The function \code{\link{get.numAtRisk}} allows to obtain the number-at-risk statistic
#'  in the raw and weighted data.

#'
#'
#' @references M.-L. Charpignon, B. Vakulenko-Lagun, B. Zheng, C. Magdamo, B. Su, K.E. Evans, S. Rodriguez, et al. 2022. Causal inference in medical records and complementary systems pharmacology for metformin drug repurposing towards dementia. Nature Communications 13:7652.
#' @references F. Li, K.L. Morgan, and A.M. Zaslavsky. 2018. Balancing Covariates via Propensity Score Weighting. Journal of the American Statistical Association 113 (521): 390–400.
#'
#' @import survival
#' @import inline
#'
#' @importFrom stats pnorm qnorm quantile sd var rexp as.formula binomial glm predict
#' @importFrom methods is
#' @importFrom utils txtProgressBar setTxtProgressBar globalVariables
#' @importFrom parallel detectCores makeCluster stopCluster
#' @importFrom foreach foreach "%dopar%"
#' @importFrom data.table rbindlist
#' @importFrom purrr map
#' @importFrom doParallel registerDoParallel
#'
#' @docType package
#' @name causalCmprsk
NULL
