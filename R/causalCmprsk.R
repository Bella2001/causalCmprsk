#' Package for Estimation of Average Treatment Effects on Time-to-Event Outcomes with Competing Events.
#'
#' The package accompanies the paper of Charpignon et al. (2020).
#' It can be applied to data with any number of competing events, including the case of only one type of event.
#' The method uses propensity scores weighting for emulation of baseline randomization.
#' The package implements different types of weights: ATE, stabilized ATE,
#' ATT, ATC and overlap weights, as described in Li et al. (2018).
#'
#' The \pkg{causalCmprsk} package provides two main functions:
#' \code{\link{fit.cox}} for Cox-based estimation
#'  and \code{\link{fit.nonpar}} that does not assume any model for potential outcomes.
#'  The illustrative examples in these functions include analysis of the
#'  Right Heart Catheterization data (Connors et al., 1996) which is available online (\url{http://biostat.mc.vanderbilt.edu/DataSets}).

#'
#'
#' @references A.F. Connors, T. Speroff, N.V. Dawson, C. Thomas, F.E. Harrell, D. Wagner, N. Desbiens, et al. 1996. The Effectiveness of Right Heart Catheterization in the Initial Care of Critically Ill Patients. Journal of the American Medical Association 276: 889–97.
#' @references M.-L. Charpignon, B. Vakulenko-Lagun, B. Zheng, C. Magdamo, B. Su, K.E. Evans, S. Rodriguez, et al. 2020. Uncovering the Links Between Metformin, Dementia and Aging Using Emulated Trials in EHR and Systems Pharmacology.
#' @references F. Li, K.L. Morgan, and A.M. Zaslavsky. 2018. Balancing Covariates via Propensity Score Weighting. Journal of the American Statistical Association 113 (521): 390–400.
#'
#' @import survival
#' @import tidyverse
#' @import inline
#' @import data.table
#' @import ggalt
#' @import cobalt
#' @import ggsci
#' @import modEvA
#' @import naniar
#' @import readxl
#' @import hrbrthemes
#' @import magrittr
#' @import summarytools
#'
#'
#' @importFrom stats pnorm qnorm quantile sd var
#' @useDynLib causalCmprsk, .registration = TRUE
#'
#' @docType package
#' @name causalCmprsk
NULL