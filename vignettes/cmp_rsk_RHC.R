## ----label= fig1, fig.cap="A single-world competing risks model with two events *discharge* and *in-hospital death*.", echo=FALSE, fig.align="center", out.width = "60%"----
    knitr::include_graphics('RHC_cmp_rsk.png')

## ----cran_packages, echo=TRUE, include=TRUE, message=FALSE, warning=FALSE-----
library(survival)
library(tidyverse)
library(ggalt) # for stepribbon
library(cobalt)
library(ggsci) # Nature Color palette
library(modEvA) # for H-L GOF test
library(naniar) # for exploring missingness 
library(knitr)
library(kableExtra)
library(summarytools)  
library(DT) 
library(Hmisc)

library(causalCmprsk)

opts_chunk$set(results = 'asis',      # Can also be set at the chunk-level
                comment = NA,
                prompt  = FALSE,
                cache   = FALSE)
st_options(plain.ascii = FALSE,        # Always use this option in Rmd documents
          style        = "rmarkdown",  # Always use this option in Rmd documents
          footnote     = NA,           # Makes html-rendered results more concise
          subtitle.emphasis = FALSE)   # Improves layout with some rmardown themes

## ----code1, echo=TRUE---------------------------------------------------------
getHdata(rhc) # loading data using the Hmisc package
rhc_raw <- rhc

## -----------------------------------------------------------------------------
miss.report <- rhc_raw %>% miss_var_summary()
kable(miss.report[1:10,]) %>%
 kable_styling(bootstrap_options = "striped", full_width = F)

## -----------------------------------------------------------------------------
rhc_cleaning1 <- rhc_raw %>%
  mutate(RHC = as.numeric(swang1 == "RHC"), 
         trt=swang1,
         E = ifelse(is.na(dthdte), 1, ifelse(dschdte==dthdte, 2, 1)), 
         T = dschdte - sadmdte,
         T.death = ifelse(is.na(dthdte), lstctdte - sadmdte, dthdte - sadmdte), # censored time to death
         D = ifelse(death =="No", 0, 1), # death indicator
         sex_Male = as.numeric(sex == "Male"),
         race_black = as.numeric(race == "black"),
         race_other = as.numeric(race == "other"),
         income_11_25K = as.numeric(income == "$11-$25k"),
         income_25_50K = as.numeric(income == "$25-$50k"),
         income_50K = as.numeric(income == "> $50k"),
         ninsclas_Private_Medicare = as.numeric(ninsclas == "Private & Medicare"),
         ninsclas_Medicare = as.numeric(ninsclas == "Medicare"),
         ninsclas_Medicare_Medicaid = as.numeric(ninsclas == "Medicare & Medicaid"),
         ninsclas_Medicaid = as.numeric(ninsclas == "Medicaid"),
         ninsclas_No_Insurance = as.numeric(ninsclas == "No insurance"),
# combine cat1 with cat2, i.e the primary disease category with the secondary disease category:         
         cat_CHF = as.numeric(cat1 == "CHF" | (!is.na(cat2))&(cat2 == "CHF") ),
         cat_Cirrhosis = as.numeric(cat1 == "Cirrhosis" | (!is.na(cat2))&(cat2 == "Cirrhosis")),
         cat_Colon_Cancer = as.numeric(cat1 == "Colon Cancer" | (!is.na(cat2))&(cat2 == "Colon Cancer")),
         cat_Coma = as.numeric(cat1 == "Coma" | (!is.na(cat2))&(cat2 == "Coma")),
         cat_COPD = as.numeric(cat1 == "COPD" | (!is.na(cat2))&(cat2 == "COPD")),
         cat_Lung_Cancer = as.numeric(cat1 == "Lung Cancer" | (!is.na(cat2))&(cat2 == "Lung Cancer")),
         cat_MOSF_Malignancy = as.numeric(cat1 == "MOSF w/Malignancy" | (!is.na(cat2))&(cat2 == "MOSF w/Malignancy")),
         cat_MOSF_Sepsis = as.numeric(cat1 == "MOSF w/Sepsis" | (!is.na(cat2))&(cat2 == "MOSF w/Sepsis")),
         dnr1_Yes = as.numeric(dnr1 == "Yes"),
         card_Yes = as.numeric(card == "Yes"),
         gastr_Yes = as.numeric(gastr == "Yes"),
         hema_Yes = as.numeric(hema == "Yes"),
         meta_Yes = as.numeric(meta == "Yes"),
         neuro_Yes = as.numeric(neuro == "Yes"),
         ortho_Yes = as.numeric(ortho == "Yes"),
         renal_Yes = as.numeric(renal == "Yes"),
         resp_Yes = as.numeric(resp == "Yes"),
         seps_Yes = as.numeric(seps == "Yes"),
         trauma_Yes = as.numeric(trauma == "Yes"),
         ca_Yes = as.numeric(ca == "Yes"),
         ca_Metastatic = as.numeric(ca == "Metastatic")
  )

# variables selection and data reordering:
rhc_full <- rhc_cleaning1 %>% 
  select(ptid, RHC, trt, T, T.death, E, D, sex_Male, 
         age, edu, race_black, race_other, income_11_25K, income_25_50K, income_50K,
         ninsclas_Private_Medicare, ninsclas_Medicare, ninsclas_Medicare_Medicaid,
         ninsclas_Medicaid, ninsclas_No_Insurance, 
         cat_CHF, cat_Cirrhosis, cat_Colon_Cancer, cat_Coma, cat_COPD, cat_Lung_Cancer, 
         cat_MOSF_Malignancy, cat_MOSF_Sepsis,
         # diagnoses:
         dnr1_Yes, card_Yes, gastr_Yes, hema_Yes, meta_Yes, neuro_Yes, ortho_Yes, renal_Yes, 
         resp_Yes, seps_Yes, trauma_Yes,
         ca_Yes, ca_Metastatic,
         # lab tests:
         wtkilo1, hrt1, meanbp1, resp1, temp1,
         aps1, das2d3pc, scoma1, 
         surv2md1, alb1, bili1, crea1, hema1, paco21, 
         pafi1, ph1, pot1, sod1, wblc1,
         # all variables with "hx" are preexisting conditions: 
         amihx, cardiohx, chfhx, chrpulhx, dementhx, 
         gibledhx, immunhx, liverhx, malighx, psychhx, 
         renalhx, transhx, 
         death, sadmdte, dschdte, dthdte, lstctdte)

## -----------------------------------------------------------------------------
# omit 1 obs with missing discharge date for the length-of-stay analysis:
rhc <- rhc_full[!is.na(rhc_raw$dschdte),]
rhc_full$T.death.30 <- pmin(30, rhc_full$T.death)
rhc_full$D.30 <- ifelse(rhc_full$T.death <=30, 1, 0)

## ---- echo=TRUE---------------------------------------------------------------
E <- as.factor(rhc$E)
levels(E) <- c("discharge", "in-hospital death")
t <- addmargins(table(E, useNA="no"))
kable(t) %>% kable_styling(bootstrap_options = "striped", full_width = F)

## ---- echo=TRUE---------------------------------------------------------------
D <- as.factor(rhc_full$D)
levels(D) <- c("censored", "died")
t <- addmargins(table(D, useNA="no"))
kable(t) %>% kable_styling(bootstrap_options = "striped", full_width = F)

## ---- echo=TRUE---------------------------------------------------------------
D.30 <- as.factor(rhc_full$D.30)
levels(D.30) <- c("censored", "died")
t <- addmargins(table(D.30, useNA="no"))
kable(t) %>% kable_styling(bootstrap_options = "striped", full_width = F)

## ---- echo=TRUE---------------------------------------------------------------
t <- addmargins(table(rhc_full$trt, useNA="no"))
kable(t) %>% kable_styling(bootstrap_options = "striped", full_width = F)

## -----------------------------------------------------------------------------
covs.names <- c("age", "sex_Male", "edu", "race_black", "race_other",
                "income_11_25K", "income_25_50K", "income_50K", 
                "ninsclas_Private_Medicare", "ninsclas_Medicare", "ninsclas_Medicare_Medicaid",
                "ninsclas_Medicaid", "ninsclas_No_Insurance",
                "cat_CHF", "cat_Cirrhosis", "cat_Colon_Cancer", "cat_Coma", "cat_COPD", "cat_Lung_Cancer", 
                "cat_MOSF_Malignancy", "cat_MOSF_Sepsis", 
                "dnr1_Yes", "wtkilo1", "hrt1", "meanbp1", 
                "resp1", "temp1", 
                "card_Yes", "gastr_Yes", "hema_Yes", "meta_Yes", "neuro_Yes", "ortho_Yes", 
                "renal_Yes", "resp_Yes", "seps_Yes", "trauma_Yes", 
                "ca_Yes", "ca_Metastatic",
                "amihx", "cardiohx", "chfhx", "chrpulhx",
                "dementhx", "gibledhx", "immunhx", "liverhx", 
                "malighx", "psychhx", "renalhx", "transhx",
                "aps1", "das2d3pc", "scoma1", "surv2md1",
                "alb1", "bili1", "crea1", "hema1", "paco21", "pafi1", 
                "ph1", "pot1", "sod1", "wblc1")


## -----------------------------------------------------------------------------
ps.stab.ATE <- get.weights(df=rhc, A="RHC", C=covs.names, wtype = "stab.ATE") 

## ---- fig.align="center", fig.cap="Visual test of goodness-of-fit of a PS model and Hosmer-Lemeshow test.", fig.height=5, fig.width=5----
gof <- HLfit(obs=rhc$RHC, pred=ps.stab.ATE$ps, n.bins=50, bin.method = "quantiles", plot.bin.size=FALSE)

## ----label=better, fig.align="center", fig.cap="Checking goodness-of-fit of a PS model with selection of significant covariates.", fig.height=5, fig.width=5----
ind <- ps.stab.ATE$summary.glm$coefficients[,4]<0.05
ps.stab.ATE.2 <- get.weights(df=rhc, A="RHC", C=covs.names[ind], wtype = "stab.ATE")
gof <- HLfit(obs=rhc$RHC, pred=ps.stab.ATE.2$ps, n.bins=50, bin.method = "quantiles", plot.bin.size=FALSE)

## ---- fig.align="center", label=fig9, fig.cap="Checking overlap in PS distributions.", fig.height=5, fig.width=5----
rhc.ps <- rhc %>% mutate(ps = ps.stab.ATE$ps, stab.ATE.w = ps.stab.ATE$w)
ggplot(rhc.ps, aes(x = ps, fill = trt, color=trt)) +  geom_density(alpha = 0.5) + scale_fill_npg() + scale_color_npg()+ theme(plot.title = element_text(hjust = 0.5), legend.position = c(0.8, 0.8),
          legend.background=element_rect(fill="transparent"))

## -----------------------------------------------------------------------------
ps.overlap <- get.weights(df=rhc, A="RHC", C=covs.names, wtype = "overlap")
rhc.ps <- rhc.ps %>%   mutate(overlap.w = ps.overlap$w)

## ---- fig.align="center", label=fig10, fig.width=6.5, fig.height=8, message=FALSE, warning=FALSE, fig.cap="Testing covariate balance of the PS model with all covariates."----
#fig.asp=1.1, 
covs <- subset(rhc.ps, select = covs.names)
suppressWarnings(love.plot(RHC ~ covs, data=rhc.ps,
              weights = list(overlap=rhc.ps$overlap.w, stab.ATE=rhc.ps$stab.ATE.w), 
              threshold = .1, abs=TRUE, binary = "std", stars="raw", s.d.denom="pooled",
          var.order = "unadjusted",   
          title = "Comparing balancing weights") + scale_fill_npg() + scale_color_npg())

## ----figcompare, label=figcompare, fig.align="center", fig.width=6.5, fig.height=8,  message=FALSE, warning=FALSE, fig.cap="Testing covariate balance of the PS model with selection of significant covariates."----
ps.overlap.2 <- get.weights(df=rhc, A="RHC", C=covs.names[ind], wtype = "overlap")
rhc.ps <- rhc.ps %>% mutate(ps.2 = ps.stab.ATE.2$ps, stab.ATE.w.2 = ps.stab.ATE.2$w,
                         overlap.w.2 = ps.overlap.2$w)
covs <- subset(rhc.ps, select = covs.names)
suppressWarnings(love.plot(RHC ~ covs, data=rhc.ps,
              weights = list(overlap=rhc.ps$overlap.w.2, stab.ATE=rhc.ps$stab.ATE.w.2), 
              threshold = .1, abs=TRUE, binary = "std", stars="raw", s.d.denom="pooled",
          var.order = "unadjusted",   
          title = "Comparing balancing weights") + scale_fill_npg() + scale_color_npg())

## ---- fig.align="center", out.width=c('33%', '33%', '34%'), fig.show='hold',  message=FALSE, warning=FALSE, fig.cap="Comparing distributions of `aps1` across two treatment arms for: (A) unadjusted, (B) stab.ATE, (C) overlap weights"----
rhc.ps$stab.ATE.w.1 <- rhc.ps$overlap.w.1<- rep(NA, nrow(rhc.ps))
rhc.ps$stab.ATE.w.1[rhc.ps$trt=="RHC"] <-
  rhc.ps$stab.ATE.w[rhc.ps$trt=="RHC"]/sum(rhc.ps$stab.ATE.w[rhc.ps$trt=="RHC"])
rhc.ps$stab.ATE.w.1[rhc.ps$trt=="No RHC"] <-
  rhc.ps$stab.ATE.w[rhc.ps$trt=="No RHC"]/sum(rhc.ps$stab.ATE.w[rhc.ps$trt=="No RHC"])
rhc.ps$overlap.w.1[rhc.ps$trt=="RHC"] <-
  rhc.ps$overlap.w[rhc.ps$trt=="RHC"]/sum(rhc.ps$overlap.w[rhc.ps$trt=="RHC"])
rhc.ps$overlap.w.1[rhc.ps$trt=="No RHC"] <-
  rhc.ps$overlap.w[rhc.ps$trt=="No RHC"]/sum(rhc.ps$overlap.w[rhc.ps$trt=="No RHC"])
ggplot(rhc.ps, aes(x = aps1, fill = trt, color=trt)) + ggtitle("A")+
  geom_density(alpha = 0.5) + scale_fill_npg() + scale_color_npg()+
  theme(plot.title = element_text(hjust = 0.5), legend.position = c(0.7, 0.7),
          legend.background=element_rect(fill="transparent"))
ggplot(rhc.ps, aes(x = aps1, fill = trt, color=trt, weight=stab.ATE.w.1)) + ggtitle("B")+ geom_density(alpha = 0.5) + scale_fill_npg() + scale_color_npg()+
  theme(plot.title = element_text(hjust = 0.5), legend.position = c(0.7, 0.7),
          legend.background=element_rect(fill="transparent"))
ggplot(rhc.ps, aes(x = aps1, fill = trt, color=trt, weight=overlap.w.1)) + ggtitle("C")+  geom_density(alpha = 0.5) + scale_fill_npg() + scale_color_npg()+
  theme(plot.title = element_text(hjust = 0.5), legend.position = c(0.7, 0.7),
          legend.background=element_rect(fill="transparent"))

## ---- fig.align="center", out.width=c('33%', '33%', '34%'), message=FALSE, warning=FALSE, fig.cap="Comparing distributions of `meanbp1` across two treatment arms for: (A) unadjusted, (B) stab.ATE, (C) overlap weights", fig.show='hold'----
ggplot(rhc.ps, aes(x = meanbp1, fill = trt, color=trt)) + ggtitle("A")+
  geom_density(alpha = 0.5) + scale_fill_npg() + scale_color_npg()+
  theme(plot.title = element_text(hjust = 0.5), legend.position = c(0.7, 0.7),
          legend.background=element_rect(fill="transparent"))
ggplot(rhc.ps, aes(x = meanbp1, fill = trt, color=trt, weight=stab.ATE.w.1)) + ggtitle("B")+ geom_density(alpha = 0.5) + scale_fill_npg() + scale_color_npg()+
  theme(plot.title = element_text(hjust = 0.5), legend.position = c(0.7, 0.7),
          legend.background=element_rect(fill="transparent"))
ggplot(rhc.ps, aes(x = meanbp1, fill = trt, color=trt, weight=overlap.w.1)) + ggtitle("C")+  geom_density(alpha = 0.5) + scale_fill_npg() + scale_color_npg()+
  theme(plot.title = element_text(hjust = 0.5), legend.position = c(0.7, 0.7),
          legend.background=element_rect(fill="transparent"))

## ---- label=table1------------------------------------------------------------
d <- rhc.ps[, covs.names]
add.names <- function(x, a) { paste(x,a, sep="")}
# original population:
descr.orig <- round( cbind( 
  apply(d, 2, wtd.mean, weights=rep(1, nrow(d))),
  sqrt(apply(d, 2, wtd.var, weights=rep(1, nrow(d)))), 
  t(apply(d, 2, wtd.quantile, weights=rep(1, nrow(d)), probs=c(.25, .5, .75)))
), dig=2)
descr.ATE <- round( cbind( 
  apply(d, 2, wtd.mean, weights=rhc.ps$stab.ATE.w),
  sqrt(apply(d, 2, wtd.var, weights=rhc.ps$stab.ATE.w)), 
  t(apply(d, 2, wtd.quantile, weights=rhc.ps$stab.ATE.w, probs=c(.25, .5, .75)))
), dig=2)
descr.overlap <- round( cbind( 
  apply(d, 2, wtd.mean, weights=rhc.ps$overlap.w),
  sqrt(apply(d, 2, wtd.var, weights=rhc.ps$overlap.w)), 
  t(apply(d, 2, wtd.quantile, weights=rhc.ps$overlap.w, probs=c(.25, .5, .75)))
), dig=2)
df3 <- data.frame(cbind(descr.orig, descr.ATE, descr.overlap))
row.names(df3) <- covs.names
names(df3) <- rep( c("mean", "sd", "q1", "med", "q3"), times=3)

kable(df3, caption="Comparing marginal distributions of covariates for three populations. q1, med and q3 correspond to 1st, 2nd and 3rd quartiles of a distribution.") %>%
  kable_styling(bootstrap_options = c("striped", "hover", "condensed", "responsive"), full_width = FALSE) %>%
   add_header_above(c("covariate "= 1, "original population" = 5, "stab.ATE population" = 5, "overlap population"=5)) %>%
  column_spec(2, background = "lightgrey") %>%
  column_spec(7, background = "lightgrey") %>%
  column_spec(12, background = "lightgrey") %>%  
     scroll_box(width = "100%", height = "500px")

## ----label=ATEwdistr, fig.height=5, fig.width=5, fig.align="center" ,  message=FALSE, warning=FALSE, fig.cap="Distribution of stab.ATE weights."----
ggplot(rhc.ps, aes(x = stab.ATE.w, fill = trt, color=trt)) +  
  geom_histogram(alpha = 0.4, bins=30) + scale_fill_npg() + scale_color_npg()+ theme(plot.title = element_text(hjust = 0.5), legend.position = c(0.8, 0.8),
          legend.background=element_rect(fill="transparent")) + ggtitle("Distribution of stab.ATE weights")

## ---- label=influential-------------------------------------------------------
kable(rhc.ps[rhc.ps$stab.ATE.w>10, c("trt", "stab.ATE.w", "overlap.w", "ps")],
      caption="The list of most influential observations in the stab.ATE weighting.")

## ---- label=overlapdistr, fig.height=5, fig.width=5, fig.align="center" ,  message=FALSE, warning=FALSE, fig.cap="Distribution of overlap weights."----
ggplot(rhc.ps, aes(x = overlap.w, fill = trt, color=trt)) +  
  geom_histogram(alpha = 0.4, bins=30) + scale_fill_npg() + scale_color_npg()+ theme(plot.title = element_text(hjust = 0.5), legend.position = c(0.8, 0.8),
          legend.background=element_rect(fill="transparent")) + ggtitle("Distribution of overlap weights")

## -----------------------------------------------------------------------------
# overlap weights========================:
# Nonparametric estimation:
res.overlap <- fit.nonpar(df=rhc, X="T", E="E", A="RHC", C=covs.names, 
                          wtype="overlap", cens=0, conf.level=0.95, bs=TRUE, nbs.rep=80, seed=17, parallel = FALSE) # The small number of bootstrap replications (nbs.rep=80) was chosen for illustration purposes.  
# Cox-based estimation:
res.cox.overlap <- fit.cox(df=rhc, X="T", E="E", A="RHC", C=covs.names, 
                          wtype="overlap", cens=0, conf.level=0.95, bs=TRUE, nbs.rep=80, seed=17, parallel = FALSE)
# stab.ATE weights==========================:
# Nonparametric estimation:
res.stab.ATE <- fit.nonpar(df=rhc, X="T", E="E", A="RHC", C=covs.names, 
                          wtype="stab.ATE", cens=0, conf.level=0.95, bs=TRUE, nbs.rep=80,
                          seed=17, parallel = FALSE)
# Cox-based estimation:
res.cox.stab.ATE <- fit.cox(df=rhc, X="T", E="E", A="RHC", C=covs.names, 
                          wtype="stab.ATE", cens=0, conf.level=0.95, bs=TRUE, nbs.rep=120,
                          seed=17, parallel = FALSE)
# unadjusted analysis========================:
# Nonparametric estimation:
res.un <- fit.nonpar(df=rhc, X="T", E="E", A="RHC", C=covs.names, 
                          wtype="unadj", cens=0, conf.level=0.95, bs=TRUE, nbs.rep=80,
                     seed=17, parallel = FALSE)
# Cox-based estimation:
res.cox.un <- fit.cox(df=rhc, X="T", E="E", A="RHC", C=covs.names, 
                          wtype="unadj", cens=0, conf.level=0.95, bs=TRUE, nbs.rep=80,
                      seed=17, parallel = FALSE)

## ----echo=TRUE, message=FALSE, warning=FALSE----------------------------------
get.logCumHR <- function(res, res.cox)
{
  df1 <- rbind( data.frame(time=res$time, Outcome=1,
                      CumHazR=res$trt.eff[[paste("Ev=", 1, sep="")]]$log.CumHazR, 
                      CIL.CumHazR=res$trt.eff[[paste("Ev=", 1, sep="")]]$log.CumHazR.CI.L, 
                      CIU.CumHazR=res$trt.eff[[paste("Ev=", 1, sep="")]]$log.CumHazR.CI.U),
            data.frame(time=res$time, Outcome=2,
                      CumHazR=res$trt.eff[[paste("Ev=", 2, sep="")]]$log.CumHazR, 
                      CIL.CumHazR=res$trt.eff[[paste("Ev=", 2, sep="")]]$log.CumHazR.CI.L, 
                      CIU.CumHazR=res$trt.eff[[paste("Ev=", 2, sep="")]]$log.CumHazR.CI.U))
  df2 <- data.frame(time=c(min(res.cox$time), max(res.cox$time)), Outcome=rep(3, 2),
                         CumHazR=rep(res.cox$trt.eff[[paste("Ev=", 1, sep="")]]$log.CumHazR, 2),    
                         CIL.CumHazR=rep(NA, 2),
                         CIU.CumHazR=rep(NA, 2))
  df3 <- data.frame(time=c(min(res.cox$time), max(res.cox$time)), Outcome=rep(4, 2),
                         CumHazR=rep(res.cox$trt.eff[[paste("Ev=", 2, sep="")]]$log.CumHazR, 2),  
                         CIL.CumHazR=rep(NA, 2),
                         CIU.CumHazR=rep(NA, 2))
  df <- rbind(df1, df2, df3)
  df$Outcome <- factor(df$Outcome)
  levels(df$Outcome) <- c("Discharge-nonpar", "Death-nonpar", 
                          paste("Discharge-Cox:", round(res.cox$trt.eff[[paste("Ev=", 1, sep="")]]$log.CumHazR,dig=2),
                                ", 95% CI=[",
                                round(res.cox$trt.eff[[paste("Ev=", 1, sep="")]]$log.CumHazR.CI.L,dig=2),
                                ",",
                                round(res.cox$trt.eff[[paste("Ev=", 1, sep="")]]$log.CumHazR.CI.U,dig=2),
                                "]", sep=""), 
                          paste("Death-Cox:", round(res.cox$trt.eff[[paste("Ev=", 2, sep="")]]$log.CumHazR,dig=2), ", 95% CI=[",
                                round(res.cox$trt.eff[[paste("Ev=", 2, sep="")]]$log.CumHazR.CI.L,dig=2),
                                ",",
                                round(res.cox$trt.eff[[paste("Ev=", 2, sep="")]]$log.CumHazR.CI.U,dig=2),
                                "]", sep=""))
  return(df)
}
logCumHR.overlap <- get.logCumHR(res.overlap, res.cox.overlap)
logCumHR.stab.ATE <- get.logCumHR(res.stab.ATE, res.cox.stab.ATE)
logCumHR.unadj <- get.logCumHR(res.un, res.cox.un)

## ----echo=TRUE, fig.align="center", fig.show='hold', fig.width=4.5, fig.height=4.5,  message=FALSE, warning=FALSE, fig.cap="log(Cum.HR(t)) of RHC vs No-RHC for discharge and in-hospital death."----

plot.logCumHazR <- function(df, title)
{
  p <- ggplot(df, aes(x=time, y=CumHazR, color=Outcome, fill=Outcome, shape=Outcome)) + ggtitle(title) +
    geom_step(size=1.1) + 
      geom_ribbon(aes(ymin=CIL.CumHazR, ymax=CIU.CumHazR), alpha=0.2, stat="stepribbon") +
    scale_fill_npg() + scale_color_npg()
  p <- p + xlab("time from admission to ICU (days)") + ylab("log.CumHR(t)") + 
    theme(axis.text.x = element_text(face="bold", angle=45),
          axis.text.y = element_text(face="bold"), plot.title = element_text(hjust = 0.5))+
    theme(legend.position = c(0.5, 0.3),
          legend.background=element_rect(fill="transparent"),
          panel.grid.minor = element_line(size = .5,colour = "gray92"), 
          panel.grid.major = element_line(size = .5,colour = "#C0C0C0")) +
    geom_vline(xintercept=30, linetype="dashed")
  p
}

plot.logCumHazR(df=logCumHR.overlap, title="Overlap weights") 
plot.logCumHazR(df=logCumHR.stab.ATE, title="Stabilized ATE weights") 
plot.logCumHazR(df=logCumHR.unadj, title="Unadjusted") 

## ----echo=TRUE, message=FALSE, warning=FALSE----------------------------------
get.CumHaz <- function(res1, res2, res3, Ev)
{
  df <- rbind( data.frame(time=res1$time, population=1, TRT=1,
                          CumHaz=res1$trt.1[[paste("Ev=", Ev, sep="")]]$CumHaz, 
                          CIL.CumHaz=res1$trt.1[[paste("Ev=", Ev, sep="")]]$CumHaz.CI.L, 
                          CIU.CumHaz=res1$trt.1[[paste("Ev=", Ev, sep="")]]$CumHaz.CI.U),
              data.frame(time=res2$time, population=2, TRT=1,
                          CumHaz=res2$trt.1[[paste("Ev=", Ev, sep="")]]$CumHaz, 
                          CIL.CumHaz=res2$trt.1[[paste("Ev=", Ev, sep="")]]$CumHaz.CI.L, 
                          CIU.CumHaz=res2$trt.1[[paste("Ev=", Ev, sep="")]]$CumHaz.CI.U),
             data.frame(time=res3$time, population=3, TRT=1,
                          CumHaz=res3$trt.1[[paste("Ev=", Ev, sep="")]]$CumHaz, 
                          CIL.CumHaz=res3$trt.1[[paste("Ev=", Ev, sep="")]]$CumHaz.CI.L, 
                          CIU.CumHaz=res3$trt.1[[paste("Ev=", Ev, sep="")]]$CumHaz.CI.U),
             data.frame(time=res1$time, population=1, TRT=0,
                          CumHaz=res1$trt.0[[paste("Ev=", Ev, sep="")]]$CumHaz, 
                          CIL.CumHaz=res1$trt.0[[paste("Ev=", Ev, sep="")]]$CumHaz.CI.L, 
                          CIU.CumHaz=res1$trt.0[[paste("Ev=", Ev, sep="")]]$CumHaz.CI.U),
              data.frame(time=res2$time, population=2, TRT=0,
                          CumHaz=res2$trt.0[[paste("Ev=", Ev, sep="")]]$CumHaz, 
                          CIL.CumHaz=res2$trt.0[[paste("Ev=", Ev, sep="")]]$CumHaz.CI.L, 
                          CIU.CumHaz=res2$trt.0[[paste("Ev=", Ev, sep="")]]$CumHaz.CI.U),
             data.frame(time=res3$time, population=3, TRT=0,
                          CumHaz=res3$trt.0[[paste("Ev=", Ev, sep="")]]$CumHaz, 
                          CIL.CumHaz=res3$trt.0[[paste("Ev=", Ev, sep="")]]$CumHaz.CI.L, 
                          CIU.CumHaz=res3$trt.0[[paste("Ev=", Ev, sep="")]]$CumHaz.CI.U)
             )
  df$population <- factor(df$population)
  levels(df$population) <- c("overlap", "stab.ATE", "unadjusted")
  df
}

df.CumHaz <- list()
E.set <- sort(unique(rhc$E))
E.set <- E.set[E.set!=0] # 0 is for censoring status
for (k in E.set)
  df.CumHaz[[k]] <- get.CumHaz(res.overlap, res.stab.ATE, res.un, k)

plot.CumHaz <- function(df, title, ymax)
{
  p <- ggplot(df, aes(x=time, y=CumHaz, color=population, fill=population, shape=population)) + ggtitle(title) +
    geom_step(size=1.1) + 
 #     geom_ribbon(aes(ymin=CIL.CumHaz, ymax=CIU.CumHaz), alpha=0.2, stat="stepribbon") +
    scale_fill_npg() + scale_color_npg()
  p <- p + xlab("time from admission to ICU (days)") + ylab("Cumulative Hazard (t)") + ylim(0, ymax)+
    theme(axis.text.x = element_text(face="bold", angle=45),
          axis.text.y = element_text(face="bold"), plot.title = element_text(hjust = 0.5))+
    theme(legend.position = c(0.7, 0.3),
          legend.background=element_rect(fill="transparent"),
          panel.grid.minor = element_line(size = .5,colour = "gray92"), 
          panel.grid.major = element_line(size = .5,colour = "#C0C0C0")) +
    geom_vline(xintercept=30, linetype="dashed")
  p
}

## ----echo=TRUE, out.width=c('45%', '45%'), fig.width=4.5, fig.height=4.5, fig.align="center" ,  message=FALSE, warning=FALSE, fig.cap="Cumulative hazards for discharge and in-hospital death. Comparison of 3 target populations.", fig.show='hold'----
plot.CumHaz(df=df.CumHaz[[1]][df.CumHaz[[1]]$TRT==1,], title="RHC: CumHaz of Discharge", ymax=6)
plot.CumHaz(df=df.CumHaz[[1]][df.CumHaz[[1]]$TRT==0,], title="No-RHC: CumHaz of Discharge", ymax=6)
plot.CumHaz(df=df.CumHaz[[2]][df.CumHaz[[2]]$TRT==1,], title="RHC: CumHaz of Death", ymax=2.5)
plot.CumHaz(df=df.CumHaz[[2]][df.CumHaz[[2]]$TRT==0,], title="No-RHC: CumHaz of Death", ymax=2.5)

## ----echo=TRUE, fig.align="center",  message=FALSE, warning=FALSE-------------
get.CIF.RMT <- function(res)
{
  df.CIF.RMT <- rbind( data.frame(time=res$time, Event_TRT=1,
                        CIF=res$trt.1[[paste("Ev=", 1, sep="")]]$CIF, 
                        RMT=res$trt.1[[paste("Ev=", 1, sep="")]]$RMT,
                        CumHaz=res$trt.1[[paste("Ev=", 1, sep="")]]$CumHaz,
                        CIL.CIF=res$trt.1[[paste("Ev=", 1, sep="")]]$CIF.CI.L,
                        CIU.CIF=res$trt.1[[paste("Ev=", 1, sep="")]]$CIF.CI.U,
                        CIL.RMT=res$trt.1[[paste("Ev=", 1, sep="")]]$RMT.CI.L,
                        CIU.RMT=res$trt.1[[paste("Ev=", 1, sep="")]]$RMT.CI.U,
                        CIL.CumHaz=res$trt.1[[paste("Ev=", 1, sep="")]]$CumHaz.CI.L,
                        CIU.CumHaz=res$trt.1[[paste("Ev=", 1, sep="")]]$CumHaz.CI.U),
            data.frame(time=res$time, Event_TRT=2,
                       CIF=res$trt.0[[paste("Ev=", 1, sep="")]]$CIF, 
                       RMT=res$trt.0[[paste("Ev=", 1, sep="")]]$RMT,
                       CumHaz=res$trt.0[[paste("Ev=", 1, sep="")]]$CumHaz,
                      CIL.CIF=res$trt.0[[paste("Ev=", 1, sep="")]]$CIF.CI.L,
                      CIU.CIF=res$trt.0[[paste("Ev=", 1, sep="")]]$CIF.CI.U,
                      CIL.RMT=res$trt.0[[paste("Ev=", 1, sep="")]]$RMT.CI.L,
                      CIU.RMT=res$trt.0[[paste("Ev=", 1, sep="")]]$RMT.CI.U,
                      CIL.CumHaz=res$trt.0[[paste("Ev=", 1, sep="")]]$CumHaz.CI.L,
                      CIU.CumHaz=res$trt.0[[paste("Ev=", 1, sep="")]]$CumHaz.CI.U),
             data.frame(time=res$time, Event_TRT=3,
                        CIF=res$trt.1[[paste("Ev=", 2, sep="")]]$CIF, 
                        RMT=res$trt.1[[paste("Ev=", 2, sep="")]]$RMT,
                        CumHaz=res$trt.1[[paste("Ev=", 2, sep="")]]$CumHaz,
                        CIL.CIF=res$trt.1[[paste("Ev=", 2, sep="")]]$CIF.CI.L,
                        CIU.CIF=res$trt.1[[paste("Ev=", 2, sep="")]]$CIF.CI.U,
                        CIL.RMT=res$trt.1[[paste("Ev=", 2, sep="")]]$RMT.CI.L,
                        CIU.RMT=res$trt.1[[paste("Ev=", 2, sep="")]]$RMT.CI.U,
                        CIL.CumHaz=res$trt.1[[paste("Ev=", 2, sep="")]]$CumHaz.CI.L,
                        CIU.CumHaz=res$trt.1[[paste("Ev=", 2, sep="")]]$CumHaz.CI.U),
            data.frame(time=res$time,  Event_TRT=4,
                       CIF=res$trt.0[[paste("Ev=", 2, sep="")]]$CIF,
                       RMT=res$trt.0[[paste("Ev=", 2, sep="")]]$RMT,
                       CumHaz=res$trt.0[[paste("Ev=", 2, sep="")]]$CumHaz,
                      CIL.CIF=res$trt.0[[paste("Ev=", 2, sep="")]]$CIF.CI.L,
                      CIU.CIF=res$trt.0[[paste("Ev=", 2, sep="")]]$CIF.CI.U,
                      CIL.RMT=res$trt.0[[paste("Ev=", 2, sep="")]]$RMT.CI.L,
                      CIU.RMT=res$trt.0[[paste("Ev=", 2, sep="")]]$RMT.CI.U,
                      CIL.CumHaz=res$trt.0[[paste("Ev=", 2, sep="")]]$CumHaz.CI.L,
                      CIU.CumHaz=res$trt.0[[paste("Ev=", 2, sep="")]]$CumHaz.CI.U)
            )
  df.CIF.RMT$Event_TRT <- factor(df.CIF.RMT$Event_TRT)
  levels(df.CIF.RMT$Event_TRT) <- c("Discharge-RHC", "Discharge-No RHC", "In-hospital death-RHC", "In-hospital death-No RHC")
  return(df.CIF.RMT)
}

df.CIF.RMT <- get.CIF.RMT(res.overlap)
pest.overlap.30 <- get.pointEst(res.overlap, 30) 
pest.overlap.100 <- get.pointEst(res.overlap, 100) 

## ----label=CumHaz4, echo=TRUE, fig.align="center", fig.width=4.5, fig.height=4.5,  message=FALSE, warning=FALSE, fig.cap="Cause-specific cumulative hazards.", fig.show='hold'----
p <- ggplot(df.CIF.RMT, aes(x=time, y=CumHaz, color=Event_TRT, fill=Event_TRT, shape=Event_TRT)) +
  geom_step(size=1.1) + #ggtitle("Cumulative hazards") + 
    geom_ribbon(aes(ymin=CIL.CumHaz, ymax=CIU.CumHaz), alpha=0.2, stat="stepribbon") +
   scale_fill_npg() + scale_color_npg() #+xlim(0, 150)
p <- p +  geom_vline(xintercept=30, linetype="dashed")+
  xlab("time from admission to ICU (days)") + ylab("Cumulative hazard of event by time t") +
  theme(axis.text.x = element_text(face="bold", angle=45),
        axis.text.y = element_text(face="bold"), plot.title = element_text(hjust = 0.5))+
  theme(legend.position = c(0.3, 0.8),
        legend.background=element_rect(fill="transparent"),
        panel.grid.minor = element_line(size = .5,colour = "gray92"),
        panel.grid.major = element_line(size = .5,colour = "#C0C0C0"))
p

## ----label=Risk4, echo=TRUE, fig.align="center", fig.width=4.5, fig.height=4.5,  message=FALSE, warning=FALSE, fig.cap="Cause-specific risk curves", fig.show='hold'----
p <- ggplot(df.CIF.RMT, aes(x=time, y=CIF, color=Event_TRT, fill=Event_TRT, shape=Event_TRT)) +
  geom_step(size=1.1) + #ggtitle("Risk curves") + 
      geom_ribbon(aes(ymin=CIL.CIF, ymax=CIU.CIF), alpha=0.2, stat="stepribbon") +
   scale_fill_npg() + scale_color_npg() +xlim(0, 150)
p <- p +  geom_vline(xintercept=30, linetype="dashed")+
  xlab("time from admission to ICU (days)") + ylab("Probability of event by time t") +
  theme(axis.text.x = element_text(face="bold", angle=45),
        axis.text.y = element_text(face="bold"), plot.title = element_text(hjust = 0.5))+
  theme(legend.position = c(0.7, 0.25),
        legend.background=element_rect(fill="transparent"),
        panel.grid.minor = element_line(size = .5,colour = "gray92"),
        panel.grid.major = element_line(size = .5,colour = "#C0C0C0"))
p

## ----label=RMT4, echo=TRUE, fig.align="center", fig.width=4.5, fig.height=4.5, message=FALSE, warning=FALSE, fig.cap="Restricted-mean-time-lost/gained due to specific event.", fig.show='hold'----
p <- ggplot(df.CIF.RMT, aes(x=time, y=RMT, color=Event_TRT, 
                            fill=Event_TRT, shape=Event_TRT)) +
  geom_line(size=1.1) + #ggtitle("Restricted-mean-time-lost/gained due to specific event") +  
  geom_ribbon(aes(ymin=CIL.RMT, ymax=CIU.RMT), alpha=0.2) + #, stat="stepribbon") +
  scale_fill_npg() + scale_color_npg()
 p <- p +  geom_vline(xintercept=30, linetype="dashed")+
  xlab("time from admission to ICU (days)") + ylab("average time lost due to death (gained by recovery)") + 
  theme(axis.text.x = element_text(face="bold", angle=45),
        axis.text.y = element_text(face="bold"), plot.title = element_text(hjust = 0.5))+
  theme(legend.position = c(0.3, 0.8),
        legend.background=element_rect(fill="transparent"),
        panel.grid.minor = element_line(size = .5,colour = "gray92"),
        panel.grid.major = element_line(size = .5,colour = "#C0C0C0"))
p

## ----echo=TRUE, fig.align="center",  message=FALSE, warning=FALSE-------------
get.RD.RMT <- function(res)
{
  df <- rbind( data.frame(time=res$time, Outcome=1,
                      RD=res$trt.eff[[paste("Ev=", 1, sep="")]]$RD, 
                      CIL.RD=res$trt.eff[[paste("Ev=", 1, sep="")]]$RD.CI.L, 
                      CIU.RD=res$trt.eff[[paste("Ev=", 1, sep="")]]$RD.CI.U,
                      ATE.RMT=res$trt.eff[[paste("Ev=", 1, sep="")]]$ATE.RMT, 
                      CIL.ATE.RMT=res$trt.eff[[paste("Ev=", 1, sep="")]]$ATE.RMT.CI.L, 
                      CIU.ATE.RMT=res$trt.eff[[paste("Ev=", 1, sep="")]]$ATE.RMT.CI.U),
            data.frame(time=res$time, Outcome=2,
                      RD=res$trt.eff[[paste("Ev=", 2, sep="")]]$RD, 
                      CIL.RD=res$trt.eff[[paste("Ev=", 2, sep="")]]$RD.CI.L, 
                      CIU.RD=res$trt.eff[[paste("Ev=", 2, sep="")]]$RD.CI.U,
                      ATE.RMT=res$trt.eff[[paste("Ev=", 2, sep="")]]$ATE.RMT, 
                      CIL.ATE.RMT=res$trt.eff[[paste("Ev=", 2, sep="")]]$ATE.RMT.CI.L, 
                      CIU.ATE.RMT=res$trt.eff[[paste("Ev=", 2, sep="")]]$ATE.RMT.CI.U))
  
  df$Outcome <- factor(df$Outcome)
  levels(df$Outcome) <- c("Discharge", "In-hospital death")
  return(df)
}
df.RD.RMT <- get.RD.RMT(res.overlap)

## ----label=fig16, echo=TRUE, fig.align="center", fig.width=4.5, fig.height=4.5,   message=FALSE, warning=FALSE, fig.cap="Risk differences for discharge and in-hospital death.", fig.show='hold'----
p <- ggplot(df.RD.RMT, aes(x=time, y=RD, color=Outcome, 
                            fill=Outcome, shape=Outcome)) +
  geom_step(size=1.1) + #ggtitle("Risk differences for discharge and in-hospital death") + 
    geom_ribbon(aes(ymin=CIL.RD, ymax=CIU.RD), alpha=0.2, stat="stepribbon") +
  scale_fill_npg() + scale_color_npg()
p <- p +  geom_vline(xintercept=30, linetype="dashed")+
  xlab("time from admission to ICU (days)") + 
  ylab("Risk difference") + xlim(0, 150) +
  theme(axis.text.x = element_text(face="bold", angle=45),
        axis.text.y = element_text(face="bold"), plot.title = element_text(hjust = 0.5))+
  theme(legend.position = c(0.8, 0.5),
        legend.background=element_rect(fill="transparent"),
        panel.grid.minor = element_line(size = .5,colour = "gray92"),
        panel.grid.major = element_line(size = .5,colour = "#C0C0C0"))
p

## ----echo=TRUE, fig.align="center", fig.width=4.5, fig.height=4.5,  message=FALSE, warning=FALSE, fig.cap="Difference in average time lost/gained due to treatment.", fig.show='hold'----
p <- ggplot(df.RD.RMT, aes(x=time, y=ATE.RMT, color=Outcome, 
                            fill=Outcome, shape=Outcome)) +
  geom_line(size=1.1) + #ggtitle("Difference in average time lost/gained due to treatment") + 
    geom_ribbon(aes(ymin=CIL.ATE.RMT, ymax=CIU.ATE.RMT), alpha=0.2) + #, stat="stepribbon") +
  scale_fill_npg() + scale_color_npg()
p <- p +  geom_vline(xintercept=30, linetype="dashed")+
  xlab("time from admission to ICU (days)") + 
  ylab("time (days)") + 
  theme(axis.text.x = element_text(face="bold", angle=45),
        axis.text.y = element_text(face="bold"), plot.title = element_text(hjust = 0.5))+
  theme(legend.position = c(0.8, 0.5),
        legend.background=element_rect(fill="transparent"),
        panel.grid.minor = element_line(size = .5,colour = "gray92"),
        panel.grid.major = element_line(size = .5,colour = "#C0C0C0"))
p

## -----------------------------------------------------------------------------
get.cumPr <- function(res) # 1-KM(t)
{
  df <- rbind( data.frame(time=res$time, TRT=1,
                          CumPr=1-exp(-res$trt.1[[paste("Ev=", 1, sep="")]]$CumHaz - 
                                   res$trt.1[[paste("Ev=", 2, sep="")]]$CumHaz) ),
               data.frame(time=res$time, TRT=2, 
                          CumPr=1-exp(-res$trt.0[[paste("Ev=", 1, sep="")]]$CumHaz - 
                                   res$trt.0[[paste("Ev=", 2, sep="")]]$CumHaz))
  )
  df$TRT <- factor(df$TRT)
  levels(df$TRT) <- c("RHC", "No-RHC")
  df
}
df.cumPr <- get.cumPr(res.overlap)

## ----echo=TRUE, fig.align="center", fig.width=4.5, fig.height=4.5, message=FALSE, warning=FALSE, fig.cap="Cumulative distribution of length-of-stay in RHC vs No-RHC groups.", fig.show='hold'----
p <- ggplot(df.cumPr, aes(x=time, y=CumPr, color=TRT, fill=TRT, shape=TRT)) + 
     geom_step(size=1.1) + #ggtitle("Cumulative distribution of length-of-stay")+
    scale_fill_npg() + scale_color_npg()
  p <- p + xlab("time (years)") + ylab("Probability of discharge") + ylim(0, 1)+xlim(0,100)+
    theme(axis.text.x = element_text(face="bold", angle=45),
          axis.text.y = element_text(face="bold"), plot.title = element_text(hjust = 0.5))+
    theme(legend.position = c(0.6, 0.3),
          legend.background=element_rect(fill="transparent"),
          panel.grid.minor = element_line(size = .5,colour = "gray92"), 
          panel.grid.major = element_line(size = .5,colour = "#C0C0C0")) +
    geom_vline(xintercept=30, linetype="dashed")
  p

## -----------------------------------------------------------------------------
# nonparametric:
# acm=all-cause mortality
res.acm.30 <- fit.nonpar(df=rhc_full, X="T.death.30", E="D.30", A="RHC", C=covs.names, wtype="overlap", cens=0, conf.level=0.95, bs=TRUE, nbs.rep=80, seed=17, parallel = FALSE)
# Cox-based estimation:
res.cox.acm.30 <- fit.cox(df=rhc_full, X="T.death.30", E="D.30", A="RHC", C=covs.names, wtype="overlap", cens=0, conf.level=0.95, bs=TRUE, nbs.rep=80, seed=17, parallel = FALSE)

## ----echo=TRUE, fig.align="center",  fig.width=5, fig.height=5, message=FALSE, warning=FALSE, fig.cap="log-ratio of cumulative hazards for 30-day mortality, using overlap weighting."----
df1.acm.30 <- rbind( data.frame(time=res.acm.30$time, 
                                CumHazR=res.acm.30$trt.eff[[paste("Ev=", 1, sep="")]]$log.CumHazR, Outcome=1,
                        CIL.CumHazR=res.acm.30$trt.eff[[paste("Ev=", 1, sep="")]]$log.CumHazR.CI.L, 
                        CIU.CumHazR=res.acm.30$trt.eff[[paste("Ev=", 1, sep="")]]$log.CumHazR.CI.U))
df2.acm.30 <- data.frame(time=c(min(res.cox.acm.30$time), max(res.cox.acm.30$time)),
                       CumHazR=rep(res.cox.acm.30$trt.eff[[paste("Ev=", 1, sep="")]]$log.CumHazR, 2), 
                       Outcome=rep(3, 2),
                       CIL.CumHazR=rep(NA, 2),
                       CIU.CumHazR=rep(NA, 2))

df.acm.30 <- rbind(df1.acm.30, df2.acm.30)
df.acm.30$Outcome <- factor(df.acm.30$Outcome)
levels(df.acm.30$Outcome) <- c("Death-nonpar",  
                        paste("Death-Cox:", round(res.cox.acm.30$trt.eff[[paste("Ev=", 1, sep="")]]$log.CumHazR,dig=2),
                                ", 95% CI=[",
                                round(res.cox.acm.30$trt.eff[[paste("Ev=", 1, sep="")]]$log.CumHazR.CI.L,dig=2),
                                ",",
                                round(res.cox.acm.30$trt.eff[[paste("Ev=", 1, sep="")]]$log.CumHazR.CI.U,dig=2),
                                "]", sep=""))

p <- ggplot(df.acm.30, aes(x=time, y=CumHazR, color=Outcome, fill=Outcome, shape=Outcome)) +
  #geom_line(size=1) +
 # ggtitle("log.CumHR(t) of 30-day mortality") +
  geom_step(size=1.1) + 
    geom_ribbon(aes(ymin=CIL.CumHazR, ymax=CIU.CumHazR), alpha=0.2, stat="stepribbon") +
  scale_fill_npg() + scale_color_npg()
p <- p +
  xlab("time from admission to ICU (days)") + ylim(0, 0.75)+
  ylab("log-ratio of CumHaz(t)") + #geom_hline(yintercept=0)+ #ylim(-0.35, 0.05)+
  theme(axis.text.x = element_text(face="bold", angle=45),
        axis.text.y = element_text(face="bold"), plot.title = element_text(hjust = 0.5))+
  theme(legend.position = c(0.6, 0.8),
        legend.background=element_rect(fill="transparent"),
        panel.grid.minor = element_line(size = .5,colour = "gray92"), 
        panel.grid.major = element_line(size = .5,colour = "#C0C0C0")) 
p

# getting an estimate at a time point = 30 days: 
est30 <- get.pointEst(res.acm.30, 30) 
#cat("HR=", round(exp(est30$trt.eff[[1]]$log.CumHazR), dig=2),
#    ", 95% CI=[", round(exp(est30$trt.eff[[1]]$log.CumHazR.CI.L), dig=2), ",", 
#                       round(exp(est30$trt.eff[[1]]$log.CumHazR.CI.U), dig=2), "]\n")

## ----echo=TRUE, fig.align="center",  message=FALSE, warning=FALSE-------------
get.RD.RMT <- function(res)
{
  df <- rbind( data.frame(time=res$time, 
                      RD=res$trt.eff[[paste("Ev=", 1, sep="")]]$RD, 
                      CIL.RD=res$trt.eff[[paste("Ev=", 1, sep="")]]$RD.CI.L, 
                      CIU.RD=res$trt.eff[[paste("Ev=", 1, sep="")]]$RD.CI.U,
                      ATE.RMT=res$trt.eff[[paste("Ev=", 1, sep="")]]$ATE.RMT, 
                      CIL.ATE.RMT=res$trt.eff[[paste("Ev=", 1, sep="")]]$ATE.RMT.CI.L, 
                      CIU.ATE.RMT=res$trt.eff[[paste("Ev=", 1, sep="")]]$ATE.RMT.CI.U))
  
  return(df)
}
df.acm.RD.RMT <- get.RD.RMT(res.acm.30)

## ----label=fig162, echo=TRUE, fig.align="center", fig.width=4.5, fig.height=4.5,   message=FALSE, warning=FALSE, fig.cap="Risk differences for 30-day mortality.", fig.show='hold'----
p <- ggplot(df.acm.RD.RMT, aes(x=time, y=RD)) +
  geom_step(size=1.1) +  
    geom_ribbon(aes(ymin=CIL.RD, ymax=CIU.RD), alpha=0.2, stat="stepribbon") +
  scale_fill_npg() + scale_color_npg()
p <- p +  #geom_vline(xintercept=30, linetype="dashed")+
  xlab("time from admission to ICU (days)") + 
  ylab("Risk difference") + #xlim(0, 150) +
  theme(axis.text.x = element_text(face="bold", angle=45),
        axis.text.y = element_text(face="bold"), plot.title = element_text(hjust = 0.5))+
  theme(legend.position = c(0.8, 0.5),
        legend.background=element_rect(fill="transparent"),
        panel.grid.minor = element_line(size = .5,colour = "gray92"),
        panel.grid.major = element_line(size = .5,colour = "#C0C0C0"))
p

#cat("Risk difference=", round((est30$trt.eff[[1]]$RD), dig=4),
#    ", 95% CI=[", round((est30$trt.eff[[1]]$RD.CI.L), dig=4), ",", 
#                       round((est30$trt.eff[[1]]$RD.CI.U), dig=4), "]\n")

## ---- label=figRMT30, echo=TRUE, fig.align="center", fig.width=4.5, fig.height=4.5,  message=FALSE, warning=FALSE, fig.cap="Difference in average time lost due to RHC in the 30-day mortality.", fig.show='hold'----
p <- ggplot(df.acm.RD.RMT, aes(x=time, y=ATE.RMT)) +
  geom_line(size=1.1) + #ggtitle("Difference in average time lost/gained due to treatment") + 
    geom_ribbon(aes(ymin=CIL.ATE.RMT, ymax=CIU.ATE.RMT), alpha=0.2) + #, stat="stepribbon") +
  scale_fill_npg() + scale_color_npg()
p <- p +  #geom_vline(xintercept=30, linetype="dashed")+
  xlab("time from admission to ICU (days)") + 
  ylab("RMT differences (days)") + 
  theme(axis.text.x = element_text(face="bold", angle=45),
        axis.text.y = element_text(face="bold"), plot.title = element_text(hjust = 0.5))+
  theme(legend.position = c(0.8, 0.5),
        legend.background=element_rect(fill="transparent"),
        panel.grid.minor = element_line(size = .5,colour = "gray92"),
        panel.grid.major = element_line(size = .5,colour = "#C0C0C0"))
p
#cat("RMTL difference=", round((est30$trt.eff[[1]]$ATE.RMT), dig=2),
#    ", 95% CI=[", round((est30$trt.eff[[1]]$ATE.RMT.CI.L), dig=2), ",", 
#                       round((est30$trt.eff[[1]]$ATE.RMT.CI.U), dig=2), "]\n")

