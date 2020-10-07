# October 07, 2020
# compare sequential and parallel vesions in the analysis of the RHC data set

library(survival)
library(dplyr)
library(data.table) # check if I use it!!!
library(ggplot2)
library(ggalt) # for stepribbon
library(tidyverse)
library(cobalt)
library(ggsci) # Nature Color palette
library(modEvA) # for H-L GOF test
library(naniar) # for exploring missingness
library(readxl)
library(knitr)
library(kableExtra)
library(magrittr)
library(summarytools) # for summary distributions of covariates - check if I use it!!!
library(DT) # check if I use it!!!!
library(Hmisc)

library(causalCmprsk)

column_types_rhc <-
  cols(urin1 = "d", meanbp1 = "d", resp1 = "d",
       swang1 = col_factor(c("No RHC", "RHC")),
       death = col_factor(c("No", "Yes")),
       sex = col_factor(c("Male", "Female")),
       cat1 = col_factor(c("ARF", "CHF", "Cirrhosis", "Colon Cancer", "Coma", "COPD",
                           "Lung Cancer", "MOSF w/Malignancy", "MOSF w/Sepsis")),
       # Cat1 - Primary disease category
       cat2 = col_factor(c("ARF", "CHF", "Cirrhosis", "Colon Cancer", "Coma", "COPD",
                           "Lung Cancer", "MOSF w/Malignancy", "MOSF w/Sepsis")),
       # Cat2 - Secondary disease category
       dnr1 = col_factor(c("No", "Yes")),
       card = col_factor(c("No", "Yes")),
       gastr = col_factor(c("No", "Yes")),
       hema = col_factor(c("No", "Yes")),
       meta = col_factor(c("No", "Yes")),
       neuro = col_factor(c("No", "Yes")),
       ortho = col_factor(c("No", "Yes")),
       renal = col_factor(c("No", "Yes")),
       resp = col_factor(c("No", "Yes")),
       seps = col_factor(c("No", "Yes")),
       trauma = col_factor(c("No", "Yes")),
       income = col_factor(c("Under $11k", "$11-$25k", "$25-$50k", "> $50k")),
       ninsclas = col_factor(c("Private", "Private & Medicare", "Medicare",
                               "Medicare & Medicaid", "Medicaid", "No insurance")),
       race = col_factor(c("white", "black", "other")),
       ca = col_factor(c("No", "Yes", "Metastatic"))
  )
suppressWarnings(rhc_raw <- read_csv("http://biostat.mc.vanderbilt.edu/wiki/pub/Main/DataSets/rhc.csv",
                                     col_types = column_types_rhc))

miss.report <- rhc_raw %>% miss_var_summary()
kable(miss.report[1:10,]) %>%
  kable_styling(bootstrap_options = "striped", full_width = F)

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


# omit 1 obs with missing discharge date for the length-of-stay analysis:
rhc <- rhc_full[!is.na(rhc_raw$dschdte),]
rhc_full$T.death.30 <- pmin(30, rhc_full$T.death)
rhc_full$D.30 <- ifelse(rhc_full$T.death <=30, 1, 0)


E <- as.factor(rhc$E)
levels(E) <- c("discharge", "in-hospital death")
t <- addmargins(table(E, useNA="no"))
kable(t) %>% kable_styling(bootstrap_options = "striped", full_width = F)

D <- as.factor(rhc_full$D)
levels(D) <- c("censored", "died")
t <- addmargins(table(D, useNA="no"))
kable(t) %>% kable_styling(bootstrap_options = "striped", full_width = F)

D.30 <- as.factor(rhc_full$D.30)
levels(D.30) <- c("censored", "died")
t <- addmargins(table(D.30, useNA="no"))
kable(t) %>% kable_styling(bootstrap_options = "striped", full_width = F)

t <- addmargins(table(rhc_full$trt, useNA="no"))
kable(t) %>% kable_styling(bootstrap_options = "striped", full_width = F)


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

# Nonparametric estimation:
res.overlap <- fit.nonpar(df=rhc, T="T", E="E", A="RHC", C=covs.names, wtype="overlap", bs=TRUE, nbs.rep=10, seed=17, cens=0, conf.level=0.95)
par.res.overlap <- parallel.fit.nonpar(df=rhc, T="T", E="E", A="RHC", C=covs.names, wtype="overlap", bs=TRUE, nbs.rep=10, seed=17, cens=0, conf.level=0.95)

# Cox-based estimation:
res.cox.overlap <- fit.cox(df=rhc, T="T", E="E", A="RHC", C=covs.names, wtype="overlap", bs=TRUE, nbs.rep=10, seed=17, cens=0, conf.level=0.95)
par.res.cox.overlap <- parallel.fit.cox(df=rhc, T="T", E="E", A="RHC", C=covs.names, wtype="overlap", bs=TRUE, nbs.rep=10, seed=17, cens=0, conf.level=0.95)


# # ------------- check point estimates in sequential and parallel versions of code -------------
# get.CumHaz <- function(res1, res2, Ev)
# {
#   df <- rbind( data.frame(time=res1$time, population=1, TRT=1,
#                           CumHaz=res1$trt.1[[paste("Ev=", Ev, sep="")]]$CumHaz,
#                           CIL.CumHaz=res1$trt.1[[paste("Ev=", Ev, sep="")]]$CumHaz.CI.L,
#                           CIU.CumHaz=res1$trt.1[[paste("Ev=", Ev, sep="")]]$CumHaz.CI.U),
#                data.frame(time=res2$time, population=2, TRT=1,
#                           CumHaz=res2$trt.1[[paste("Ev=", Ev, sep="")]]$CumHaz,
#                           CIL.CumHaz=res2$trt.1[[paste("Ev=", Ev, sep="")]]$CumHaz.CI.L,
#                           CIU.CumHaz=res2$trt.1[[paste("Ev=", Ev, sep="")]]$CumHaz.CI.U),
#                # data.frame(time=res3$time, population=3, TRT=1,
#                #            CumHaz=res3$trt.1[[paste("Ev=", Ev, sep="")]]$CumHaz,
#                #            CIL.CumHaz=res3$trt.1[[paste("Ev=", Ev, sep="")]]$CumHaz.CI.L,
#                #            CIU.CumHaz=res3$trt.1[[paste("Ev=", Ev, sep="")]]$CumHaz.CI.U),
#                data.frame(time=res1$time, population=1, TRT=0,
#                           CumHaz=res1$trt.0[[paste("Ev=", Ev, sep="")]]$CumHaz,
#                           CIL.CumHaz=res1$trt.0[[paste("Ev=", Ev, sep="")]]$CumHaz.CI.L,
#                           CIU.CumHaz=res1$trt.0[[paste("Ev=", Ev, sep="")]]$CumHaz.CI.U),
#                data.frame(time=res2$time, population=2, TRT=0,
#                           CumHaz=res2$trt.0[[paste("Ev=", Ev, sep="")]]$CumHaz,
#                           CIL.CumHaz=res2$trt.0[[paste("Ev=", Ev, sep="")]]$CumHaz.CI.L,
#                           CIU.CumHaz=res2$trt.0[[paste("Ev=", Ev, sep="")]]$CumHaz.CI.U)
#                # data.frame(time=res3$time, population=3, TRT=0,
#                #            CumHaz=res3$trt.0[[paste("Ev=", Ev, sep="")]]$CumHaz,
#                #            CIL.CumHaz=res3$trt.0[[paste("Ev=", Ev, sep="")]]$CumHaz.CI.L,
#                #            CIU.CumHaz=res3$trt.0[[paste("Ev=", Ev, sep="")]]$CumHaz.CI.U)
#   )
#   df$population <- factor(df$population)
#   levels(df$population) <- c("sequential", "parallel")
#   df
# }
#
# df.CumHaz <- list()
# E.set <- sort(unique(rhc$E))
# E.set <- E.set[E.set!=0] # 0 is for censoring status
# for (k in E.set)
#   df.CumHaz[[k]] <- get.CumHaz(res.overlap, par.res.overlap, k)
#
#
# plot.CumHaz <- function(df, title, ymax)
# {
#   p <- ggplot(df, aes(x=time, y=CumHaz, color=population, fill=population, shape=population)) + ggtitle(title) +
#     geom_step(size=1.1) +
#     #     geom_ribbon(aes(ymin=CIL.CumHaz, ymax=CIU.CumHaz), alpha=0.2, stat="stepribbon") +
#     scale_fill_npg() + scale_color_npg()
#   p <- p + xlab("time from admission to ICU (days)") + ylab("Cumulative Hazard (t)") + ylim(0, ymax)+
#     theme(axis.text.x = element_text(face="bold", angle=45),
#           axis.text.y = element_text(face="bold"), plot.title = element_text(hjust = 0.5))+
#     theme(legend.position = c(0.7, 0.3),
#           legend.background=element_rect(fill="transparent"),
#           panel.grid.minor = element_line(size = .5,colour = "gray92"),
#           panel.grid.major = element_line(size = .5,colour = "#C0C0C0")) +
#     geom_vline(xintercept=30, linetype="dashed")
#   p
# }
#
# plot.CumHaz(df=df.CumHaz[[1]][df.CumHaz[[1]]$TRT==1,], title="RHC: CumHaz of Discharge", ymax=6)
# plot.CumHaz(df=df.CumHaz[[1]][df.CumHaz[[1]]$TRT==0,], title="No-RHC: CumHaz of Discharge", ymax=6)
# plot.CumHaz(df=df.CumHaz[[2]][df.CumHaz[[2]]$TRT==1,], title="RHC: CumHaz of Death", ymax=2.5)
# plot.CumHaz(df=df.CumHaz[[2]][df.CumHaz[[2]]$TRT==0,], title="No-RHC: CumHaz of Death", ymax=2.5)



