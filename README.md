# causalCmprsk - Nonparametric and Cox-based Estimation of ATE in Competing Risks

The `causalCmprsk` package is designed for estimation of average treatment effects (ATE) of two static treatment regimes on time-to-event outcomes with K competing events (K can be 1). The method uses propensity scores weighting for emulation of baseline randomization. The package accompanies the paper of Charpignon, Vakulenko-Lagun, Zheng, Magdamo et al., *Uncovering the links between metformin, dementia and aging using emulated trials in EHR and systems pharmacology* (submitted to Nature Medicine, 2020).

The `causalCmprsk` package provides two main functions: `fit.cox` which assumes the Cox proportional hazards regression for potential outcomes, and `fit.nonpar` that does not make any modeling assumptions for potential outcomes. 

The  `causalCmprsk` package can be installed by
```{r}
devtools::install_github("Bella2001/causalCmprsk")
```
The examples of how to use `causalCmprsk` package are [here]( https://htmlpreview.github.io/?https://github.com/Bella2001/causalCmprsk/blob/causalCmprsk/index.html).
 
 
 
