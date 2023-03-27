# causalCmprsk - Nonparametric and Cox-based Estimation of Average Treatment Effects in Competing Risks

The `causalCmprsk` package is designed for estimation of average treatment effects (ATE) of point interventions/treatments on time-to-event outcomes with K competing events (K can be 1). The method assumes that there is no unmeasured confounding and uses propensity scores weighting for emulation of baseline randomization. 

The `causalCmprsk` package provides two main functions: `fit.cox` which assumes the Cox proportional hazards regression for potential outcomes, and `fit.nonpar` that does not make any modeling assumptions for potential outcomes. 

# Installation 
The  `causalCmprsk` package can be installed by
```{r}
devtools::install_github("Bella2001/causalCmprsk")
```
# Examples
The examples of how to use `causalCmprsk` package on real data can be found [here](https://htmlpreview.github.io/?https://github.com/Bella2001/causalCmprsk/blob/dev/index.html).

# References
 
- M.-L. Charpignon, B. Vakulenko-Lagun, B. Zheng, C. Magdamo  et al., *Causal inference in medical records and complementary systems pharmacology for metformin drug repurposing towards dementia*, 2022, Nature Communications.



