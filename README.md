# tclboot an R package for
## Conditional inference in small sample scenarios using a resampling approach


This R package provides an implementation of a non-parametric resampling technique in the context of multidimensional hypothesis testing of assumptions of the Rasch model based on conditional distributions. It is suggested in small sample size scenarios, in which asymptotic theory is not applicable, to approximate the exact sampling distribution of various well-known Chi^2 test statistics like Wald, likelihood ratio, score and gradient tests as well as others. A procedure to compute the power function of the tests is also presented. A number of examples of scenarios are discussed in which the power function does not converge to one with an increasing deviation from the hypothesis to be tested, i.e. the respective assumption of the model. Finally, an attempt to modify the critical region of the test is made aiming at improving the power.


## Installation


The package can be installed using the `devtools`-package:

```
install.packages("devtools")
devtools::install_github("akurz1/tclboot")
```

## Figures 3 and 4 

Draxler, C., & Kurz, A. (submitted). Conditional inference in small sample scenarios. Stats(Special Issue "Re-sampling Methods for Statistical Inference of the 2020s"). Abstract Retrieved from https://www.mdpi.com/journal/stats/special_issues/resampling_statistics

#### Figure 3: 

Local power curve for the test U_1 and its modiﬁed version referring to an extreme item in a scenario with sample size n = 30 and number of items k = 10. The nominal level of α is set to .05 (dotted horizontal line).
```
res3 <- MCplot_simdata_crit(N=30, k=10, nstats = "sqs", alpha = 0.05,
+                                xlim = 5, legend.position = "top", legend.direction = "horizontal",
+                                ctr_alpha = 1,
+                                ctr_mod = "high")

Fig3 <- res3$power_item$I2  # ggplot object
```

#### Figure 4: 

Local power curve for the test U_1 and its modiﬁed version referring to all (free) items in a scenario with sample size n = 30 and number of items k = 10. The nominal level of α is set to .05 (dotted horizontal line).
```
res4 <- MCplot_simdata_crit(N=30, k=10, nstats = "sqs", alpha = 0.05,
                               xlim = 5, legend.position = c(.5, .85), legend.direction = "horizontal",
                               ctr_alpha = .2,
                               ctr_mod = "both")

ggpubr::ggarrange(plotlist = res3$power_item[-10], ncol=3, nrow = 3, common.legend = TRUE, legend = "top")
```

