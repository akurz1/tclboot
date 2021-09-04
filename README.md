# tclboot an R package
## Conditional inference in small sample scenarios using a resampling approach


This R package provides an implementation of a non-parametric resampling technique in the context of multidimensional hypothesis testing of assumptions of the Rasch model based on conditional distributions. It is suggested in small sample size scenarios, in which asymptotic theory is not applicable, to approximate the exact sampling distribution of various well-known Chi^2 test statistics like Wald, likelihood ratio, score and gradient tests as well as others. A procedure to compute the power function of the tests is also presented. A number of examples of scenarios are discussed in which the power function does not converge to one with an increasing deviation from the hypothesis to be tested, i.e. the respective assumption of the model. Finally, an attempt to modify the critical region of the test is made aiming at improving the power.


## Installation


The package can be installed using the `devtools`-package:

```
install.packages("devtools")
devtools::install_github("akurz1/tclboot")
```

## Figures 1 and 2

#### Figure 1:
Figure 1. Local power curve for tests U_1 and U_2 and sample sizes n = 10, 20, 30, 50. Power curve displayed as a function of δ_2 (item 2, an extreme item) and δ_4 (item 4, a moderate item). The nominal level α is set to 0.05 (dotted horizontal line).

```
res_U1 <- MCplot_simdata_ext(npers = c(10,20,30,50), k=6, nstats = "U1", xlim = 4,
+                            legend.position = "top",
+                            legend.direction = "horizontal", tcolor = TRUE) #tcolor = TRUE)

res_U2 <- MCplot_simdata_ext(npers = c(10,20,30,50), k=6, nstats = "U2", xlim = 4,
+                            legend.position = "top",
+                            legend.direction = "horizontal", tcolor = TRUE) #tcolor = TRUE)

plotlist <- list(
+   res_U1$plotlist$I2,
+   res_U1$plotlist$I4,
+   res_U2$plotlist$I2,
+   res_U2$plotlist$I4)

ggpubr:: ggarrange(plotlist = plotlist, ncol=2, nrow = 2, common.legend = TRUE, legend = "top")
```

#### Figure 2:
Local power curve for tests W, LR, RS and G and sample sizes n = 10, 20, 30, 50. Power curve displayed as a function of δ_2 (item 2, an extreme item) and δ_4 (item 4, a moderate item). The nominal level α is set to 0.05 (dotted horizontal line).

```
res1 <- MCplot_simdata_ext(npers = c(10,20,30,50), k=6, nstats = "W", xlim = 4,
                           legend.position = "top",
                           legend.direction = "horizontal", tcolor = TRUE)

res2 <- MCplot_simdata_ext(npers = c(10,20,30,50), k=6, nstats = "LR", xlim = 4,
                           legend.position = "top",
                           legend.direction = "horizontal", tcolor = TRUE)
                           
res3 <- MCplot_simdata_ext(npers = c(10,20,30,50), k=6, nstats = "RS", xlim = 4,
                           legend.position = "top",
                           legend.direction = "horizontal", tcolor = TRUE)

res4 <- MCplot_simdata_ext(npers = c(10,20,30,50), k=6, nstats = "G", xlim = 4,
                           legend.position = "top",
                           legend.direction = "horizontal", tcolor = TRUE)

l1 <- expression(italic("W")) # Wald
l2 <- expression(italic("LR"))  # abs
l3 <- expression(italic("RS"))  # RS
l4 <- expression(italic("G"))  # G

plotlist <- list(
  res1$plotlist$I2 + ggplot2::labs(tag = l1),
  res1$plotlist$I4,
  NULL,
  res2$plotlist$I2 + ggplot2::labs(tag = l2),
  res2$plotlist$I4,
  res3$plotlist$I2 + ggplot2::labs(tag = l3),
  res3$plotlist$I4,
  NULL,
  res4$plotlist$I2 + ggplot2::labs(tag = l4),
  res4$plotlist$I4)

ggpubr:: ggarrange(plotlist = plotlist,widths = c(1,1, 0.1, 1, 1),  ncol=5, nrow = 2, common.legend = TRUE, legend = "top")

```



## Figures 3 and 4 

Draxler, C., & Kurz, A. (submitted). Conditional inference in small sample scenarios. Stats(Special Issue "Re-sampling Methods for Statistical Inference of the 2020s"). Abstract Retrieved from https://www.mdpi.com/journal/stats/special_issues/resampling_statistics

#### Figure 3: 

Local power curve for the test U_1 and its modiﬁed version referring to an extreme item in a scenario with sample size n = 30 and number of items k = 10. The nominal level of α is set to .05 (dotted horizontal line).
```
res3 <- MCplot_simdata_crit(N=30, k=10, nstats = "U1", alpha = 0.05,
+                                xlim = 5, legend.position = "top", legend.direction = "horizontal",
+                                ctr_alpha = 1,
+                                ctr_mod = "high")

Fig3 <- res3$power_item$I2  # ggplot object
```

#### Figure 4: 

Local power curve for the test U_1 and its modiﬁed version referring to all (free) items in a scenario with sample size n = 30 and number of items k = 10. The nominal level of α is set to .05 (dotted horizontal line).
```
res4 <- MCplot_simdata_crit(N=30, k=10, nstats = "U1", alpha = 0.05,
                               xlim = 5, legend.position = c(.5, .85), legend.direction = "horizontal",
                               ctr_alpha = .2,
                               ctr_mod = "both")

ggpubr::ggarrange(plotlist = res3$power_item[-10], ncol=3, nrow = 3, common.legend = TRUE, legend = "top")
```

