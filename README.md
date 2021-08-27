# tclboot
## Conditional inference in small sample scenarios


This R package provides an implementation of a non-parametric resampling technique in the context of multidimensional hypothesis testing of assumptions of the Rasch model based on conditional distributions. It is suggested in small sample size scenarios, in which asymptotic theory is not applicable, to approximate the exact sampling distribution of various well-known Chi^2 test statistics like Wald, likelihood ratio, score and gradient tests as well as others. A procedure to compute the power function of the tests is also presented. A number of examples of scenarios are discussed in which the power function does not converge to one with an increasing deviation from the hypothesis to be tested, i.e. the respective assumption of the model. Finally, an attempt to modify the critical region of the test is made aiming at improving the power.


## Installation


The package can be installed using the `devtools`-package:

```
install.packages("devtools")
devtools::install_github("akurz1/tclboot")
```
