## Posterior predictive simulations for protein sequence analyses in MrBayes using dirichlet priors and model averaging

To conduct these simulations you need:

- [R](http://www.r-project.com) with the following packages:

      - Phangorn


The script is called *pps_dpp_ma.R*. It contains a set of funcitons so simulate protein data for MrBayes analyses that use model averaging or informative priors on the substitution matrix. The latter approach is basically a GTR matrix with informative priors, as described by Huelsenbeck et al. (2008).

Use the example files in the R prompt as follows:

First load phangorn and the script with the simulation functions.

```r
library(phangorn)
source('pps_dpp_ma.R')
```

Load the data and the posterior samples of trees and parameter estimates for MrBayes. This first example is for Dirichlet priors on the rate matrix. To run a similar analysis refer to the MrBayes manual.


```r
trees_ma<- read.nexus('denv2_dpp.run1.t')
log_ma <- read.table('denv2_dpp.run1.p', head = T, skip = 1)
denv_data <- read.phyDat('example_aa_data.fasta', format = 'fasta', type = 'AA')
```

To obtain the posterior predicive simulations use the function pps_dpp. Then use multlik and homogen_test to calculate the multinomial likelihood and the chi-squared test statistics for every posterior predictive data set. Note that n_samples is the number of posterior predictive data sets to simulate. In practice, about 1000 data sets is useful, but for the purpose of this example we will use 100 only.

```r
sims_dpp <- pps_dpp(log_out = log_ma, tree_out = trees_ma, s_len = 3391, n_samples = 100)
```

```
## [1] "simulating with alpha =  0.170504"
## 110 sequences with 3388 character and 639 different site patterns.
## The states are a r n d c q e g h i l k m f p s t w y v 
## [1] "simulating with alpha =  0.139211"
```

```r
sims_multlik <- sapply(sims_dpp, function(x) multlik(x))
sims_chisq <- sapply(sims_dpp, function(x) homogen_test(x))
```
The P values can be calculated as follows:

```r
# For the multinomial likelihood
pval_multlik <- sum(sims_multlik > multlik(denv_data)) / length(sims_dpp)
print(pval_multlik)
```

```
## [1] 0.64
```

```r
# For the chi-squared
pval_chisq <- sum(sims_chisq > homogen_test(denv_data)) / length(sims_dpp)
print(pval_chisq)
```

```
## [1] 0.38
```

For model averaging, we can use a similar procedure. The log file sholud be like that shown here.


```r
trees_ma <- read.nexus('denv2_MA.run1.t')
log_ma <- read.table('denv2_MA.run1.p', head = T, skip = 1)

sims_ma <- pps_ma(log_out = log_ma, tree_out = trees_ma, s_len = 3391, n_samples = 100)
```

```
## [1] "simulated with settings 1634000 -17908.83 733.5583 0.2997991 0.1027881 1"
## 110 sequences with 3392 character and 865 different site patterns.
## The states are a r n d c q e g h i l k m f p s t w y v 
## [1] "simulated with settings 2048500 -17912.11 733.6556 0.2996645 0.1025421 1"
## 110 sequences with 3392 character and 911 different site patterns.
## The states are a r n d c q e g h i l k m f p s t w y v 
## [1] "simulated with settings 1215000 -17935.58 735.7994 0.2966951 0.1100538 1"
## 110 sequences with 3392 character and 868 different site patterns.
```

```r
sims_multlik <- sapply(sims_ma, function(x) multlik(x))
sims_chisq <- sapply(sims_ma, function(x) homogen_test(x))

pval_multlik <- sum(sims_multlik > multlik(denv_data)) / length(sims_ma)
pval_chisq <- sum(sims_chisq > homogen_test(denv_data)) / length(sims_ma)
```



## References

Huelsenbeck, J. P., Joyce, P., Lakner, C., & Ronquist, F. (2008). Bayesian analysis of amino acid substitution models. Philosophical Transactions of the Royal Society B: Biological Sciences, 363(1512), 3941-3953.
