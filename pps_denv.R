library(phangorn)
source('pps_dpp_ma.R')


log_out <- read.table('DENV2_MA.run1.p', head = T, skip = 1)
tree_out <- read.nexus('DENV2_MA.run1.t')
s_len <- 3391

data <- read.phyDat('DENV2_aa.nexus', type = 'protein')



denv_ma <- pps_ma(log_out, tree_out, s_len, n_samples = 100)

ma_ml <- sapply(denv_ma, function(x) multlik(x))
emp_ml <- multlik(data)

emp_chisq <- homogen_test(data)

pps_chisq <- sapply(denv_ma, function(x) homogen_test(x))
