library(phangorn)
source('pps_dpp_ma.R')

emp_data <- read.phyDat('denv_no_outliers_AA.fasta', format = 'fasta', type = 'AA')
emp_ml <- multlik(emp_data)
emp_chisq <- homogen_test(emp_data)

log_out <- read.table('denv2_dpp.run1.p', head = T, skip = 1)
log_out$alpha <- log_out$alpha + 0.2
tree_out <- read.nexus('denv2_dpp.run1.t')

s_dpp <- pps_dpp(log_out, tree_out, s_len = 3391, n_samples = 100)

mliks <- sapply(s_dpp,  function(x) multlik(x))
chisqs <- sapply(s_dpp, function(x) homogen_test(x) )



log_out <- read.table('denv2_MA.run1.p', head = T, skip = 1)
log_out$alpha <- log_out$alpha + 0.2
tree_out <- read.nexus('denv2_MA.run1.t')

s_dpp <- pps_ma(log_out, tree_out, s_len = 3391, n_samples = 20)

mliks <- sapply(s_dpp,  function(x) multlik(x))
chisqs <- sapply(s_dpp, function(x) homogen_test(x) )

