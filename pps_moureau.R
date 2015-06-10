library(phangorn)
library(methods)
source('pps_dpp.R')

par(mfrow = c(2, 2))
par(bg = 'black')
par(fg = 'white')
par(col.lab = 'white')
par(col.axis = 'white')





hein_data <- rem_gaps(phyDat(read.nexus.data('Heinze_dpp.nexus'), type = 'AA'))
hein_slen <- ncol(as.character(hein_data))

# Model averaging of the Q matrix
hein_ma_log <- read.table('../mb_runs/Heinze_mixed.p', head = T, skip = 1)
hein_ma_tre <- read.nexus('../mb_runs/Heinze_mixed.t')

s_ma <- pps_ma(hein_ma_log, hein_ma_tre, hein_slen, n_samples = 100)
ml_ma <- sapply(s_ma, function(x) multlik(x))
x2_ma <- sapply(s_ma, function(x) homogen_test(x))
hist(ml_ma)
hist(x2_ma)

# Dirichlet prior on Q matrix
hein_dpp_log <- read.table('Heinze_dpp.run1.p', head = T, skip = 1)
hein_dpp_tre <- read.nexus('Heinze_dpp.run1.t')

s_dpp <- pps_dpp(hein_dpp_log, hein_dpp_tre, hein_slen, n_samples = 100)
ml_dpp <- sapply(s_dpp, function(x) multlik(x))
x2_dpp <- sapply(s_dpp, function(x) homogen_test(x))
hist(ml_dpp)
hist(x2_dpp)


stop('analysing heinze')







mour_data <- rem_gaps(phyDat(read.nexus.data('Moureau_dpp.nexus'), type = 'AA'))
mour_slen <- ncol(as.character(mour_data))

# Model averaging of the Q matrix
mour_ma_log <- read.table('../mb_runs/Moureau_mixed.p', head = T, skip = 1)
mour_ma_tre <- read.nexus('../mb_runs/Moureau_mixed.t')

s_ma <- pps_ma(mour_ma_log, mour_ma_tre, mour_slen, n_samples = 100)
ml_ma <- sapply(s_ma, function(x) multlik(x))
x2_ma <- sapply(s_ma, function(x) homogen_test(x))
hist(ml_ma)
hist(x2_ma)

# Dirichlet prior on Q matrix
mour_dpp_log <- read.table('Moureau_dpp.run1.p', head = T, skip = 1)
mour_dpp_tre <- read.nexus('Moureau_dpp.run1.t')

s_dpp <- pps_dpp(mour_dpp_log, mour_dpp_tre, mour_slen, n_samples = 100)
ml_dpp <- sapply(s_dpp, function(x) multlik(x))
x2_dpp <- sapply(s_dpp, function(x) homogen_test(x))
hist(ml_dpp)
hist(x2_dpp)

stop('analysing Moureau')






pet_data <- rem_gaps(phyDat(read.nexus.data('Petterson_dpp.nexus'), type = 'AA'))
pet_slen <- ncol(as.character(pet_data))

# Model averaging of the Q matrix
pet_ma_log <- read.table('../mb_runs/Petterson_mixed.p', head = T, skip = 1)
pet_ma_tre <- read.nexus('../mb_runs/Petterson_mixed.t')

s_ma <- pps_ma(pet_ma_log, pet_ma_tre, pet_slen, n_samples = 100)
ml_ma <- sapply(s_ma, function(x) multlik(x))
x2_ma <- sapply(s_ma, function(x) homogen_test(x))
hist(ml_ma)
hist(x2_ma)

# Dirichlet prior on Q matrix
pet_dpp_log <- read.table('Pettersson_dpp.run1.p', head = T, skip = 1)
pet_dpp_tre <- read.nexus('Pettersson_dpp.run1.t')

s_dpp <- pps_dpp(pet_dpp_log, pet_dpp_tre, pet_slen, n_samples = 100)
ml_dpp <- sapply(s_dpp, function(x) multlik(x))
x2_dpp <- sapply(s_dpp, function(x) homogen_test(x))
hist(ml_dpp)
hist(x2_dpp)


stop('analysing Pettersson')







