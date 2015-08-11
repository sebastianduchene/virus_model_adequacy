library(ape)

denv <- read.dna('DENV2_NT.fasta', format = 'fasta')
pet <- read.dna('Petterson_et_al_NT.fasta', format = 'fasta')
mour <- read.dna('Moureau_et_al_NT.fasta', format = 'fasta')
hein <- read.dna('Heinze_et_al_NT.fasta', format = 'fasta')

plot_sites <- function(data_mat){
    bf_sites <- sapply(1:nrow(data_mat), function(x) base.freq(data_mat[x, ]))
    cols <- c('red', 'blue', 'orange', 'green')
    plot(bf_sites[1, ], type = 'l', col = cols[1], ylim = c(0.1, 0.35))
    for(i in 2:4){lines(bf_sites[i, ], col = cols[i] )}
}

par(mfrow = c(1, 4))
plot_sites(denv)
plot_sites(pet)
plot_sites(mour)
plot_sites(hein)
