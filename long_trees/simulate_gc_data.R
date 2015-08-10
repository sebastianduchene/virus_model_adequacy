library(phangorn)

# Simulate data sets with long trees

qmat <- c(0.1, 0.3, 0.4, 0.5, 0.8, 1)
freqs <- c(0.49, 0.1, 0.1, 0.49)
g_cats <- phangorn:::discrete.gamma(alpha = 1, k = 4)



sim_data <- function(tr, q, bf, rates, l){
    s <- list()

    for(r in 1:length(rates)){
        s[[r]] <-  as.DNAbin(simSeq(tr, l = l/4, Q = q, bf = bf, model = 'GTR', rate = rates[r]))
    }
    concat_list <- function(c_list){
        if(length(c_list) == 1){
            return(c_list)
        }else if(length(c_list) == 2){
            return(cbind(c_list[[1]], c_list[[2]]))
        }else{
            return(cbind(c_list[[1]], concat_list(c_list[-1])))
        }
    }

    return(concat_list(s))
}


# Simulate gtr+g long trees, all the same tree

for(i in 1:10){
    tr <- rtree(50)
    tr$edge.length <- rlnorm(length(tr$edge.length), meanlog = -1.0, sdlog = 0.5)
    sum(tr$edge.length)
    s <- sim_data(tr, qmat, freqs, g_cats, 200)
    write.tree(tr, file = paste('sim_', i, '.tree', sep = ''))
    write.dna(s, file  = paste('sim_', i, '.fasta', sep = ''), format = 'fasta')
}

