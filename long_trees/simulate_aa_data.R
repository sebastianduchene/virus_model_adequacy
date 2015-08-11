library(phangorn)

# Simulate data sets with long trees


sim_data <- function(tr, rates, l){
    s <- list()
    for(r in 1:length(rates)){
        s[[r]] <-  simSeq(tr, l = l/4, model = 'JTT', rate = rates[r], type = 'AA')
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


# Simulate JTT+g long trees, all the same tree

g_cats <- phangorn:::discrete.gamma(alpha = 1, k = 4)

for(i in 1:10){
    tr <- rtree(50)
    tr$edge.length <- rlnorm(length(tr$edge.length), meanlog = -0.8, sdlog = 0.5)
    print(sum(tr$edge.length))
    s <- sim_data(tr, g_cats, 200)
    write.tree(tr, file = paste('sim_AA_', i, '.tree', sep = ''))
    write.phyDat(s, file  = paste('sim_AA_', i, '.fasta', sep = ''), format = 'fasta')
}

