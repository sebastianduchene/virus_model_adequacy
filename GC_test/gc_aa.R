library(phangorn)

simulate_hetero <- function(ntax = 50, slen = 1000){

    # 5 taxa have different base frequencies
    tr1 <- rtree(ntax - 4)
    tr1$edge.length <- rlnorm(length(tr1$edge.length), meanlog = -1.5, sd = 0.3)
    tr1$tip.label <- paste0('A', 1:length(tr1$tip.label))
    s1 <- simSeq(tr1, l = slen, model = 'JTT', type = 'AA')

    tr2 <- rtree(5)
    tr2$edge.length <- rlnorm(length(tr2$edge.length), meanlog = -1.5, sd = 0.3)
    tr2$tip.label <- paste0('B', 1:length(tr2$tip.label))
    s2 <- simSeq(tr2, l = slen, model = 'LG', type = 'AA', rootseq = as.character(s1)[1, ] )

    tr3 <- bind.tree(tr1, tr2, where = 1)

    s3 <- rbind(as.character(s1)[-1, ], as.character(s2))
    s3 <- phyDat(s3, type = 'AA')

    return(list(tr3, s3))
}


gc_test <- function(aa_data, parallel = F, nsims = 10, rm_gaps = TRUE, model = NULL){
    require(phangorn)

    print(aa_data)

    concat_list <- function(c_list){
        if(length(c_list) == 2 ){
            return(c(c_list[[1]], c_list[[2]]))
        }else if(length(c_list) > 2){
            return(c(c_list[[1]], concat_list(c_list[-1])))
        }else{
            return(c_list)
        }
    }

    multlik <- function(al){
        if(!is.matrix(al)) al <- as.character(al)
        nsites <- ncol(al)
        al_patterns <- table(sapply(1:nsites, function(x) paste(al[, x], collapse = '')))
        return(sum(sapply(al_patterns, function(x) (log(x) * x))) - (nsites*log(nsites)))
    }

    rem_gaps <- function(aa_data){
        if(!is.matrix(aa_data)){
            aa_data <- as.character(aa_data)
        }
        has_gap <- function(x){
            aas <- c( 'a', 'r', 'n', 'd', 'c', 'q', 'e', 'g', 'h', 'i', 'l', 'k', 'm', 'f', 'p', 's', 't', 'w', 'y', 'v')
            return(!(any(aas %in% x) & !any(c('-', '?') %in% x) & length(unique(x)) <= 20))
        }
        gap_sites <- sapply(1:ncol(aa_data), function(d) has_gap(as.character(aa_data[, d])))
        return(phyDat(aa_data[, !gap_sites], type = 'AA'))
    }

    start_tree <- nj(dist.ml(aa_data))

    if(rm_gaps){
        aa_data <- rem_gaps(aa_data)
    }

    if(is.null(model)) stop('please specify a substitution model. It should be one of LG, WAG, JTT, LG+G, WAG+G, or JTT+G')

    if(model == 'LG'){
        opt_data <- function(aa_data) optim.pml(pml(tree = start_tree, data = aa_data, k = 1, model = 'LG'), optNni = T, model = 'LG')
    }else if(model == 'WAG'){
        opt_data <- function(aa_data) optim.pml(pml(tree = start_tree, data = aa_data, k = 1, model = 'WAG'), optNni = T, model = 'WAG')
    }else if(model == 'JTT'){
        opt_data <- function(aa_data) optim.pml(pml(tree = start_tree, data = aa_data, k = 1, model = 'JTT'), optNni = T, model = 'JTT')
    }else  if(model == 'LG+G'){
        opt_data <- function(aa_data) optim.pml(pml(tree = start_tree, data = aa_data, k = 4, model = 'LG'), optNni = T, optGamma = T, model = 'LG')
    }else if(model == 'WAG+G'){
        opt_data <- function(aa_data) optim.pml(pml(tree = start_tree, data = aa_data, k = 4, model = 'WAG'), optNni = T, optGamma = T, model = 'WAG')
    }else if(model == 'JTT+G'){
        opt_data <- function(aa_data) optim.pml(pml(tree = start_tree, data = aa_data, k = 4, model = 'JTT'), optNni = T, optGamma = T, model = 'JTT')
    }else{
        stop('please specify a substitution model. It should be one of LG, WAG, JTT, LG+G, WAG+G, or JTT+G')
    }

    s_len <- length(as.character(aa_data)[1, ])

    get_sim_rep <- function(mle){
        if(model == 'LG'){
            sim_dat <- simSeq(mle$tree, l =  s_len, model = 'LG', type = 'AA')
            sim_opt <- opt_data(sim_dat)
        }else if(model == 'WAG'){
            sim_dat <- simSeq(mle$tree, l = s_len, model = 'WAG', type = 'AA')
            sim_opt <- opt_data(sim_dat)
        }else if(model == 'JTT'){
            sim_dat <- simSeq(mle$tree, l = s_len, model = 'JTT', type = 'AA')
            sim_opt <- opt_data(sim_dat)
        }else if(model == 'LG+G'){
                    rates = phangorn:::discrete.gamma(mle$shape, k = 4)
                    sim_dat_all<- lapply(rates, function(r) simSeq(mle$tree, l = round(s_len/4, 0), rate = r, model = 'LG', type = 'AA'))
                    sim_dat <- concat_list(sim_dat_all)
                    sim_opt <- opt_data(sim_dat)
                }else if(model == 'WAG+G'){
                    rates = phangorn:::discrete.gamma(mle$shape, k = 4)
                    sim_dat_all<- lapply(rates, function(r) simSeq(mle$tree, l = round(s_len/4, 0), rate = r, model = 'WAG', type = 'AA'))
                    sim_dat <- concat_list(sim_dat_all)
                    sim_opt <- opt_data(sim_dat)
                }else if(model == 'JTT+G'){
                    rates = phangorn:::discrete.gamma(mle$shape, k = 4)
                    sim_dat_all<- lapply(rates, function(r) simSeq(mle$tree, l = round(s_len/4, 0), rate = r, model = 'JTT', type = 'AA'))
                    sim_dat <- concat_list(sim_dat_all)
                    sim_opt <- opt_data(sim_dat)
                }else{
                    stop('please specify a substitution model. It should be one of LG, WAG, JTT, LG+G, WAG+G, or JTT+G')
                }
 #
        print(sim_dat)
        print(sim_opt)
#
        return(multlik(sim_dat) - sim_opt$logLik)
    }

    mle_aa_data <- opt_data(aa_data)
    unc_lik_aa_data <- multlik(aa_data)
    con_lik_aa_data <- mle_aa_data$logLik
    delta_stat <- unc_lik_aa_data - con_lik_aa_data

#
    print(mle_aa_data)
#

    if(parallel){
        require(foreach)
        require(doParallel)
        cl <- makeCluster(6)
        registerDoParallel(cl)
        delta_sims <- foreach(x = 1:nsims, .packages = c('phangorn', 'ape'), .combine = c) %dopar% get_sim_rep(mle_aa_data)
        stopCluster(cl)
    }else{
        delta_sims<- vector()
        for(i in 1:nsims){
            delta_sims[i] <- get_sim_rep(mle_aa_data)
        }
    }

    return(list(delta_stat, delta_sims, sum(delta_stat > delta_sims) / nsims, mle_aa_data))
}




    concat_list <- function(c_list){
        if(length(c_list) == 2 ){
            return(c(c_list[[1]], c_list[[2]]))
        }else if(length(c_list) > 2){
            return(c(c_list[[1]], concat_list(c_list[-1])))
        }else{
            return(c_list)
        }
    }






library(phangorn)
library(methods)

args <- commandArgs(trailingOnly = TRUE)
file_name <- args[1]
model_run <- args[2]

data_file <- read.phyDat(file_name, format = 'fasta', type = 'AA')

run_1 <- gc_test(data_file, parallel = T, nsims = 100, model = model_run)

out_name <- gsub('[+]', '', paste0(model_run, '_',  gsub('[.].*$', '.Rdata',  file_name), collapse = ''))

save(run_1, file= out_name)

