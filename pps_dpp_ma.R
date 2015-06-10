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

concat_list <- function(c_list){
    if(length(c_list) == 2 ){
        return(cbind(c_list[[1]], c_list[[2]]))
    }else if(length(c_list) > 2){
        return(cbind(c_list[[1]], concat_list(c_list[-1])))
    }else{
        return(c_list)
    }
}

multlik <- function(al){
    if(class(al) != 'DNAbin') al <- as.character(al)
    if(!is.matrix(al)) stop('Please supply the sequences as a matrix')
    nsites <- ncol(al)
    al_patterns <- table(sapply(1:nsites, function(x) paste(al[, x], collapse = '')))
    return(sum(sapply(al_patterns, function(x) (log(x) * x))) - (nsites*log(nsites)))
}

homogen_test <- function(seq_data){
    if(!is.matrix(seq_data)) seq_data <- as.character(seq_data)
    sites <- unique(as.character(seq_data))
    count_sites <- function(seq) sapply(sites, function(x) sum(x == seq))
    contingency_table <- t(sapply(1:nrow(seq_data), function(x) count_sites(seq_data[x, ])))
    chi_test <- chisq.test(contingency_table)
#    print(chi_test)
    return(chi_test$statistic)
}

## Function to simulate PPS using gtr+g
pps_dpp <- function(log_out, tree_out, s_len, n_samples = 10){
    if(length(tree_out) != nrow(log_out)) stop('The trees and log file have different number of samples')
    p_samples <- sample(1:length(tree_out), n_samples)    

    q_params <- grep('^r.+', colnames(log_out))
    bf_params <- grep('^pi.+', colnames(log_out))
    pps <- list()
    for(i in 1:length(p_samples)){
        post_sample <- log_out[p_samples[i], ]
        post_tree <- tree_out[[p_samples[i]]]
        rates <- phangorn:::discrete.gamma(alpha = post_sample$alpha, k = 4)
        q <- as.numeric(post_sample[q_params])
        bf <- as.numeric(post_sample[bf_params])
        print(paste('simulating with alpha = ', post_sample$alpha))
        pps[[i]] <-  phyDat(concat_list(lapply(rates, function(x) as.character(simSeq(post_tree, Q = q, bf = bf, type = 'AA', l = s_len / 4, rate = x)))), type = 'AA')
	print(pps[[i]])
    }
    return(pps)
}

# Function to simulate PPS using model AA model averaging
pps_ma <- function(log_out, tree_out, s_len, n_samples = 10){
    models <- c('PP', 'JTT', 'Dayhoff', 'mtREV24', 'mtmam', 'WAG', 'RtREV', 'cpREV', 'VT', 'Blosum62')
    model_list <- list()
    model_list[['PP']][[1]] <- rep(1, 190)
    model_list[['PP']][[2]] <- rep(1/20, 20)

    for(m in models[-1]){
        model_params <- eval(parse(text = paste0('phangorn:::.', m)))
        model_list[[m]][[1]] <- model_params$Q
        model_list[[m]][[2]] <- model_params$bf
    }
    
    simulate_sample <- function(log_data, tree_data, sampled_index, al_len){
        log_sample <- log_data[sampled_index, ]
        sampled_model <- models[log_sample$aamodel]
        sampled_alpha <- log_sample$alpha
        sampled_tree <- tree_data[[sampled_index]]

        gamma_rates <- phangorn:::discrete.gamma(alpha = sampled_alpha, k = 4)
        simulated_data <- phyDat(concat_list(lapply(gamma_rates, function(x) simSeq(sampled_tree, l = round(al_len / 4, 0), Q = model_list[[sampled_model]][[1]], bf = model_list[[sampled_model]][[2]], type = 'AA'))), type = 'AA')
#        simulated_data <- phyDat(concat_list(lapply(gamma_rates, function(x) simSeq(sampled_tree, l = 3000, Q = model_list[[sampled_model]][[1]], bf = model_list[[sampled_model]][[2]], type = 'AA'))),type = 'AA' )
	return(simulated_data)
    }

    p_samples <- sample(1:length(tree_out), n_samples)
    pps <- list()
    for(i in 1:length(p_samples)){
    	  pps[[i]] <- simulate_sample(log_out, tree_out, p_samples[i], al_len = s_len)
	  print(paste(c('simulated with settings', log_out[p_samples[i], ]), collapse = ' '))
	  print(pps[[i]])
    }
    return(pps)
}


#library(phangorn)

#tr1 <- rtree(30)
#tr1$edge.length <- rlnorm(n = length(tr1$edge.length), meanlog = -3.5, sdlog = 0.3)
#Q <- c(1, 1, 1, 1, 1, rep(0.05, 185))
#Q <-  rep(1, 190)
#bf <- c(0.25, 0.25, 0.25, 0.25, rep(0, 16))

#aa1 <- simSeq(tr1, type = 'AA', l = 1000, Q = Q, bf = bf)
#write.phyDat(aa1, file= 'poisson_aa.fasta', format = 'fasta')
#write.nexus.data(aa1, file = 'poisson_aa.nexus')
#opt1 <- optim.pml(pml(tr1, data = aa1), optQ = T, optEdge = F)


#log_out <-  read.table('test_aa.log', head = T, skip = 1)
#tree_out <- read.nexus('test_aa.trees')
#s_len <- 1000

#pps_test <- simulate_pps_aa(log_out, tree_out, 1000, 5)

