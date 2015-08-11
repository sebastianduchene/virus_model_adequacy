library(phangorn)

# JC

res_mat <- array(NA, dim = c(3, 10, 3))

dimnames(res_mat) <- list(c('JC', 'GTR', 'GTRG'), 1:10, c('P', 'Tdist', 'TLdiff')   )

JC_runs <- dir(pattern = '^JC.+')
for(r in 1:length(JC_runs)){
    load(JC_runs[r])
    res_mat[1, r,  1] <- run_1[[3]]
    tree <- read.tree(    gsub('Rdata', 'tree', gsub('JC_', '', JC_runs[r])))
    res_mat[1, r, 2] <- dist.topo(tree, run_1[[4]]$tree)
    res_mat[1, r, 3] <- abs(round((sum(tree$edge.length) - sum(run_1[[4]]$tree$edge.length)) / sum(tree$edge.length) , 2) * 100)
}

#GTR

GTR_runs <- dir(pattern = '^GTR_.+')
for(r in 1:length(GTR_runs)){
    load(GTR_runs[r])
    res_mat[2, r,  1] <- run_1[[3]]
    tree <- read.tree(    gsub('Rdata', 'tree', gsub('GTR_', '', GTR_runs[r])))
    res_mat[2, r, 2] <- dist.topo(tree, run_1[[4]]$tree)
    res_mat[2, r, 3] <- abs(round((sum(tree$edge.length) - sum(run_1[[4]]$tree$edge.length)) / sum(tree$edge.length) , 2) * 100)
}

#GTRG
GTRG_runs <- dir(pattern = '^GTRG.+')
for(r in 1:length(GTRG_runs)){
    load(GTRG_runs[r])
    res_mat[3, r,  1] <- run_1[[3]]
    tree <- read.tree(    gsub('Rdata', 'tree', gsub('GTRG_', '', GTRG_runs[r])))
    res_mat[3, r, 2] <- dist.topo(tree, run_1[[4]]$tree)
    res_mat[3, r, 3] <- abs(round((sum(tree$edge.length) - sum(run_1[[4]]$tree$edge.length)) / sum(tree$edge.length) , 2) * 100)
}

par(mfrow = c(3, 1))
for(i in 1:3){
    image(   t(res_mat[nrow(res_mat):1, ,i ])   )
}

res_means <- matrix(NA, 3, 3)
dimnames(res_means) <- list(c('JC', 'GTR', 'GTRG'), c('P', 'Tdist', 'TLdiff'))
for(z in 1:dim(res_mat)[3]){
    res_means[, z] <- rowMeans(res_mat[, , z])
}


write.table(res_means, file = 'long_trees_nt_results.csv', sep = ',')
