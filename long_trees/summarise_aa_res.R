library(phangorn)

# JC

res_mat <- array(NA, dim = c(4, 10, 3))

dimnames(res_mat) <- list(c('LG', 'LGG', 'JTT', 'JTTG'), 1:10, c('P', 'Tdist', 'TLdiff')   )

LG_runs <- dir(pattern = '^LG_.+')
for(r in 1:length(LG_runs)){
    load(LG_runs[r])
    res_mat[1, r,  1] <- run_1[[3]]
    tree <- read.tree(    gsub('Rdata', 'tree', gsub('LG_', '', LG_runs[r])))
    res_mat[1, r, 2] <- dist.topo(tree, run_1[[4]]$tree)
    res_mat[1, r, 3] <- abs(round((sum(tree$edge.length) - sum(run_1[[4]]$tree$edge.length)) / sum(tree$edge.length) , 2) * 100)
}

#GTR

LGG_runs <- dir(pattern = '^LGG_.+')
for(r in 1:length(LGG_runs)){
    load(LGG_runs[r])
    res_mat[2, r,  1] <- run_1[[3]]
    tree <- read.tree(    gsub('Rdata', 'tree', gsub('LGG_', '', LGG_runs[r])))
    res_mat[2, r, 2] <- dist.topo(tree, run_1[[4]]$tree)
    res_mat[2, r, 3] <- abs(round((sum(tree$edge.length) - sum(run_1[[4]]$tree$edge.length)) / sum(tree$edge.length) , 2) * 100)
}

#JTT
JTT_runs <- dir(pattern = '^JTT_.+')
for(r in 1:length(JTT_runs)){
    load(JTT_runs[r])
    res_mat[3, r,  1] <- run_1[[3]]
    tree <- read.tree(    gsub('Rdata', 'tree', gsub('JTT_', '', JTT_runs[r])))
    res_mat[3, r, 2] <- dist.topo(tree, run_1[[4]]$tree)
    res_mat[3, r, 3] <- abs(round((sum(tree$edge.length) - sum(run_1[[4]]$tree$edge.length)) / sum(tree$edge.length) , 2) * 100)
}

#JTTG
JTTG_runs <- dir(pattern = '^JTTG_.+')
for(r in 1:length(JTTG_runs)){
    load(JTTG_runs[r])
    res_mat[4, r,  1] <- run_1[[3]]
    tree <- read.tree(    gsub('Rdata', 'tree', gsub('JTTG_', '', JTTG_runs[r])))
    res_mat[4, r, 2] <- dist.topo(tree, run_1[[4]]$tree)
    res_mat[4, r, 3] <- abs(round((sum(tree$edge.length) - sum(run_1[[4]]$tree$edge.length)) / sum(tree$edge.length) , 2) * 100)
}



par(mfrow = c(3, 1))
for(i in 1:3){
    image(   t(res_mat[nrow(res_mat):1, ,i ])   )
}

res_means <- matrix(NA, 4, 3)
dimnames(res_means) <- list(c('LG', 'LGG', 'JTT', 'JTTG'), c('P', 'Tdist', 'TLdiff'))
for(z in 1:dim(res_mat)[3]){
    res_means[, z] <- rowMeans(res_mat[, , z])
}


write.table(res_means, file = 'long_trees_aa_results.csv', sep = ',')
