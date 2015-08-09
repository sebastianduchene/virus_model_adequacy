#
def fasta_to_matrix(file_name):
    import numpy as np
    import re
    raw_lines = open(file_name, 'r').readlines()
    seqs = [re.sub('\n', '', i) for i in raw_lines if len(re.findall('^>', i)) == 0]
    seq_mat = np.empty(shape = (len(seqs), len(seqs[0])), dtype = str)
    
    for s in range(len(seqs)):
        seq_mat[s, :] = list(seqs[s])
    return seq_mat

#
def dict_to_matrix(seq_dict):
    import numpy as np
    import re
    
    seq_mat = np.empty(shape = (len(seq_dict), len(seq_dict[seq_dict.keys()[0]])), dtype = str)
    for s in range(len(seq_dict)):
        seq_mat[s, :] = list(seq_dict[seq_dict.keys()[s]])
    return seq_mat

#
def multlik(seq_mat):
    from scipy import stats
    nsites = seq_mat.shape[1]
    patterns = stats.itemfreq([''.join(seq_mat[:, i]) for i in range(nsites)])
    patterns_al = [float(i) for i in patterns[:, 1]]
    return sum([(np.log(i) * i) for i in patterns_al]) - nsites*np.log(nsites)

#
def chisq(seq_mat):
    from scipy import stats
    sites = list(set(seq_mat.ravel()))
    
    def count_sites(seq, sites):
        counts = list()
        for s in sites:
            counts.append(sum([i == s for i in seq]))
        return counts

    contingency_table = np.empty(shape = (seq_mat.shape[0], len(sites)))

    for seq in range(seq_mat.shape[0]):
        contingency_table[seq, :] = count_sites(seq_mat[seq, :], sites)

    return stats.chi2_contingency(contingency_table)[0]

#
def get_params_gtrg(sampled_params, line):
    import numpy
    rates = list(sampled_params.ix[line, [i for i in sampled_params.columns if len(re.findall('^r[(].*', i)) > 0  ]])
    freqs = np.array((sampled_params.ix[line, [i for i in sampled_params.columns if len(re.findall('^pi[(].*', i)) > 0 ]  ]))
    freqs = list(freqs / freqs.sum())
    alpha = sampled_params.ix[line, 'alpha']
    return (rates, freqs, alpha)

#
def get_params_codons(sampled_params, line):    
    freqs = np.array(sampled_params.ix[line, [i for i in sampled_params.columns if len(re.findall('pi[(][A-Z]', i)) > 0]])
    freqs = list(freqs / freqs.sum())
    omegas = list(sampled_params.ix[line, [i for i in sampled_params.columns if len(re.findall('omega.*', i)) > 0]])
    return (freqs, omegas)

#
def sim_gtr(rates, freqs, alpha, tree, nsites):
    import pyvolve
    custom_mu = {}
    
    for r, val in zip(['AC', 'AG', 'AT', 'CG', 'CT', 'CG'], rates):
        custom_mu[r]= val
    
    gtr_model = pyvolve.Model('nucleotide', {'mu':custom_mu, 'state_freqs':freqs}, alpha = alpha, num_categories = 4)
    gtr_partition = pyvolve.Partition(models = gtr_model, size = nsites)
    tr = pyvolve.read_tree(tree = tree)
    gtr_evolver = pyvolve.Evolver(partitions = gtr_partition, tree = tr)
    gtr_evolver()
    return dict_to_matrix(gtr_evolver.get_sequences())

#
def sim_codon(freqs, omegas, tree, nsites):
    import pyvolve
    
    tr = pyvolve.read_tree(tree = tree)

    #temporary convert tree
    
    gy_model = pyvolve.Model('MG', {'omega':omegas, 'state_freqs':freqs})
    # Note that the number of sites should be divided by 3!
# test using three partitions
    gy_partition = pyvolve.Partition(models = gy_model, size = nsites/3)
    gy_evolver = pyvolve.Evolver(partitions = gy_partition, tree = tr)
    gy_evolver()
    return dict_to_matrix(gy_evolver.get_sequences())

#
def var_sites(mat_dat):
    return sum([len(set(''.join(mat_dat[:, i]))) > 1 for i in range(mat_dat.shape[1])]) / float(mat_dat.shape[1])

#
def fix_tree(tree):
    import dendropy
    tr = dendropy.Tree.get(data = tree, schema = 'newick')
    for edge in tr.postorder_edge_iter():
        if not edge.length is None:
            edge.length = edge.length * 3
    return tr.as_string(schema = 'newick')

#
def find_tree(target_tree, tree_list, sch1 = 'newick', sch2 = 'nexus'):
    '''
    Note that target_tree and tree_list must
    '''
    import dendropy as dp
    tr = dp.Tree.get_from_path(target_tree, schema = sch1, rooting = 'force-unrooted')
    tr.update_bipartitions()    
    trs = dp.TreeList.get_from_path(tree_list, schema = sch2, taxon_namespace = tr.taxon_namespace, rooting = 'force-unrooted')

    for t in trs:
        t.update_bipartitions()
    
    tree_dists = list()
    
    for t_post in trs:
        tree_dists.append(dp.calculate.treecompare.unweighted_robinson_foulds_distance(tr, t_post))
        
    tree_lens = [t.length() for t in trs]
    true_len = tr.length()

#    return [any([d == 0 for d in tree_dists]), true_len < max(tree_lens) and true_len > min(tree_lens)]
    return [any([d == 0 for d in tree_dists]), abs(true_len - mean(tree_lens)) / true_len ]


