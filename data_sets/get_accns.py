import re, os, sys

print sys.argv

lines = open(sys.argv[1], 'r').readlines()

taxon_names = [i for i in lines if '>' in i]

accns = list()

for n in taxon_names:
    n = re.sub('NC_', 'NC', n)
    n_temp = re.split('_|[|]', n)
    accns.append(re.sub('>', '', n_temp[0])+', ')
    print n_temp[0]

open(re.sub('[.].*', '_accns.txt', sys.argv[1]), 'w').writelines(accns)


