To conduct the GC test using these scripts you need:

-[R](http://www.r-project.com] with the following packages:

  - APE

  - Phangorn

  - doParallel

  - foreach

The latter two packages allow parallel computation of likelihood optimisation, which makes the analysis much faster. The default is to use 6 processors, but this can be changed according to the system.

The sequence data sholud be in fasta format. Put the data in the same directory as the scripts gc_aa.R and gc_dna.R, and follow these instructions:

- Open a terminal window (in OSX or a Unix-like system), navigate to the folder with the sequence data and the scripts.

- Type the commands below using the example alignments provided, which consist of simulated sequences:

The basic command is:
Rscript gc_datatype.R alignment.fasta Model. For example, for DNA sequences:

```
Rscript gc_dna.R sim_jc.fasta GTR
```

GTR can be substituted by JC or by GTR+G. 


For amino acid data:

```
Rscript gc_aa.R sim_codon.fasta LG
```

LG can be substituted by LG+G, JTT, or JTT+G. 

The output is an Rdata file with the name of the model used. To see the results, open R and type the code below. For this example I will use the output from the nucleotide model.





```r
library(phangorn)
```

```
## Loading required package: ape
```

```r
load('GTR_sim_jc.Rdata')
run_1
```

```
## [[1]]
## [1] 5817.369
## 
## [[2]]
##   [1] 5818.149 5780.858 5884.129 5827.654 5805.620 5763.238 5787.836
##   [8] 5785.555 5822.102 5820.445 5813.145 5745.985 5760.055 5790.661
##  [15] 5873.559 5873.329 5742.587 5786.859 5709.623 5815.669 5807.970
##  [22] 5812.859 5725.112 5806.480 5783.277 5786.290 5824.713 5822.715
##  [29] 5858.824 5751.376 5816.987 5859.807 5724.665 5897.526 5775.786
##  [36] 5836.121 5869.086 5787.608 5825.904 5746.229 5796.300 5798.980
##  [43] 5840.050 5901.179 5741.204 5799.923 5837.900 5771.140 5857.100
##  [50] 5899.198 5749.023 5818.192 5803.883 5786.883 5848.533 5836.865
##  [57] 5898.650 5837.564 5891.216 5829.746 5746.061 5808.906 5776.649
##  [64] 5846.440 5764.435 5807.752 5713.025 5741.671 5824.918 5832.002
##  [71] 5854.719 5740.061 5783.044 5840.532 5846.704 5908.720 5857.089
##  [78] 5859.481 5816.694 5813.063 5865.848 5729.712 5761.930 5842.252
##  [85] 5722.080 5779.345 5763.580 5873.798 5750.314 5763.476 5812.707
##  [92] 5814.357 5758.394 5777.568 5874.901 5714.406 5887.876 5798.069
##  [99] 5810.548 5710.517
## 
## [[3]]
## [1] 0.59
## 
## [[4]]
## 
##  loglikelihood: -12719.58 
## 
## unconstrained loglikelihood: -6902.21 
## 
## Rate matrix:
##          a        c        g        t
## a 0.000000 1.243602 1.230352 1.087466
## c 1.243602 0.000000 1.132066 1.101655
## g 1.230352 1.132066 0.000000 1.000000
## t 1.087466 1.101655 1.000000 0.000000
## 
## Base frequencies:  
## 0.2466399 0.2631685 0.2522741 0.2379175
```


**run_1** is a list with the following items.

  - The delta statistic of the empirical data.

  - The delta statistics for the bootstrap replicates.

  - The **P** value of the delta statistic.

  - The optimised tree in pml format. For details, please see the help file for pml, using **?pml** at the R prompt.


## Changing the default settings

To conduct more bootstrap replicates, as used in our paper, open any of the scripts in a text editor and change the n_sims argument of the following line:

run_1 <- gc_test(data_file, parallel = T, nsims = 100, model = model_run)

It is also possible to turn off parallelisation by selecting parallel = F.

