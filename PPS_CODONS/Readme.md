# Posterior predictive simulations for codon models 

To replicate these analyses you need:

- Python (>2.7) with the following modules

    - [numpy](http://www.numpy.org/)
    
    - [pandas](http://pandas.pydata.org)
    
    - [dendropy](https://pythonhosted.org/DendroPy/)
    
    - [pyvolve](http://sjspielman.org/pyvolve/)
    
There are a few steps to these analyses, which I have explained in [this Ipython notebook](). The first requirement, however, is to analyse the data using a codon model in MrBayes, in this paper we have used the M3 model, but there are others available. An example input file for MrBayes in nexus format can be found [here](https://raw.githubusercontent.com/sebastianduchene/virus_model_adequacy/master/PPS_CODONS/example_M3_run.nexus).


