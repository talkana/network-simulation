## Pipeline for simulation of tree-child species networks and corresponding gene trees
### author: Natalia Rutecka
### Introduction
The pipeline has been developed to conduct a simulation study
described in Wawerka et al. 2021. It simulates tree-child
phylogenetic networks and corresponding gene trees.
The folder "parameters" contains a recommended set of 
biological parameters initially proposed in 
Molloy and Warnow 2020.
### Requirements 
The pipeline uses some functions from embretnet package, 
available at: https://bitbucket.org/pgor17/embretnet.

Before running the pipeline you need to download 
[SimPhy](https://github.com/adamallo/SimPhy) and add the main executable (simphy) and
INDELIble_wrapper.pl to your PATH. You should also download 
[INDELible](http://abacus.gene.ucl.ac.uk/software/indelible/download.php).

Depending on the choice of tree inference method you might need to download
[PhyML](http://www.atgc-montpellier.fr/phyml/download.php) (for ML method) and/or [ninja](https://wheelerlab.org/software/ninja/) (for NJ method).

### Usage
To perform simulations, run *paralsim.py*. 

Example usage: Simulate phylogenetic networks with numbers of leaves from {20, 40, 60} and numbers of reticulations from {4, 8}.
Simulate 100 networks for each parameter set and 1 gene tree per network. Use neighbour-joining method (recommended for large datasets) and the recommended biological parameters. Use 5 cores.

`python3 paralsim.py -l 20 40 60 -r 4 8 -n 100 -d 1 -ml False -p parameters -c 5 -o output_folder`

For more options, run:
`python3 paralsim.py -h`

### Questions
If you have any questions, feel free to contact me at n.rutecka@student.uw.edu.pl
