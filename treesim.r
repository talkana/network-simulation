library("TreeSim")
library(argparser, quietly=TRUE)
p = arg_parser("Simulate a species tree with a given number of leaves and other parameters from Molloy and Warnow 2020")
p = add_argument(p, "--l", help="number of leaves", type="integer")
argv = parse_args(p)
br = 1.8*10**(-9) #birth rate
dr = 0 #death rate
height = 1800000337.5
tree = sim.bd.taxa.age(argv$l, 1, br, dr, 1, height, T)
newick = write.tree(tree[[1]])
print(newick)
