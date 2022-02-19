import subprocess
import argparse
import os
from retnetsim import simulate
import random

def parse_input():
    parser = argparse.ArgumentParser(
        description='Builds random reticulation network and generates corresponding ML gene trees with parameters from Molloy and Warnow 2020')
    parser.add_argument('-g', required=True, help='configuration file for gene tree sampling (simphy)')
    parser.add_argument('-a', required=True, help='configuration file for building alignment (indelible)')
    parser.add_argument('-o', required=True, help='name of the output folder')
    parser.add_argument('-r', type=int, required=True, help='number of reticulations')
    parser.add_argument('-t', type=int, required=True, help='number of tips in species network')
    parser.add_argument('-n', type=int, required=True, help='number of species networks to generate')
    parser.add_argument('-d', type=int, required=True, help='number of gene trees per displayed tree')
    parser.add_argument('-dn', type=int, required=True, help='number of displayed trees to use per each network')
    parser.add_argument('-ml', default=False, type=bool, required=False, help='should i use maximum likelihood method to infer trees from multiple sequence alignment (if false neighbour-joing method is used)')
    parser.add_argument('-s', default=414, type=int, required=False, help='Random seed for sequence simulation')
    parsed = parser.parse_args()
    return parsed



def change_taxa_names(file):
    sequences = open(file).read().split("\n")
    newpath = file + "_renamed"
    new = open(newpath, "w")
    new.write(sequences[0] + "\n")
    for seq in sequences[1:]:
        if seq and seq[0] == "t":
            id = seq.split()[0]
            new_id = id.split("_")[0]
            seq = seq.replace(id, new_id)
            new.write(seq + "\n")
    new.close()
    os.remove(file)
    os.rename(newpath, file)

def run_ML(path):
    for file in os.listdir(path):
        if ".phy" in file:
            filepath = path + "/" + file
            change_taxa_names(filepath)
            subprocess.run(["phyml", "-i", filepath, "-b", "-4", "-m", "GTR", "-a", "e", "--quiet"])

def run_NJ(path):
    for file in os.listdir(path):
        if ".phy" in file:
            filepath = path + "/" + file
            change_taxa_names(filepath)
            new_path = filepath.replace("_TRUE.phy", ".fasta")
            tree_path = path + "/tree_NJ"
            subprocess.run(["ninja", "--in", new_path, "--out", tree_path, "--alph_type", "d"])


def save_nexus(newlist, path):
    f = open(path, "w")
    f.write("#NEXUS\nbegin trees;\n")
    for i in range(len(newlist)):
        f.write("\ttree " + str(i) + "=" + newlist[i] + ";\n")
    f.write("end;")
    f.close()


def is_good_alignment(reppath):
    for folder in os.listdir(reppath):
        if folder[-1] in "0123456789":
            for file in os.listdir(reppath + "/" + folder):
                if ".phy" in file:
                    filepath = reppath + "/" + folder + "/" + file
                    sequences = open(filepath).read().split()[3::2]
                    if len(sequences) != len(set(sequences)):
                        return False
    return True


def has_positive_lengths(path):
    i = 0
    newick = open(path).read()
    while i < len(newick):
        if newick[i] == ":":
            i += 1
            st = i
            while newick[i].isalnum() or newick[i] in ".-+":
                i += 1
            bl = float(newick[st:i])
            if bl < 0: return False
        else:
            i += 1
    return True


def is_good(reppath):
    for folder in os.listdir(reppath):
        if folder[-1] in "0123456789":
            good = False
            for file in os.listdir(reppath + "/" + folder):
                if "phyml_tree" in file:
                    good = True
                elif "tree_NJ" in file and has_positive_lengths("/".join([reppath, folder, file])):
                    good = True
            if not good:
                return False
    return True

def main():
    args = parse_input()
    os.mkdir(args.o)
    for n in range(1, args.n + 1):
        reppath = args.o + "/" + str(n)
        nw = simulate(args.t, args.r)
        os.mkdir(reppath)

        # save network
        nwpath = reppath + "/" + "network"
        nwf = open(nwpath, "w")
        nwf.write(nw)
        nwf.close()

        # generate and save displayed trees
        process = subprocess.Popen(['python3 embretnet/embnet.py', "-n", nw, "-pnd"], stdout=subprocess.PIPE)
        (output, error) = process.communicate()
        output = output.decode()
        dtr = []
        if error:
            print(error)
            raise ChildProcessError
        else:
            sp = output.splitlines()
            for elt in sp[1:]:
                if elt[0] == "(":
                    dtr.append(elt)
            dtr = random.sample(dtr, args.dn)
            dpath = reppath + "/" + "displayed_trees"
            save_nexus(dtr, dpath)

            all_done = False
            while not all_done:
                # run gene tree and sequence simulation
                simphy_args = open(args.g).read().split()
                simphy_com = ['simphy', "-o", reppath, "-sr", dpath, "-rl",
                              "f:" + str(args.d), "-rs", str(len(dtr))]
                simphy_com.extend(simphy_args)
                subprocess.call(simphy_com)
                subprocess.call(['INDELIble_wrapper.pl', reppath, args.a, str(args.s), "1"])
                # filter files with good score and run ML on them
                for folder in os.listdir(reppath):
                    if folder[-1] in "0123456789":
                        if args.ml:
                            run_ML(reppath + "/" + folder)
                        else:
                            run_NJ(reppath + "/" + folder)
                all_done = is_good(reppath)


if __name__ == "__main__":
    main()
