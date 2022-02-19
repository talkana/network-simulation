import os
import multiprocessing as mp
import argparse
import subprocess
import threading

parser = argparse.ArgumentParser(
    description='Run singlesim.py in parallel for different numbers of leaves and reticulations')
parser.add_argument('-c', required=True, type=int, help='Number of cores to use')
parser.add_argument('-o', required=True, help='Name of output folder')
parser.add_argument('-n', required=True, type=int, help='Number of networks per category')
parser.add_argument('-p', required=True, help='Folder with parameters')
parser.add_argument('-l', required=True, nargs='+', help='Numbers of leaves (space separated)')
parser.add_argument('-r', required=True, nargs='+', help='Numbers of reticulations (space separated)')
parser.add_argument('-d', required=True, type=int, help='Numbers of displayed trees per network')
parser.add_argument('-g', required=False, type=int, default=1, help='Numbers of gene trees per displayed tree')
parser.add_argument('-ml', default=False, type=bool, required=False,
                    help='should i use maximum likelihood method to infer trees from multiple sequence alignment (if false neighbour-joing method is used)')
parser.add_argument('-s', type=int, required=False, default=414, help='Random seed for sequence simulation')
args = parser.parse_args()
os.mkdir(args.o)

# make queue with all parameters
q = mp.Queue()
for p in os.listdir(args.p):
    if p != "indelible.txt":
        ppath = args.o + "/" + p
        os.mkdir(ppath)
        for n in args.l:
            for r in args.r:
                if int(n) > int(r):  # TC networks
                    dn = str(args.d)
                    g = str(args.g)
                    outfile = ppath + "/" + "n" + n + "r" + r
                    params = ["-t", n, "-r", r, "-o", outfile, "-g", args.p + "/" + p, "-dn", dn, "-d", g, "-n", str(args.n), "-ml", str(args.ml), "-s", str(args.s)]
                    print(params)
                    q.put(params)


# run in parallel
def worker():
    while True:
        if q.empty():
            return
        par = q.get()
        command = ["python3", "singlesim.py", "-a", "parameters/indelible.txt"]
        command.extend(par)
        subprocess.run(command)


for i in range(args.c):
    t = threading.Thread(target=worker)
    t.start()
