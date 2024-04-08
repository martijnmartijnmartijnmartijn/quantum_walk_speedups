"""
@author: Martijn Brehm (m.a.brehm@uva.nl)
@date: 15/01/2024

This script takes as input CSV files containing experimental data created by
the "run_experiment.py" script and plots the results. For results of the
classical backtracking solver, plots the computed time complexity of the
quantum backtracking algorithm(s). For results of a modern classical solver,
plots its own time complexity.

Saving some runs of the script here for convenience.
Compare quantum and modern SAT solver on random satisfiable instances:
    python3 plot_queries_over_n.py data/BT-random_sat-3-10-40-50.csv    data/CDCL-random_sat-3-10-70-200.csv
Compare quantum and modern SAT solver on random unsatisfiable instances:
    python3 plot_queries_over_n.py data/BT-random_sat-3-10-40-50.csv data/BT-random_unsat-3-10-40-50.csv
Compare quantum and modern SAT solver on community satisfiable instances:
    python3 plot_queries_over_n.py data/BT-community_sat-3-10-30-50.csv
"""
from sys import argv
import numpy as np
import matplotlib.pyplot as plt
from math import log, exp, ceil
from gate_complexity.satcomplexity import bt
from gate_complexity.satcomplexity_grover import grover

# Algorithms included in this list will be plotted.
algs = {
    "m22" : 1,
    "Detection" : 8,
    "Binary search" : 9,
    "Grover" : 10,
    # "Detect opt $W$, $R=1$" : 12,
    # "B. search opt $W$, $R=n$" : 13,

}

# The ratio clauses/variables where a phase transition occurs for uniformly random 3-SAT, 4-SAT, ..., 13-SAT instances.
ratios = (4.267,9.931,21.117,43.37,87.79,176.54,354.01,708.92,1418.71,2838.28,
          5677.41,11355.67,22712.20)

# Indices into data from .csv files.
N_VARS = 0
REP = 1
TIME = 2
SIZE = 3
DEPTH = 4
N_SOLUTIONS = 5
N_SOLUTIONS_TOTAL = 6
FIRST_SOL_SIZE = 7
FIRST_SOL_DEPTH = 8

# Function for number of quantum queries of detection alg, and for rounding.
round_e = lambda number : "{0:.3e}".format(number).replace("e+0", "\cdot 10^") if len(str(number)) > 5 else str(number) # TODO Move to library?

# Define linestyles and colors that are used, and readable names for data types.
styles = ['-', '--', '-.', ':', 'None', ' ', '', 'solid', 'dashed', 'dashdot', 'dotted']

colors = ['#1f77b4', '#ff7f0e', '#2ca02c', '#d62728', '#9467bd', '#8c564b', '#e377c2', '#7f7f7f', '#bcbd22', '#17becf']
DELTA = 0.001

# Set up plot.
plt.figure(figsize=[13.5, 7])
plt.subplots_adjust(left=0.07, right=0.75)
plt.rcParams.update({"text.usetex":True, "font.family":"serif", "font.size":13})

plt.title("Running time of classical and quantum SAT solver for $\delta={}$".format(DELTA))
plt.ylabel("Number of seconds")
plt.xlabel("Number of variables in 3SAT instance")
plt.yscale("log")

# Load in the data from the given file.
print(argv[1:])
for i, file in enumerate(argv[1:]):
    print(file)
    name = file.split('/')[1].split('-')
    solver = name[0]
    k = name[1]
    sat = "Satisfiable" if name[2] == "sat" else "Unsatisfiable"
    mode = "Community $\\beta={},T={}$".format(name[3].split('_')[1], name[3].split('_')[2]) if "community" in name[3] else "Random"
    n1 = int(name[4])
    n2 = int(name[5])
    stepsize = int(name[6])
    reps = int(name[7].split('.')[0])
    dims = ((n2 - n1) // stepsize + 1, reps)
    data = np.loadtxt(open(file, "rb"), delimiter=",", skiprows=1).T

    # For CDCL, plot only classical, for BT plot only quantum algs.
    xs = range(n1, n2 + 1, stepsize)
    for j,alg in enumerate(algs):
        print(alg)
        if solver == "CDCL" and j > 0:
            break
        if solver == "BT" and j == 0:
            continue
        print("okay")
        # Extract complexity, compute (least squares fit of) median and std dev.
        d = np.reshape(data[algs[alg]], dims)
        ys = np.median(d, axis=1)
        ys_std_dev = np.std(d, axis=1)
        slope, intercept = np.polyfit(xs, np.log(ys), 1)
        label = "{} {} {}-SAT\n{}\n${}\cdot 2^".format(sat, mode, k, alg, round_e(exp(intercept))) + "{" + str(round(slope / np.log(2),3)) + "n}$"

        # Plot median complexity and fit over n.
        plt.scatter(xs, ys, marker="o", color=colors[j], linewidth=0.01)
        # plt.fill_between(xs, ys - ys_std_dev, ys + ys_std_dev, alpha=0.1, color=colors[j])
        # if (k > 6 and not "community" in title) or ("community" in title and T > 2.5):
        style = '-' if sat == "Satisfiable" else '--'
        print(style, sat)
        plt.plot(xs, np.exp(slope * xs + intercept), linestyle=style, label=label, color=colors[j])
        # else:
        #     plt.plot(xs, np.exp(slope * xs + intercept), linestyle='-', color=colors[j])

# Add ticks on x-axis.
xint = []
locs, labels = plt.xticks()
for each in locs:
    xint.append(int(each))
plt.xticks(range(min(xint),max(xint)+1,2))

# Render legend and show plot.
plt.legend(bbox_to_anchor=(1,1), loc="upper left")
plt.show()
