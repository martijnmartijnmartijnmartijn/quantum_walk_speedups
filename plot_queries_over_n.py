"""
@author: Martijn Brehm (m.a.brehm@uva.nl)
@date: 15/01/2024

This script takes as input CSV files containing experimental data created by
the "run_experiment.py" script and plots the results.

Specifically, for each file, the median runtime over n is computed for each algorithm (classical, detection, search, Grover) and an exponential fit to these medians is made. The medians and these exponential fits are then plotted over n.

Next,
Next, for varyinh

For results of the
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

round_e = lambda number : "{0:.3e}".format(number).replace("e-", "\cdot 10^{") + "}"

# Algorithms included in this list will be plotted.
algs = {
    "m22" : 1,
    "Detection" : 8, # 8 = query complexity 9 = t-depth 10 = t-coun
    "Binary search" : 11, # 11 = query complexity 12 = t-depth 13 = t-coun
    "Grover" : 14, # 14 = query complexity 15 = t-depth 16 = t-coun
}

# The ratio clauses/variables where a phase transition occurs for uniformly random 3-SAT, 4-SAT, ..., 15-SAT instances.
ratios = (4.267, 9.931, 21.117, 43.37, 87.79, 176.54, 354.01, 708.92, 1418.71, 2838.28, 5677.41, 11355.67, 22712.20)
S_PER_DAY = 86400

# Set up plot.
colors = ['#1f77b4', '#ff7f0e', '#2ca02c', '#d62728', '#9467bd', '#8c564b', '#e377c2', '#7f7f7f', '#bcbd22', '#17becf']
plt.figure(figsize=[13.5, 7])
plt.subplots_adjust(left=0.07, right=0.75)
plt.rcParams.update({"text.usetex":True, "font.family":"serif", "font.size":13})

DELTA = 0.001
plt.title("Running time of classical and quantum SAT solver for $\delta={}$".format(DELTA))
plt.ylabel("Number of seconds")
plt.xlabel("Number of variables in 3SAT instance")
plt.yscale("log")

# Save exp. fits by clause size k / temperature T per alg, for sat and unsat.
random_per_k = {alg : {"sat" : {}, "unsat" : {}} for alg in algs }
community_per_T = {alg : {"sat" : {}, "unsat" : {}} for alg in algs }

# Load in the data from the given file.
for i, file in enumerate(argv[1:]):
    print(file)
    name = file.split('/')[1].split('-')
    solver, k, n1, n2 = name[0], int(name[1]), int(name[4]), int(name[5])
    stepsize = int(name[6])
    reps = int(name[7].split('.')[0])
    dims = ((n2 - n1) // stepsize + 1, reps)
    data = np.loadtxt(open(file, "rb"), delimiter=",", skiprows=1).T

    sat = "Satisfiable" if name[2] == "sat" else "Unsatisfiable"
    mode = "Community $\\beta={},T={}$".format(name[3].split('_')[1], name[3].split('_')[2]) if "community" in name[3] else "Random"

    # For CDCL, plot only classical, for BT plot only quantum algs.
    xs = range(n1, n2 + 1, stepsize)
    for j,alg in enumerate(algs):
        if (solver == "CDCL" and j > 0) or solver == "BT" and j == 0:
            continue

        # Extract complexity, compute (least squares fit of) median.
        d = np.reshape(d, dims)
        ys = np.median(d, axis=1)
        ys_std_dev = np.std(d, axis=1)
        slope, intercept = np.polyfit(xs, np.log(ys), 1)
        print(slope, intercept)
        # label = "{} {} {}-SAT\n{}\n${}\cdot 2^".format(sat, mode, k, alg, round_e(exp(intercept))) + "{" + str(round(slope / np.log(2),3)) + "n}$"
        label = ""

        # Save exponential fit for later plots.
        random_per_k[alg][name[2]][k] = (slope, intercept)
        if "Community" in mode:
            community_per_T[alg][name[2]][name[3].split('_')[2]] = (slope, intercept)

        # Plot median complexity and fit over n.
        plt.scatter(xs, ys, marker="o", color=colors[j], linewidth=0.01)
        # if (k > 6 and not "community" in title) or ("community" in title and T > 2.5):
        style = '-' if sat == "Satisfiable" else '--'
        plt.plot(xs, np.exp(slope * xs + intercept), linestyle=style, label=label, color=colors[j])
        plt.annotate("$k={}$".format(k), xy=(xs[-1] + 2, ys[-1]), color=colors[j])

# Add ticks on x-axis, add legend, and show plot.
xint = []
locs, labels = plt.xticks()
for each in locs:
    xint.append(int(each))
plt.xticks(range(min(xint),max(xint)+1,10))
plt.legend(bbox_to_anchor=(1,1), loc="upper left")
plt.show()

# Plot max-n over k.
plt.title("Largest instance size solvable by each algorithm in one day for $\delta={}$".format(DELTA))
plt.ylabel("Number of variables $n$")
plt.xlabel("Clause size $k$")
plt.yscale("linear")
for j,alg in enumerate(algs):
    max_ns = []
    ks = list(random_per_k[alg]["sat"].keys())
    ks.sort()
    for k in ks:
        slope, intercept = random_per_k[alg]["sat"][k]
        max_ns.append((log(S_PER_DAY) - intercept) / slope)
    plt.plot(ks, max_ns, color=colors[j], label="Satisfiable\n" + str(alg))

    max_ns = []
    ks = list(random_per_k[alg]["unsat"].keys())
    ks.sort()
    for k in ks:
        slope, intercept = random_per_k[alg]["unsat"][k]
        max_ns.append((log(S_PER_DAY) - intercept) / slope)
        print(ks, max_ns)
    plt.plot(ks, max_ns, color=colors[j], label="Unsatisfiable\n" + str(alg), linestyle="--")
plt.legend()
plt.show()

# Plot speed-up over k.
plt.title("Speed-up of quantum algorithms over classical on largest instance solvable by quantum algorithm in one day for $\delta={}$".format(DELTA))
plt.ylabel("Speed-up factor")
plt.xlabel("Clause size $k$")
plt.yscale("log")
for j,alg in enumerate(algs):
    speed_ups = []
    ks = list(random_per_k[alg]["sat"].keys())
    ks.sort()
    print(ks)
    for k in ks:
        slope, intercept = random_per_k[alg]["sat"][k]
        max_n = (log(S_PER_DAY) - intercept) / slope
        slope, intercept = random_per_k["m22"]["sat"][k]
        speed_ups.append(exp(slope * max_n + intercept) / S_PER_DAY)
    plt.plot(ks, speed_ups, color=colors[j], label="Satisfiable\n" + str(alg))
    speed_ups = []
    ks = list(random_per_k[alg]["unsat"].keys())
    ks.sort()
    for k in ks:
        slope, intercept = random_per_k[alg]["unsat"][k]
        max_n = (log(S_PER_DAY) - intercept) / slope
        slope, intercept = random_per_k["m22"]["unsat"][k]
        speed_ups.append(exp(slope * max_n + intercept) / S_PER_DAY)
    plt.plot(ks, speed_ups, color=colors[j], label="Unsatisfiable\n" + str(alg), linestyle="--")
plt.legend()
plt.show()
