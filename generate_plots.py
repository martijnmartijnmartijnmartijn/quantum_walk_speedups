from os import system, listdir

files = ["data/" + f for f in listdir("data/")]

# Parameter ranges s.t. we can solve 30 reps of the largests instances in 5m.
# ranges_BT_random = {
#     3 : (10, 49, 30, 3), # n=46 1m36s   n=49 4m0s
#     4 : (10, 40, 30, 3), # n=37 3m8s    n=40 7m40s
#     5 : (10, 34, 30, 3), # n=31         n=34 17m5s
#     6 : (10, 28, 30, 3), # n=25 2m58s   n=28 9m50s
#     7 : (10, 22, 30, 3), # n=19 44s     n=22 ...
#     8 : (10, 22, 30, 3), # n=19 3m45s   n=22 12m37s
#     9 : (13, 19, 30, 3), # n=16 1m45s   n=19 9m37s
#     10 : (13, 16, 30, 3),# n=13 1m6s    n=16 6m34s
#     11 : (13, 16, 30, 3),# n=13 2m8s    n=16 15m2s
#     12 : (16, 40, 10, 3),# n=16 8m57s               15, 25, 7, 1
#     13 : (16, 40, 10, 3),# n=16 115m                15, 20, 5, 1
#     14 : (16, 40, 10, 3),# n=16 269m
#     # 15 : (16, 40, 10, 3), # 16, 20, 5, 1
# }
# ranges_CDCL_random = {
#     3 : (10, 253, 100, 3),# n=250 4m47s n=253 5m37s
#     4 : (10, 91, 100, 3), # n=88 4m27s  n=91 6m36s
#     5 : (10, 58, 100, 3), # n=55 4m33s  n=58 8m56s
#     6 : (10, 43, 100, 3), # n=40 3m1s   n=43 8m1s
#     7 : (10, 31, 100, 3), # n=31 3m15s  n=34 5m59s
#     8 : (10, 31, 100, 3), # n=28 2m8s   n=31 9m41s
#     9 : (13, 100, 100, 3), # n=19 24s ... stopped this accidentally
#     10 : (13, 100, 100, 3),# n=22 1m59s     n=25 6m18s
#     11 : (13, 100, 100, 3),# n=19 2m9s      n=22 37m38s
#     12 : (16, 100, 100, 3),# n=16 33m40s
#     13 : (16, 100, 100, 3),# n=16 21m30s
#     14 : (16, 100, 100, 3),# n=16 56m44s
#     # 15 : (16, 100, 100, 3),
# }

# 1. Compare algs on 8-SAT (where speed-ups exist on random instances) for data with increasing community structure T=0,0.25,0.75,
# for k in ranges_CDCL_random:
#     print(k)
#     n1, n2, reps, stepsize = ranges_BT_random[k]
#     if not any([f for f in files if "BT-{}-".format(k) in f]):
#         system("python3 run_experiment.py {} {} {} {} {} BT".format(k, n1, n2, stepsize, reps))
    # n1, n2, reps, stepsize = ranges_CDCL_random[k]
    # if not any([f for f in files if "CDCL-{}-".format(k) in f]):
    #     system("python3 run_experiment.py {} {} {} {} {} CDCL".format(k, n1, n2, stepsize, reps))
# system("python3 plot_queries_over_n.py " + " ".join([f for f in files if "random" in f and "unsat" in f]))
# system("python3 plot_queries_over_n.py " + " ".join([f for f in files if "random" in f and "-sat-" in f]))
system("python3 plot_queries_over_n.py " + " ".join([f for f in files if "random" in f]))

# 2. Compare algs on

# random_sat = [f for f in files if "random" in f and "-sat-" in f]
# random_unsat = [f for f in files if "random" in f and "-unsat-" in f]
# system("python3 plot_queries_over_n.py " + " ".join(random_sat))
# system("python3 plot_queries_over_n.py " + " ".join(random_unsat))

# # Compare algs on random k-SAT for k=3,4,...,14.
# random_sat = [f for f in files if "random" in f and "-sat-" in f]
# random_unsat = [f for f in files if "random" in f and "-unsat-" in f]
# system("python3 plot_queries_over_n.py " + " ".join(random_sat))
# system("python3 plot_queries_over_n.py " + " ".join(random_unsat))


# params = [(0, 0), (0, 1), (0, 2), (0, 4), (0, 100), (0.5, 0), (0.5, 1), (0.5, 2), (0.5, 4), (0, 100)]

# for beta, T in params:
#     print(beta, T)
#     # system("python3 run_experiment.py 3 10 30 community_sat BT 0 0 {} {}".format(beta, T))


# for T in [0, 0.25, 0.5, 0.75, 1, 1.25, 1.5, 2, 2.5, 3, 4]:
#     system("python3 run_experiment.py 7 10 22 20 community_sat BT 0 0.8 {}".format(T))
#     system("python3 run_experiment.py 7 10 22 20 community_sat CDCL 0 0.8 {}".format(T))



# # system("python3 run_experiment.py 3 10 30 30 community_sat BT 0 0.8 0")
# # system("python3 run_experiment.py 3 10 30 30 random_sat CDCL")
# # system("python3 run_experiment.py 4 10 30 30 community_sat BT 0 0.8 0")
# # system("python3 run_experiment.py 4 10 30 30 random_sat CDCL")
# # system("python3 run_experiment.py 5 10 30 30 community_sat BT 0 0.8 0")
# # system("python3 run_experiment.py 5 10 30 30 random_sat CDCL")
# # system("python3 run_experiment.py 6 10 30 30 community_sat BT 0 0.8 0")
# # system("python3 run_experiment.py 6 10 30 30 random_sat CDCL")
# # system("python3 run_experiment.py 7 10 30 30 community_sat BT 0 0.8 0")
# # system("python3 run_experiment.py 7 10 30 30 random_sat CDCL")
# # system("python3 run_experiment.py 8 10 30 30 community_sat BT 0 0.8 0")
# # system("python3 run_experiment.py 8 10 30 30 random_sat CDCL 0 0.8 0")
# # system("python3 run_experiment.py 9 18 30 30 community_sat BT 0 0.8 0")
# # system("python3 run_experiment.py 9 10 30 30 random_sat CDCL")
# # system("python3 run_experiment.py 10 10 30 30 community_sat BT 0 0.8 0")
# # system("python3 run_experiment.py 10 10 30 30 random_sat CDCL")


# Plot random instances (detect, search, grover, m22)
#     python3 plot_queries_over_n.py data/BT-random_sat-3-10-30-30.csv data/BT-random_sat-4-10-29-30.csv data/BT-random_sat-5-10-28-20.csv data/BT-random_sat-6-10-22-20.csv data/BT-random_sat-7-10-18-15.csv data/BT-random_sat-8-10-16-15.csv data/BT-random_sat-9-15-19-15.csv data/CDCL-random_sat-3-10-27-30.csv data/CDCL-random_sat-4-10-30-30.csv data/CDCL-random_sat-5-10-30-20.csv data/CDCL-random_sat-6-10-30-20.csv data/CDCL-random_sat-7-10-30-15.csv data/CDCL-random_sat-8-10-30-15.csv data/CDCL-random_sat-9-10-25-30.csv

# Plot community structure instances (detect, search, grover, m22)
# python3 plot_queries_over_n.py data/BT-community_sat_0.8_0-7-10-22-20.csv data/BT-community_sat_0.8_0.25-7-10-22-20.csv data/BT-community_sat_0.8_0.5-7-10-22-20.csv data/BT-community_sat_0.8_0.75-7-10-22-20.csv data/BT-community_sat_0.8_1-7-10-22-20.csv data/BT-community_sat_0.8_1.25-7-10-22-20.csv data/BT-community_sat_0.8_1.5-7-10-22-20.csv data/BT-community_sat_0.8_2-7-10-22-20.csv data/BT-community_sat_0.8_2.5-7-10-22-20.csv data/BT-community_sat_0.8_3-7-10-22-20.csv data/BT-community_sat_0.8_4-7-10-22-20.csv data/CDCL-community_sat_0.8_0-7-10-22-20.csv data/CDCL-community_sat_0.8_0.25-7-10-22-20.csv data/CDCL-community_sat_0.8_0.5-7-10-22-20.csv data/CDCL-community_sat_0.8_0.75-7-10-22-20.csv data/CDCL-community_sat_0.8_1-7-10-22-20.csv data/CDCL-community_sat_0.8_1.25-7-10-22-20.csv data/CDCL-community_sat_0.8_1.5-7-10-22-20.csv data/CDCL-community_sat_0.8_2-7-10-22-20.csv data/CDCL-community_sat_0.8_2.5-7-10-22-20.csv data/CDCL-community_sat_0.8_3-7-10-22-20.csv data/CDCL-community_sat_0.8_4-7-10-22-20.csv
