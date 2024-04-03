# quantum_walk_speedups

This repository allows you to generate several types of random SAT instances,
solve them with a basic backtracking solver and a more modern sat solver,
and plot the resulting time complexity. Most importantly, the data
gathered from the backtracking solver can be used to determine the query
complexity of a quantum walk algorithm that solves SAT instances by walking
on the classical backtracking tree (an algorithm due to Montanaro 2018,
`https://theoryofcomputing.org/articles/v014a015`), and subsequently the
time complexity of these queries is determined by computing the depth of
a circuit to implement the queries.

- `run_experiment.py` can be used to generate SAT instances of a given type,
  which are then solved, saving the relevant data to files into a .csv file
  in the `data/` folder. See the documentation in the file itself for precise
  usage instructions.
- `plot_queries_over_n.py` can be called using the .csv files output by
  `run_experiment.py` to plot the measured time complexity of the classical
  solvers, and the computed time complexity of the quantum algorithm(s).
- `library.py` contains some functions used to determine the optimal
  configuration of the quantum walk algorithm (i.e. determine C and number
  of repetitions). These functions are not used elsewhere in the code.
- `data/` contains all the .csv files generated by `run_experiment.py`.
- `solvers/` contains the used backtracking solver in `bt.c` (modified from the
  solver used in `https://doi.org/10.5523/bris.19va21gun3c7629f291kmd6w37`) which
  can be compiled using the given makefile.
- `date_complexity/` contains two files from the code from
  `https://doi.org/10.22331/q-2019-07-18-167` which compute the circuit
  size/depth for circuits needed to implement Grover/quantum walk queries.
- `generators-IJCAI17` contains files from
  `https://www.ijcai.org/proceedings/2017/89` that are used to generate random
  SAT instances with community structure.
