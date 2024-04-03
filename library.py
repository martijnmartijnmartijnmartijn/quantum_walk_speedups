from math import log, ceil, floor, pi, sin, sqrt
import matplotlib.pyplot as plt
from sympy.abc import x, i
from sympy import solveset, binomial, summation, S, Symbol, N, sqrt as sympysqrt
from subprocess import check_output

x = Symbol('x', positive=True)
y = Symbol('y', positive=True)

def get_C(ep, n):
    """
    For the given ep and n, solves the following equation for C, assuming C>1:
        ep = 1/(1+C)^n \sum_{i=0}^{floor(n/2)} \binom(n,i) * C^i
    """
    exp = (1/(1+x))**n * summation(binomial(n, i) * (x**i) , (i, 0, floor(n/2))) - ep
    if not x in exp.free_symbols: # For n=1, the expression is already solved.
        return [exp]
    for C in  solveset(exp, x, domain=S.Reals):
        if C > 1:
            return C
    return None

def get_optimal_C_and_n(ep, R=0, W=0):
    """
    Given a maximum error probaiblity epsilon, computes the optimal setting of
    the constants C and n for Belovs' detection algorithm, assuming we have set
    a=\sqrt{b} and b=1/(2+2C)^2.

    Note that this optimisation is independent of R and W. However, if R and W
    are set to positive integers (and the debug print statement are
    turned on) then the number of queries actually made by the optimal algorithm
    given an isntance with R and W are printed.
    """
    sqrt_RW_factors = []
    queries = []
    Cs = []
    i = 0
    while i < 2 or sqrt_RW_factors[-1] < sqrt_RW_factors[-2]:
        n = 2 * i + 1
        Cs.append(get_C(ep, n))
        sqrt_RW_factors.append(n * sqrt(Cs[-1]*(1+Cs[-1]**2)))
        i += 1
    return round(Cs[-2], 8), n - 2, round(sqrt_RW_factors[-2], 8)


# # # # # # # # # # # # # # # # # # # # # # #
# Old stuff used for Piddock's algorithm.
# # # # # # # # # # # # # # # # # # # # # # #
# def create_graph(strings, eta=1, x=1):
#     """
#     Take a list of bitstring given as output by the SAT solver, and turn it
#     into the corresponding graph. We suppose a full binary tree numbered from
#     left-to-right, top-to-bottom. So:

#                                 1           layer 0
#                               /   \
#                             2       3       layer 1
#                           /   \    /  \
#                          4     5  6    7    layer 2
#                         / \   / \  \
#                        8  9  10 11  13 ...  layer 3
#                       / \
#                      16 17  ....            layer 4

#     The number of a string a_1a_2a_3...a_n is then
#         2^n + a_1 2^{n-1} + a_2 2^{n-2} + ... + a_i 2^{n-i} + ... + a_n 2^{n-n}
#     For example, the string 011, i.e. go left, right, right, gives us
#         2^3 + 0 2^2 + 2^1 + 2^0 = 8 + 2 + 1 = 11
#     and the string 101 gives us
#         2^3 + 2^2 + 0 2^1 + 2^0 = 8 + 4 + 1 = 13
#     """
#     get_number = lambda s : 2**len(s) + sum([2**(len(s) - i - 1) for i in range(len(s)) if s[i] == '1'])
#     edges = []
#     for j,string in enumerate(strings):
#         prev_node = 1 # Root vertex.
#         node = 0
#         for i in range(1, len(string) + 1): # Iterate over each prefix of the string.
#             node = get_number(string[:i]) # Find the node for this prefix.
#             edges.append((prev_node, node, {"weight" : 1})) # Add an edge to the previous vertex.
#             prev_node = node
#         edges.append((node, FINAL_VERTEX, {"weight" : FINAL_EDGE_WEIGHT_INVERTED}))
#     G = nx.Graph(edges)
#     return G

# def add_extra_edges(G, eta, x):
#     """
#     After creating the graph, we add an additional edge after each marked vertex
#     with weight x, and we connect the end-points of these marked edges to a
#     unique final vertex. We give these final edges weight 1000000000, so that
#     the extra energy at these edges does not contribute significant energy.

#     Root node = 1
#     Starting node = 0 with edge (0,1).
#     Final node = 1 with edges (k, -1)
#     """
#     G.add_edge(STARTING_VERTEX,1,weight=eta) # Add the new source vertex 0.
#     marked_vertices = list(G[FINAL_VERTEX].keys())
#     for i, marked_vertex in enumerate(marked_vertices): # For each marked vertex.
#         G.remove_edge(marked_vertex, FINAL_VERTEX) # Remove the edge from marked vertex to sink vertex
#         G.add_edge(marked_vertex, -2 - i, weight=x) # Add a new edge with weight x in between marked vertex and sink.
#         G.add_edge(-2 - i, FINAL_VERTEX, weight=FINAL_EDGE_WEIGHT_INVERTED)
#     return G

# def graph_set_x(G, x):
#     """
#     Takes a graph that already has extra edges in between the marked edges and
#     the final edge -1. Changes the weight of these in-between edges to x.
#     """
#     # marked_vertices = list(G[FINAL_VERTEX].keys())
#     for marked_vertex in G[FINAL_VERTEX]:
#         for vertex in G[marked_vertex]:
#             if vertex != FINAL_VERTEX:
#                 G[vertex][marked_vertex]["weight"] = x
#     return G
#     # eff_res = nx.resistance_distance(G, 0, final_node, weight="weight")

# def get_effective_resistance(G, s=STARTING_VERTEX, t=FINAL_VERTEX):
#     """
#     STARTING_VERTEX is the newly adding starting vertex, which connects to the
#     root of the tree with weight 1/eta. You can also use the root, which
#     should have a value of 1.
#     """
#     return nx.resistance_distance(G, s, t, weight="weight")

# def test_graph_stuff():
#     G = create_graph(["0001110010","000111001110000","000111001110001","000111001110100","000111001110101"])
#     G = create_graph(["000","111","011"])
#     print(get_effective_resistance(G, s=1))
#     nx.draw(G, with_labels=True, font_weight='bold')
#     plt.show()
#     add_extra_edges(G, eta=4, x=2)
#     print(get_effective_resistance(G))
#     nx.draw(G, with_labels=True, font_weight='bold')
#     plt.show()
#     graph_set_x(G, x=0.01)
#     nx.draw(G, with_labels=True, font_weight='bold')
#     plt.show()

# def round_estimate(s, a):
#     """
#     If a is the amplitude that amplitude estimation returns, round it to the
#     output it would have given if s bits would have been used.

#     Specifically, iterate from k=0 to 2^s-1. Compute sin^2(pi k/2^s) and compare
#     to a. Return k and sin^2(pi k/2^s) corresponding to the closest estimate.
#     """
#     min = 1
#     min_estimate = 1
#     min_k = -1
#     # print("Starting rounding of {}".format(a))
#     for k in range(2**s):
#         estimate = sin(pi * k / (2**s))**2
#         d = abs(a-estimate)
#         if d < min:
#             min = d
#             min_estimate = estimate
#             min_k = k
#     return min_estimate

# def simulate_algorithm(solutions, W, s=7, Piddock=True):
#     """
#     Given the set of bitstrings describing the paths to marked vertices,
#     computes, eta ...
#     """
#     # Step 1b. We compute r_1 and eta.
#     G = create_graph(solutions)
#     R_sigma_M = get_effective_resistance(G, s=1)
#     if not Piddock:
#         return [(0, 0, 0), 0, 0, 0, R_sigma_M, [0]]
#     r_1 = ceil(log(W * (R_sigma_M * epsilon_a(s) + R_sigma_M) / (1 - 2 * epsilon_a(s)), 2))
#     eta = 2**(r_1)/W

#     # Step 1c. We compute b, r_2, r_1(x).
#     x_max = 2**ceil(log(len(solutions),2))
#     x = x_max * eta
#     r_2 = 1
#     add_extra_edges(G, eta=eta, x=eta)
#     R_a = get_effective_resistance(G)
#     tilde_a = round_estimate(s, eta/R_a)
#     graph_set_x(G, x)
#     R_b = get_effective_resistance(G)
#     tilde_b = round_estimate(s, eta/R_b)

#     # Save all effective resistances. We can then round with a different s.
#     R_xs = [R_a] # TODO think about saving this differently.

#     # Already hit the stopping condition, go downwards.
#     downwards = 2 * tilde_b < tilde_a
#     while True:
#         if downwards:
#             x /= 2
#         else:
#             x *= 2
#         graph_set_x(G, x)
#         R_b = get_effective_resistance(G)
#         tilde_b = round_estimate(s, eta/R_b)
#         if downwards and 2 * tilde_b > tilde_a:
#             x *= 2
#             break
#         elif not downwards and 2 * tilde_b < tilde_a:
#             break
#     r_2 = log(x/eta,2)

#     # Step 1d.
#     a = eta
#     b = x

#     # Step 2.
#     d = R_b - 2 * R_a
#     p = (R_a + d) / (log(b/a, 2) * (2 * R_a  + d))
#     return [(s, r_2, p), eta, x, r_1, R_sigma_M, R_xs]
