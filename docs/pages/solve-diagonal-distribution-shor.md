# The <code>solve_diagonal_distribution_shor</code> executable

## Synopsis
```console
Synopsis: mpirun solve_diagonal_distribution_shor \
   [ -search-bound <bound> ] \
      <distribution> { <distribution> }
```

Simulates the quantum algorithm by sampling the distribution, and solves the simulated outputs for d given r using Shor's original post-processing algorithm modified to search over t as explained in [[E19p]](https://arxiv.org/pdf/1905.09084.pdf).

In total 10^3 problem instances are consider to gather statistics.

The results are written to the console and to <code>logs/solve-diagonal-shor.txt</code>.

> <b>Note:</b> This is an MPI program. The node with rank zero acts as server. All other nodes are clients, requesting jobs from and reporting back to the server node. A minimum of two nodes is hence required.

### Mandatory command line arguments
Entries <code>\<distribution\></code> where
- <code>\<distribution\></code> is the path to the distribution

### Optional command line arguments
Flags specifying the search bound (defaults to zero):
- <code>-search-bound \<bound\></code> sets the search bound on |t| to <code>\<bound\></code>

Note that alpha_d = round(alpha_r d/r) + Delta for some small Delta.
The search bound on t is related to Delta, as is explained in detail in Section 5.1 of [[E19p]](https://arxiv.org/pdf/1905.09084.pdf).
Searching over t as proposed in [[E19p]](https://arxiv.org/pdf/1905.09084.pdf) increases the success probability.

## Interpreting the output
The log file <code>logs/solve-diagonal-shor.txt</code> is on the format
```
# Processing: diagonal-distribution-det-dim-2048-m-2048-sigma-0-s-1.txt
m: 2048 sigma: 0 s: 1 n: 1 -- success: 640 -- fail: 360 (198) -- prepare:     3.350 ms solve:     0.144 ms [    0.067,     0.732]

# Processing: diagonal-distribution-det-dim-2048-m-2048-sigma-1-s-1.txt
m: 2048 sigma: 1 s: 1 n: 1 -- success: 833 -- fail: 167 (82) -- prepare:     3.310 ms solve:     0.170 ms [    0.070,     1.568]

# Processing: diagonal-distribution-det-dim-2048-m-2048-sigma-2-s-1.txt
m: 2048 sigma: 2 s: 1 n: 1 -- success: 896 -- fail: 104 (43) -- prepare:     4.097 ms solve:     0.171 ms [    0.073,     0.278]

# Processing: diagonal-distribution-det-dim-2048-m-2048-sigma-5-s-1.txt
m: 2048 sigma: 5 s: 1 n: 1 -- success: 991 -- fail: 9 (6) -- prepare:     3.713 ms solve:     0.178 ms [    0.107,     5.578]
```
where we find m, l, n -- #success -- #fail -- prep-time -- solve-time, and
- m is the bit length of the order r
- sigma is the padding length
- s is the tradeoff factor
- n is the number of runs
- #success is the number of problem instances that were successfully solved
- #fail is the number of problem instances not solved, where the count within parenthesis is the number of problem instances that failed due to sampling errors
- prep-time is the average time in ms required to setup the problem instances
- solve-time is the average [min, max] time in ms required to solve the problem instances
