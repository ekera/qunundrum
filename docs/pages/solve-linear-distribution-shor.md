# The <code>solve_linear_distribution_shor</code> executable

## Synopsis
```console
Synopsis: mpirun solve_linear_distribution_shor
   [ -search-bound-j <bound> ] [ -search-bound-cofactors <bound> ]
      <distribution> { <distribution> }
```

Simulates the quantum algorithm by sampling the distribution, and solves the simulated outputs for the order $r$ using Shor's original post-processing algorithm based on continued fraction expansion.

In total $10^3$ problem instances are consider to gather statistics.

The results are written to the console and to <code>logs/solve-linear-shor.txt</code>.

> <b>Note:</b> This is an MPI program. The node with rank zero acts as server. All other nodes are clients, requesting jobs from and reporting back to the server node. A minimum of two nodes is hence required.

### Mandatory command line arguments
Arguments <code>\<distribution\></code> where
- <code>\<distribution\></code> is the path to the distribution

### Optional command line arguments
Flag specifying the search bound (defaults to 2^8):
- <code>-search-bound-j \<bound\></code> sets the search bound for $j$ to <code>\<bound\></code>

Flag specifying the search bound (defaults to 2^16):
- <code>-search-bound-cofactors \<bound\></code> sets the search bound for cofactors to <code>\<bound\></code>

## Interpreting the output
The log file <code>logs/solve-linear-shor.txt</code> is on the format
```
# Processing: linear-distribution-det-dim-2048-r-m-2048-s-1.txt
# Bounds: (t: 256, cofactor: 65536)
# Timestamp: 2024-02-29 01:46:32 CET
m: 2048 s: 1 n: 1 -- success: 999 -- fail: 1 (0) -- prepare:    20.329 ms solve:    20.699 ms [    9.937,  5111.765] 

# Processing: linear-distribution-det-dim-2048-r-m-2048-s-1.txt
# Bounds: (t: 0, cofactor: 1)
# Timestamp: 2024-02-29 01:48:36 CET
m: 2048 s: 1 n: 1 -- success: 220 -- fail: 780 (0) -- prepare:     0.949 ms solve:    11.039 ms [    9.886,    12.118] 
```
where we find $m$, $s$ or $\ell$, $n$ — #success — #fail — prep-time — solve-time, and where
- $m$ is the bit length of the order $r$,
- $s$ is the tradeoff factor such that $\ell = \lceil m / s \rceil$, if $s$ was specified when the distribution was generated, otherwise $\ell$ is explicitly stated instead,
- $n$ is the number of runs,
- #success is the number of problem instances that were successfully solved,
- #fail is the number of problem instances not solved, where the count within parenthesis is the number of problem instances that failed due to sampling errors,
- prep-time is the average time in ms required to setup the problem instances,
- solve-time is the average [min, max] time in ms required to solve the problem instances.
