# The <code>solve_diagonal_distribution_shor</code> executable

## Synopsis
```console
Synopsis: mpirun solve_diagonal_distribution_shor \
   [ -t-bound <t-bound> ] [ -eta-bound <eta-bound> ] \
      <distribution> { <distribution> }
```

Simulates the quantum algorithm by sampling the distribution, and solves the simulated outputs for the logarithm $d$ given the order $r$ using Shor's original post-processing algorithm modified to search over $t$ and $\eta$ as explained in [[E19p]](https://doi.org/10.48550/arXiv.1905.09084).

In total $10^3$ problem instances are consider to gather statistics.

The results are written to the console and to <code>logs/solve-diagonal-shor.txt</code>.

> <b>Note:</b> This is an MPI program. The node with rank zero acts as server. All other nodes are clients, requesting jobs from and reporting back to the server node. A minimum of two nodes is hence required.

### Mandatory command line arguments
Arguments <code>\<distribution\></code> where
- <code>\<distribution\></code> is the path to the distribution

### Optional command line arguments
Flag specifying the search bound $B_t$ in $t$ (defaults to zero):
- <code>-t-bound \<t-bound\></code> sets the search bound $B_t$ to <code>\<t-bound\></code>

   All $t \in [-B_t, B_t] \cap \mathbb Z$ are searched when when solving $(j, k)$ for $d$ given $r$.

   The search bound $B_\Delta$ in $\Delta$ is selected as a function of $B_t$.

Flag specifying the search bound $B_\eta$ in $\eta$ (defaults to zero):
- <code>-eta-bound \<eta-bound\></code> sets the search bound $B_\eta$ to <code>\<eta-bound\></code>

   All $\eta \in [-B_\eta, B_\eta] \cap \mathbb Z$ are searched when sampling $j$ and $\eta$, and when solving $(j, k)$ for $d$ given $r$.

## Interpreting the output
The log file <code>logs/solve-diagonal-shor.txt</code> is on the format
```
# Processing: diagonal-distribution-det-dim-2048-m-2048-sigma-0-s-1.txt
# Bounds: (eta = 100 (100), t = 10000)
# Timestamp: 2024-02-28 21:00:24 CET
m: 2048 sigma: 0 s: 1 n: 1 -- success: 1000 -- fail: 0 (0) -- prepare:     1.360 ms solve:    36.747 ms [    0.049,  3412.087] 

# Processing: diagonal-distribution-det-dim-2048-m-2048-sigma-0-s-1.txt
# Bounds: (eta = 0 (100), t = 0)
# Timestamp: 2024-02-29 11:45:30 CET
m: 2048 sigma: 0 s: 1 n: 1 -- success: 612 -- fail: 388 (44) -- prepare:     1.362 ms solve:     0.071 ms [    0.047,     0.320] 

# Processing: diagonal-distribution-det-dim-2048-m-2048-sigma-0-s-1.txt
# Bounds: (eta = 0 (100), t = 10000)
# Timestamp: 2024-02-29 11:47:32 CET
m: 2048 sigma: 0 s: 1 n: 1 -- success: 793 -- fail: 207 (0) -- prepare:     1.822 ms solve:     8.095 ms [    0.046,    40.735] 

# Processing: diagonal-distribution-det-dim-2048-m-2048-sigma-5-s-1.txt
# Bounds: (eta = 0 (0), t = 10000)
# Timestamp: 2024-02-29 11:47:38 CET
m: 2048 sigma: 5 s: 1 n: 1 -- success: 993 -- fail: 7 (7) -- prepare:     1.782 ms solve:     0.072 ms [    0.048,     1.863] 
```
where we find $m$, $\varsigma$, $s$ or $\ell$, $n$ — #success — #fail — prep-time — solve-time, and where
- $m$ is the bit length of the order $r$,
- $\varsigma$ is the bit length of the padding,
- $s$ is the tradeoff factor such that $\ell = \lceil m / s \rceil$, if $s$ was specified when the distribution was generated, otherwise $\ell$ is explicitly stated instead,
- $n$ is the number of runs,
- #success is the number of problem instances that were successfully solved,
- #fail is the number of problem instances not solved, where the count within parenthesis is the number of problem instances that failed due to sampling errors,
- prep-time is the average time in ms required to setup the problem instances, and
- solve-time is the average [min, max] time in ms required to solve the problem instances.