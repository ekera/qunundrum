# The <code>solve_linear_distribution_shor</code> executable

## Synopsis
```console
Synopsis: mpirun solve_linear_distribution_shor
   [ -search-bound-j <bound> ] [ -search-bound-cofactors <bound> ]
      <distribution> { <distribution> }
```

Simulates the quantum algorithm by sampling the distribution, and solves the simulated outputs for r using Shor's original post-processing algorithm based on continued fraction expansion.

The results are written to the console and to <code>logs/solve-linear-shor.txt</code>.

> <b>Note:</b> This is an MPI program. The node with rank zero acts as server. All other nodes are clients, requesting jobs from and reporting back to the server node. A minimum of two nodes is hence required.

### Mandatory command line arguments
Entries <code>\<distribution\></code> where
- <code>\<distribution\></code> is the path to the distribution

### Optional command line arguments
Flags specifying the search bound (defaults to 2^8):
- <code>-search-bound-j \<bound\></code> sets the search bound for j to <code>\<bound\></code>

Flags specifying the search bound (defaults to 2^16):
- <code>-search-bound-cofactors \<bound\></code> sets the search bound for cofactors to <code>\<bound\></code>

## Interpreting the output
The log file <code>logs/solve-linear-shor.txt</code> is on the format
```
# Processing: linear-distribution-max-dim-2048-r-m-2048-s-1.txt
m: 2048 s: 1 n: 1 -- success: 1000 -- fail: 0 (0) -- prepare:     1.789 ms solve:    36.435 ms [   14.173,   650.876] 
```
where we find m, s, n -- #success -- #fail -- prep-time -- solve-time, and
- m is the bit length of the order r
- s is the tradeoff factor
- n is the number of runs
- #success is the number of problem instances that were successfully solved
- #fail is the number of problem instances not solved, where the count within parenthesis is the number of problem instanced that failed due to sampling errors
- prep-time is the average time required to setup the problem instance
- solve-time is the average [min, max] time in ms required to solve the problem instance
