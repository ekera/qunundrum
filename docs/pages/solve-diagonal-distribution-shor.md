# The <code>solve_diagonal_distribution_shor</code> executable

## Synopsis
```console
Synopsis: mpirun solve_diagonal_distribution_shor \
   [ -search-bound <bound> ] \
      <distribution> { <distribution> }
```

Simulates the quantum algorithm by sampling the distribution, and solves the simulated outputs for d given r using Shor's original post-processing algorithm improved with our searching over t.

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
The bound on t is related to Delta, as is explained in detail in Section 3 in the paper on revisiting Shor's algorithm for computing general discrete logarithms [[E19]](https://arxiv.org/pdf/1905.09084.pdf).
Searching over t as proposed in the aforementioned paper increases the success probability.

## Interpreting the output
The log file <code>logs/solve-diagonal-shor.txt</code> is on the format
```
# Processing: diagonal-distribution-det-dim-2048-m-2048-l-0.txt
m: 2048 l: 0 n: 1 -- success: 799 -- fail: 201 (201) -- prepare:     1.688 ms solve:     0.159 ms [    0.068,    12.158] 

# Processing: diagonal-distribution-det-dim-2048-m-2048-l-1.txt
m: 2048 l: 1 n: 1 -- success: 900 -- fail: 100 (100) -- prepare:     1.244 ms solve:     0.139 ms [    0.068,     0.209] 

# Processing: diagonal-distribution-det-dim-2048-m-2048-l-2.txt
m: 2048 l: 2 n: 1 -- success: 943 -- fail: 57 (57) -- prepare:     1.333 ms solve:     0.149 ms [    0.070,     6.230] 

# Processing: diagonal-distribution-det-dim-2048-m-2048-l-5.txt
m: 2048 l: 5 n: 1 -- success: 990 -- fail: 10 (10) -- prepare:     1.607 ms solve:     0.183 ms [    0.067,    13.168] 
```
where we find m, l, n -- #success -- #fail -- prep-time -- solve-time, and
- m is the bit length of the order r
- l is the padding length
- n is the number of runs
- #success is the number of problem instances that were successfully solved
- #fail is the number of problem instances not solved, where the count within parenthesis is the number of problem instances that failed due to sampling errors
- prep-time is the average time in ms required to setup the problem instances
- solve-time is the average [min, max] time in ms required to solve the problem instances
