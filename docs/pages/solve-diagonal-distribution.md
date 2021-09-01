# The <code>solve_diagonal_distribution</code> distribution

## Synopsis
```console
Synopsis: mpirun solve_diagonal_distribution \
   [ -adaptive | -non-adaptive | -non-adaptive-early-abort ] \
      [ -closest | -enumerate ] [ -timeout <timeout> ] \
         [ -lll | -lll-then-bkz | -bkz | -hkz ] \
            <distribution> <n> { <distribution> <n> }
```

Simulates the quantum algorithm by sampling the distribution, and solves the simulated outputs of n runs for d given r.

In total 10^3 problem instances are consider to gather statistics.

The results are written to the console and to <code>logs/solve-diagonal.txt</code>.

> <b>Note:</b> This is an MPI program. The node with rank zero acts as server. All other nodes are clients, requesting jobs from and reporting back to the server node. A minimum of two nodes is hence required.

### Mandatory command line arguments
Tuples <code>\<distribution\> <n></code> where
- <code>\<distribution\></code> is the path to the distribution
- <code>\<n\></code> is the number of runs

### Optional command line arguments
Flags specifying the search strategy (defaults to <code>-adaptive</code>):
- <code>-adaptive</code> increment or decrement n to find the minimum
- <code>-non-adaptive</code> attempt to solve only for the n specified
- <code>-non-adaptive-early-abort</code> abort immediately if there are too many failures

Flags specifying whether an enumeration is performed (defaults to <code>-closest</code>):
- <code>-closest</code> solves by finding the closest vector to a given lattice vector
- <code>-enumerate</code> solves by enumerating the lattice

Flags specifying the enumeration timeout (defaults to 300 s):
- <code>-timeout \<timeout\></code> sets the enumeration timeout to <code>\<timeout\></code> seconds

Once this timeout is elapsed the enumeration is aborted. This is reported as a failure to solve.

Flags specifying the lattice reduction algorithm (defaults to <code>-lll-then-bkz</code>):
- <code>-lll</code> use Lenstra-Lenstra-Lovász (LLL)
- <code>-bkz</code> use block Korkin-Zolotarev (BKZ)
- <code>-hkz</code> use Hermite Korkin-Zolotarev (HKZ)
- <code>-lll-then-bkz</code> use LLL and then BKZ if solving the LLL-reduced basis fails

## Interpreting the output
The log file <code>logs/solve-diagonal.txt</code> is on the format
```
# Processing: diagonal-distribution-det-dim-2048-m-2048-sigma-12-s-30.txt
m: 2048 sigma: 12 s: 30 n: 35 -- success: 1000 -- fail: 0 (0) -- prepare:   118.375 ms solve:  4348.937 ms [ 3635.028,  4540.994] ** decrementing
m: 2048 sigma: 12 s: 30 n: 34 -- success: 986 -- fail: 11 (9) -- prepare:   116.460 ms solve:  4136.399 ms [ 2199.815, 124892.005] ** stopping

# Processing: diagonal-distribution-det-dim-2048-m-2048-sigma-5-s-1.txt
m: 2048 sigma: 5 s: 1 n: 1 -- success: 992 -- fail: 8 (2) -- prepare:     4.694 ms solve:     1.079 ms [    0.497,    20.561] ** stopping
```
where we find m, sigma, s, n -- #success -- #fail -- prep-time -- solve-time, and
- m is the bit length of the order r
- sigma is the padding length
- s is the tradeoff factor
- n is the number of runs
- #success is the number of problem instances that were successfully solved
- #fail is the number of problem instances not solved, where the count within parenthesis is the number of problem instances that failed due to sampling errors
- prep-time is the average time in ms required to setup the problem instances
- solve-time is the average [min, max] time in ms required to solve the problem instances