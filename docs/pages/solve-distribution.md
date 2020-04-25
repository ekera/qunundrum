# The <code>solve_distribution</code> executable

## Synopsis
```console
Synopsis: mpirun solve_distribution \
   [ -adaptive | -non-adaptive | -non-adaptive-early-abort ] \
      [ -closest | -enumerate ] [ -timeout <timeout> ] \
         [ -lll | -lll-then-bkz | -bkz | -hkz ] \
            <distribution> <n> { <distribution> <n> }
```

Simulates the quantum algorithm by sampling the distribution, and solves the simulated outputs of n runs for d and r.

The results are written to the console and to <code>logs/solve.txt</code>.

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
- <code>-closest</code> solves by finding the closest vectors to two given lattice vectors
- <code>-enumerate</code> solves by enumerating the lattice

Flags specifying the enumeration timeout (defaults to 300s):
- <code>-timeout \<timeout\></code> sets the enumeration timeout to <code>\<timeout\></code> seconds

Once this timeout is elapsed the enumeration is aborted. This is reported as a failure to solve.
There are separate timeouts when enumerating for d and r respectively.

Flags specifying the lattice reduction algorithm (defaults to <code>-lll-then-bkz</code>):
- <code>-lll</code> use Lenstra-Lenstra-Lovász (LLL)
- <code>-bkz</code> use block Korkin-Zolotarev (BKZ)
- <code>-hkz</code> use Hermite Korkin-Zolotarev (HKZ)
- <code>-lll-then-bkz</code> use LLL and then BKZ if solving the LLL-reduced basis fails

## Interpreting the output
The log file <code>logs/solve.txt</code> is on the format
```
# Processing: distribution-det-dim-heuristic-sigma-optimal-m-2048-s-30.txt
m: 2048 s: 30 n: 35 -- success: 993 -- fail: 7 (5) -- prepare:    98.479 ms solve:  3052.485 ms [ 1695.249, 61805.845] ** decrementing
m: 2048 s: 30 n: 34 -- success: 984 -- fail: 11 (1) -- prepare:    94.748 ms solve:  3114.389 ms [ 1604.829, 56063.432] ** stopping
```
where we find m, s, n -- #success -- #fail -- prep-time -- solve-time, and
- m is the bit length of the order r
- s is the tradeoff factor
- n is the number of runs
- #success is the number of problem instances that were successfully solved
- #fail is the number of problem instances not solved, where the count within parenthesis is the number of problem instanced that failed due to sampling errors
- prep-time is the average time required to setup the problem instance
- solve-time is the average [min, max] time in ms required to solve the problem instance
