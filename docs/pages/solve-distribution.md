# The <code>solve_distribution</code> executable

## Synopsis
```console
Synopsis: mpirun solve_distribution \
   [ -adaptive | -non-adaptive | -non-adaptive-early-abort ] \
      [ -closest | -enumerate ] [ -timeout <timeout> ] \
         [ -detect-smooth-order ] \
            [ -lll | -lll-then-bkz | -bkz | -hkz ] \
               <distribution> <n> { <distribution> <n> }
```

Simulates the quantum algorithm by sampling the distribution, and solves the simulated outputs of $n$ runs for the logarithm $d$ and order $r$.

In total $10^3$ problem instances are consider to gather statistics.

The results are written to the console and to <code>logs/solve.txt</code>.

> <b>Note:</b> This is an MPI program. The node with rank zero acts as server. All other nodes are clients, requesting jobs from and reporting back to the server node. A minimum of two nodes is hence required.

### Mandatory command line arguments
Tuples <code>\<distribution\> \<n\></code> where
- <code>\<distribution\></code> is the path to the distribution
- <code>\<n\></code> is the number of runs

### Optional command line arguments
Flags specifying the search strategy (defaults to <code>-adaptive</code>):
- <code>-adaptive</code> increment or decrement $n$ to find the minimum
- <code>-non-adaptive</code> attempt to solve only for the $n$ specified
- <code>-non-adaptive-early-abort</code> abort immediately if there are too many failures

Flags specifying whether an enumeration is performed (defaults to <code>-closest</code>):
- <code>-closest</code> solves by finding the closest vector to a given lattice vector
- <code>-enumerate</code> solves by enumerating the lattice

Flag specifying whether smooth orders are detected (defaults to no detection):
- <code>-detect-smooth-order</code> detects and leverages smooth orders when solving

Flag specifying the enumeration timeout (defaults to 300 s):
- <code>-timeout \<timeout\></code> sets the enumeration timeout to <code>\<timeout\></code> seconds

   Once this timeout is elapsed the enumeration is aborted. This is reported as a failure to solve. There are separate timeouts when enumerating for $d$ and $r$ respectively.

Flags specifying the lattice reduction algorithm (defaults to <code>-lll-then-bkz</code>):
- <code>-lll</code> use Lenstra-Lenstra-Lovász (LLL)
- <code>-bkz</code> use block Korkin-Zolotarev (BKZ)
- <code>-hkz</code> use Hermite Korkin-Zolotarev (HKZ)
- <code>-lll-then-bkz</code> use LLL and then BKZ if solving the LLL-reduced basis fails

## Interpreting the output
The log file <code>logs/solve.txt</code> is on the format
```
# Processing: distribution-det-dim-heuristic-sigma-optimal-m-2048-s-30.txt
# Search strategy: Adaptive
# Solution method: Closest
# Reduction algorithm: LLL then BKZ
# Detect smooth order: False
# Timestamp: 2024-02-29 01:36:00 CET
m: 2048 s: 30 n: 35 -- success: 992 -- fail: 8 (6) -- prepare:    37.095 ms solve:  1268.404 ms [ 1174.381, 33361.921] ** decrementing
m: 2048 s: 30 n: 34 -- success: 978 -- fail: 11 (2) -- prepare:    33.243 ms solve:  1356.987 ms [ 1064.196, 30762.823] ** stopping
```
where we find $m$, $s$ or $\ell$, $n$ — #success — #fail — prep-time — solve-time, and where
- $m$ is the bit length of the order $r$,
- $s$ is the tradeoff factor such that $\ell = \lceil m / s \rceil$, if $s$ was specified when the distribution was generated, otherwise $\ell$ is explicitly stated instead,
- $n$ is the number of runs,
- #success is the number of problem instances that were successfully solved,
- #fail is the number of problem instances not solved, where the count within parenthesis is the number of problem instances that failed due to sampling errors,
- prep-time is the average time in ms required to setup the problem instances, and
- solve-time is the average [min, max] time in ms required to solve the problem instances.
