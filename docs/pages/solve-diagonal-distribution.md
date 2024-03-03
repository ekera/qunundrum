# The <code>solve_diagonal_distribution</code> distribution

## Synopsis
```console
Synopsis: mpirun solve_diagonal_distribution \
   [ -delta-bound <delta-bound> ] [ -eta-bound <eta-bound> ] \
      [ -adaptive | -non-adaptive | -non-adaptive-early-abort ] \
         [ -closest | -enumerate ] [ -timeout <timeout> ] \
            [ -lll | -lll-then-bkz | -bkz | -hkz ] \
               <distribution> <n> { <distribution> <n> }
```

Simulates the quantum algorithm by sampling the distribution, and solves the simulated outputs of $n$ runs for the logarithm $d$ given the order $r$.

This when accepting to search all combinations of peak indices $\eta_1, \ldots, \eta_n$ such that $\eta_i \in [-B_\eta, B_\eta] \cap \mathbb Z$.

In total $10^3$ problem instances are consider to gather statistics.

The results are written to the console and to <code>logs/solve-diagonal.txt</code>.

> <b>Note:</b> This is an MPI program. The node with rank zero acts as server. All other nodes are clients, requesting jobs from and reporting back to the server node. A minimum of two nodes is hence required.

### Mandatory command line arguments
Tuples <code>\<distribution\> <n></code> where
- <code>\<distribution\></code> is the path to the distribution
- <code>\<n\></code> is the number of runs

### Optional command line arguments
Flag specifying the search bound $B_\Delta$ in $\Delta$ (defaults to <code>DELTA_BOUND</code>):
- <code>-delta-bound \<eta-bound\></code> sets $B_\Delta$ to <code>\<eta-bound\></code>

   All $\Delta \in [-B_\Delta, B_\Delta] \cap \mathbb Z$ are searched when sampling $k_i$ given $j_i$ and $\eta_i$.

Flag specifying the search bound $B_\eta$ in $\eta$ (defaults to same as for the distribution):
- <code>-eta-bound \<eta-bound\></code> sets $B_\eta$ to <code>\<eta-bound\></code>

   All $\eta_i \in [-B_\eta, B_\eta] \cap \mathbb Z$ are searched when sampling $j_i$ and $\eta_i$, and when solving the pairs $(j_i, k_i)$ for $d$ given $r$.

Flags specifying the search strategy (defaults to <code>-adaptive</code>):
- <code>-adaptive</code> increment or decrement $n$ to find the minimum
- <code>-non-adaptive</code> attempt to solve only for the $n$ specified
- <code>-non-adaptive-early-abort</code> abort immediately if there are too many failures

Flags specifying whether an enumeration is performed (defaults to <code>-closest</code>):
- <code>-closest</code> solves by finding the closest vector to a given lattice vector
- <code>-enumerate</code> solves by enumerating the lattice

Flag specifying the enumeration timeout (defaults to 300 s):
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
# Bounds: (eta = <all> (0), delta = 1000000)
# Search strategy: Adaptive
# Solution method: Closest
# Reduction algorithm: LLL then BKZ
# Timestamp: 2024-02-29 11:18:15 CET
m: 2048 sigma: 12 s: 30 n: 34 -- success: 987 -- fail: 11 (5) -- prepare:    69.406 ms solve:  2412.413 ms [ 1817.261, 78989.007] ** incrementing
m: 2048 sigma: 12 s: 30 n: 35 -- success: 997 -- fail: 3 (2) -- prepare:    42.009 ms solve:  2244.208 ms [ 2022.457, 85395.739] ** stopping
```
where we find $m$, $\varsigma$, $s$ or $\ell$, $n$ — #success — #fail — prep-time — solve-time, and
where
- $m$ is the bit length of the order $r$,
- $\varsigma$ is the bit length of the padding,
- $s$ is the tradeoff factor such that $\ell = \lceil m / s \rceil$, if $s$ was specified when the distribution was generated, otherwise $\ell$ is explicitly stated instead,
- $n$ is the number of runs,
- #success is the number of problem instances that were successfully solved,
- #fail is the number of problem instances not solved, where the count within parenthesis is the number of problem instances that failed due to sampling errors,
- prep-time is the average time in ms required to setup the problem instances, and
- solve-time is the average [min, max] time in ms required to solve the problem instances.
