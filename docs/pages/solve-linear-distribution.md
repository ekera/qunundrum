# The <code>solve_linear_distribution</code> distribution

## Synopsis
```console
Synopsis: mpirun solve_linear_distribution \
   [ -adaptive | -non-adaptive | -non-adaptive-early-abort ] \
      [ -closest | -enumerate ] [ -timeout <timeout> ] \
         [ -detect-smooth-order ] \
            [ -lll | -lll-then-bkz | -bkz | -hkz ] \
               <distribution> <n> { <distribution> <n> }
```

Simulates the quantum algorithm by sampling the distribution, and solves the simulated outputs of $n$ runs for the logarithm $d$ or order $r$.

In total $10^3$ problem instances are consider to gather statistics.

The results are written to the console and to <code>logs/solve-linear.txt</code>.

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

   Once this timeout is elapsed the enumeration is aborted. This is reported as a failure to solve.

Flags specifying the lattice reduction algorithm (defaults to <code>-lll-then-bkz</code>):
- <code>-lll</code> use Lenstra-Lenstra-Lovász (LLL)
- <code>-bkz</code> use block Korkin-Zolotarev (BKZ)
- <code>-hkz</code> use Hermite Korkin-Zolotarev (HKZ)
- <code>-lll-then-bkz</code> use LLL and then BKZ if solving the LLL-reduced basis fails

## Interpreting the output
The log file <code>logs/solve-linear.txt</code> is on the format
```
# Processing: linear-distribution-det-dim-2048-d-m-2048-s-30.txt
# Search strategy: Adaptive
# Solution method: Closest
# Reduction algorithm: LLL then BKZ
# Detect smooth order: False
# Timestamp: 2024-02-29 01:36:22 CET
m: 2048 s: 30 n: 35 -- success: 994 -- fail: 6 (4) -- prepare:    32.841 ms solve:  1542.360 ms [ 1360.476, 34093.709] ** decrementing
m: 2048 s: 30 n: 34 -- success: 986 -- fail: 11 (2) -- prepare:    19.435 ms solve:  1601.184 ms [ 1083.961, 31584.762] ** stopping

# Processing: linear-distribution-det-dim-2048-r-m-2048-s-30.txt
# Search strategy: Adaptive
# Solution method: Closest
# Reduction algorithm: LLL then BKZ
# Detect smooth order: False
# Timestamp: 2024-02-29 01:37:35 CET
m: 2048 s: 30 n: 34 -- success: 992 -- fail: 8 (3) -- prepare:    19.220 ms solve:  1439.850 ms [ 1224.399, 31690.782] ** decrementing
m: 2048 s: 30 n: 33 -- success: 976 -- fail: 11 (1) -- prepare:    19.029 ms solve:  1438.394 ms [  961.454, 28802.589] ** stopping

# Processing: linear-distribution-det-dim-2048-d-m-2048-s-1.txt
# Search strategy: Adaptive
# Solution method: Enumerate (timeout: 300 s)
# Reduction algorithm: LLL then BKZ
# Detect smooth order: False
# Timestamp: 2024-02-29 01:41:04 CET
m: 2048 s: 1 n: 1 -- success: 1000 -- fail: 0 (0) -- prepare:     2.167 ms solve:     5.964 ms [    0.329,    26.722] ** stopping

# Processing: linear-distribution-det-dim-2048-r-m-2048-s-1.txt
# Search strategy: Adaptive
# Solution method: Enumerate (timeout: 300 s)
# Reduction algorithm: LLL then BKZ
# Detect smooth order: False
# Timestamp: 2024-02-29 01:41:05 CET
m: 2048 s: 1 n: 1 -- success: 1000 -- fail: 0 (0) -- prepare:     1.357 ms solve:     0.795 ms [    0.323,     5.122] ** stopping
```
where we find $m$, $s$ or $\ell$, $n$ — #success — #fail — prep-time — solve-time, and where
- $m$ is the bit length of the logarithm $d$ or order $r$,
- $s$ is the tradeoff factor such that $\ell = \lceil m / s \rceil$, if $s$ was specified when the distribution was generated, otherwise $\ell$ is explicitly stated instead,
- $n$ is the number of runs,
- #success is the number of problem instances that were successfully solved,
- #fail is the number of problem instances not solved, where the count within parenthesis is the number of problem instances that failed due to sampling errors,
- prep-time is the average time in ms required to setup the problem instances,
- solve-time is the average [min, max] time in ms required to solve the problem instances.
