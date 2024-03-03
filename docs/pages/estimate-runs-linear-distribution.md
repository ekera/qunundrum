# The <code>estimate_runs_linear_distribution</code> executable

## Synopsis
```console
Synopsis: mpirun estimate_runs_linear_distribution \
   [ -v-bound <v-bound> ] <distribution> { <distribution> }
```
Estimates the number of runs required to solve a linear distribution for a logarithm $d$ or order $r$ depending on the type of distribution input with at least 99\% success probability when accepting to enumerate at most $B_v$ vectors in the lattice.

The results are written to the console and to <code>logs/estimate-runs-linear.txt</code>.

> <b>Note:</b> This is an MPI program. The node with rank zero acts as server. All other nodes are clients, requesting jobs from and reporting back to the server node. A minimum of two nodes is hence required.

### Mandatory command line arguments
Arguments <code>\<distribution\></code> where
- <code>\<distribution\></code> is the path to a distribution

### Optional command line arguments
Flag specifying the bound $B_v$ in $v$ (defaults to 2):
- <code>-v-bound <v-bound\></code> sets the bound $B_v$ to <code>\<v-bound\></code>

## Interpreting the output
The log file <code>logs/estimate-runs-linear.txt</code> is on the format
```
# Processing: linear-distribution-det-dim-2048-d-m-2048-s-30.txt
# Bounds: (v = 2)
# Timestamp: 2024-02-29 01:22:16 CET
m: 2048 s: 30 n: 30 -- tau: 6.746676 v: 1.67072E+74 <2839>
m: 2048 s: 30 n: 31 -- tau: 6.789410 v: 3.19131E+56 <3041>
m: 2048 s: 30 n: 35 -- tau: 6.862463 v: 6.70879E-16 <3373>
m: 2048 s: 30 n: 34 -- tau: 6.821922 v: 312.246 <3244>

# Processing: linear-distribution-det-dim-2048-r-m-2048-s-30.txt
# Bounds: (v = 2)
# Timestamp: 2024-02-29 01:22:22 CET
m: 2048 s: 30 n: 30 -- tau: 5.735868 v: 6.16829E+64 <1394>
m: 2048 s: 30 n: 31 -- tau: 5.753717 v: 3.36692E+46 <1447>
m: 2048 s: 30 n: 34 -- tau: 5.839342 v: 1.38688E-08 <1675>
m: 2048 s: 30 n: 33 -- tau: 5.798844 v: 1.35113E+10 <1533>

# Processing: linear-distribution-det-dim-2048-d-m-2048-s-1.txt
# Bounds: (v = 2)
# Timestamp: 2024-02-29 01:22:40 CET
m: 2048 s: 1 n: 1 -- tau: 4.288488 v: 1202.5 <96>
m: 2048 s: 1 n: 2 -- tau: 4.776080 v: 7.54776E-612 <191>

# Processing: linear-distribution-det-dim-2048-r-m-2048-s-1.txt
# Bounds: (v = 2)
# Timestamp: 2024-02-29 01:22:41 CET
m: 2048 s: 1 n: 1 -- tau: 3.306757 v: 310.504 <46>
m: 2048 s: 1 n: 2 -- tau: 3.745119 v: 8.87243E-613 <109>
```
where we find $m$, $s$ or $\ell$, $n$ — $\tau$, $v$, \<#errors\> on each line, and where
- $m$ is the bit length of the logarithm $d$ or order $r$,
- $s$ is the tradeoff factor such that $\ell = \lceil m / s \rceil$, if $s$ was specified when the distribution was generated, otherwise $\ell$ is explicitly stated instead,
- $n$ is the number of runs,
- $\tau$ is a value related to the radius $R$,
- $v$ is the number of vectors that we at most expect to have to enumerate, and
- #errors is the number of sets of $n$ simulated outputs for which $R$ was set to infinity due to sampling errors.

Note that the executable first computes estimates for $n = s$ and $n = s + 1$. It then interpolates a guess for the correct $n$ by using that $v$ is reduced by approximately a constant factor for every increment of $n$ by one. Depending on whether this guess is above or below the bound, $n$ is then increased or decreased until the smallest $n$ that satisfies the bound on $v$ has been found.
