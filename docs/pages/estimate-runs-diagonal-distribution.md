# The <code>estimate_runs_diagonal_distribution</code> executable

## Synopsis
```console
Synopsis: mpirun estimate_runs_diagonal_distribution \
   [ -v-bound <v-bound> ] [ -delta-bound <delta-bound> ] \
      [ -eta-bound <eta-bound> ] \
         <distribution> { <distribution> }
```

Estimates the number of runs required to solve a diagonal distribution for the logarithm $d$ given the order $r$ with at least 99\% success probability when accepting to enumerate at most $B_v$ vectors in the lattice.

This when accepting to search all combinations of peak indices $\eta_1, \ldots, \eta_n$ such that $\eta_i \in [-B_\eta, B_\eta] \cap \mathbb Z$.

The results are written to the console and to <code>logs/estimate-runs-diagonal.txt</code>.

> <b>Note:</b> This is an MPI program. The node with rank zero acts as server. All other nodes are clients, requesting jobs from and reporting back to the server node. A minimum of two nodes is hence required.

### Mandatory command line arguments
Arguments <code>\<distribution\></code> where
- <code>\<distribution\></code> is the path to a distribution

### Optional command line arguments
Flag specifying the bound $B_v$ in $v$ (defaults to $2$):
- <code>-v-bound \<v-bound\></code> sets the bound $B_v$ to <code>\<v-bound\></code>

Flag specifying the search bound $B_\Delta$ in $\Delta$ (defaults to <code>DELTA_BOUND</code>):
- <code>-delta-bound \<eta-bound\></code> sets $B_\Delta$ to <code>\<eta-bound\></code>

   All $\Delta \in [-B_\Delta, B_\Delta] \cap \mathbb Z$ are searched when sampling $k$ given $j$ and $\eta$.

Flag specifying the search bound $B_\eta$ in $\eta$ (defaults to same as for the distribution):
- <code>-eta-bound \<eta-bound\></code> sets $B_\eta$ to <code>\<eta-bound\></code>

   All $\eta_i \in [-B_\eta, B_\eta] \cap \mathbb Z$ are searched when sampling $j_i$ and $\eta_i$.

## Interpreting the output
The log file <code>logs/estimate-runs-diagonal.txt</code> is on the format
```
# Processing: diagonal-distribution-det-dim-2048-m-2048-sigma-12-s-30.txt
# Bounds: (eta = <all> (0), delta = 1000000, v = 2)
# Timestamp: 2024-02-29 08:47:24 CET
m: 2048 sigma: 12 s: 30 n: 30 -- tau: 6.049514 v: 5.21266E+67 <1461>
m: 2048 sigma: 12 s: 30 n: 31 -- tau: 6.028474 v: 1.49257E+49 <1467>
m: 2048 sigma: 12 s: 30 n: 34 -- tau: 6.141849 v: 2.13424E-05 <1583>
m: 2048 sigma: 12 s: 30 n: 33 -- tau: 6.103036 v: 1.75424E+13 <1557>
```
where we find $m$, $\varsigma$, $s$ or $\ell$, $n$ — $\tau$ $v$ \<#errors\> on each line, and where
- $m$ is the bit length of the order $r$,
- $\varsigma$ is the bit length of the padding,
- $s$ is the tradeoff factor such that $\ell = \lceil m / s \rceil$, if $s$ was specified when the distribution was generated, otherwise $\ell$ is explicitly stated instead,
- $n$ is the number of runs,
- $\tau$ is a value related to the radius $R$,
- $v$ is the number of vectors that we at most expect to have to enumerate, and
- #errors is the number of sets of $n$ simulated outputs for which $R$ was set to infinity due to sampling errors.

Note that the executable first computes estimates for $n = s$ and $n = s + 1$. It then interpolates a guess for the correct $n$ by using that $v$ is reduced by approximately a constant factor for every increment of $n$ by one. Depending on whether this guess is above or below the bound, $n$ is then increased or decreased until the smallest $n$ that satisfies the bound on $v$ has been found.
