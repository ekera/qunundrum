# The <code>estimate_runs_distribution</code> executable

## Synopsis
```console
Synopsis: mpirun estimate_runs_distribution \
   [ -v-bound <v-bound> ] <distribution> { <distribution> }
```

Estimates the number of runs required to solve a distribution for both a logarithm $d$ and order $r$ with at least 99\% success probability when accepting to enumerate at most $B_v$ vectors in the lattice.

The results are written to the console and to <code>logs/estimate-runs.txt</code>.

> <b>Note:</b> This is an MPI program. The node with rank zero acts as server. All other nodes are clients, requesting jobs from and reporting back to the server node. A minimum of two nodes is hence required.

### Mandatory command line arguments
Arguments <code>\<distribution\></code> where
- <code>\<distribution\></code> is the path to a distribution

### Optional command line arguments
Flag specifying the bound $B_v$ in $v$ (defaults to 2):
- <code>-v-bound <v-bound\></code> sets the bound $B_v$ to <code>\<v-bound\></code>

## Interpreting the output
The log file <code>logs/estimate-runs.txt</code> is on the format
```
# Processing: distribution-det-dim-heuristic-sigma-optimal-m-2048-s-30.txt
# Bounds: (v = 2)
# Timestamp: 2024-02-29 01:27:36 CET
m: 2048 s: 30 n: 30 -- tau_d 6.725372 v_d: 1.05704E+74 <3135> -- tau_r: 5.949782 v_r: 6.11467E+66 <3135>
m: 2048 s: 30 n: 31 -- tau_d 6.777422 v_d: 2.44615E+56 <3356> -- tau_r: 6.002945 v_r: 8.47257E+48 <3356>
m: 2048 s: 30 n: 35 -- tau_d 6.890261 v_d: 1.34245E-15 <3968> -- tau_r: 6.138421 v_r: 9.55378E-24 <3968>
m: 2048 s: 30 n: 34 -- tau_d 6.834617 v_d: 424.872 <3729> -- tau_r: 6.079734 v_r: 4.72936E-06 <3729>
```
where we find $m$, $s$ or $\ell$, $n$ — $\tau_d$, $v_d$, \<#errors\> — $\tau_r$, $v_r$, \<#errors\> on each line, and where
- $m$ is bit length of the order $r$,
- $s$ is the tradeoff factor such that $\ell = \lceil m / s \rceil$, if $s$ was specified when the distribution was generated, otherwise $\ell$ is explicitly stated instead,
- $n$ is the number of runs,
- $\tau_d$ is a value related to the radius $R_d$ for $d$,
- $v_d$ is the number of vectors that we at most expect to have to enumerate for $d$,
- $\tau_r$ is a value related to the radius $R_r$ for $r$,
- $v_r$ is the number of vectors that we at most expect to have to enumerate for $r$, and
- #errors is the number of sets of $n$ simulated outputs for which $R_d$ or $R_r$ was set to infinity due to sampling errors.

Note that the executable first computes estimates for $n = s$ and $n = s + 1$. It then interpolates a guess for the correct $n$ by using that $v_d$ is reduced by approximately a constant factor for every increment of $n$ by one, and that $v_d$ is greater in general than $v_r$. Depending on whether this guess is above or below the bound, $n$ is then increased or decreased until the smallest $n$ that satisfies the bound on both $v_d$ and $v_r$ has been found.