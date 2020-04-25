# The <code>estimate_runs_distribution</code> executable

## Synopsis
```console
Synopsis: mpirun estimate_runs_distribution \
   [ -bound <bound> ] <distribution> { <distribution> }
```

Estimates the number of runs required to solve a distribution for both a logarithm d and order r with at least 99% success probability when accepting to enumerate at most B vectors in the lattice.

The results are written to the console and to <code>logs/estimate-runs.txt</code>.

> <b>Note:</b> This is an MPI program. The node with rank zero acts as server. All other nodes are clients, requesting jobs from and reporting back to the server node. A minimum of two nodes is hence required.

### Mandatory command line arguments
Entries <code>\<distribution\></code> where
- <code>\<distribution\></code> is the path to a distribution

### Optional command line arguments
Flags controlling the bound B (defaults to 2):
- <code>-bound <bound\></code> sets the bound B to <code>\<bound\></code>

## Interpreting the output
The log file <code>logs/estimate-runs.txt</code> is on the format
```
# Processing: distribution-det-dim-heuristic-sigma-optimal-m-2048-s-30.txt
2048 30 30 -- 6.760466 2.2469E+74 <3180> -- 5.946987 5.75828E+66 <3180>
2048 30 31 -- 6.786973 3.02337E+56 <3426> -- 5.995831 7.23577E+48 <3426>
2048 30 35 -- 6.879788 1.03371E-15 <3863> -- 6.120637 6.12975E-24 <3863>
2048 30 34 -- 6.818136 284.849 <3662> -- 6.051322 2.37386E-06 <3662>
```
where we find m, s, n -- tau_d v_d \<#errors\> -- tau_r v_r \<#errors\> on each line, and 
- m is the bit length of the order r
- s is the tradeoff factor
- n is the number of runs
- tau_d is a constant related to the radius R_d for d
- v_d is the number of vectors that we at most expect to have to enumerate for d
- tau_r is a constant related to the radius R_r for r
- v_r is the number of vectors that we at most expect to have to enumerate for r
- #errors is the number of sets of n simulated outputs for which R_d or R_r was set to infinity due to sampling errors

Note that the executable first computes estimates for n = s and n = s + 1. It then interpolates a guess for the correct n by using that v_d is reduced by approximately a constant factor for every increment of n by one, and that v_d is greater in general than v_r. Depending on whether this guess is above or below the bound, n is then increased or decreased until the smallest n that satisfies the bound on both v_d and v_r has been found.
