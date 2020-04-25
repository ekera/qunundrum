# The <code>estimate_runs_linear_distribution</code> executable

## Synopsis
```console
Synopsis: mpirun estimate_runs_linear_distribution \
   [ -bound <bound> ] <distribution> { <distribution> }
```

Estimates the number of runs required to solve a linear distribution for a logarithm d or order r depending on the type of distribution input with at least 99% success probability when accepting to enumerate at most B vectors in the lattice.

The results is written to the console and to <code>logs/estimate-runs-linear.txt</code>.

> <b>Note:</b> This is an MPI program. The node with rank zero acts as server. All other nodes are clients, requesting jobs from and reporting back to the server node. A minimum of two nodes is hence required.

### Mandatory command line arguments
Entries <code>\<distribution\></code> where
- <code>\<distribution\></code> is the path to a distribution

### Optional command line arguments
Flags controlling the bound B (defaults to 2):
- <code>-bound <bound\></code> sets the bound B to <code>\<bound\></code>

## Interpreting the output
The log file <code>logs/estimate-runs-linear.txt</code> is on the format
```
# Processing: linear-distribution-det-dim-2048-d-m-2048-s-30.txt
2048 30 30 -- 6.762732 2.35905E+74 <2850>
2048 30 31 -- 6.771063 2.12435E+56 <2928>
2048 30 35 -- 6.884955 1.17595E-15 <3406>
2048 30 34 -- 6.834919 427.997 <3275>

# Processing: linear-distribution-det-dim-2048-r-m-2048-s-30.txt
2048 30 30 -- 5.741114 6.90434E+64 <1324>
2048 30 31 -- 5.752430 3.27218E+46 <1430>
2048 30 34 -- 5.837407 1.32327E-08 <1649>
2048 30 33 -- 5.789114 1.07427E+10 <1562>

# Processing: linear-distribution-max-dim-2048-rsa-n-2048-delta-20-m-1023-l-1003.txt
1023 0 1 -- 4.355016 1.38282E+09 <94>
1023 0 2 -- 4.849680 3.47726E-291 <175>

# Processing: linear-distribution-max-dim-2048-r-m-2048-s-1.txt
2048 1 1 -- 3.372572 340.15 <51>
2048 1 2 -- 3.856419 1.11801E-612 <94>
```
where we find m, s, n -- tau v \<#errors\> on each line, and 
- m is the bit length of the logarithm d or order r
- s is the tradeoff factor (set to zero if l was explicitly specified)
- n is the number of runs
- tau is a constant related to the radius R
- v is the number of vectors that we at most expect to have to enumerate
- #errors is the number of sets of n simulated outputs for which R was set to infinity due to sampling errors

Note that the executable first computes estimates for n = s and n = s + 1. It then interpolates a guess for the correct n by using that v is reduced by approximately a constant factor for every increment of n by one. Depending on whether this guess is above or below the bound, n is then increased or decreased until the smallest n that satisfies the bound on v has been found.
