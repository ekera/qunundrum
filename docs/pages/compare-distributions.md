# The <code>compare_distributions</code> executable

## Synopsis
```console
Synopsis: compare_distributions <distribution1> <distribution2>
```

Compares the slices in two distributions and prints the result to the console.

### Mandatory command line arguments
An entry <code>\<distribution1\></code> <code><distribution2\></code> where
- <code>\<distribution1\></code> is the paths to the first distribution
- <code>\<distribution2\></code> is the paths to the second distribution

## Interpreting the output
The console output contains information on the format
```console
$ ./compare_distributions distributions/distribution-det-dim-heuristic-sigma-heuristic-m-2048-s-1.txt distributions/distribution-det-dim-heuristic-sigma-optimal-m-2048-s-30.txt
Loading distribution "distributions/distribution-det-dim-heuristic-sigma-heuristic-m-2048-s-1.txt"...
Loading distribution "distributions/distribution-det-dim-heuristic-sigma-optimal-m-2048-s-30.txt"...

(2018, 2025) (2018, 2025): 1.15892e-16 1.15892e-16 -- 0
(2018, -2025) (2018, -2025): 1.15892e-16 1.15892e-16 -- 0
(2018, 2026) (2018, 2026): 2.31783e-16 2.31783e-16 -- 0

(..)

(-2058, -2057) (-2058, -2057): 1.01496e-08 1.01496e-08 -- 4.5826e-16
(-2058, 2058) (-2058, 2058): 1.441e-10 1.47732e-10 -- 3.63189e-12
(-2058, -2058) (-2058, -2058): 1.74753e-05 1.74753e-05 -- 1.29956e-16

0.999891 0.999891 -- 2.99359e-11 (max: 3.63189e-12)
```

On each line of the printout, we find (±m1+td, ±m1+tr) (±m2+td, ±m2+tr): p1 p2 -- delta where
- m1 is the order length in bits for the first distribution
- m2 is the order length in bits for the second distribution
- td, tr are integer parameters that run over a sub-interval of [-30, 30]
- p1 is the probability mass captured by the slice for which m1+td <= alpha_d < m1+td+1 and m1+tr <= alpha_r < m1+tr+1 in the first distribution
- p2 is the probability mass captured by the slice for which m2+td <= alpha_d < m2+td+1 and m2+tr <= alpha_r < m2+tr+1 in the second distribution
- delta = abs(p1 - p2) is the absolute value of the difference between p1 and p2

If a slice is missing from both distributions for a given combination (td, tr), nothing is printed.

If a slice is only present in one of the distribution for a given combination (td, tr), a warning is printed.

This is followed by the total probability mass captured by the slices in first and second distribution, respectively, the sum of delta, and the maximum delta, on the last line. Note that the totals on the last line only consider the slices that exist in both distributions for the above ranges of td, tr.
