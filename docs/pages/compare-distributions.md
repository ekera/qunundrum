# The <code>compare_distributions</code> executable

## Synopsis
```console
Synopsis: compare_distributions <distribution1> <distribution2>
```

Compares the slices in two distributions and prints the result to the console.

### Mandatory command line arguments
Arguments <code>\<distribution1\></code> <code><distribution2\></code> where
- <code>\<distribution1\></code> is the path to the first distribution
- <code>\<distribution2\></code> is the path to the second distribution

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

On each line of the printout, we find $(\pm m_1 + t_d, \pm m_1 + t_r)$, $(\pm m_2 + t_d, \pm m_2 + t_r)$: $p_1$, $p_2$ — $\delta$ where
- $m_1$ is the bit length of the order $r$ for the first distribution,
- $m_2$ is the bit length of the order $r$ for the second distribution,
- $t_d$, $t_r$ are integer parameters that run over a sub-interval of $[-30, 30] \cap \mathbb Z$,
- $p_1$ is the probability mass captured by the slice for which $m_1 + t_d \le \alpha_d < m_1 + t_d + 1$ and $m_1 + t_r \le \alpha_r < m_1 + t_r + 1$ in the first distribution,
- $p_2$ is the probability mass captured by the slice for which $m_2 + t_d \le \alpha_d < m_2 + t_d + 1$ and $m_2 + t_r \le \alpha_r < m_2 + t_r + 1$ in the second distribution, and
- $\delta = |p_1 - p_2|$ is the absolute value of the difference between $p_1$ and $p_2$.

If a slice is missing from both distributions for a given combination $(t_d, t_r)$, nothing is printed.

If a slice is only present in one of the distribution for a given combination $(t_d, t_r)$, a warning is printed.

This is followed by the total probability mass captured by the slices in first and second distribution, respectively, the sum of $\delta$, and the maximum $\delta$, on the last line. Note that the totals on the last line only consider the slices that exist in both distributions for the above ranges of $t_d$ and $t_r$.
