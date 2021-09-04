# The <code>compare_linear_distributions</code> executable

## Synopsis
```console
Synopsis: compare_linear_distributions <distribution1> <distribution2>
```

Compares the slices in two linear distributions and prints the result to the console.

### Mandatory command line arguments
An entry <code>\<distribution1\></code> <code><distribution2\></code> where
- <code>\<distribution1\></code> is the paths to the first distribution
- <code>\<distribution2\></code> is the paths to the second distribution

## Interpreting the output
The console output contains information on the format
```console
$ ./compare_linear_distributions distributions/linear-distribution-det-dim-2048-d-m-2048-s-30.txt distributions/collapsed-d-distribution-det-dim-heuristic-sigma-optimal-m-2048-s-30.txt
Loading distribution "distributions/linear-distribution-det-dim-2048-d-m-2048-s-30.txt"...
Loading distribution "distributions/collapsed-d-distribution-det-dim-heuristic-sigma-optimal-m-2048-s-30.txt"...

2018 2018: 1.27743e-09 1.27743e-09 -- 5.69146e-16
2019 2019: 2.55487e-09 2.55487e-09 -- 6.74726e-16
2020 2020: 5.10974e-09 5.10974e-09 -- 8.85885e-16
2021 2021: 1.02195e-08 1.02195e-08 -- 5.64611e-16
2022 2022: 2.0439e-08 2.0439e-08 -- 6.65654e-16
2023 2023: 4.08779e-08 4.08779e-08 -- 8.67741e-16
2024 2024: 8.17558e-08 8.17558e-08 -- 5.68057e-16
2025 2025: 1.63512e-07 1.63512e-07 -- 6.72548e-16
2026 2026: 3.27023e-07 3.27023e-07 -- 1.3451e-15
2027 2027: 6.54047e-07 6.54047e-07 -- 2.04573e-15
2028 2028: 1.30809e-06 1.30809e-06 -- 4.09146e-15
2029 2029: 2.61619e-06 2.61619e-06 -- 7.78012e-15
2030 2030: 5.23237e-06 5.23237e-06 -- 1.55602e-14
2031 2031: 1.04647e-05 1.04647e-05 -- 3.11205e-14
2032 2032: 2.09295e-05 2.09295e-05 -- 6.2241e-14
2033 2033: 4.1859e-05 4.1859e-05 -- 1.23861e-13
2034 2034: 8.3718e-05 8.3718e-05 -- 2.47722e-13
2035 2035: 0.000167436 0.000167436 -- 4.95443e-13
2036 2036: 0.000334872 0.000334872 -- 9.90155e-13
2037 2037: 0.000669743 0.000669743 -- 1.98031e-12
2038 2038: 0.00133948 0.00133948 -- 3.95994e-12
2039 2039: 0.00267893 0.00267893 -- 7.91969e-12
2040 2040: 0.00535755 0.00535755 -- 1.58378e-11
2041 2041: 0.0107127 0.0107127 -- 3.16626e-11
2042 2042: 0.0214064 0.0214064 -- 6.32044e-11
2043 2043: 0.0426606 0.0426606 -- 1.25586e-10
2044 2044: 0.0841153 0.0841153 -- 2.44702e-10
2045 2045: 0.158948 0.158948 -- 4.40719e-10
2046 2046: 0.254137 0.254137 -- 5.80204e-10
2047 2047: 0.224628 0.224628 -- 2.57175e-10
2048 2048: 0.094866 0.094866 -- 1.17928e-10
2049 2049: 0.0486841 0.0486841 -- 9.94612e-11
2050 2050: 0.0245339 0.0245339 -- 1.20292e-10
2051 2051: 0.0122926 0.0122926 -- 2.10369e-10
2052 2052: 0.00614957 0.00614957 -- 4.20905e-10
2053 2053: 0.00307519 0.00307519 -- 1.10974e-09
2054 2054: 0.00153765 0.00153765 -- 1.04106e-09
2055 2055: 0.00076883 0.000768809 -- 2.10207e-08
2056 2056: 0.000384416 0.000384262 -- 1.54048e-07
2057 2057: 0.000192208 0.00019125 -- 9.58245e-07
2058 2058: 9.61041e-05 8.46107e-05 -- 1.14934e-05

0.999904 0.999891 -- 1.26316e-05 (max: 1.14934e-05)
```

On each line of the printout, we find (m1+t) (m2+t): p1 p2 -- delta where
- m1 is the order or logarithm length in bits for the first distribution
- m2 is the order or logarithm length in bits for the second distribution
- t is an integer parameter that runs over a sub-interval of [-30, 30]
- p1 is the probability mass captured by the slice for which m1+t <= abs(alpha) < m1+t+1 in the first distribution
- p2 is the probability mass captured by the slice for which m2+t <= abs(alpha) < m2+t+1 in the second distribution
- delta = abs(p1 - p2) is the absolute value of the difference between p1 and p2

If a slice is missing from both distributions for a given t, nothing is printed.

If a slice is only present in one of the distribution for a given t, a warning is printed.

This is followed by the total probability mass captured by the first and second distribution, respectively, the sum of delta, and the maximum delta, on the last line. Note that the totals on the last line only consider the slices that exist in both distributions for the above range of t.

Note furthermore that the above comparison assumes that the two distributions are symmetric, in that only slices for positive m1, m2 are compared and p1, p2 doubled accordingly.
