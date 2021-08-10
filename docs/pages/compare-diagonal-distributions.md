# The <code>compare_diagonal_distributions</code> executable

## Synopsis
```console
Synopsis: compare_diagonal_distributions <distribution1> <distribution2>
```

Compares the slices in two diagonal distributions and prints the result to the console.

### Mandatory command line arguments
An entry <code>\<distribution1\></code> <code><distribution2\></code> where
- <code>\<distribution1\></code> is the paths to the first distribution
- <code>\<distribution2\></code> is the paths to the second distribution

## Interpreting the output
The console output contains information on the format
```console
$ ./compare_diagonal_distributions distributions/diagonal-distribution-det-dim-2048-m-2048-sigma-5-s-1.txt distributions/diagonal-distribution-det-dim-2048-m-2048-sigma-2-s-1.txt 
Loading distribution "distributions/diagonal-distribution-det-dim-2048-m-2048-sigma-5-s-1.txt"...
Loading distribution "distributions/diagonal-distribution-det-dim-2048-m-2048-sigma-2-s-1.txt"...

2018 2018: 1.94434e-09 1.94434e-09 -- 0
2019 2019: 3.88868e-09 3.88868e-09 -- 0
2020 2020: 7.77736e-09 7.77736e-09 -- 0
2021 2021: 1.55547e-08 1.55547e-08 -- 0
2022 2022: 3.11095e-08 3.11095e-08 -- 0
2023 2023: 6.22189e-08 6.22189e-08 -- 0
2024 2024: 1.24438e-07 1.24438e-07 -- 0
2025 2025: 2.48876e-07 2.48876e-07 -- 0
2026 2026: 4.97751e-07 4.97751e-07 -- 0
2027 2027: 9.95503e-07 9.95503e-07 -- 0
2028 2028: 1.99101e-06 1.99101e-06 -- 0
2029 2029: 3.98201e-06 3.98201e-06 -- 0
2030 2030: 7.96402e-06 7.96402e-06 -- 0
2031 2031: 1.5928e-05 1.5928e-05 -- 0
2032 2032: 3.18561e-05 3.18561e-05 -- 0
2033 2033: 6.37122e-05 6.37122e-05 -- 0
2034 2034: 0.000127424 0.000127424 -- 0
2035 2035: 0.000254849 0.000254849 -- 0
2036 2036: 0.000509697 0.000509697 -- 0
2037 2037: 0.00101939 0.00101939 -- 0
2038 2038: 0.00203877 0.00203877 -- 0
2039 2039: 0.00407745 0.00407745 -- 0
2040 2040: 0.00815412 0.00815412 -- 0
2041 2041: 0.016302 0.016302 -- 0
2042 2042: 0.0325541 0.0325541 -- 0
2043 2043: 0.0647103 0.0647103 -- 0
2044 2044: 0.126282 0.126282 -- 0
2045 2045: 0.228819 0.228819 -- 0
2046 2046: 0.305721 0.305721 -- 0
2047 2047: 0.112178 0.112178 -- 0
2048 2048: 0.0471674 0.0471674 -- 0
Warning: Failed to find slice 1 in both distributions. This may be due to a difference in the number of padding bits.
Warning: Failed to find slice 2 in both distributions. This may be due to a difference in the number of padding bits.
Warning: Failed to find slice 3 in both distributions. This may be due to a difference in the number of padding bits.

0.950043 0.950043 -- 0 (max: 0)
```

On each line of the printout, we find (m1+t) (m2+t): p1 p2 -- delta where
- m1 is the order length in bits for the first distribution
- m2 is the order length in bits for the second distribution
- t is an integer parameter that runs over a sub-interval of [-30, 30]
- p1 is the probability mass captured by all slices in alpha_d for which m1+t <= abs(alpha_r) < m1+t+1 in the first distribution
- p2 is the probability mass captured by all slices in alpha_d for which m2+t <= abs(alpha_r) < m2+t+1 in the second distribution
- delta = abs(p1 - p2) is the absolute value of the difference between p1 and p2

If no slice is present in either distribution for a given t, nothing is printed.

If one or more slices are present in only one of the distribution for a given t, a warning is printed.

This is followed by the total probability mass captured by the first and second distribution, respectively, the sum of delta, and the maximum delta, on the last line. Note that the totals on the last line only consider the slices that exist in both distributions for the above range of t.

Note furthermore that the above comparison assumes that the two distributions are symmetric, in that only slices for positive m1, m2 are compared and p1, p2 doubled accordingly.
