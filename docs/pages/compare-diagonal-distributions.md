# The <code>compare_diagonal_distributions</code> executable

## Synopsis
```console
Synopsis: compare_diagonal_distributions <distribution1> <distribution2>
```

Compares the slices in two diagonal distributions and prints the result to the console.

### Mandatory command line arguments
Arguments <code>\<distribution1\></code> <code><distribution2\></code> where
- <code>\<distribution1\></code> is the path to the first distribution
- <code>\<distribution2\></code> is the path to the second distribution

## Interpreting the output
The console output contains information on the format
```console
$ ./compare_diagonal_distributions distributions/diagonal-distribution-det-dim-2048-m-2048-sigma-10-s-30.txt distributions/diagonal-distribution-det-dim-2048-m-2048-sigma-0-s-30.txt 
Loading distribution "distributions/diagonal-distribution-det-dim-2048-m-2048-sigma-10-s-30.txt"...
Loading distribution "distributions/diagonal-distribution-det-dim-2048-m-2048-sigma-0-s-30.txt"...

-- slices for eta: 0
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
Warning: Failed to find slice -1 in one distribution. This may be due to a difference in the number of padding bits.
Warning: Failed to find slice 0 in one distribution. This may be due to a difference in the number of padding bits.
Warning: Failed to find slice 1 in one distribution. This may be due to a difference in the number of padding bits.
Warning: Failed to find slice 2 in one distribution. This may be due to a difference in the number of padding bits.
Warning: Failed to find slice 3 in one distribution. This may be due to a difference in the number of padding bits.
Warning: Failed to find slice 4 in one distribution. This may be due to a difference in the number of padding bits.
Warning: Failed to find slice 5 in one distribution. This may be due to a difference in the number of padding bits.
Warning: Failed to find slice 6 in one distribution. This may be due to a difference in the number of padding bits.
Warning: Failed to find slice 7 in one distribution. This may be due to a difference in the number of padding bits.
Warning: Failed to find slice 8 in one distribution. This may be due to a difference in the number of padding bits.

-- slices for eta: 1
2018 2018: 1.26409e-17 3.41096e-12 -- 3.41095e-12
2019 2019: 2.52819e-17 6.82192e-12 -- 6.82189e-12
2020 2020: 5.05637e-17 1.36438e-11 -- 1.36438e-11
2021 2021: 1.01127e-16 2.72877e-11 -- 2.72876e-11
2022 2022: 2.02255e-16 5.45753e-11 -- 5.45751e-11
2023 2023: 4.0451e-16 1.0915e-10 -- 1.0915e-10
2024 2024: 8.09021e-16 2.183e-10 -- 2.183e-10
2025 2025: 1.61805e-15 4.36599e-10 -- 4.36598e-10
2026 2026: 3.23611e-15 8.73191e-10 -- 8.73188e-10
2027 2027: 6.47226e-15 1.74635e-09 -- 1.74635e-09
2028 2028: 1.29447e-14 3.49259e-09 -- 3.49258e-09
2029 2029: 2.58904e-14 6.98474e-09 -- 6.98471e-09
2030 2030: 5.17842e-14 1.39677e-08 -- 1.39676e-08
2031 2031: 1.03582e-13 2.79281e-08 -- 2.7928e-08
2032 2032: 2.0722e-13 5.58272e-08 -- 5.5827e-08
2033 2033: 4.1466e-13 1.11539e-07 -- 1.11538e-07
2034 2034: 8.30206e-13 2.22614e-07 -- 2.22614e-07
2035 2035: 1.66395e-12 4.4338e-07 -- 4.43378e-07
2036 2036: 3.3421e-12 8.79386e-07 -- 8.79383e-07
2037 2037: 6.74115e-12 1.72945e-06 -- 1.72944e-06
2038 2038: 1.37115e-11 3.34298e-06 -- 3.34296e-06
2039 2039: 2.83505e-11 6.23335e-06 -- 6.23332e-06
2040 2040: 6.0499e-11 1.07463e-05 -- 1.07462e-05
2041 2041: 1.36877e-10 1.53528e-05 -- 1.53526e-05
2042 2042: 3.42583e-10 1.24319e-05 -- 1.24315e-05
2043 2043: 9.99672e-10 7.6553e-06 -- 7.6543e-06
2044 2044: 3.50911e-09 0.000481955 -- 0.000481951
2045 2045: 1.38921e-08 0.00840754 -- 0.00840753
2046 2046: 4.3625e-08 0.10323 -- 0.103229
Warning: Failed to find slice -1 in one distribution. This may be due to a difference in the number of padding bits.
Warning: Failed to find slice 0 in one distribution. This may be due to a difference in the number of padding bits.
Warning: Failed to find slice 1 in one distribution. This may be due to a difference in the number of padding bits.
Warning: Failed to find slice 2 in one distribution. This may be due to a difference in the number of padding bits.
Warning: Failed to find slice 3 in one distribution. This may be due to a difference in the number of padding bits.
Warning: Failed to find slice 4 in one distribution. This may be due to a difference in the number of padding bits.
Warning: Failed to find slice 5 in one distribution. This may be due to a difference in the number of padding bits.
Warning: Failed to find slice 6 in one distribution. This may be due to a difference in the number of padding bits.
Warning: Failed to find slice 7 in one distribution. This may be due to a difference in the number of padding bits.
Warning: Failed to find slice 8 in one distribution. This may be due to a difference in the number of padding bits.

-- slices for eta: -1
2018 2018: 1.26409e-17 3.41096e-12 -- 3.41095e-12
2019 2019: 2.52819e-17 6.82192e-12 -- 6.82189e-12
2020 2020: 5.05637e-17 1.36438e-11 -- 1.36438e-11
2021 2021: 1.01127e-16 2.72877e-11 -- 2.72876e-11
2022 2022: 2.02255e-16 5.45754e-11 -- 5.45752e-11
2023 2023: 4.04509e-16 1.09151e-10 -- 1.0915e-10
2024 2024: 8.09018e-16 2.18302e-10 -- 2.18301e-10
2025 2025: 1.61803e-15 4.36606e-10 -- 4.36605e-10
2026 2026: 3.23605e-15 8.73219e-10 -- 8.73216e-10
2027 2027: 6.47205e-15 1.74647e-09 -- 1.74646e-09
2028 2028: 1.29439e-14 3.49305e-09 -- 3.49303e-09
2029 2029: 2.58869e-14 6.98655e-09 -- 6.98652e-09
2030 2030: 5.17703e-14 1.39749e-08 -- 1.39749e-08
2031 2031: 1.03527e-13 2.79571e-08 -- 2.7957e-08
2032 2032: 2.06998e-13 5.59431e-08 -- 5.59429e-08
2033 2033: 4.13776e-13 1.12002e-07 -- 1.12002e-07
2034 2034: 8.26668e-13 2.24469e-07 -- 2.24468e-07
2035 2035: 1.6498e-12 4.50796e-07 -- 4.50795e-07
2036 2036: 3.2855e-12 9.09052e-07 -- 9.09049e-07
2037 2037: 6.51475e-12 1.84811e-06 -- 1.8481e-06
2038 2038: 1.28059e-11 3.8176e-06 -- 3.81758e-06
2039 2039: 2.47285e-11 8.13144e-06 -- 8.13141e-06
2040 2040: 4.60136e-11 1.83324e-05 -- 1.83324e-05
2041 2041: 7.89827e-11 4.55977e-05 -- 4.55976e-05
2042 2042: 1.11764e-10 0.000131819 -- 0.000131819
2043 2043: 8.84641e-11 0.000459811 -- 0.000459811
2044 2044: 5.31111e-11 0.00188964 -- 0.00188964
2045 2045: 2.83224e-09 0.00795752 -- 0.00795752
2046 2046: 3.00071e-08 0.0232952 -- 0.0232952
Warning: Failed to find slice -1 in one distribution. This may be due to a difference in the number of padding bits.
Warning: Failed to find slice 0 in one distribution. This may be due to a difference in the number of padding bits.
Warning: Failed to find slice 1 in one distribution. This may be due to a difference in the number of padding bits.
Warning: Failed to find slice 2 in one distribution. This may be due to a difference in the number of padding bits.
Warning: Failed to find slice 3 in one distribution. This may be due to a difference in the number of padding bits.
Warning: Failed to find slice 4 in one distribution. This may be due to a difference in the number of padding bits.
Warning: Failed to find slice 5 in one distribution. This may be due to a difference in the number of padding bits.
Warning: Failed to find slice 6 in one distribution. This may be due to a difference in the number of padding bits.
Warning: Failed to find slice 7 in one distribution. This may be due to a difference in the number of padding bits.
Warning: Failed to find slice 8 in one distribution. This may be due to a difference in the number of padding bits.

(..)

0.790698 0.996209 -- 0.205511 (max: 0.103229)

$ ./compare_diagonal_distributions distributions/diagonal-distribution-det-dim-2048-m-2048-sigma-5-s-1.txt distributions/diagonal-distribution-det-dim-2048-m-2048-sigma-0-s-1.txt
Loading distribution "distributions/diagonal-distribution-det-dim-2048-m-2048-sigma-5-s-1.txt"...
Loading distribution "distributions/diagonal-distribution-det-dim-2048-m-2048-sigma-0-s-1.txt"...

-- slices for eta: 0
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
Warning: Failed to find slice -1 in one distribution. This may be due to a difference in the number of padding bits.
Warning: Failed to find slice 0 in one distribution. This may be due to a difference in the number of padding bits.
Warning: Failed to find slice 1 in one distribution. This may be due to a difference in the number of padding bits.
Warning: Failed to find slice 2 in one distribution. This may be due to a difference in the number of padding bits.
Warning: Failed to find slice 3 in one distribution. This may be due to a difference in the number of padding bits.

Warning: Slices for eta in [1, 100] have been ignored in the above comparison. These slices are only present in the second distribution.

0.790698 0.790698 -- 0 (max: 0)
```

On each line of the printout, we find $m_1+t$, $m_2+t$: $p_1$, $p_2$ — $\delta$ where
- $m_1$ is the bit length of the order $r$ for the first distribution,
- $m_2$ is the bit length of the order $r$ for the second distribution,
- $t$ is an integer parameter that runs over a sub-interval of $[-30, 30] \cap \mathbb Z$,
- $p_1$ is the probability mass captured by all slices which $m_1 + t$ <= $| \alpha_r | < m_1 + t + 1$ in the first distribution,
- $p_2$ is the probability mass captured by all slices which $m_2 + t$ <= $| \alpha_r | < m_2 + t + 1$ in the second distribution and
- $\delta = |p_1 - p_2|$ is the absolute value of the difference between $p_1$ and $p_2$.

If no slice is present in either distribution for a given $t$, nothing is printed.

If one or more slices are present in only one of the distribution for a given $t$, a warning is printed.

This is followed by the total probability mass captured by the first and second distribution, respectively, the sum of $\delta$, and the maximum $\delta$, on the last line. Note that the totals on the last line only consider the slices that exist in both distributions for the above range of t.

Note furthermore that the above comparison assumes that the two distributions are symmetric, in that only slices for positive $m_1$, $m_2$ are compared and $p_1$, $p_2$ doubled accordingly.
