# The <code>info_linear_distribution</code> executable

## Synopsis
```console
Synopsis: info_linear_distribution <distribution>
```

Prints information on a linear distribution to the console.

### Mandatory command line arguments
An argument <code>\<distribution\></code> where
- <code>\<distribution\></code> is the path to the distribution

## Interpreting the output
The console output contains information on the format
```console
$ ./info_linear_distribution distributions/linear-distribution-det-dim-2048-d-m-2048-s-30.txt
Importing the distribution from "distributions/linear-distribution-det-dim-2048-d-m-2048-s-30.txt"...
Done importing the distribution.

Precision: 192 bits

Parameters:
 m: 2048
 s: 30
 l: 69
 d: 30460294969952480466359944543065930847636338231876963441861263379307409111477531857159412622854886791355025809646176785576524205395222673574326839481353941451128423961795765897813675000697999260501817390040111107221404097648999923218044478853922164062675324489479683307881055802048520145093733180009378708372783039703715042069148935628065372774267860379245095268624820987327737544057352470300082619246647841885572823310298409488067137586044846960970986439672887024218352899975862961482301098046561103362064197724520937259485103053360948141025395770910393006502737792940908063003378133610093965058745405743747817032187

 Flags: 00000000

Slices:
 Probability <= 0.1: 4 slice(s) (0 - 3) (0.478765)
  Dimension 2048: 4 slice(s)
 Probability <= 0.01: 14 slice(s) (4 - 17) (0.47521416)
  Dimension 2048: 14 slice(s)
 Probability <= 0.001: 12 slice(s) (18 - 29) (0.040266555)
  Dimension 2048: 12 slice(s)
 Probability <= 0.0001: 12 slice(s) (30 - 41) (0.005034991)
  Dimension 2048: 12 slice(s)
 Probability <= 1e-05: 12 slice(s) (42 - 53) (0.00060225461)
  Dimension 2048: 12 slice(s)
 Probability <= 1e-06: 6 slice(s) (54 - 59) (1.8313308e-05)
  Dimension 2048: 6 slice(s)
 Probability <= 1e-07: 6 slice(s) (60 - 65) (2.2891635e-06)
  Dimension 2048: 6 slice(s)
 Probability <= 1e-08: 8 slice(s) (66 - 73) (3.065844e-07)
  Dimension 2048: 8 slice(s)
 Probability <= 1e-09: 6 slice(s) (74 - 79) (1.788409e-08)
  Dimension 2048: 6 slice(s)
 Probability <= 1e-10: 2 slice(s) (80 - 81) (1.277435e-09)
  Dimension 2048: 2 slice(s)

Slice: 0 (alpha: 2046, dimension: 2048, probability: 0.127068)
Slice: 1 (alpha: -2046, dimension: 2048, probability: 0.127068)
Slice: 2 (alpha: 2047, dimension: 2048, probability: 0.112314)

(..)

Slice: 79 (alpha: -2019, dimension: 2048, probability: 1.27743e-09)
Slice: 80 (alpha: 2018, dimension: 2048, probability: 6.38717e-10)
Slice: 81 (alpha: -2018, dimension: 2048, probability: 6.38717e-10)

Total count: 82
Total probability: 0.999903894643643375531473
```

Most elements of the printout are rather self-explanatory. The precision is the default precision set when generating the distribution. Note that the precision when computing the distribution may be considerably higher in some steps of the computation to ensure numeric stability, and that lower precision is used to store the distribution on file. The bit flags are used internally to keep track of how the distribution was generated.

As for the other parameters
  - $m$ is the bit length of the logarithm $d$ or order $r$,
  - $s$ is the tradeoff factor such that $\ell = \lceil m / s \rceil$, or zero if $\ell$ was explicitly specified,
  - $\ell$ is the bit length of the second register, and such that $m + \ell$ is the bit length of the first register,
  - $r$ is the order, and
  - $d$ is the logarithm.

The slices section contains a summary of the number of slices that capture a certain share of the probability mass and their dimensions. This is followed by information about each slice in the distribution. The total number of slices in the distribution, and the total probability mass captured by these slices, is last reported.
