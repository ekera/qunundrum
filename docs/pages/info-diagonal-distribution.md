# The <code>info_diagonal_distribution</code> executable

## Synopsis
```console
Synopsis: info_diagonal_distribution <distribution>
```

Prints information on a diagonal distribution to the console.

### Mandatory command line arguments
An argument <code>\<distribution\></code> where
- <code>\<distribution\></code> is the path to the distribution

## Interpreting the output
The console output contains information on the format
```console
$ ./info_diagonal_distribution distributions/diagonal-distribution-det-dim-2048-m-2048-sigma-12-s-30.txt
Importing the distribution from "distributions/diagonal-distribution-det-dim-2048-m-2048-sigma-12-s-30.txt"...
Done importing the distribution.

Precision: 192 bits

Parameters:
 m:         2048
 sigma:     12
 s:         30
 l:         69

 eta_bound: 0

 r:         30959135869724094223561946967759777687318543652613328887524049258265981260886206372938806681042204705610549935046567627727525429118430339136061372053091244615966840835393590159240524459947353686947729505644176567461420475219958791347385832110082776343412671298825110972145681729023386555945211722398867532698299006188695002392503693847145093186575996726834038856940577486661479544050711638520321533582061933278231459319002866121409661642907458459385642418227936246612499534284235479071206963509080588500807548704828061927447303633043844756810571062624394575985934858964201670695882625659870744553251307284295175096540
 d:         30460294969952480466359944543065930847636338231876963441861263379307409111477531857159412622854886791355025809646176785576524205395222673574326839481353941451128423961795765897813675000697999260501817390040111107221404097648999923218044478853922164062675324489479683307881055802048520145093733180009378708372783039703715042069148935628065372774267860379245095268624820987327737544057352470300082619246647841885572823310298409488067137586044846960970986439672887024218352899975862961482301098046561103362064197724520937259485103053360948141025395770910393006502737792940908063003378133610093965058745405743747817032187

 Flags:     00000000

Slices:
 Probability <= 0.1: 4 slice(s) (0 - 3) (0.53454005)
  Dimension 2048: 4 slice(s)
 Probability <= 0.01: 12 slice(s) (4 - 15) (0.40779586)
  Dimension 2048: 12 slice(s)
 Probability <= 0.001: 14 slice(s) (16 - 29) (0.052583708)
  Dimension 2048: 14 slice(s)
 Probability <= 0.0001: 12 slice(s) (30 - 41) (0.004446153)
  Dimension 2048: 12 slice(s)
 Probability <= 1e-05: 12 slice(s) (42 - 53) (0.0005549799)
  Dimension 2048: 12 slice(s)
 Probability <= 1e-06: 6 slice(s) (54 - 59) (2.7874071e-05)
  Dimension 2048: 6 slice(s)
 Probability <= 1e-07: 8 slice(s) (60 - 67) (3.7331346e-06)
  Dimension 2048: 8 slice(s)
 Probability <= 1e-08: 6 slice(s) (68 - 73) (2.1776618e-07)
  Dimension 2048: 6 slice(s)
 Probability <= 1e-09: 6 slice(s) (74 - 79) (2.7220773e-08)
  Dimension 2048: 6 slice(s)
 Probability <= 1e-10: 2 slice(s) (80 - 81) (1.9443409e-09)
  Dimension 2048: 2 slice(s)

Slice: 0 (alpha_r: 2046, eta: 0, dimension: 2048, probability: 0.152861)
Slice: 1 (alpha_r: -2046, eta: 0, dimension: 2048, probability: 0.152861)
Slice: 2 (alpha_r: 2045, eta: 0, dimension: 2048, probability: 0.114409)

(..)

Slice: 79 (alpha_r: -2019, eta: 0, dimension: 2048, probability: 1.94434e-09)
Slice: 80 (alpha_r: 2018, eta: 0, dimension: 2048, probability: 9.7217e-10)
Slice: 81 (alpha_r: -2018, eta: 0, dimension: 2048, probability: 9.7217e-10)

Total count: 82
Total probability: 0.999952605676525317873575
```

Most elements of the printout are rather self-explanatory. The precision is the default precision set when generating the distribution. Note that the precision when computing the distribution may be considerably higher in some steps of the computation to ensure numeric stability, and that lower precision is used to store the distribution on file. The bit flags are used internally to keep track of how the distribution was generated.

As for the other parameters
- $m$ is the bit length of the order $r$,
- $\varsigma$ is the padding length in $m$, and such that $m + \varsigma$ is the bit length of the first register,
- $s$ is the tradeoff factor such that $\ell = \lceil m / s \rceil$, or zero if $\ell$ was explicitly specified,
- $\ell$ is the bit length of the second register,
- $r$ is the order, and
- $d$ is the logarithm.

The slices section contains a summary of the number of slices that capture a certain share of the probability mass and their dimensions. This is followed by information about each slice in the distribution. The total number of slices in the distribution, and the total probability mass captured by these slices, is last reported.