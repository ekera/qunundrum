# The <code>info_diagonal_distribution</code> executable

## Synopsis
```console
Synopsis: info_diagonal_distribution <distribution>
```

Prints information on a diagonal distribution to the console.

### Mandatory command line arguments
An entry <code>\<distribution\></code> where
- <code>\<distribution\></code> is the path to the distribution

## Interpreting the output
The console output contains information on the format
```console
$ ./info_diagonal_distribution distributions/diagonal-distribution-det-dim-2048-m-2048-l-2.txt
Importing the distribution from "distributions/diagonal-distribution-det-dim-2048-m-2048-l-2.txt"...
Done importing the distribution.

Precision: 192 bits

Parameters:
 m: 2048
 l: 2
 r: 30959135869724094223561946967759777687318543652613328887524049258265981260886206372938806681042204705610549935046567627727525429118430339136061372053091244615966840835393590159240524459947353686947729505644176567461420475219958791347385832110082776343412671298825110972145681729023386555945211722398867532698299006188695002392503693847145093186575996726834038856940577486661479544050711638520321533582061933278231459319002866121409661642907458459385642418227936246612499534284235479071206963509080588500807548704828061927447303633043844756810571062624394575985934858964201670695882625659870744553251307284295175096540
 d: 30460294969952480466359944543065930847636338231876963441861263379307409111477531857159412622854886791355025809646176785576524205395222673574326839481353941451128423961795765897813675000697999260501817390040111107221404097648999923218044478853922164062675324489479683307881055802048520145093733180009378708372783039703715042069148935628065372774267860379245095268624820987327737544057352470300082619246647841885572823310298409488067137586044846960970986439672887024218352899975862961482301098046561103362064197724520937259485103053360948141025395770910393006502737792940908063003378133610093965058745405743747817032187

 Flags: 00000000

Slices:
 Probability <= 0.1: 2 slice(s) (0 - 1) (0.23472893)
  Dimension 2048: 2 slice(s)
 Probability <= 0.01: 16 slice(s) (2 - 17) (0.5248806)
  Dimension 2048: 16 slice(s)
 Probability <= 0.001: 38 slice(s) (18 - 55) (0.13183268)
  Dimension 2048: 38 slice(s)
 Probability <= 0.0001: 122 slice(s) (56 - 177) (0.040536544)
  Dimension 2048: 122 slice(s)
 Probability <= 1e-05: 344 slice(s) (178 - 521) (0.011851516)
  Dimension 2048: 344 slice(s)
 Probability <= 1e-06: 334 slice(s) (522 - 855) (0.0014162694)
  Dimension 2048: 334 slice(s)
 Probability <= 1e-07: 262 slice(s) (856 - 1117) (0.00010043237)
  Dimension 2048: 262 slice(s)
 Probability <= 1e-08: 278 slice(s) (1118 - 1395) (1.108774e-05)
  Dimension 2048: 278 slice(s)
 Probability <= 1e-09: 276 slice(s) (1396 - 1671) (1.0613734e-06)
  Dimension 2048: 276 slice(s)
 Probability <= 1e-10: 262 slice(s) (1672 - 1933) (1.0054359e-07)
  Dimension 2048: 262 slice(s)
 Probability <= 1e-11: 256 slice(s) (1934 - 2189) (1.0168105e-08)
  Dimension 2048: 256 slice(s)
 Probability <= 1e-12: 228 slice(s) (2190 - 2417) (9.2284383e-10)
  Dimension 2048: 228 slice(s)
 Probability <= 1e-13: 124 slice(s) (2418 - 2541) (5.8757003e-11)
  Dimension 2048: 124 slice(s)

Slice: 0 (alpha_r: 2046, delta: 0, dimension: 2048, probability: 0.117364)
Slice: 1 (alpha_r: -2046, delta: 0, dimension: 2048, probability: 0.117364)
Slice: 2 (alpha_r: 2045, delta: 0, dimension: 2048, probability: 0.0890348)

(..)

Slice: 2539 (alpha_r: -2018, delta: 20, dimension: 2048, probability: 1.23249e-13)
Slice: 2540 (alpha_r: 2018, delta: -20, dimension: 2048, probability: 1.23247e-13)
Slice: 2541 (alpha_r: -2018, delta: -20, dimension: 2048, probability: 1.23247e-13)

Total count: 2542
Total probability: 0.94535922248991521594369
```

Most elements of the printout are rather self-explanatory. The precision is the default precision set when generating the distribution. Note that the precision when computing the distribution may be considerably higher in some steps of the computation to ensure numeric stability, and that lower precision is used to store the distribution on file. The bit flags are used internally to keep track of how the distribution was generated.

As for the other parameters
  - m is the length in bits of the order
  - l is the padding length
  - r is the order
  - d is the logarithm

The slices section contains a summary of the number of slices that capture a certain share of the probability mass and their dimensions. This is followed by information about each slice in the distribution. The total number of slices in the distribution, and the total probability mass captured by these slices, is last reported.
