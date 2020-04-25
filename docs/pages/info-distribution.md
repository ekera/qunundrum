# The <code>info_distribution</code> executable

## Synopsis
```console
Synopsis: info_distribution <distribution>
```

Prints information on a distribution to the console.

### Mandatory command line arguments
An entry <code>\<distribution\></code> where
- <code>\<distribution\></code> is the path to the distribution

## Interpreting the output
The console output contains information on the format
```console
$ ./info_distribution distributions/distribution-det-dim-heuristic-sigma-optimal-m-2048-s-30.txt
Importing the distribution from "distributions/distribution-det-dim-heuristic-sigma-optimal-m-2048-s-30.txt"...
Done importing the distribution.

Precision: 192 bits

Parameters:
 m: 2048
 s: 30
 l: 69
 r: 30959135869724094223561946967759777687318543652613328887524049258265981260886206372938806681042204705610549935046567627727525429118430339136061372053091244615966840835393590159240524459947353686947729505644176567461420475219958791347385832110082776343412671298825110972145681729023386555945211722398867532698299006188695002392503693847145093186575996726834038856940577486661479544050711638520321533582061933278231459319002866121409661642907458459385642418227936246612499534284235479071206963509080588500807548704828061927447303633043844756810571062624394575985934858964201670695882625659870744553251307284295175096540
 d: 30460294969952480466359944543065930847636338231876963441861263379307409111477531857159412622854886791355025809646176785576524205395222673574326839481353941451128423961795765897813675000697999260501817390040111107221404097648999923218044478853922164062675324489479683307881055802048520145093733180009378708372783039703715042069148935628065372774267860379245095268624820987327737544057352470300082619246647841885572823310298409488067137586044846960970986439672887024218352899975862961482301098046561103362064197724520937259485103053360948141025395770910393006502737792940908063003378133610093965058745405743747817032187

Slices:
 Probability <= 0.01: 20 slice(s) (0 - 19) (0.43964045) (2)
  Dimension 128: 18 slice(s)
  Dimension 256: 2 slice(s)
 Probability <= 0.001: 130 slice(s) (20 - 149) (0.44469007) (52)
  Dimension 128: 90 slice(s)
  Dimension 256: 40 slice(s)
 Probability <= 0.0001: 268 slice(s) (150 - 417) (0.098847694) (140)
  Dimension 128: 156 slice(s)
  Dimension 256: 112 slice(s)
 Probability <= 1e-05: 394 slice(s) (418 - 811) (0.014560801) (224)
  Dimension 128: 176 slice(s)
  Dimension 256: 218 slice(s)
 Probability <= 1e-06: 474 slice(s) (812 - 1285) (0.0018732608) (280)
  Dimension 128: 212 slice(s)
  Dimension 256: 262 slice(s)
 Probability <= 1e-07: 608 slice(s) (1286 - 1893) (0.00025392533) (324)
  Dimension 128: 312 slice(s)
  Dimension 256: 296 slice(s)
 Probability <= 1e-08: 588 slice(s) (1894 - 2481) (2.2237962e-05) (316)
  Dimension 128: 554 slice(s)
  Dimension 256: 34 slice(s)
 Probability <= 1e-09: 602 slice(s) (2482 - 3083) (2.4792181e-06) (298)
  Dimension 128: 580 slice(s)
  Dimension 256: 22 slice(s)
 Probability <= 1e-10: 752 slice(s) (3084 - 3835) (3.2243144e-07) (308)
  Dimension 128: 728 slice(s)
  Dimension 256: 24 slice(s)
 Probability <= 1e-11: 664 slice(s) (3836 - 4499) (2.6667493e-08) (344)
  Dimension 128: 664 slice(s)
 Probability <= 1e-12: 508 slice(s) (4500 - 5007) (2.1058233e-09) (232)
  Dimension 128: 508 slice(s)
 Probability <= 1e-13: 512 slice(s) (5008 - 5519) (2.303732e-10) (200)
  Dimension 128: 512 slice(s)
 Probability <= 1e-14: 340 slice(s) (5520 - 5859) (1.3237832e-11) (148)
  Dimension 128: 340 slice(s)
 Probability <= 1e-15: 260 slice(s) (5860 - 6119) (1.0999101e-12) (104)
  Dimension 128: 260 slice(s)
 Probability <= 1e-16: 244 slice(s) (6120 - 6363) (1.0977292e-13) (92)
  Dimension 128: 244 slice(s)

Slice: 0 (alpha_d: 2047, alpha_r: 2046, dimension: 128, probability: 0.0437698, error: 6.32627e-11)
Slice: 1 (alpha_d: -2047, alpha_r: -2046, dimension: 128, probability: 0.0437698, error: 6.32627e-11)
Slice: 2 (alpha_d: 2046, alpha_r: 2046, dimension: 128, probability: 0.0369427, error: 3.43951e-11)

(..)

Slice: 6361 (alpha_d: -2029, alpha_r: -2056, dimension: 128, probability: 1.00698e-16, error: 4.27716e-18) **
Slice: 6362 (alpha_d: -2029, alpha_r: 2056, dimension: 128, probability: 1.00698e-16, error: 4.27716e-18) **
Slice: 6363 (alpha_d: 2029, alpha_r: -2056, dimension: 128, probability: 1.00698e-16, error: 4.27716e-18) **

Filtered number of slices: 5890 / 6364
Filtered probability: 0.999574766661506914917478
Filtered error: 1.06832964729267942872048e-06

Total count: 6364
Total probability: 0.999891265127786492276137
Total error: 2.276648537893891726088e-05
```

Most elements of the printout are rather self-explanatory. The precision is the default precision set when generating the distribution. Note that the precision when computing the distribution may be considerably higher in some steps of the computation to ensure numeric stability, and that lower precision is used to store the distribution on file.

As for the other parameters
  - m is the length in bits of the order
  - s is the tradeoff factor
  - l is the padding length, or zero, if l was explicitly specified
  - r is the order
  - d is the logarithm

The slices section contains a summary of the number of slices that capture a certain share of the probability mass and their dimensions. Note that slices are in general <i>computed</i> for a much large dimension than is reported in the printout, but scaled down when <i>stored</i> in the distribution to save disk space. It is the stored dimension that is reported.

Note furthermore that the values within parenthesis in the summary section is the probability captured by the slices in the group, followed by the number of slices for which at least one warning flag was raised during the course of the computation. When the slice is computed, a probability estimate and an upper bound on the error in the estimate are computed in many different points throughout the region of the argument plane covered by the slice.

A warning is raised if the error bound is close to or exceeds the probability estimate in at least one point. This needs not be an issue in itself, but it indicates that the approximation used to produce the estimate is being pushed close to its limit. This typically happens for large tradeoff factors as in the above example.

The summary section is followed by detailed information about each slice in the distribution. Slices where the error is within a factor 10^2 of the total probability mass captured by the frame are marked with \*\*.

The total number of slices in the distribution, the total probability mass captured by these slices, and the total approximation error for these slices, is last reported, along with the corresponding count, probability and error for a filtered distribution in which the slices marked with \*\* have been left out, reducing the error at the expense of loosing probability mass. The slices marked with \*\* can be filtered out is desired, see the [<code>filter_distribution</code>](filter-distribution.md) executable.