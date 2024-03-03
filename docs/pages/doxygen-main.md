# Qunundrum
This repository contains the source code of MPI programs for computing the probability distributions induced by:

* Shor's order-finding algorithm [[Shor94]](https://doi.org/10.1109/SFCS.1994.365700) [[Shor97]](https://doi.org/10.1137/S0097539795293172)
* Seifert's order-finding algorithm with tradeoffs [[Seifert01]](https://doi.org/10.1007/3-540-45353-9_24)
* Ekerå–Håstad's algorithm for computing short discrete logarithms with tradeoffs [[EH17]](https://doi.org/10.1007/978-3-319-59879-6_20) [(Slides)](https://2017.pqcrypto.org/conference/slides/quantum/ekera-hastad-complex.pdf) [[E20]](https://doi.org/10.1007/s10623-020-00783-2)
* Ekerå–Håstad's algorithm for factoring RSA integers with tradeoffs [[EH17]](https://doi.org/10.1007/978-3-319-59879-6_20) [(Slides)](https://2017.pqcrypto.org/conference/slides/quantum/ekera-hastad-complex.pdf) [[E20]](https://doi.org/10.1007/s10623-020-00783-2)
* Ekerå's algorithm for computing general discrete logarithms and orders with tradeoffs [[E21]](https://doi.org/10.1515/jmc-2020-0006)
* Shor's algorithm for computing general discrete logarithms [[Shor94]](https://doi.org/10.1109/SFCS.1994.365700) [[Shor97]](https://doi.org/10.1137/S0097539795293172) with modifications [[E19p]](https://doi.org/10.48550/arXiv.1905.09084)

Once computed the distributions may be sampled to simulate the quantum algorithms. This is possible for large cryptographically relevant problem instances. Note however that the solution to the problem (i.e. the group order, the discrete logarithm, or in some cases both) must be known.

This repository furthermore contains the source code of MPI programs that estimate the number of vectors that need to be enumerated in the classical lattice-based post-processing algorithms of Ekerå and Ekerå–Håstad, and of MPI programs that execute the post-processing algorithms with respect to simulated outputs from the quantum algorithms. For completeness, implementations of Shor's original post-processing algorithms are also provided.

### Other resources
See [[E21b]](https://doi.org/10.1007/s11128-021-03069-1), [[E22p]](https://doi.org/10.48550/arXiv.2201.07791) and the [factoritall](https://github.com/ekera/factoritall) repository for resources on efficiently factoring any integer after a single call to an order-finding algorithm. See [[E23p]](https://doi.org/10.48550/arXiv.2309.01754) for a more recent analysis of the single-run complexity on Ekerå–Håstad's algorithm. See the [Quaspy](https://github.com/ekera/quaspy) repository for resources on simulating the above algorithms in [Phython3](https://www.python.org) when not making tradeoffs.

### Remarks
Note that the source code in this repository was developed for academic research purposes. It grew out of our research project in an organic manner as research questions were posed and answered. It is distributed "as is" without warranty of any kind, either expressed or implied. For further details on the terms of use, see the [license](LICENSE.md).

It is possible to further optimize portions of the code. However, the current code performs sufficiently well for our purposes. Note furthermore that the portions of the code that pertain to Shor's original algorithm for computing general discrete logarithms are based on a heuristic that lacks an error bound. These portions, and the heuristic, are currently a work in progress.

## Installing and compiling
To compile and run these programs under e.g. [Ubuntu 22.04 LTS](https://releases.ubuntu.com/releases/22.04), first execute:

```console
$ sudo apt install libgmp-dev libmpfr-dev libfplll-dev libopenmpi-dev
$ sudo apt install gcc g++ make openmpi-bin
```

This installs libraries and header files for [GMP](https://gmplib.org), [MPFR](https://www.mpfr.org) and [fpLLL](https://github.com/fplll/fplll), as well as libraries, headers and binaries for [OpenMPI](https://www.open-mpi.org) and for compiling C and C++ sources. You may then proceed to compile the executables:
```console
$ make
```

Under other Linux and Unix distributions, ensure that the tools contained in the aforementioned packages are installed and available in your search paths prior to running the above command. Under other operating systems, you may need to setup build scripts yourself.

### Building the documentation
To build the documentation using [Doxygen](http://www.doxygen.nl), first execute:

```console
$ sudo apt install doxygen graphviz
```

You may then proceed to build the documentation in HTML format from the source files:
```console
$ make documentation
```

## About and acknowledgments
This source code was developed by [Martin Ekerå](mailto:ekera@kth.se), in part at [KTH, the Royal Institute of Technology](https://www.kth.se/en), in Stockholm, [Sweden](https://www.sweden.se). Valuable comments and advice were provided by Johan Håstad throughout the development process.

Funding and support for this work was provided by the Swedish NCSA that is a part of the [Swedish Armed Forces](https://www.mil.se).

Computations were performed on the [Beskow Cray XC40](https://www.pdc.kth.se/about/history-of-pdc/recent-systems-at-pd/beskow-1.737436) supercomputer and its pre- and post-processing cluster [Tegner](https://www.pdc.kth.se/about/history-of-pdc/recent-systems-at-pd/tegner-1.737437) at [PDC](https://www.pdc.kth.se) at [KTH](https://www.kth.se/en). Access was provided by the [Swedish National Infrastructure for Computing (SNIC)](https://www.snic.se). This version of the source code is intended to be run on generic Linux-based clusters.

More recently, computations were performed on the [Dardel HPE Cray EX](https://www.pdc.kth.se/hpc-services/computing-systems/about-the-dardel-hpc-system-1.1053338) supercomputer at [PDC](https://www.pdc.kth.se) at [KTH](https://www.kth.se/en). Access was provided by the [National Academic Infrastructure for Super­computing in Sweden (NAISS)](https://www.naiss.se).
