# The <code>generate_diagonal_distribution</code> executable

## Synopsis
```console
Synopsis: mpirun generate_diagonal_distribution \
   [ -dim <dimension> ] [ -eta-bound <eta-bound> ] \
      [ -det | -rnd | -exp <d> <r> ] \
         ( [ -s ] <m> <sigma> <s> { <m> <sigma> <s> } | \
             -l   <m> <sigma> <l> { <m> <sigma> <l> } )
```

Computes a part of the distribution induced by Shor's algorithm for computing general discrete logarithms when modified as in [[E19p]](https://doi.org/10.48550/arXiv.1905.09084).

The full distribution is two-dimensional in $(\alpha_d, \alpha_r)$. This executable computes and stores a distribution in $(\alpha_r, \eta)$ for $\eta \in [-B_\eta, B_\eta] \cap \mathbb Z$ via the expression for $f_\eta(\theta_r)$ given in [[E19p]](https://doi.org/10.48550/arXiv.1905.09084).

When the distribution is to be sampled by other executables, $(\alpha_r, \eta)$ is first sampled from the part of the distribution computed by this executable. Then $j$ is sampled uniformly at random from all values of $j$ that yield $\alpha_r$. Finally, $k$ is then sampled given $j$ and $\eta$ via the angle $\phi_\eta$ and the expression for $h(\phi_\eta)$ given in [[E19p]](https://doi.org/10.48550/arXiv.1905.09084)

The distribution is said to be diagonal since constructive interference is expected to arise on the "diagonal" in the argument plane where the angle $\phi_\eta$ is small.

The distribution generated will be assigned an appropriate name and written to the <code>distributions</code> directory. If this directory does not exist, it will be created. If the distribution already exists, an error will be reported.

> <b>Note:</b> This is an MPI program. The node with rank zero acts as server. All other nodes are clients, requesting jobs from and reporting back to the server node. A minimum of two nodes is hence required.

### Mandatory command line arguments
Tuples <code>\<m\></code> <code>\<sigma\></code> <code>\<s\></code> where

- <code>\<m\></code> is the bit length $m$ of the order $r$
- <code>\<sigma\></code> is the padding length $\varsigma$
- <code>\<s\></code> is the tradeoff factor $s$; used to set $\ell = \lceil m / s \rceil$

or, if the <code>-l</code> flag is specified, tuples <code>\<m\></code> <code>\<sigma\></code> <code>\<l\></code> where
- <code>\<m\></code> and <code>\<sigma\></code> are as above
- <code>\<l\></code> is the parameter $\ell$

### Optional command line arguments
Flags specifying the value of $d$ and $r$ (defaults to <code>-det</code>):
- <code>-det</code> selects $d$ and $r$ deterministically by reading from Catalan's constant
- <code>-rnd</code> selects $r$ uniformly at random from $(2^{m-1}, 2^m)$ and $d$ uniformly at random from $[r/2, r)$
- <code>-exp \<d\> \<r\></code> explicitly sets $d$ and $r$ to <code>\<d\></code> and <code>\<r\></code> where it is required that $0 < d < r$ and $2^{m-1} < r < 2^m$

Flag specifying the bound $B_\eta$ (defaults to 0):
- <code>-eta-bound \<eta-bound\></code> sets $B_\eta$ to <code>\<eta-bound\></code>

   The bound $B_\eta$ controls how many peak indices $\eta$ are included in the distribution: More specifically, the peaks with peak indices $\eta \in [-B_\eta, B_\eta] \cap \mathbb Z$ are included.

Flag specifying the dimension (defaults to 2048):
- <code>-dim \<dimension\></code> sets the dimension to <code>\<dimension\></code>

   The dimension specifies the resolution of the histogram. Must be a power of two.
