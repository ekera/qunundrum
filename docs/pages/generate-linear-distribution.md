# The <code>generate_linear_distribution</code> executable

## Synopsis
```console
Synopsis: mpirun generate_linear_distribution \
   [ -d | -r ] [ -dim <dimension> ] \
      [ -min | -max | -det | -rnd | -exp <value> ] \
         ( [ -s ] <m> <s> { <m> <s> } | -l <m> <l> { <m> <l> } )
```

Computes the distribution induced by Ekerå–Håstad's algorithm for computing short discrete logarithms if the <code>-d</code> flag is specified. This is the default. Computes the distribution induced by Shor's (for $s = 1$) and Seifert's (for $s > 1$) order-finding algorithms if the <code>-r</code> flag is specified.

All of the aforementioned distributions are distributions in a argument $\alpha$, or equivalently, angle $\theta$. In this software such distributions are called linear. The executable is named accordingly.

The distribution generated will be assigned an appropriate name and written to the <code>distributions</code> directory. If this directory does not exist, it will be created. If the distribution already exists, an error will be reported.

> <b>Note:</b> This is an MPI program. The node with rank zero acts as server. All other nodes are clients, requesting jobs from and reporting back to the server node. A minimum of two nodes is hence required.

### Mandatory command line arguments
Tuples <code>\<m\></code> <code>\<s\></code> where
- <code>\<m\></code> is the bit length $m$ of the logarithm $d$ or order $r$
- <code>\<s\></code> is the tradeoff factor $s$; used to set $\ell = \lceil m / s \rceil$

or, if the <code>-l</code> flag is specified, tuples <code>\<m\></code> <code>\<l\></code> where
- <code>\<m\></code> is as above
- <code>\<l\></code> is the parameter $\ell$

### Optional command line arguments
Flags specifying the value of the logarithm $d$ or order $r$ (defaults to <code>-max</code>):
- <code>-min</code> selects minimal $d$ or $r$ equal to $2^{m - 1} + 1$
- <code>-max</code> selects maximal $d$ or $r$ equal to $2^m - 1$
- <code>-det</code> selects $d$ or $r$ on $(2^{m-1}, 2^m)$ deterministically by reading from Catalan's constant
- <code>-rnd</code> selects $d$ or $r$ uniformly at random from $(2^{m-1}, 2^m)$
- <code>-exp \<value\></code> explicitly sets $d$ or $r$ to <code>\<value\></code> where $0 < d < 2^m$ and $2^{m-1} < r < 2^m$

Flag specifying the dimension (defaults to 2048):
- <code>-dim \<dimension\></code> sets the dimension to <code>\<dimension\></code>

   The dimension specifies the resolution of the histogram. Must be a power of two.
