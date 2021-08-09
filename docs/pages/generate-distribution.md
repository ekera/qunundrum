# The <code>generate_distribution</code> executable

## Synopsis
```console
Synopsis: mpirun generate_distribution \
   [ -det | -rnd | -exp <d> <r> ] [ -dim-heuristic | -dim <dimension> ] \
      [ -approx-quick | [ -sigma-heuristic | -sigma-optimal ] ] \
         ( [ -s ] <m> <s> { <m> <s> } | -l <m> <l> { <m> <l> } )
```

Computes the distribution induced by Ekerå's algorithm for general discrete logarithms and orders with tradeoffs.

The distribution generated will be assigned an appropriate name and written to the <code>distributions</code> directory. If this directory does not exist, it will be created. If the distribution already exists, an error will be reported. Information for the distribution will be stored along with the distribution, as will collapsed marginal distributions and, if applicable, filtered distributions.

> <b>Note:</b> This is an MPI program. The node with rank zero acts as server. All other nodes are clients, requesting jobs from and reporting back to the server node. A minimum of two nodes is hence required.

### Mandatory command line arguments
Tuples <code>\<m\></code> <code>\<s\></code> where
- <code>\<m\></code> is the length m in bits of the order r
- <code>\<s\></code> is the tradeoff factor s; used to set l = ceil(m / s)

or, if the <code>-l</code> flag is specified, tuples <code>\<m\></code> <code>\<l\></code> where
- <code>\<m\></code> is as above
- <code>\<l\></code> is the parameter l

### Optional command line arguments
Flags controlling the values of d and r (defaults to <code>-det</code>):
- <code>-det</code> selects d and r on (2^(m-1), 2^m) deterministically by reading from Catalan's constant
- <code>-rnd</code> selects r uniformly at random from (2^(m-1), 2^m) and d uniformly at random from [r/2, r)
- <code>-exp \<d\> \<r\></code> explicitly sets d and r to <code>\<d\></code> and <code>\<r\></code> where 0 < <code>\<d\></code> < <code>\<r\></code> and 2^(m-1) < <code>\<r\></code> < 2^m

Flags controlling the approximation method (defaults to <code>-approx-with-error-bound</code>):
- <code>-approx-with-error-bound</code> uses the error-bounded approximation in the paper on computing general discrete logarithms and orders with tradeoffs
- <code>-approx-quick</code> uses the quick and dirty approximation in the introduction to the paper on computing general discrete logarithms and orders with tradeoffs

Flags controlling the selection method for sigma (defaults to <code>-sigma-heuristic</code>):
- <code>-sigma-heuristic</code> adaptively selects sigma according to a heuristic
- <code>-sigma-optimal</code> adaptively exhaust sigma and picks the optimal value

Flags controlling the dimension (defaults to <code>-dim-heuristic</code>):
- <code>-dim-heuristic</code> adaptively sets the dimension according to a heuristic
- <code>-dim \<dimension\></code> sets the dimension to <code>\<dimension\></code>

The dimension specifies the resolution of the histogram. Must be a power of two. The dimension is 2^nu in the paper on computing general discrete logarithms with tradeoffs. See this reference for a description of how the histogram is constructed.
