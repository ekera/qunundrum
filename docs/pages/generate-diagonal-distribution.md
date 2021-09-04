# The <code>generate_diagonal_distribution</code> executable

## Synopsis
```console
Synopsis: mpirun generate_diagonal_distribution \
   [ -dim <dimension> ] [ -det | -rnd | -exp <d> <r> ] \
      ( [ -s ] <m> <sigma> <s> { <m> <sigma> <s> } | \
          -l   <m> <sigma> <l> { <m> <sigma> <l> } )
```

Computes a part of the distribution induced by Shor's algorithm for computing general discrete logarithms when modified as in [[E19p]](https://arxiv.org/pdf/1905.09084.pdf).

The full distribution is two-dimensional in (alpha_d, alpha_r), or equivalently, in (theta_d, theta_r).

This executable computes and stores a one-dimensional distribution in alpha_r using the expression for f(theta_r) given in [[E19p]](https://arxiv.org/pdf/1905.09084.pdf). The distribution in phi = theta_d - d/r theta_r is computed on the fly using the expression for h(phi) given in [[E19p]](https://arxiv.org/pdf/1905.09084.pdf) by the executables that use the distribution.

The distribution is said to be diagonal since constructive interference is expected to arise on the diagonal in the argument plane where the angle phi is small.

The distribution generated will be assigned an appropriate name and written to the <code>distributions</code> directory. If this directory does not exist, it will be created. If the distribution already exists, an error will be reported.

> <b>Note:</b> This is an MPI program. The node with rank zero acts as server. All other nodes are clients, requesting jobs from and reporting back to the server node. A minimum of two nodes is hence required.

### Mandatory command line arguments
Tuples <code>\<m\></code> <code>\<sigma\></code> <code>\<s\></code> where

- <code>\<m\></code> is the length m in bits of the order r
- <code>\<sigma\></code> is the padding length sigma
- <code>\<s\></code> is the tradeoff factor s; used to set l = ceil((m + sigma) / s)

or, if the <code>-l</code> flag is specified, tuples <code>\<m\></code> <code>\<sigma\></code> <code>\<l\></code> where
- <code>\<m\></code> and <code>\<sigma\></code> are as above
- <code>\<l\></code> is the parameter l

### Optional command line arguments
Flags controlling the value of d and r (defaults to <code>-det</code>):
- <code>-det</code> selects d and r on (2^(m-1), 2^m) deterministically by reading from Catalan's constant
- <code>-rnd</code> selects r uniformly at random from (2^(m-1), 2^m) and d uniformly at random from [r/2, r)
- <code>-exp \<d\> \<r\></code> explicitly sets d and r to <code>\<d\></code> and <code>\<r\></code> where 0 < <code>\<d\></code> < <code>\<r\></code> and 2^(m-1) < <code>\<r\></code> < 2^m

Flags controlling the dimension (defaults to 2048):
- <code>-dim \<dimension\></code> sets the dimension to <code>\<dimension\></code>

The dimension specifies the resolution of the histogram. Must be a power of two. The dimension is 2^nu in [[E19p]](https://arxiv.org/pdf/1905.09084.pdf). See [[E19p]](https://arxiv.org/pdf/1905.09084.pdf) for a complete description of how the histogram is constructed in two steps.
