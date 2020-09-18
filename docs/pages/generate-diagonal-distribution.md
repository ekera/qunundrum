# The <code>generate_diagonal_distribution</code> executable

## Synopsis
```console
Synopsis: mpirun generate_diagonal_distribution \
   [ -dim <dimension> ] [ -bound-delta <bound> ] \
      [ -det | -rnd | -exp <d> <r> ] <m> <l> { <m> <l> }
```

Computes the distribution induced by Shor's algorithm for computing general discrete logarithms when tweaked as in [[E19]](https://arxiv.org/pdf/1905.09084.pdf).

The aforementioned distribution is a distribution in two arguments (alpha_d, alpha_r), or equivalently, in two angle pairs (theta_d, theta_r).
The argument pairs are related in that alpha_d = round(alpha_r d/r) + Delta for some small Delta.
In this software such distributions are called diagonal as all probability mass in on the diagonal were alpha_d is approximated round(alpha_r d/r).
The executable is named accordingly.

The distribution generated will be assigned an appropriate name and written to the <code>distributions</code> directory. If this directory does not exist, it will be created. If the distribution already exists, an error will be reported.

> <b>Note:</b> This is an MPI program. The node with rank zero acts as server. All other nodes are clients, requesting jobs from and reporting back to the server node. A minimum of two nodes is hence required.

### Mandatory command line arguments
Tuples <code>\<m\></code> <code>\<l\></code> where
   
- <code>\<m\></code> is the bit length of the order r
- <code>\<l\></code> is the number of padding bits l

### Optional command line arguments
Flags controlling the value of d and r (defaults to <code>-det</code>):
- <code>-det</code> selects d and r on [2^(m-1), 2^m) deterministically by reading from Catalan's constant
- <code>-rnd</code> selects r uniformly at random from [2^(m-1), 2^m) and d uniformly at random from [r/2, r)
- <code>-exp \<d\> \<r\></code> explicitly sets d and r to <code>\<d\></code> and <code>\<r\></code> where 0 < <code>\<d\></code> < <code>\<r\></code> and 2^(m-1) <= <code>\<r\></code> < 2^m

Flags controlling the bound on |Delta| (defaults to: 20):
- <code>-bound-delta \<bound\></code> sets the bound on |Delta| to <code>\<bound\></code>

Flags controlling the dimension (defaults to <code>-dim-heuristic</code>):
- <code>-dim-heuristic</code> adaptively sets the dimension according to a heuristic
- <code>-dim \<dimension\></code> sets the dimension to <code>\<dimension\></code>

The dimension specifies the resolution of the histogram. Must be a power of two. The dimension is 2^nu in the paper on revisiting Shor's algorithm. See this reference for a description of how the histogram is constructed.
