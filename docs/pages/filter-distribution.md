# The <code>filter_distribution</code> executable

## Synopsis
```console
Synopsis: filter_distribution <distribution> { <distribution> }
```

Filters out slices from a distribution. A slice is filtered out if the quotient between the total probability captured by the slice, and the upper bound on the total error in the slice, is less than $10^2$.

The filtered distribution is automatically assigned an appropriate name and written to the <code>distributions</code> directory. If this directory does not exist, it will be created.

### Mandatory command line arguments
Arguments <code>\<distribution\></code> where
- <code>\<distribution\></code> is the path to a distribution
