# Normalised Mutual Information for a pair of categorical variables

Computes the symmetric Normalised Mutual Information (NMI) between two
categorical variables as an association measure suitable for use as an
edge weight in an undirected `catgraph`.

## Usage

``` r
nmi_assoc(x, y, adjusted = FALSE)
```

## Arguments

- x:

  A factor, character, or logical vector.

- y:

  A factor, character, or logical vector of the same length as `x`.

- adjusted:

  Logical. If `TRUE`, returns Adjusted Mutual Information (AMI), which
  corrects for chance. Default `FALSE`.

## Value

A named list with:

- `effect_size`:

  NMI or AMI value, numeric in \\\[0, 1\]\\ (AMI may be slightly
  negative due to numerical noise; clamped to 0).

- `metric`:

  Character: `"nmi"` or `"ami"`.

- `adjusted`:

  Logical indicating whether AMI was computed.

- `type`:

  Contingency-table type: `"2x2"` or `"RxC"` or `"degenerate"`.

- `statistic`:

  The chi-square statistic (from the underlying
  [`compute_assoc`](https://atinakosta.github.io/catgraph/reference/compute_assoc.md)
  call), for reference.

- `p_value`:

  Chi-square p-value, for reference.

- `df`:

  Degrees of freedom.

- `n`:

  Number of pairwise-complete observations.

## Details

NMI is defined as:

\$\$NMI(X, Y) = \frac{I(X; Y)}{\sqrt{H(X) \cdot H(Y)}}\$\$

where \\I(X;Y) = H(X) + H(Y) - H(X,Y)\\ is the mutual information and
\\H(\cdot)\\ denotes Shannon entropy (in nats, using natural logarithm).
The geometric-mean normalisation ensures the result lies in \\\[0, 1\]\\
and is symmetric: \\NMI(X, Y) = NMI(Y, X)\\.

**Small-sample behaviour.** Raw mutual information is upward-biased when
contingency tables are sparse. `nmi_assoc()` offers two corrections via
the `adjusted` argument:

- `adjusted = FALSE` (default): returns plain NMI computed from observed
  frequencies. Suitable when tables are well-populated.

- `adjusted = TRUE`: returns Adjusted Mutual Information (AMI), which
  subtracts the expected MI under random permutation and rescales, so
  that the expected value under independence is 0. Recommended for
  sparse tables or variables with many categories.

**Relationship to Cramer's V.** Both NMI and Cramer's V are symmetric,
bounded in \\\[0, 1\]\\, and equal 0 under independence. They measure
different things: Cramer's V quantifies departure from independence
relative to a chi-square null; NMI quantifies the fraction of
uncertainty in one variable explained by the other. The two will
generally agree on strong vs. weak associations but can differ on
variables with unequal numbers of categories or highly skewed marginals.

Entropy is computed using natural logarithms (nats). The choice of base
does not affect NMI because it cancels in the ratio. Zero-count cells
contribute 0 to entropy sums (the standard convention: 0 \* log(0) = 0).

The AMI formula follows Vinh, Epps & Bailey (2010), adapted for
two-variable contingency tables rather than clustering partitions.
Expected MI is computed analytically from the marginal counts using the
hypergeometric model.

## References

Cover, T. M., & Thomas, J. A. (2006). *Elements of Information Theory*
(2nd ed.). Wiley.
[doi:10.1002/047174882X](https://doi.org/10.1002/047174882X)

Vinh, N. X., Epps, J., & Bailey, J. (2010). Information theoretic
measures for clusterings comparison: Variants, properties, normalisation
and correction for chance. *Journal of Machine Learning Research*, 11,
2837–2854. <https://jmlr.org/papers/v11/vinh10a.html>

## See also

[`effect_size`](https://atinakosta.github.io/catgraph/reference/effect_size.md),
[`compute_assoc`](https://atinakosta.github.io/catgraph/reference/compute_assoc.md)

## Examples

``` r
set.seed(1)
x <- sample(c("A", "B", "C"), 120, replace = TRUE)
y <- sample(c("yes", "no"),    120, replace = TRUE)
nmi_assoc(x, y)
#> $effect_size
#> [1] 0.01694453
#> 
#> $metric
#> [1] "nmi"
#> 
#> $adjusted
#> [1] FALSE
#> 
#> $type
#> [1] "RxC"
#> 
#> $statistic
#> [1] 3.507813
#> 
#> $p_value
#> [1] 0.1730965
#> 
#> $df
#> [1] 2
#> 
#> $n
#> [1] 120
#> 
nmi_assoc(x, y, adjusted = TRUE)
#> $effect_size
#> [1] 0.007281291
#> 
#> $metric
#> [1] "ami"
#> 
#> $adjusted
#> [1] TRUE
#> 
#> $type
#> [1] "RxC"
#> 
#> $statistic
#> [1] 3.507813
#> 
#> $p_value
#> [1] 0.1730965
#> 
#> $df
#> [1] 2
#> 
#> $n
#> [1] 120
#> 
```
