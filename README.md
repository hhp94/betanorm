
<!-- README.md is generated from README.Rmd. Please edit that file -->

# bmiqpp

<!-- badges: start -->

<!-- badges: end -->

**bmiqpp** provides fast beta-mixture quantile (BMIQ) calibration of DNA
methylation beta values against a gold-standard profile. Mixture
component parameters are estimated with a Newton–Raphson procedure and
step-halving, implemented in C++ via
[Rcpp](https://cran.r-project.org/package=Rcpp) and
[RcppArmadillo](https://cran.r-project.org/package=RcppArmadillo).

The package also includes a small quantile-normalization helper that
maps each sample row onto a target distribution (no missing values).

## Installation

Install the development version from
[GitHub](https://github.com/hhp94/bmiqpp):

``` r
# install.packages("pak")
pak::pak("hhp94/bmiqpp")
```

## Quick start

Load the package and the built-in example data: a gold-standard
annotation table and a small beta matrix (6 samples × ~21k probes).

``` r
library(bmiqpp)

data(horvath_goldstandard)
data(GPL21145_sample)

dim(GPL21145_sample)
#> [1]     6 21368
head(horvath_goldstandard$goldstandard2, 3)
#> [1] 0.8070658 0.7452180 0.0625570
```

Calibrate samples to the gold-standard mean betas. By default this uses
a **three-component** mixture (unmethylated / intermediate / methylated)
with outer EM iterations `niter = 5`.

``` r
fit <- bmiq_calibration(
  datM = GPL21145_sample,
  goldstandard.beta = horvath_goldstandard$goldstandard2,
  verbose = FALSE
)

fit
#> BMIQ calibration result
#>   Samples:   6
#>   Succeeded: 6
#>   Failed:    0
#>   H skipped: 0 (U/M-only)
calibrated <- as.matrix(fit)
dim(calibrated)
#> [1]     6 21368
```

### Two-component mixtures

Set `nL = 2` for unmethylated / methylated only (no intermediate H
step). In that mode `doH` defaults to `FALSE`.

``` r
fit2 <- bmiq_calibration(
  datM = GPL21145_sample,
  goldstandard.beta = horvath_goldstandard$goldstandard2,
  nL = 2L,
  verbose = FALSE
)
fit2
#> BMIQ calibration result
#>   Samples:   6
#>   Succeeded: 6
#>   Failed:    0
```

### Sample-level error handling

By default a failed sample stops the call. To keep going and record
failures:

``` r
fit <- bmiq_calibration(
  datM = datM,
  goldstandard.beta = gold,
  on.sample.error = "continue",
  failed.sample = "NA",      # or "original"
  h.policy = "optional",     # keep U/M if H fails
  fit.policy = "usable"      # accept finite (not fully converged) fits
)
```

## Quantile normalization

`quantile_norm()` is a minimal, no-missing-data port of
`preprocessCore::normalize.quantiles.use.target()`. Each row is ranked
and mapped onto a target vector.

``` r
set.seed(1)
X <- matrix(rnorm(5 * 100), nrow = 5)
target <- sort(rnorm(100))
Y <- quantile_norm(X, target)
dim(Y)
#> [1]   5 100
```

## Main functions

| Function | Role |
|----|----|
| `bmiq_calibration()` | BMIQ calibration to a gold-standard beta profile |
| `quantile_norm()` | Quantile-normalize matrix rows to a target distribution |
| `horvath_goldstandard` | Gold-standard probe annotation + mean betas |
| `GPL21145_sample` | Example beta matrix for demos and tests |

## References

Teschendorff AE, Marabita F, Lechner M, Bartlett T, Tegner J,
Gomez-Cabrero D, Beck S (2013). A beta-mixture quantile normalization
method for correcting probe design bias in Illumina Infinium 450 k DNA
methylation data. *Bioinformatics* 29(2), 189–196.
<https://doi.org/10.1093/bioinformatics/bts680>

Bolstad BM, Irizarry RA, Åstrand M, Speed TP (2003). A comparison of
normalization methods for high density oligonucleotide array data based
on variance and bias. *Bioinformatics* 19(2), 185–193.
<https://doi.org/10.1093/bioinformatics/19.2.185>
