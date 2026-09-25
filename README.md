
<!-- README.md is generated from README.Rmd. Please edit that file -->

# betanorm

<!-- badges: start -->

<!-- badges: end -->

**betanorm** is a small toolkit of fast C++ array-normalization
utilities ([Rcpp](https://cran.r-project.org/package=Rcpp) /
[RcppArmadillo](https://cran.r-project.org/package=RcppArmadillo)):

- **Quantile normalization** - map each sample row onto a target
  distribution (no missing values; a focused port of
  `preprocessCore::normalize.quantiles.use.target()`).
- **BMIQ-style gold-standard calibration** - for DNA methylation beta
  matrices when a pipeline expects beta-mixture quantile calibration to
  a reference profile. The default path (`nL = 3`, `niter = 5`) is
  intentionally legacy-compatible, not a claim that five outer EM steps
  are statistical convergence.

## Installation

Install the development version from
[GitHub](https://github.com/hhp94/betanorm):

``` r
# install.packages("pak")
pak::pak("hhp94/betanorm")
```

## Quantile normalization

`quantile_norm()` maps each matrix row onto a target vector by rank:
equal length uses Bolstad-style indexing; unequal length uses linear
quantile interpolation. Inputs must be finite with no missing values.

``` r
library(betanorm)
```

``` r
set.seed(1)
X <- matrix(rnorm(5 * 100), nrow = 5)
target <- sort(rnorm(100))
Y <- quantile_norm(X, target)
dim(Y)
#> [1]   5 100
```

## Gold-standard calibration (BMIQ-style)

For methylation betas, `bmiq_calibration()` fits a two- or
three-component beta mixture per sample and maps components toward a
gold-standard profile. Use this when you need that procedure
(for example, clock-style workflows); prefer `quantile_norm()` when a plain
target distribution is enough.

Built-in example data: gold-standard annotation and a small beta matrix
(6 samples x ~21k probes).

``` r
data(horvath_goldstandard)
data(GPL21145_sample)

dim(GPL21145_sample)
#> [1]     6 21368
head(horvath_goldstandard$goldstandard2, 3)
#> [1] 0.8070658 0.7452180 0.0625570
```

### Legacy defaults (`nL = 3`, `niter = 5`)

Three-state BMIQ with a hard outer EM cap of **five** iterations is the
intentional default. That pair reproduces common legacy implementations
(Horvath-style / MEAT-style clocks and similar pipelines). **Do not
raise `niter` on the `nL = 3` path** unless you deliberately want
results that diverge from those implementations; `bmiq_calibration()`
warns if you do.

If you want the mixture EM to run longer or toward convergence, switch
to **`nL = 2`** (U/M only; no intermediate H step) and raise `niter` /
tighten `fit.policy` there instead.

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
#>   H skipped: 0 (U plus upper-M)
calibrated <- as.matrix(fit)
dim(calibrated)
#> [1]     6 21368
```

### Two-component path (`nL = 2`)

Use `nL = 2` for unmethylated / methylated only (`doH` defaults to
`FALSE`) and whenever you prefer longer EM rather than the legacy
five-iteration three-state cap.

``` r
fit2 <- bmiq_calibration(
  datM = GPL21145_sample,
  goldstandard.beta = horvath_goldstandard$goldstandard2,
  nL = 2L,
  niter = 25L,
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
  failed.sample = "NA", # or "original"
  h.policy = "optional", # keep U/M if H fails
  fit.policy = "usable" # accept finite (not fully converged) fits
)
```

## Main functions

| Function | Role |
|----|----|
| `quantile_norm()` | Quantile-normalize matrix rows to a target distribution |
| `bmiq_calibration()` | BMIQ-style calibration to a gold-standard beta profile |
| `horvath_goldstandard` | Gold-standard probe annotation + mean betas |
| `GPL21145_sample` | Example beta matrix for demos and tests |

## Licensing

betanorm is licensed under GPL (\>= 3) because parts of it are derived
from GPL-family code. Its quantile normalization is derived from
preprocessCore (LGPL \>= 2), and its BMIQ calibration and beta-mixture
EM are derived from Steve Horvath's BMIQcalibration, which builds on
Andrew Teschendorff's BMIQ and on the RPMM package (GPL \>= 2). The
original authors are credited as copyright holders in `DESCRIPTION`, and
`inst/COPYRIGHTS` records which file derives from which source.

## References

Bolstad BM, Irizarry RA, Astrand M, Speed TP (2003). A comparison of
normalization methods for high density oligonucleotide array data based
on variance and bias. *Bioinformatics* 19(2), 185-193.
<https://doi.org/10.1093/bioinformatics/19.2.185>

Bolstad BM. preprocessCore: A collection of pre-processing functions. R
package, Bioconductor. <https://doi.org/10.18129/B9.bioc.preprocessCore>

Teschendorff AE, Marabita F, Lechner M, Bartlett T, Tegner J,
Gomez-Cabrero D, Beck S (2013). A beta-mixture quantile normalization
method for correcting probe design bias in Illumina Infinium 450 k DNA
methylation data. *Bioinformatics* 29(2), 189-196.
<https://doi.org/10.1093/bioinformatics/bts680>

Horvath S (2013). DNA methylation age of human tissues and cell types.
*Genome Biology* 14(10), R115.
<https://doi.org/10.1186/gb-2013-14-10-r115>

Houseman EA, Koestler DC. RPMM: Recursively Partitioned Mixture Model. R
package, CRAN. <https://cran.r-project.org/package=RPMM>
