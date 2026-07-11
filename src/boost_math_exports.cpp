// -----------------------------------------------------------------------------
// Rcpp exports of Boost.Math wrappers for numerical agreement tests vs base R.
// Production BMIQ path still uses R::digamma / R::trigamma and R pbeta/qbeta;
// Equality-test entry points vs base R (also a dry run for any future
// worker path that cannot call R special functions; see dev/parallel.md).
// -----------------------------------------------------------------------------

// [[Rcpp::depends(BH)]]
// [[Rcpp::depends(Rcpp)]]

#include <Rcpp.h>

#include <bmiqpp/boost_math.hpp>

// [[Rcpp::export]]
Rcpp::NumericVector boost_digamma_cpp(Rcpp::NumericVector x)
{
    const R_xlen_t n = x.size();
    Rcpp::NumericVector out(n);
    for (R_xlen_t i = 0; i < n; ++i)
    {
        if (Rcpp::NumericVector::is_na(x[i]))
        {
            out[i] = NA_REAL;
            continue;
        }
        out[i] = bmiqpp::boost_math::digamma(x[i]);
    }
    return out;
}

// [[Rcpp::export]]
Rcpp::NumericVector boost_trigamma_cpp(Rcpp::NumericVector x)
{
    const R_xlen_t n = x.size();
    Rcpp::NumericVector out(n);
    for (R_xlen_t i = 0; i < n; ++i)
    {
        if (Rcpp::NumericVector::is_na(x[i]))
        {
            out[i] = NA_REAL;
            continue;
        }
        out[i] = bmiqpp::boost_math::trigamma(x[i]);
    }
    return out;
}

// [[Rcpp::export]]
Rcpp::NumericVector boost_pbeta_cpp(
    Rcpp::NumericVector q,
    double shape1,
    double shape2,
    bool lower_tail = true)
{
    const R_xlen_t n = q.size();
    Rcpp::NumericVector out(n);
    for (R_xlen_t i = 0; i < n; ++i)
    {
        if (Rcpp::NumericVector::is_na(q[i]))
        {
            out[i] = NA_REAL;
            continue;
        }
        out[i] = bmiqpp::boost_math::pbeta(
            q[i], shape1, shape2, lower_tail);
    }
    return out;
}

// [[Rcpp::export]]
Rcpp::NumericVector boost_qbeta_cpp(
    Rcpp::NumericVector p,
    double shape1,
    double shape2,
    bool lower_tail = true)
{
    const R_xlen_t n = p.size();
    Rcpp::NumericVector out(n);
    for (R_xlen_t i = 0; i < n; ++i)
    {
        if (Rcpp::NumericVector::is_na(p[i]))
        {
            out[i] = NA_REAL;
            continue;
        }
        out[i] = bmiqpp::boost_math::qbeta(
            p[i], shape1, shape2, lower_tail);
    }
    return out;
}
