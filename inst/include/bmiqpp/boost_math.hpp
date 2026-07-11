// -----------------------------------------------------------------------------
// Thread-safe Boost.Math wrappers (digamma / trigamma / pbeta / qbeta).
//
// Used today for unit tests against base R. Intended later for per-sample
// worker threads (see dev/parallel.md); production BMIQ still uses R-side
// special functions for quantile maps / modes.
//
// Pure C++ only — no R / Rcpp types here.
// -----------------------------------------------------------------------------
#ifndef BMIQPP_BOOST_MATH_HPP_
#define BMIQPP_BOOST_MATH_HPP_

#include <boost/math/special_functions/beta.hpp>
#include <boost/math/special_functions/digamma.hpp>
#include <boost/math/special_functions/trigamma.hpp>

namespace bmiqpp
{
namespace boost_math
{

// -----------------------------------------------------------------------------
// Digamma / trigamma  (R: digamma(), trigamma(); C API: R::digamma, R::trigamma)
// -----------------------------------------------------------------------------
inline double digamma(double x)
{
    return boost::math::digamma(x);
}

inline double trigamma(double x)
{
    return boost::math::trigamma(x);
}

// -----------------------------------------------------------------------------
// Regularized incomplete beta and its inverse.
//
// R correspondence (shape1 = a, shape2 = b):
//   pbeta(q, a, b, lower.tail = TRUE)  == ibeta(a, b, q)
//   pbeta(q, a, b, lower.tail = FALSE) == ibetac(a, b, q)
//   qbeta(p, a, b, lower.tail = TRUE)  == ibeta_inv(a, b, p)
//   qbeta(p, a, b, lower.tail = FALSE) == ibetac_inv(a, b, p)
// -----------------------------------------------------------------------------
inline double pbeta(
    double q,
    double shape1,
    double shape2,
    bool lower_tail = true)
{
    if (lower_tail)
    {
        return boost::math::ibeta(shape1, shape2, q);
    }
    return boost::math::ibetac(shape1, shape2, q);
}

inline double qbeta(
    double p,
    double shape1,
    double shape2,
    bool lower_tail = true)
{
    if (lower_tail)
    {
        return boost::math::ibeta_inv(shape1, shape2, p);
    }
    return boost::math::ibetac_inv(shape1, shape2, p);
}

} // namespace boost_math
} // namespace bmiqpp

#endif // BMIQPP_BOOST_MATH_HPP_
