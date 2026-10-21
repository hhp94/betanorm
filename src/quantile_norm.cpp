// Derived from preprocessCore (normalize.quantiles.use.target, get_ranks),
// copyright Ben Bolstad, LGPL (>= 2). Distributed here under GPL (>= 3).
// See inst/COPYRIGHTS.

#include <Rcpp.h>

#include <algorithm>
#include <cfloat>
#include <cmath>
#include <vector>

namespace
{

    // A strict weak ordering with NaN last: every NaN sorts after every
    // number and is equivalent to every other NaN (as -0 is to 0). A bare
    // `a < b` makes NaN equivalent to everything, which is not a strict weak
    // ordering and violates std::sort's precondition (UB: measured to
    // silently drop the NaN and return plausible wrong numbers, and to leave
    // a target unsorted). quantile_norm_fixed() and betanorm's
    // quantile_norm() reject NaN in obj and in target before this point;
    // this keeps the kernel defined if that ever moves. For NaN-free input
    // every comparison is unchanged, so both sorts are too.
    inline bool less_nan_last(double a, double b)
    {
        if (std::isnan(a))
        {
            return false;
        }
        if (std::isnan(b))
        {
            return true;
        }
        return a < b;
    }

    struct DataItem
    {
        double value;
        int index;

        bool operator<(const DataItem &other) const
        {
            return less_nan_last(value, other.value);
        }
    };

    // equal lengths (n == m): integer rank -> target[k-1]; half-rank (ties)
    // -> mean of adjacent target values. Bolstad using_target() path.
    inline double target_equal_rank(
        double rank,
        const double *target)
    {
        const double rank_floor = std::floor(rank);
        const std::size_t k = static_cast<std::size_t>(rank_floor);

        if (rank - rank_floor > 0.4)
            return 0.5 * (target[k - 1] + target[k]);

        return target[k - 1];
    }

    // 1-based target rank, clamped into [1, m]. The clamp is done in double:
    // casting a negative double to size_t is UB, not a wrap.
    inline std::size_t clamp_rank(double rank, std::size_t m)
    {
        double k = std::floor(rank);
        if (!(k >= 1.0))
        {
            k = 1.0;
        }
        if (k > static_cast<double>(m))
        {
            k = static_cast<double>(m);
        }
        return static_cast<std::size_t>(k);
    }

    // unequal lengths: linear interpolation of the matching target quantile.
    // inv_nm1 = 1/(n-1), m1 = m-1 are hoisted by the caller.
    inline double target_unequal_rank(
        double rank,
        double inv_nm1,
        std::size_t m,
        double m1,
        const double *target)
    {
        const double percentile = (rank - 1.0) * inv_nm1;

        double target_index = 1.0 + m1 * percentile;

        const double target_floor =
            std::floor(target_index + 4.0 * DBL_EPSILON);

        double fraction = target_index - target_floor;

        if (std::fabs(fraction) <= 4.0 * DBL_EPSILON)
        {
            fraction = 0.0;
        }
        // Both early returns clamp, matching the interpolating branch below.
        // With rank in [1, n] the index already lies in [1, m], so the clamp
        // never fires on valid input and the result is unchanged.
        if (fraction == 0.0)
        {
            return target[clamp_rank(target_floor + 0.5, m) - 1];
        }

        if (fraction == 1.0)
        {
            return target[clamp_rank(target_floor + 1.5, m) - 1];
        }

        const std::size_t k = static_cast<std::size_t>(
            std::floor(target_floor + 0.5));

        if (k < m && k > 0)
        {
            return (1.0 - fraction) * target[k - 1] +
                   fraction * target[k];
        }

        if (k >= m)
            return target[m - 1];

        return target[0];
    }

    // Ranks one sample and writes its normalized values. The sample is a row
    // of obj, its values `stride` apart (R is column-major): they are
    // gathered into items, sorted, and each tie group's value is scattered to
    // the same row of the result. map_rank is applied once per tie group.
    template <typename MapRank>
    void normalize_row(
        const double *row,
        double *result_row,
        std::size_t stride,
        std::size_t n,
        std::vector<DataItem> &items,
        MapRank map_rank)
    {
        for (std::size_t variable = 0; variable < n; ++variable)
        {
            items[variable].value = row[variable * stride];
            items[variable].index = static_cast<int>(variable);
        }

        std::sort(items.begin(), items.end());

        std::size_t first = 0;
        while (first < n)
        {
            std::size_t last = first;
            while (
                last + 1 < n &&
                items[last].value == items[last + 1].value)
            {
                ++last;
            }

            // arithmetic-mean rank, same as preprocessCore::get_ranks().
            const double rank =
                (static_cast<double>(first) +
                 static_cast<double>(last) + 2.0) /
                2.0;

            const double normalized = map_rank(rank);

            for (std::size_t k = first; k <= last; ++k)
                result_row[static_cast<std::size_t>(items[k].index) * stride] =
                    normalized;

            first = last + 1;
        }
    }

} // namespace

// [[Rcpp::export(rng = false)]]
Rcpp::NumericMatrix qnorm_target_rows_cpp(
    const Rcpp::NumericMatrix &obj,
    const Rcpp::NumericVector &target)
{
    if (target.size() == 0)
        Rcpp::stop("target has length 0");

    std::vector<double> sorted_target(target.begin(), target.end());
    std::sort(sorted_target.begin(), sorted_target.end(), less_nan_last);

    const int n_samples = obj.nrow();
    const int n_variables = obj.ncol();
    const std::size_t stride = static_cast<std::size_t>(n_samples);
    const std::size_t n = static_cast<std::size_t>(n_variables);
    const std::size_t m = sorted_target.size();
    const double *target_ptr = sorted_target.data();
    const double *obj_data = obj.begin();

    // Each sample is gathered from its row of obj and scattered to the same
    // row of the result, with no transposed copy, so the result is the only
    // allocation the size of obj. Every cell is written below, so it is not
    // zero-filled.
    Rcpp::NumericMatrix result(Rcpp::no_init(n_samples, n_variables));
    double *result_data = result.begin();

    // One variable: a sample has nothing to rank against, so every value
    // maps to the smallest target value (the unequal path would divide by
    // n - 1 = 0).
    if (n == 1)
    {
        const double v = target_ptr[0];
        for (std::size_t sample = 0; sample < stride; ++sample)
            result_data[sample] = v;
    }
    else if (n == m)
    {
        std::vector<DataItem> items(n);
        for (std::size_t sample = 0; sample < stride; ++sample)
        {
            normalize_row(
                obj_data + sample,
                result_data + sample,
                stride,
                n,
                items,
                [target_ptr](double rank)
                { return target_equal_rank(rank, target_ptr); });
        }
    }
    else
    {
        // hoisted unequal-path constants (were re-derived per rank group).
        // preprocessCore divides by n - 1 at each rank instead. The product
        // with the reciprocal differs from that quotient in its last bit at
        // some ranks, which moves about one value in seven on this path, by
        // at most 7e-13 in the cases measured. Kept, so betanorm's results
        // stay as they are (betanorm's dev/DECISIONS.md, 2026-10-21).
        const double inv_nm1 = 1.0 / static_cast<double>(n - 1);
        const double m1 = static_cast<double>(m - 1);

        std::vector<DataItem> items(n);
        for (std::size_t sample = 0; sample < stride; ++sample)
        {
            normalize_row(
                obj_data + sample,
                result_data + sample,
                stride,
                n,
                items,
                [inv_nm1, m, m1, target_ptr](double rank)
                {
                    return target_unequal_rank(
                        rank, inv_nm1, m, m1, target_ptr);
                });
        }
    }

    if (obj.hasAttribute("dimnames"))
    {
        result.attr("dimnames") = obj.attr("dimnames");
    }
    return result;
}
