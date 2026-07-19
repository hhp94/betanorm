#include <Rcpp.h>

#include <algorithm>
#include <cfloat>
#include <cmath>
#include <vector>

namespace
{

    struct DataItem
    {
        double value;
        int index;

        bool operator<(const DataItem &other) const
        {
            return value < other.value;
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
        if (fraction == 0.0)
        {
            const std::size_t k = static_cast<std::size_t>(
                std::floor(target_floor + 0.5));
            return target[k - 1];
        }

        if (fraction == 1.0)
        {
            const std::size_t k = static_cast<std::size_t>(
                std::floor(target_floor + 1.5));
            return target[k - 1];
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

    // fill items from a contiguous sample column, sort, map ranks -> target.
    // map_rank is applied once per tie-group.
    template <typename MapRank>
    void normalize_column(
        double *column,
        std::size_t n,
        std::vector<DataItem> &items,
        MapRank map_rank)
    {
        for (std::size_t variable = 0; variable < n; ++variable)
        {
            items[variable].value = column[variable];
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
                column[items[k].index] = normalized;

            first = last + 1;
        }
    }

} // namespace

// [[Rcpp::export]]
Rcpp::NumericMatrix qnorm_target_rows_cpp(
    const Rcpp::NumericMatrix &obj,
    const Rcpp::NumericVector &target)
{
    if (target.size() == 0)
        Rcpp::stop("target has length 0");

    std::vector<double> sorted_target(target.begin(), target.end());
    std::sort(sorted_target.begin(), sorted_target.end());

    // variables x samples: each sample becomes a contiguous column so that
    // ranking is a linear scan over packed memory (R is column-major).
    Rcpp::NumericMatrix xt = Rcpp::transpose(obj);

    const int n_variables = xt.nrow();
    const int n_samples = xt.ncol();
    const std::size_t n = static_cast<std::size_t>(n_variables);
    const std::size_t m = sorted_target.size();
    const double *target_ptr = sorted_target.data();
    double *xt_data = xt.begin();

    if (n == 1)
    {
        const double v = target_ptr[0];
        for (int sample = 0; sample < n_samples; ++sample)
            xt_data[static_cast<std::size_t>(sample)] = v;
    }
    else if (n == m)
    {
        std::vector<DataItem> items(n);
        for (int sample = 0; sample < n_samples; ++sample)
        {
            double *column =
                xt_data + static_cast<std::size_t>(sample) * n;
            normalize_column(column, n, items, [target_ptr](double rank)
                             { return target_equal_rank(rank, target_ptr); });
        }
    }
    else
    {
        // hoisted unequal-path constants (were re-derived per rank group).
        const double inv_nm1 = 1.0 / static_cast<double>(n - 1);
        const double m1 = static_cast<double>(m - 1);

        std::vector<DataItem> items(n);
        for (int sample = 0; sample < n_samples; ++sample)
        {
            double *column =
                xt_data + static_cast<std::size_t>(sample) * n;
            normalize_column(
                column,
                n,
                items,
                [inv_nm1, m, m1, target_ptr](double rank)
                {
                    return target_unequal_rank(
                        rank, inv_nm1, m, m1, target_ptr);
                });
        }
    }

    Rcpp::NumericMatrix result = Rcpp::transpose(xt);

    if (obj.hasAttribute("dimnames"))
    {
        result.attr("dimnames") = obj.attr("dimnames");
    }
    return result;
}
