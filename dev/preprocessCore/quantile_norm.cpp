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

    // map an arithmetic-mean rank onto the sorted target.
    double target_at_rank(
        double rank,
        std::size_t n,
        const std::vector<double> &target)
    {
        const std::size_t m = target.size();

        // minimal extension for a sample with one observed value.
        if (n == 1)
        {
            return target[0];
        }
        // Bolstad's direct equal-length handling:
        // integer rank -> target value
        // half rank -> mean of adjacent target values
        if (n == m)
        {
            const double rank_floor = std::floor(rank);
            const std::size_t k =
                static_cast<std::size_t>(rank_floor);

            if (rank - rank_floor > 0.4)
                return 0.5 * (target[k - 1] + target[k]);

            return target[k - 1];
        }

        // unequal lengths: interpolate the corresponding target quantile.
        const double percentile =
            (rank - 1.0) / static_cast<double>(n - 1);

        double target_index =
            1.0 + static_cast<double>(m - 1) * percentile;

        const double target_floor =
            std::floor(target_index + 4.0 * DBL_EPSILON);

        double fraction = target_index - target_floor;

        if (std::fabs(fraction) <= 4.0 * DBL_EPSILON)
            fraction = 0.0;

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

} // namespace

// x is samples x variables.
// [[Rcpp::export]]
Rcpp::NumericMatrix qnorm_target_rows_cpp(
    const Rcpp::NumericMatrix &obj,
    const Rcpp::NumericVector &target)
{
    std::vector<double> sorted_target;
    sorted_target.reserve(target.size());

    for (R_xlen_t i = 0; i < target.size(); ++i)
    {
        if (!R_IsNA(target[i]))
        {
            sorted_target.push_back(target[i]);
        }
    }

    if (sorted_target.empty())
        Rcpp::stop("target has no non-missing values");

    std::sort(sorted_target.begin(), sorted_target.end());

    // variables x samples: each sample is now a contiguous column.
    Rcpp::NumericMatrix xt = Rcpp::transpose(x);

    const int n_variables = xt.nrow();
    const int n_samples = xt.ncol();

    std::vector<DataItem> items;
    items.reserve(n_variables);

    for (int sample = 0; sample < n_samples; ++sample)
    {
        items.clear();

        double *column =
            xt.begin() +
            static_cast<std::size_t>(sample) *
                static_cast<std::size_t>(n_variables);

        for (int variable = 0; variable < n_variables; ++variable)
        {
            const double value = column[variable];

            if (!R_IsNA(value))
            {
                DataItem item;
                item.value = value;
                item.index = variable;
                items.push_back(item);
            }
        }

        std::sort(items.begin(), items.end());

        const std::size_t n = items.size();
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

            // same arithmetic-mean rank as preprocessCore::get_ranks().
            const double rank =
                (static_cast<double>(first) +
                 static_cast<double>(last) + 2.0) /
                2.0;

            const double normalized =
                target_at_rank(rank, n, sorted_target);

            for (std::size_t k = first; k <= last; ++k)
                column[items[k].index] = normalized;

            first = last + 1;
        }
    }

    Rcpp::NumericMatrix result = Rcpp::transpose(xt);

    if (x.hasAttribute("dimnames"))
        result.attr("dimnames") = x.attr("dimnames");

    return result;
}
