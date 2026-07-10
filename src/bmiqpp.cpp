#include <RcppArmadillo.h>

#include <algorithm>
#include <cmath>
#include <limits>
#include <string>
#include <vector>

// [[Rcpp::depends(RcppArmadillo)]]

namespace
{

    // -----------------------------------------------------------------------------
    // Weighted sufficient statistics for a Beta distribution.
    //
    // The weighted mean and variance accumulator is numerically more stable than
    // calculating E[Y^2] - E[Y]^2. During EM, statistics for all nL components are
    // collected in one pass through the observations.
    // -----------------------------------------------------------------------------
    struct BetaStats
    {
        double weight = 0.0;
        double mean = 0.0;
        double m2 = 0.0;
        double sum_log_y = 0.0;
        double sum_log1m_y = 0.0;
        arma::uword positive_count = 0;

        void add(
            double y,
            double log_y,
            double log1m_y,
            double effective_weight)
        {
            if (!(effective_weight > 0.0))
            {
                return;
            }

            ++positive_count;

            const double new_weight = weight + effective_weight;
            const double delta = y - mean;

            mean += (effective_weight / new_weight) * delta;
            m2 += effective_weight * delta * (y - mean);
            weight = new_weight;

            sum_log_y += effective_weight * log_y;
            sum_log1m_y += effective_weight * log1m_y;
        }
    };

    // -----------------------------------------------------------------------------
    // Internal result from the Beta shape optimizer.
    // -----------------------------------------------------------------------------
    struct BetaFit
    {
        double a = 1.0;
        double b = 1.0;
        double loglik = NA_REAL;
        int iterations = 0;
        bool usable = false;
        bool converged = false;
        std::string status = "failed";
        std::string reason = "unknown failure";
    };

    // -----------------------------------------------------------------------------
    // Weighted Beta log-likelihood based on sufficient statistics.
    // -----------------------------------------------------------------------------
    inline double beta_loglik_stats(
        double a,
        double b,
        const BetaStats &stats)
    {
        return (a - 1.0) * stats.sum_log_y +
               (b - 1.0) * stats.sum_log1m_y +
               stats.weight *
                   (std::lgamma(a + b) -
                    std::lgamma(a) -
                    std::lgamma(b));
    }

    // -----------------------------------------------------------------------------
    // Newton optimization with step halving.
    //
    // Intentional divergence from legacy R:
    //
    //   * Legacy blc() used BFGS.
    //   * Legacy blc2() used Nelder-Mead.
    //   * This implementation uses the same Newton optimizer for both.
    //
    // A finite estimate at which Newton stalls, or reaches maxit, remains usable.
    // This is not a silent Beta(1,1) fallback: the best finite estimate found is
    // retained and its status is returned. R can apply a stricter policy by
    // requiring complete optimizer and EM convergence.
    //
    // Empty components, zero weighted variance, and invalid initial estimates are
    // hard failures because they do not provide a defensible finite Beta MLE.
    // -----------------------------------------------------------------------------
    BetaFit fit_beta_from_stats(
        const BetaStats &stats,
        int maxit,
        int max_halving,
        double score_tol,
        double min_shape,
        double armijo)
    {
        BetaFit out;

        if (!(stats.weight > 0.0) ||
            !std::isfinite(stats.weight) ||
            stats.positive_count <= 1)
        {
            out.reason = "zero or insufficient effective weight";
            return out;
        }

        const double p = stats.mean;
        const double variance = stats.m2 / stats.weight;

        if (!(p > 0.0 && p < 1.0) || !std::isfinite(p))
        {
            out.reason = "invalid weighted mean";
            return out;
        }

        if (!(variance > 0.0) || !std::isfinite(variance))
        {
            out.reason =
                "zero weighted variance; no finite stable Beta MLE";
            return out;
        }

        double concentration =
            p * (1.0 - p) / variance - 1.0;

        if (!std::isfinite(concentration))
        {
            concentration = 1.0;
        }
        else
        {
            concentration = std::max(1e-6, concentration);
        }

        double a = std::max(min_shape, p * concentration);
        double b = std::max(min_shape, (1.0 - p) * concentration);
        double loglik = beta_loglik_stats(a, b, stats);

        if (!std::isfinite(loglik))
        {
            out.reason = "non-finite initial log-likelihood";
            return out;
        }

        for (int iter = 0; iter < maxit; ++iter)
        {
            const double ab = a + b;
            const double digamma_ab = R::digamma(ab);

            // Score divided by total effective weight.
            const double score_a =
                stats.sum_log_y / stats.weight -
                R::digamma(a) +
                digamma_ab;

            const double score_b =
                stats.sum_log1m_y / stats.weight -
                R::digamma(b) +
                digamma_ab;

            // Scale-aware score in log-shape coordinates.
            const double scaled_score =
                std::max(
                    std::abs(a * score_a),
                    std::abs(b * score_b));

            if (scaled_score <= score_tol)
            {
                out.a = a;
                out.b = b;
                out.loglik = loglik;
                out.iterations = iter;
                out.usable = true;
                out.converged = true;
                out.status = "converged";
                out.reason = "score tolerance reached";
                return out;
            }

            const double trigamma_ab = R::trigamma(ab);

            // Per-unit-effective-weight information matrix.
            const double i11 = R::trigamma(a) - trigamma_ab;
            const double i22 = R::trigamma(b) - trigamma_ab;
            const double i12 = -trigamma_ab;
            const double determinant = i11 * i22 - i12 * i12;

            if (!(determinant > 0.0) ||
                !std::isfinite(determinant))
            {
                out.a = a;
                out.b = b;
                out.loglik = loglik;
                out.iterations = iter;
                out.usable = true;
                out.converged = false;
                out.status = "stalled";
                out.reason =
                    "singular or ill-conditioned information matrix";
                return out;
            }

            // Solve I * delta = score.
            const double delta_a =
                (i22 * score_a - i12 * score_b) / determinant;

            const double delta_b =
                (i11 * score_b - i12 * score_a) / determinant;

            // Total directional derivative because loglik is total rather than
            // average weighted log-likelihood.
            const double directional_derivative =
                stats.weight *
                (score_a * delta_a + score_b * delta_b);

            if (!(directional_derivative > 0.0) ||
                !std::isfinite(directional_derivative))
            {
                out.a = a;
                out.b = b;
                out.loglik = loglik;
                out.iterations = iter;
                out.usable = true;
                out.converged = false;
                out.status = "stalled";
                out.reason =
                    "Newton direction is not an ascent direction";
                return out;
            }

            double step = 1.0;
            double new_a = a;
            double new_b = b;
            double new_loglik = loglik;
            bool accepted = false;

            for (int h = 0; h <= max_halving; ++h)
            {
                new_a = a + step * delta_a;
                new_b = b + step * delta_b;

                if (new_a > min_shape &&
                    new_b > min_shape &&
                    std::isfinite(new_a) &&
                    std::isfinite(new_b))
                {
                    new_loglik =
                        beta_loglik_stats(new_a, new_b, stats);

                    if (std::isfinite(new_loglik) &&
                        new_loglik >=
                            loglik +
                                armijo *
                                    step *
                                    directional_derivative)
                    {
                        accepted = true;
                        break;
                    }
                }

                step *= 0.5;
            }

            if (!accepted)
            {
                out.a = a;
                out.b = b;
                out.loglik = loglik;
                out.iterations = iter;
                out.usable = true;
                out.converged = false;
                out.status = "stalled";
                out.reason = "step halving failed";
                return out;
            }

            a = new_a;
            b = new_b;
            loglik = new_loglik;
        }

        out.a = a;
        out.b = b;
        out.loglik = loglik;
        out.iterations = maxit;
        out.usable = true;
        out.converged = false;
        out.status = "max_iter";
        out.reason = "maximum iterations reached";

        return out;
    }

    Rcpp::List beta_fit_to_list(const BetaFit &fit)
    {
        return Rcpp::List::create(
            Rcpp::_["par"] =
                Rcpp::NumericVector::create(fit.a, fit.b),
            Rcpp::_["status"] = fit.status,
            Rcpp::_["iterations"] = fit.iterations,
            Rcpp::_["logLik"] = fit.loglik,
            Rcpp::_["reason"] = fit.reason,
            Rcpp::_["usable"] = fit.usable,
            Rcpp::_["converged"] = fit.converged);
    }

} // anonymous namespace

// -----------------------------------------------------------------------------
// Validate a complete numeric matrix.
//
// require_open = false:
//   Check raw input values in [0, 1].
//
// require_open = true:
//   Check endpoint-clipped values in (0, 1).
//
// This scanner only verifies the numeric domain. It cannot guarantee that
// mixture classes, mode windows, MAP classes, or normalization anchors are
// nonempty.
// -----------------------------------------------------------------------------
// [[Rcpp::export]]
void scan_finite_unit_interval_cpp(
    const arma::mat &x,
    std::string name = "x",
    bool require_open = false)
{
    for (arma::uword i = 0; i < x.n_elem; ++i)
    {
        const double value = x[i];

        if (!std::isfinite(value))
        {
            Rcpp::stop(
                name +
                " must contain only finite values "
                "(no NA, NaN, or +/-Inf).");
        }

        if (require_open)
        {
            if (!(value > 0.0 && value < 1.0))
            {
                Rcpp::stop(
                    name +
                    " must lie strictly inside (0, 1) after clipping.");
            }
        }
        else if (value < 0.0 || value > 1.0)
        {
            Rcpp::stop(
                name + " must have all values in [0, 1].");
        }
    }
}

// -----------------------------------------------------------------------------
// Standalone Beta optimizer for testing and diagnostics.
//
// beta_mixture_em_cpp() does not call this wrapper because it collects
// sufficient statistics for all components simultaneously.
// -----------------------------------------------------------------------------
// [[Rcpp::export]]
Rcpp::List beta_est_newton_cpp(
    const arma::vec &y,
    const arma::vec &responsibility,
    const arma::vec &observation_weight,
    int maxit = 50,
    int max_halving = 30,
    double score_tol = 1e-10,
    double min_shape = 1e-10,
    double armijo = 1e-4)
{
    const arma::uword n = y.n_elem;

    if (responsibility.n_elem != n ||
        observation_weight.n_elem != n)
    {
        Rcpp::stop(
            "y, responsibility, and observation_weight "
            "must have equal lengths");
    }

    if (maxit < 1 || max_halving < 0)
    {
        Rcpp::stop(
            "maxit must be positive and max_halving non-negative");
    }

    if (!(score_tol > 0.0) ||
        !(min_shape > 0.0) ||
        !(armijo > 0.0 && armijo < 1.0))
    {
        Rcpp::stop(
            "Invalid score_tol, min_shape, or armijo setting");
    }

    BetaStats stats;

    for (arma::uword i = 0; i < n; ++i)
    {
        const double yi = y[i];
        const double ri = responsibility[i];
        const double wi = observation_weight[i];

        if (!std::isfinite(yi) ||
            !(yi > 0.0 && yi < 1.0))
        {
            Rcpp::stop(
                "All y values must be finite and strictly inside (0, 1)");
        }

        if (!std::isfinite(ri) || ri < 0.0)
        {
            Rcpp::stop(
                "responsibility must be finite and non-negative");
        }

        if (!std::isfinite(wi) || wi < 0.0)
        {
            Rcpp::stop(
                "observation_weight must be finite and non-negative");
        }

        stats.add(
            yi,
            std::log(yi),
            std::log1p(-yi),
            ri * wi);
    }

    return beta_fit_to_list(
        fit_beta_from_stats(
            stats,
            maxit,
            max_halving,
            score_tol,
            min_shape,
            armijo));
}

// -----------------------------------------------------------------------------
// Complete Beta-mixture EM.
//
// Intentional divergences from legacy blc()/blc2():
//
//   * Newton optimization is used for all components.
//   * No silent Beta(1,1) fallback.
//   * Log-sum-exp is used during the E-step.
//   * All component sufficient statistics are collected in one pass.
//   * An unusable component is a hard error.
//   * A finite stalled or max-iteration Beta estimate is returned with status.
//   * nL is explicit and checked against the responsibility matrix.
//
// Detailed criterion history is returned only when debug = true.
// -----------------------------------------------------------------------------
// [[Rcpp::export]]
Rcpp::List beta_mixture_em_cpp(
    const arma::vec &y,
    const arma::mat &initial_responsibility,
    int nL = 3,
    Rcpp::Nullable<Rcpp::NumericVector> weights = R_NilValue,
    int maxiter = 25,
    double tol = 1e-6,
    int beta_maxit = 50,
    int beta_max_halving = 30,
    double beta_score_tol = 1e-10,
    double min_shape = 1e-10,
    double armijo = 1e-4,
    bool debug = false)
{
    const arma::uword n = y.n_elem;

    if (nL < 2)
    {
        Rcpp::stop("nL must be at least 2");
    }

    const arma::uword K =
        static_cast<arma::uword>(nL);

    if (n == 0)
    {
        Rcpp::stop("y must be non-empty");
    }

    if (initial_responsibility.n_rows != n)
    {
        Rcpp::stop(
            "initial_responsibility must have one row "
            "for every element of y");
    }

    if (initial_responsibility.n_cols != K)
    {
        Rcpp::stop(
            "ncol(initial_responsibility) must equal nL");
    }

    if (maxiter < 1 ||
        beta_maxit < 1 ||
        beta_max_halving < 0)
    {
        Rcpp::stop(
            "maxiter and beta_maxit must be positive; "
            "beta_max_halving must be non-negative");
    }

    if (!(tol > 0.0) ||
        !(beta_score_tol > 0.0) ||
        !(min_shape > 0.0) ||
        !(armijo > 0.0 && armijo < 1.0))
    {
        Rcpp::stop(
            "Invalid tolerance, minimum shape, or Armijo setting");
    }

    arma::vec log_y(n);
    arma::vec log1m_y(n);

    for (arma::uword i = 0; i < n; ++i)
    {
        const double yi = y[i];

        if (!std::isfinite(yi) ||
            !(yi > 0.0 && yi < 1.0))
        {
            Rcpp::stop(
                "All y values must be finite and strictly inside (0, 1)");
        }

        log_y[i] = std::log(yi);
        log1m_y[i] = std::log1p(-yi);
    }

    arma::vec observation_weight(n, arma::fill::ones);

    if (weights.isNotNull())
    {
        Rcpp::NumericVector r_weights(weights.get());

        if (static_cast<arma::uword>(r_weights.size()) != n)
        {
            Rcpp::stop("weights and y must have equal lengths");
        }

        for (arma::uword i = 0; i < n; ++i)
        {
            observation_weight[i] = r_weights[i];
        }
    }

    double total_observation_weight = 0.0;

    for (arma::uword i = 0; i < n; ++i)
    {
        const double wi = observation_weight[i];

        if (!std::isfinite(wi) || wi < 0.0)
        {
            Rcpp::stop(
                "weights must be finite and non-negative");
        }

        total_observation_weight += wi;
    }

    if (!(total_observation_weight > 0.0))
    {
        Rcpp::stop(
            "At least one observation weight must be positive");
    }

    arma::mat responsibility = initial_responsibility;

    for (arma::uword i = 0; i < n; ++i)
    {
        double row_sum = 0.0;

        for (arma::uword k = 0; k < K; ++k)
        {
            const double value = responsibility(i, k);

            if (!std::isfinite(value) || value < 0.0)
            {
                Rcpp::stop(
                    "initial_responsibility must be finite "
                    "and non-negative");
            }

            row_sum += value;
        }

        if (!(row_sum > 0.0))
        {
            Rcpp::stop(
                "Every initial responsibility row "
                "must have positive sum");
        }

        responsibility.row(i) /= row_sum;
    }

    arma::vec a(K, arma::fill::ones);
    arma::vec b(K, arma::fill::ones);
    arma::vec eta(K, arma::fill::zeros);
    arma::vec mu(K);
    arma::vec old_mu(K);

    mu.fill(arma::datum::inf);
    old_mu.fill(arma::datum::inf);

    arma::vec log_component(K);
    arma::vec log_norm(K);

    Rcpp::CharacterVector fit_status(nL);
    Rcpp::CharacterVector fit_reason(nL);

    std::vector<double> criterion_trace;

    double loglikelihood = NA_REAL;
    bool converged = false;
    int completed_iterations = 0;

    for (int iter = 0; iter < maxiter; ++iter)
    {
        Rcpp::checkUserInterrupt();

        old_mu = mu;
        eta.zeros();

        for (arma::uword i = 0; i < n; ++i)
        {
            const double wi = observation_weight[i];

            for (arma::uword k = 0; k < K; ++k)
            {
                eta[k] += wi * responsibility(i, k);
            }
        }

        eta /= total_observation_weight;

        std::vector<BetaStats> stats(K);

        for (arma::uword i = 0; i < n; ++i)
        {
            const double wi = observation_weight[i];

            for (arma::uword k = 0; k < K; ++k)
            {
                stats[k].add(
                    y[i],
                    log_y[i],
                    log1m_y[i],
                    wi * responsibility(i, k));
            }
        }

        for (arma::uword k = 0; k < K; ++k)
        {
            if (!(eta[k] > 0.0) ||
                !std::isfinite(eta[k]))
            {
                Rcpp::stop(
                    "Mixture component " +
                    std::to_string(k + 1) +
                    " has zero or invalid mixture weight");
            }

            const BetaFit fit = fit_beta_from_stats(
                stats[k],
                beta_maxit,
                beta_max_halving,
                beta_score_tol,
                min_shape,
                armijo);

            if (!fit.usable)
            {
                Rcpp::stop(
                    "Beta fit failed for mixture component " +
                    std::to_string(k + 1) +
                    ": " +
                    fit.reason);
            }

            a[k] = fit.a;
            b[k] = fit.b;
            mu[k] = fit.a / (fit.a + fit.b);

            fit_status[k] = fit.status;
            fit_reason[k] = fit.reason;

            log_norm[k] =
                std::lgamma(a[k] + b[k]) -
                std::lgamma(a[k]) -
                std::lgamma(b[k]);
        }

        loglikelihood = 0.0;

        for (arma::uword i = 0; i < n; ++i)
        {
            double maximum =
                -std::numeric_limits<double>::infinity();

            for (arma::uword k = 0; k < K; ++k)
            {
                log_component[k] =
                    std::log(eta[k]) +
                    log_norm[k] +
                    (a[k] - 1.0) * log_y[i] +
                    (b[k] - 1.0) * log1m_y[i];

                maximum =
                    std::max(maximum, log_component[k]);
            }

            double sum_exp = 0.0;

            for (arma::uword k = 0; k < K; ++k)
            {
                sum_exp +=
                    std::exp(log_component[k] - maximum);
            }

            if (!(sum_exp > 0.0) ||
                !std::isfinite(sum_exp))
            {
                Rcpp::stop(
                    "Non-finite mixture likelihood in E-step");
            }

            const double log_mixture =
                maximum + std::log(sum_exp);

            loglikelihood +=
                observation_weight[i] * log_mixture;

            for (arma::uword k = 0; k < K; ++k)
            {
                responsibility(i, k) =
                    std::exp(
                        log_component[k] - log_mixture);
            }
        }

        const double criterion =
            arma::max(arma::abs(mu - old_mu));

        if (debug)
        {
            criterion_trace.push_back(criterion);
        }

        completed_iterations = iter + 1;

        if (criterion < tol)
        {
            converged = true;
            break;
        }
    }

    arma::mat a_matrix(K, 1);
    arma::mat b_matrix(K, 1);
    arma::mat mu_matrix(K, 1);

    a_matrix.col(0) = a;
    b_matrix.col(0) = b;
    mu_matrix.col(0) = mu;

    Rcpp::List result = Rcpp::List::create(
        Rcpp::_["a"] = a_matrix,
        Rcpp::_["b"] = b_matrix,
        Rcpp::_["eta"] = eta,
        Rcpp::_["mu"] = mu_matrix,
        Rcpp::_["w"] = responsibility,
        Rcpp::_["llike"] = loglikelihood,
        Rcpp::_["iterations"] = completed_iterations,
        Rcpp::_["converged"] = converged,
        Rcpp::_["fit_status"] = fit_status,
        Rcpp::_["fit_reason"] = fit_reason,
        Rcpp::_["nL"] = nL);

    if (debug)
    {
        result["criterion"] =
            Rcpp::wrap(criterion_trace);
    }

    return result;
}