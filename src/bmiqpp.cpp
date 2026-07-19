#include <RcppArmadillo.h>

#include <algorithm>
#include <cmath>
#include <limits>
#include <string>
#include <vector>

// [[Rcpp::depends(RcppArmadillo)]]

namespace
{

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

    inline void set_beta_fit(
        BetaFit &out,
        double a,
        double b,
        double loglik,
        int iterations,
        bool usable,
        bool converged,
        const char *status,
        const char *reason)
    {
        out.a = a;
        out.b = b;
        out.loglik = loglik;
        out.iterations = iterations;
        out.usable = usable;
        out.converged = converged;
        out.status = status;
        out.reason = reason;
    }

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

        double concentration = p * (1.0 - p) / variance - 1.0;
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

        // Sufficient statistics are fixed for this MLE; only digamma/trigamma
        // of the shapes change across Newton steps.
        const double mean_log_y = stats.sum_log_y / stats.weight;
        const double mean_log1m_y = stats.sum_log1m_y / stats.weight;

        for (int iter = 0; iter < maxit; ++iter)
        {
            const double ab = a + b;
            const double digamma_ab = R::digamma(ab);

            const double score_a =
                mean_log_y - R::digamma(a) + digamma_ab;
            const double score_b =
                mean_log1m_y - R::digamma(b) + digamma_ab;

            const double scaled_score =
                std::max(std::abs(a * score_a), std::abs(b * score_b));

            if (scaled_score <= score_tol)
            {
                set_beta_fit(
                    out, a, b, loglik, iter, true, true,
                    "converged", "score tolerance reached");
                return out;
            }

            const double trigamma_ab = R::trigamma(ab);
            const double i11 = R::trigamma(a) - trigamma_ab;
            const double i22 = R::trigamma(b) - trigamma_ab;
            const double i12 = -trigamma_ab;
            const double determinant = i11 * i22 - i12 * i12;

            if (!(determinant > 0.0) || !std::isfinite(determinant))
            {
                set_beta_fit(
                    out, a, b, loglik, iter, true, false, "stalled",
                    "singular or ill-conditioned information matrix");
                return out;
            }

            const double delta_a =
                (i22 * score_a - i12 * score_b) / determinant;
            const double delta_b =
                (i11 * score_b - i12 * score_a) / determinant;
            const double directional_derivative =
                stats.weight *
                (score_a * delta_a + score_b * delta_b);

            if (!(directional_derivative > 0.0) ||
                !std::isfinite(directional_derivative))
            {
                set_beta_fit(
                    out, a, b, loglik, iter, true, false, "stalled",
                    "Newton direction is not an ascent direction");
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
                    new_loglik = beta_loglik_stats(new_a, new_b, stats);
                    if (std::isfinite(new_loglik) &&
                        new_loglik >=
                            loglik +
                                armijo * step * directional_derivative)
                    {
                        accepted = true;
                        break;
                    }
                }
                step *= 0.5;
            }

            if (!accepted)
            {
                set_beta_fit(
                    out, a, b, loglik, iter, true, false, "stalled",
                    "step halving failed");
                return out;
            }

            a = new_a;
            b = new_b;
            loglik = new_loglik;
        }

        set_beta_fit(
            out, a, b, loglik, maxit, true, false, "max_iter",
            "maximum iterations reached");
        return out;
    }

} // anonymous namespace

// [[Rcpp::export]]
void scan_finite_unit_interval_cpp(
    const arma::mat &x,
    std::string name = "x",
    bool require_open = false)
{
    const double *data = x.memptr();
    const arma::uword n = x.n_elem;

    for (arma::uword i = 0; i < n; ++i)
    {
        const double value = data[i];

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
                    " must lie strictly inside (0, 1).");
            }
        }
        else if (value < 0.0 || value > 1.0)
        {
            Rcpp::stop(
                name +
                " must have all values in [0, 1].");
        }
    }
}

// [[Rcpp::export]]
Rcpp::List beta_mixture_em_cpp(
    const arma::vec &y,
    const arma::mat &initial_responsibility,
    int nL = 3,
    int maxiter = 25,
    double tol = 1e-6,
    int beta_maxit = 50,
    int beta_max_halving = 30,
    double beta_score_tol = 1e-10,
    double min_shape = 1e-10,
    double armijo = 1e-4,
    bool debug = false)
{
    if (nL < 2 || nL > 3)
    {
        Rcpp::stop("nL must be 2 or 3");
    }

    const arma::uword K = static_cast<arma::uword>(nL);
    const arma::uword n = y.n_elem;
    const arma::uword n_param = 3 * K;
    const double n_obs = static_cast<double>(n);

    if (initial_responsibility.n_rows != n ||
        initial_responsibility.n_cols != K)
    {
        Rcpp::stop(
            "initial_responsibility must be n x nL "
            "(n = length(y))");
    }

    arma::vec log_y(n);
    arma::vec log1m_y(n);
    for (arma::uword i = 0; i < n; ++i)
    {
        log_y[i] = std::log(y[i]);
        log1m_y[i] = std::log1p(-y[i]);
    }

    arma::mat responsibility = initial_responsibility;

    arma::vec a(K, arma::fill::ones);
    arma::vec b(K, arma::fill::ones);
    arma::vec eta(K, arma::fill::zeros);
    arma::vec mu(K, arma::fill::zeros);

    arma::vec param_state(n_param);
    arma::vec old_param_state(n_param);
    old_param_state.fill(arma::datum::inf);

    arma::vec log_component(K);
    arma::vec log_norm(K);

    Rcpp::CharacterVector fit_status(nL);
    Rcpp::CharacterVector fit_reason(nL);

    std::vector<double> parameter_criterion_trace;
    std::vector<double> loglik_criterion_trace;

    double loglikelihood = NA_REAL;
    double previous_loglikelihood = NA_REAL;
    double parameter_criterion =
        std::numeric_limits<double>::infinity();
    double loglik_criterion =
        std::numeric_limits<double>::infinity();
    bool converged = false;
    int completed_iterations = 0;

    for (int iter = 0; iter < maxiter; ++iter)
    {
        Rcpp::checkUserInterrupt();

        // One pass: component weights (eta) and Beta sufficient stats.
        // stats[k].weight == sum_i responsibility(i, k).
        std::vector<BetaStats> stats(K);

        for (arma::uword i = 0; i < n; ++i)
        {
            for (arma::uword k = 0; k < K; ++k)
            {
                stats[k].add(
                    y[i],
                    log_y[i],
                    log1m_y[i],
                    responsibility(i, k));
            }
        }

        for (arma::uword k = 0; k < K; ++k)
        {
            eta[k] = stats[k].weight / n_obs;

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

        // E-step terms independent of observation index i.
        arma::vec log_prior(K);
        arma::vec am1(K);
        arma::vec bm1(K);
        for (arma::uword k = 0; k < K; ++k)
        {
            log_prior[k] = std::log(eta[k]) + log_norm[k];
            am1[k] = a[k] - 1.0;
            bm1[k] = b[k] - 1.0;
        }

        loglikelihood = 0.0;

        for (arma::uword i = 0; i < n; ++i)
        {
            double maximum =
                -std::numeric_limits<double>::infinity();

            for (arma::uword k = 0; k < K; ++k)
            {
                log_component[k] =
                    log_prior[k] +
                    am1[k] * log_y[i] +
                    bm1[k] * log1m_y[i];

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

            loglikelihood += log_mixture;

            for (arma::uword k = 0; k < K; ++k)
            {
                responsibility(i, k) =
                    std::exp(
                        log_component[k] - log_mixture);
            }
        }

        for (arma::uword k = 0; k < K; ++k)
        {
            param_state[3 * k] = std::log(a[k]);
            param_state[3 * k + 1] = std::log(b[k]);
            param_state[3 * k + 2] = eta[k];
        }

        parameter_criterion =
            arma::max(arma::abs(param_state - old_param_state));

        if (std::isfinite(previous_loglikelihood))
        {
            loglik_criterion =
                std::abs(loglikelihood - previous_loglikelihood) /
                (1.0 + std::abs(previous_loglikelihood));
        }
        else
        {
            loglik_criterion =
                std::numeric_limits<double>::infinity();
        }

        if (debug)
        {
            parameter_criterion_trace.push_back(
                parameter_criterion);
            loglik_criterion_trace.push_back(loglik_criterion);
        }

        completed_iterations = iter + 1;
        old_param_state = param_state;
        previous_loglikelihood = loglikelihood;

        if (completed_iterations >= 2 &&
            parameter_criterion < tol &&
            loglik_criterion < tol)
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
        Rcpp::_["parameter_criterion"] = parameter_criterion,
        Rcpp::_["loglik_criterion"] = loglik_criterion,
        Rcpp::_["fit_status"] = fit_status,
        Rcpp::_["fit_reason"] = fit_reason,
        Rcpp::_["nL"] = nL);

    if (debug)
    {
        result["parameter_criterion_trace"] =
            Rcpp::wrap(parameter_criterion_trace);
        result["loglik_criterion_trace"] =
            Rcpp::wrap(loglik_criterion_trace);
    }

    return result;
}
