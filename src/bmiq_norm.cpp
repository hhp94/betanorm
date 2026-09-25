// Beta-mixture EM derived from blc2/betaEst2 in Steve Horvath's BMIQcalibration,
// modified from RPMM (blc, betaEst), copyright E. Andres Houseman and
// Devin C. Koestler, GPL (>= 2). GPL (>= 3) here; see inst/COPYRIGHTS.

#include <Rcpp.h>

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <limits>
#include <string>
#include <vector>

namespace
{

    struct BetaStats
    {
        double weight = 0.0;
        double sum_w2 = 0.0;
        double mean = 0.0;
        double m2 = 0.0;
        double sum_log_y = 0.0;
        double sum_log1m_y = 0.0;

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

            const double new_weight = weight + effective_weight;
            const double delta = y - mean;

            mean += (effective_weight / new_weight) * delta;
            m2 += effective_weight * delta * (y - mean);
            weight = new_weight;
            sum_w2 += effective_weight * effective_weight;

            sum_log_y += effective_weight * log_y;
            sum_log1m_y += effective_weight * log1m_y;
        }

        // Kish effective sample size: (sum w)^2 / sum w^2. A soft E-step gives
        // almost every observation a positive but negligible responsibility,
        // so a raw count of contributing observations is a poor size guard;
        // n_eff is not.
        double effective_size() const
        {
            if (!(sum_w2 > 0.0))
            {
                return 0.0;
            }
            return (weight * weight) / sum_w2;
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

    // Log of the Beta normalizing constant, -log B(a, b). Shared by the
    // component log-likelihood and the E-step prior so the two sites cannot
    // drift numerically.
    inline double log_beta_norm(double a, double b)
    {
        return std::lgamma(a + b) -
               std::lgamma(a) -
               std::lgamma(b);
    }

    inline double beta_loglik_stats(
        double a,
        double b,
        const BetaStats &stats)
    {
        return (a - 1.0) * stats.sum_log_y +
               (b - 1.0) * stats.sum_log1m_y +
               stats.weight * log_beta_norm(a, b);
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
            !(stats.effective_size() > 1.0))
        {
            out.reason = "zero or insufficient effective sample size";
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
        const double mean_log1m_y =
            stats.sum_log1m_y / stats.weight;

        for (int iter = 0; iter < maxit; ++iter)
        {
            const double ab = a + b;
            const double digamma_ab = R::digamma(ab);

            const double score_a =
                mean_log_y - R::digamma(a) + digamma_ab;
            const double score_b =
                mean_log1m_y - R::digamma(b) + digamma_ab;

            const double scaled_score =
                std::max(std::abs(a * score_a),
                         std::abs(b * score_b));

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

            if (!(determinant > 0.0) ||
                !std::isfinite(determinant))
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

// NumericVector deliberately accepts both vectors and numeric matrices.
// [[Rcpp::export]]
void scan_finite_unit_interval_cpp(
    const Rcpp::NumericVector &x,
    std::string name = "x",
    bool require_open = false)
{
    const double *data = REAL(x);
    const R_xlen_t n = x.size();

    for (R_xlen_t i = 0; i < n; ++i)
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

// Gather a contiguous run of sample rows out of a samples x probes matrix into
// a probes x sample_count block. R matrices are column-major, so a sample
// (row) is strided; the returned block stores each sample contiguously, which
// makes the per-sample work and the scatter back cache-friendly without paying
// for a full transpose.
// [[Rcpp::export]]
Rcpp::NumericMatrix gather_sample_block_cpp(
    const Rcpp::NumericMatrix &x,
    int first_sample,
    int sample_count)
{
    const int n_samples = x.nrow();
    const int n_probes = x.ncol();
    const int first0 = first_sample - 1;

    if (first0 < 0 ||
        sample_count < 1 ||
        first0 + sample_count > n_samples)
    {
        Rcpp::stop("Invalid sample block");
    }

    Rcpp::NumericMatrix block(n_probes, sample_count);

    const double *x_ptr = REAL(x);
    double *block_ptr = REAL(block);

    for (int probe = 0; probe < n_probes; ++probe)
    {
        const R_xlen_t x_offset =
            static_cast<R_xlen_t>(n_samples) * probe;

        for (int local_sample = 0;
             local_sample < sample_count;
             ++local_sample)
        {
            block_ptr[
                probe +
                static_cast<R_xlen_t>(n_probes) * local_sample] =
                x_ptr[x_offset + first0 + local_sample];
        }
    }

    return block;
}

// Scatter a probes x sample_count block back into a samples x probes matrix.
// This deliberately mutates `destination` in place, which is safe only when it
// is a freshly allocated private matrix (never an alias of the caller's input);
// bmiq_calibration() allocates `calibrated` for exactly this reason.
// [[Rcpp::export]]
void scatter_sample_block_cpp(
    Rcpp::NumericMatrix destination,
    const Rcpp::NumericMatrix &block,
    int first_sample)
{
    const int n_samples = destination.nrow();
    const int n_probes = destination.ncol();
    const int sample_count = block.ncol();
    const int first0 = first_sample - 1;

    if (block.nrow() != n_probes ||
        first0 < 0 ||
        first0 + sample_count > n_samples)
    {
        Rcpp::stop("Invalid destination or sample block");
    }

    double *destination_ptr = REAL(destination);
    const double *block_ptr = REAL(block);

    for (int probe = 0; probe < n_probes; ++probe)
    {
        const R_xlen_t destination_offset =
            static_cast<R_xlen_t>(n_samples) * probe;

        for (int local_sample = 0;
             local_sample < sample_count;
             ++local_sample)
        {
            destination_ptr[
                destination_offset + first0 + local_sample] =
                block_ptr[
                    probe +
                    static_cast<R_xlen_t>(n_probes) * local_sample];
        }
    }
}

// [[Rcpp::export]]
Rcpp::List beta_mixture_em_cpp(
    const Rcpp::NumericVector &y,
    const Rcpp::NumericMatrix &initial_responsibility,
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

    const R_xlen_t K = static_cast<R_xlen_t>(nL);
    const R_xlen_t n = y.size();
    const R_xlen_t n_param = 3 * K;
    const double n_obs = static_cast<double>(n);

    if (static_cast<R_xlen_t>(
            initial_responsibility.nrow()) != n ||
        initial_responsibility.ncol() != nL)
    {
        Rcpp::stop(
            "initial_responsibility must be n x nL "
            "(n = length(y))");
    }

    const double *y_ptr = REAL(y);

    std::vector<double> log_y(
        static_cast<std::size_t>(n));
    std::vector<double> log1m_y(
        static_cast<std::size_t>(n));

    for (R_xlen_t i = 0; i < n; ++i)
    {
        const double value = y_ptr[i];

        // Own the EM's precondition here rather than relying on callers'
        // clipping: a boundary value would otherwise surface much later as an
        // obscure non-finite likelihood error.
        if (!(value > 0.0 && value < 1.0))
        {
            Rcpp::stop(
                "y must lie strictly inside (0, 1) "
                "for the Beta mixture EM");
        }

        log_y[i] = std::log(value);
        log1m_y[i] = std::log1p(-value);
    }

    // Allocate and copy rather than assigning the NumericMatrix directly.
    // Direct assignment aliases the caller's SEXP and the E-step would mutate
    // initial_responsibility. A fresh allocation also mirrors Armadillo's
    // attribute-dropping deep copy more closely than Rcpp::clone().
    Rcpp::NumericMatrix responsibility(
        initial_responsibility.nrow(),
        initial_responsibility.ncol());

    std::copy(
        initial_responsibility.begin(),
        initial_responsibility.end(),
        responsibility.begin());

    double *responsibility_ptr = REAL(responsibility);

    std::vector<double> a(
        static_cast<std::size_t>(K), 1.0);
    std::vector<double> b(
        static_cast<std::size_t>(K), 1.0);
    std::vector<double> eta(
        static_cast<std::size_t>(K), 0.0);
    std::vector<double> mu(
        static_cast<std::size_t>(K), 0.0);

    std::vector<double> param_state(
        static_cast<std::size_t>(n_param), 0.0);
    std::vector<double> old_param_state(
        static_cast<std::size_t>(n_param),
        std::numeric_limits<double>::infinity());

    std::vector<double> log_component(
        static_cast<std::size_t>(K), 0.0);
    std::vector<double> exp_shifted(
        static_cast<std::size_t>(K), 0.0);
    std::vector<double> log_prior(
        static_cast<std::size_t>(K), 0.0);
    std::vector<double> am1(
        static_cast<std::size_t>(K), 0.0);
    std::vector<double> bm1(
        static_cast<std::size_t>(K), 0.0);
    std::vector<double> prev_a(
        static_cast<std::size_t>(K), 0.0);
    std::vector<double> prev_b(
        static_cast<std::size_t>(K), 0.0);

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

        // Retain the previous iteration's component shapes so the M-step below
        // can refuse a candidate that would lower a component's weighted
        // log-likelihood (generalized-EM ascent guard).
        prev_a = a;
        prev_b = b;

        // One pass: component weights (eta) and Beta sufficient stats.
        // stats[k].weight == sum_i responsibility(i, k).
        std::vector<BetaStats> stats(
            static_cast<std::size_t>(K));

        for (R_xlen_t i = 0; i < n; ++i)
        {
            for (R_xlen_t k = 0; k < K; ++k)
            {
                stats[k].add(
                    y_ptr[i],
                    log_y[i],
                    log1m_y[i],
                    responsibility_ptr[i + n * k]);
            }
        }

        for (R_xlen_t k = 0; k < K; ++k)
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

            // Generalized-EM ascent guard: accept the candidate only if it does
            // not lower this component's weighted log-likelihood relative to the
            // previous iteration's shapes under the current responsibilities.
            // Each M-step re-initializes from method-of-moments, so a stalled
            // Newton fit could otherwise decrease the objective and break the
            // monotonicity that "usable" fits are assumed to preserve. eta is
            // always the closed-form maximizer, so guarding the Beta shapes is
            // enough to keep the outer log-likelihood non-decreasing.
            bool retained = false;

            if (iter > 0)
            {
                const double prev_loglik =
                    beta_loglik_stats(
                        prev_a[k],
                        prev_b[k],
                        stats[k]);

                if (std::isfinite(prev_loglik) &&
                    prev_loglik > fit.loglik)
                {
                    a[k] = prev_a[k];
                    b[k] = prev_b[k];
                    fit_status[k] = "converged";
                    fit_reason[k] =
                        "retained previous shapes; M-step candidate did not "
                        "improve the weighted log-likelihood";
                    retained = true;
                }
            }

            if (!retained)
            {
                a[k] = fit.a;
                b[k] = fit.b;
                fit_status[k] = fit.status;
                fit_reason[k] = fit.reason;
            }

            mu[k] = a[k] / (a[k] + b[k]);

            // E-step terms independent of observation index i; eta[k], a[k],
            // and b[k] are final for this iteration at this point.
            log_prior[k] =
                std::log(eta[k]) + log_beta_norm(a[k], b[k]);
            am1[k] = a[k] - 1.0;
            bm1[k] = b[k] - 1.0;
        }

        loglikelihood = 0.0;

        for (R_xlen_t i = 0; i < n; ++i)
        {
            double maximum =
                -std::numeric_limits<double>::infinity();

            for (R_xlen_t k = 0; k < K; ++k)
            {
                log_component[k] =
                    log_prior[k] +
                    am1[k] * log_y[i] +
                    bm1[k] * log1m_y[i];

                maximum =
                    std::max(maximum, log_component[k]);
            }

            double sum_exp = 0.0;

            for (R_xlen_t k = 0; k < K; ++k)
            {
                exp_shifted[k] =
                    std::exp(log_component[k] - maximum);
                sum_exp += exp_shifted[k];
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

            // Single-exp E-step: responsibilities reuse the max-shifted
            // exponentials rather than re-exponentiating against log_mixture.
            // Not bit-identical to the log-space form (see DECISIONS.md,
            // single-exp adoption); snapshots pin the current numerics.
            const double inv_sum = 1.0 / sum_exp;

            for (R_xlen_t k = 0; k < K; ++k)
            {
                responsibility_ptr[i + n * k] =
                    exp_shifted[k] * inv_sum;
            }
        }

        for (R_xlen_t k = 0; k < K; ++k)
        {
            param_state[3 * k] = std::log(a[k]);
            param_state[3 * k + 1] = std::log(b[k]);
            param_state[3 * k + 2] = eta[k];
        }

        parameter_criterion = 0.0;

        for (R_xlen_t j = 0; j < n_param; ++j)
        {
            parameter_criterion =
                std::max(
                    parameter_criterion,
                    std::abs(
                        param_state[j] -
                        old_param_state[j]));
        }

        if (std::isfinite(previous_loglikelihood))
        {
            loglik_criterion =
                std::abs(
                    loglikelihood -
                    previous_loglikelihood) /
                (1.0 +
                 std::abs(previous_loglikelihood));
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
            loglik_criterion_trace.push_back(
                loglik_criterion);
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

    // Preserve the old RcppArmadillo return shapes:
    // a, b, and mu are K x 1 matrices; eta is a bare vector.
    Rcpp::NumericMatrix a_matrix(nL, 1);
    Rcpp::NumericMatrix b_matrix(nL, 1);
    Rcpp::NumericMatrix mu_matrix(nL, 1);
    Rcpp::NumericVector eta_vector(nL);

    double *a_matrix_ptr = REAL(a_matrix);
    double *b_matrix_ptr = REAL(b_matrix);
    double *mu_matrix_ptr = REAL(mu_matrix);
    double *eta_vector_ptr = REAL(eta_vector);

    for (R_xlen_t k = 0; k < K; ++k)
    {
        a_matrix_ptr[k] = a[k];
        b_matrix_ptr[k] = b[k];
        mu_matrix_ptr[k] = mu[k];
        eta_vector_ptr[k] = eta[k];
    }

    Rcpp::List result = Rcpp::List::create(
        Rcpp::_["a"] = a_matrix,
        Rcpp::_["b"] = b_matrix,
        Rcpp::_["eta"] = eta_vector,
        Rcpp::_["mu"] = mu_matrix,
        Rcpp::_["w"] = responsibility,
        Rcpp::_["llike"] = loglikelihood,
        Rcpp::_["iterations"] = completed_iterations,
        Rcpp::_["converged"] = converged,
        Rcpp::_["parameter_criterion"] =
            parameter_criterion,
        Rcpp::_["loglik_criterion"] =
            loglik_criterion,
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