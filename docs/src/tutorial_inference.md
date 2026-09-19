# Inference

## Exact telegraph likelihood

For snapshot counts of a gene without division effects, the stationary law is
Beta-Poisson and can be fitted by maximum likelihood:

```julia
fit = fit_telegraph(counts)          # (k_on, k_off, k_tx, loglik, converged), rates in units of the degradation rate
telegraph_loglik(counts, fit.k_on, fit.k_off, fit.k_tx)
```

## ABC-SMC

For any model, including growing and dividing populations, use sequential
Monte Carlo approximate Bayesian computation. Supply a function that
simulates with a parameter vector and returns a distance to the observed
summaries:

```julia
sobs = moment_summaries(observed_counts)
function dist(θ, rng)
    p = set_params(model; k_on = θ[1], k_off = θ[2], k_tx = θ[3])
    E = ensemble_final(model, x0, 30.0, 500; p = p, rng = rng, threads = false)
    summary_distance(moment_summaries(E[:, 3]), sobs)
end
prior = [LogUniform(1e-3, 10.0), LogUniform(1e-3, 10.0), LogUniform(1.0, 200.0)]
post  = abc_smc(dist, prior; n_particles = 200, generations = 6, names = [:k_on, :k_off, :k_tx])
posterior_mean(post); credible_interval(post, :k_tx)
```

The distance function must be thread-safe (use the `rng` it receives).
