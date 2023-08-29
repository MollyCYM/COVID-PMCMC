####################################################Seb example code####################################################
rm(list=ls())
library(rbi)
library(rbi.helpers)

set.seed(1234)
end_time <- 101

ncores <- 8
minParticles <- max(ncores, 16)
model_str <- "
model PZ {
  const c = 0.25   // zooplankton clearance rate
  const e = 0.3    // zooplankton growth efficiency
  const m_l = 0.1  // zooplankton linear mortality
  const m_q = 0.1  // zooplankton quadratic mortality

  param mu, sigma  // mean and standard deviation of phytoplankton growth
  state P, Z       // phytoplankton, zooplankton
  noise alpha      // stochastic phytoplankton growth rate
  obs P_obs        // observations of phytoplankton

  sub parameter {
    mu ~ uniform(0.0, 1.0)
    sigma ~ uniform(0.0, 0.5)
  }

  sub proposal_parameter {
    mu ~ truncated_gaussian(mu, 0.02, 0.0, 1.0);
    sigma ~ truncated_gaussian(sigma, 0.01, 0.0, 0.5);
  }

  sub initial {
    P <- log(2.0)  // Change P_0 starts from a fixed number
    Z ~ log_normal(log(2.0), 0.1)
  }

  sub transition(delta = 1.0) {
    alpha ~ gaussian(mu, sigma)
    ode {
      dP/dt = alpha*P - c*P*Z
      dZ/dt = e*c*P*Z - m_l*Z - m_q*Z*Z
    }
  }

  sub observation {
    P_obs ~ log_normal(log(P), 0.2)
  }
}"
model <- bi_model(lines = stringi::stri_split_lines(model_str)[[1]])
bi_model <- libbi(model)
# create simulated data set
sim <-generate_dataset(bi_model, end_time = 100)

obs_lst <- bi_read(sim)["P_obs"]

bi <- sample(
  bi_model, end_time = end_time, obs = obs_lst, nsamples = 200,
  nparticles = minParticles, nthreads = ncores, proposal = 'model'
) |>
  adapt_particles(min = minParticles, max = minParticles * 500) |>
  adapt_proposal(min = 0.1, max = 0.4) |>
  sample(nsamples = 100, thin = 1)
#> Mon Aug 28 10:10:05 2023 Adapting the proposal distribution
#> Mon Aug 28 10:10:05 2023 Adapting the number of particles
#> Mon Aug 28 10:10:24 2023 16 particles, loglikelihod variance: 9.83112526189252
#> Mon Aug 28 10:10:27 2023 32 particles, loglikelihod variance: 3.29774414840902
#> Mon Aug 28 10:10:28 2023 64 particles, loglikelihod variance: 2.82750059586941
#> Mon Aug 28 10:10:30 2023 128 particles, loglikelihod variance: 2.79694072306243
#> Mon Aug 28 10:10:33 2023 256 particles, loglikelihod variance: 4.27145834376161
#> Mon Aug 28 10:10:37 2023 512 particles, loglikelihod variance: 1.5736447523017
#> Mon Aug 28 10:10:43 2023 1024 particles, loglikelihod variance: 1.26738088354077
#> Mon Aug 28 10:10:54 2023 2048 particles, loglikelihod variance: 0.502907884993449
#> Mon Aug 28 10:10:54 2023 Choosing 2048 particles.
#> Mon Aug 28 10:10:55 2023 Initial trial run
#> Mon Aug 28 10:11:13 2023 Acceptance rate 0.507537688442211, adapting size with scale 1
#> Mon Aug 28 10:11:25 2023 Acceptance rate: 0.396984924623116

bi_lst <- bi_read(bi |> sample_obs())

packageVersion("rbi")
#> [1] '1.0.0'
packageVersion("rbi.helpers")
#> [1] '0.4.0'