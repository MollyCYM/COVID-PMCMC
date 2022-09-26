rm(list=ls())
library(rbi)
library(rbi.helpers)

# Load the data
v <- read.csv("london60w.csv", header=FALSE, stringsAsFactors=FALSE) %>%
  rowSums()
y <- data.frame(value = v) %>%
  mutate(time = seq(7, by = 7, length.out = n())) %>%
  dplyr::select(time, value)
#y<- y[1:37,]
ncores <- 8
minParticles <- max(ncores, 16)

model_str <- "
model dureau {
  obs y

  state S
  state E
  state I
  state R
  state x

  input N
  param k
  param gamma
  param sigma // Noise driver
  param E0
  param I0
  param R0
  param x0

  sub parameter {
    k ~ truncated_gaussian(5.6, 6, lower = 0) // k is the period here, not the rate, i.e. 1/k is the rate
    gamma ~ truncated_gaussian(11.5, 2.9, lower = 0) // gamma is the period, not the rate
    sigma ~ uniform(0,1)
    x0 ~ uniform(-5,2)
    I0 ~ uniform(-16, -9)
    E0 ~ uniform(-16, -9)
    R0 ~ truncated_gaussian(0.15, 0.15, lower = 0, upper = 1)
  }

  sub initial {
    S <- N
    R <- R0*S
    S <- S - R

    E <- exp(E0 + log(S))
    S <- S - E
    I <- exp(I0 + log(S))
    S <- S - I
    x <- x0
  }

  sub transition(delta = 1) {
    noise e
    e ~ wiener()
    ode(alg = 'RK4(3)', h = 1.0, atoler = 1.0e-3, rtoler = 1.0e-8) {
      dx/dt = sigma*e
      dS/dt = -exp(x)*S*(0.1*I+E)/N
      dE/dt = exp(x)*S*(0.1*I+E)/N - E*(1/k+1/gamma)
      dI/dt = E/k-I*(1/gamma+0.0087)
      dR/dt = (I+E)/gamma
    }
  }

  sub observation {
    y ~ poisson(rate=E/k)
  }

  sub proposal_parameter {
    k ~ gaussian(k, 0.005)
    sigma ~ gaussian(sigma, 0.01)
    gamma ~ gaussian(gamma, 0.01)
    x0 ~ gaussian(x0, 0.05)
    E0 ~ gaussian(E0, 0.05)
    I0 ~ gaussian(I0, 0.05)
    R0 ~ gaussian(R0, 0.05)
  }
}"
model <- bi_model(lines = stringi::stri_split_lines(model_str)[[1]])
bi_model <- libbi(model)
input_lst <- list(N = 9541000)
end_time <- max(y$time)
obs_lst <- list(y = y %>% dplyr::filter(time <= end_time))

bi <- sample(bi_model, end_time = end_time, input = input_lst, obs = obs_lst, nsamples = 1000, nparticles = minParticles, nthreads = ncores, proposal = 'prior') %>% 
  adapt_particles(min = minParticles, max = minParticles*200) %>%
  adapt_proposal(min = 0.05, max = 0.4) %>%
  sample(nsamples = 5000, thin = 5) %>% # burn in 
  sample(nsamples = 5000, thin = 5)

bi_lst <- bi_read(bi %>% sample_obs)

# fitY <- bi_lst$y %>% 
#   group_by(time) %>%
#  mutate(
#     q025 = quantile(value, 0.025),
#     q25 = quantile(value, 0.25),
#     q50 = quantile(value, 0.5),
#     q75 = quantile(value, 0.75),
#     q975 = quantile(value, 0.975)
#   ) %>% ungroup() %>%
#   left_join(y %>% rename(Y = value))
# 
# 
# g1 <- ggplot(data = fitY) +
#   geom_ribbon(aes(x = time, ymin = q25, ymax = q75),alpha=0.3) +
#   geom_ribbon(aes(x = time, ymin = q025, ymax = q975),alpha=0.3) +
#   geom_line(aes(x = time, y = q50)) +
#   geom_point(aes(x = time, y = Y), colour = "Red") +
#   ylab("Daily new confirmed cases") +
#   xlab("Time-Day")
# 
# plot_df <- bi_lst$x %>% mutate(value = exp(value)) %>%
#   group_by(time) %>%
#   mutate(
#     q025 = quantile(value, 0.025),
#     q25 = quantile(value, 0.25),
#     q50 = quantile(value, 0.5),
#     q75 = quantile(value, 0.75),
#     q975 = quantile(value, 0.975)
#   ) %>% ungroup()
# 
# g2 <- ggplot(data = plot_df) +
#   geom_ribbon(aes(x = time, ymin = q25, ymax = q75),alpha = 0.3) +
#   geom_ribbon(aes(x = time, ymin = q025, ymax = q975),alpha = 0.3) +
#   geom_line(aes(x = time, y = q50)) +
#   ylab(TeX("Transmission rate ($\\beta(t)$)")) +
#   xlab("Time-Day")
# 
# plot_df <- bi_lst$x %>% mutate(value = exp(value)) %>% 
#   group_by(np) %>% mutate(value = value - value[1]) %>%
#   group_by(time) %>%
#   mutate(
#     q025 = quantile(value, 0.025),
#     q25 = quantile(value, 0.25),
#     q50 = quantile(value, 0.5),
#     q75 = quantile(value, 0.75),
#     q975 = quantile(value, 0.975)
#   ) %>% ungroup()
# 
# g3 <- ggplot(data = plot_df) +
#   geom_ribbon(aes(x = time, ymin = q25, ymax = q75), alpha = 0.3) +
#   geom_ribbon(aes(x = time, ymin = q025, ymax = q975), alpha = 0.3) +
#   geom_line(aes(x = time, y = q50)) +
#   ylab(TeX("Relative transmission rate. ($\\beta(t)-\\beta(0)$)")) +
#   xlab("Time-Day")
# 
# 
# ggarrange(g1, g2, g3, ncol = 1, nrow = 3, align = "v")