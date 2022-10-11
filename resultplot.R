# fitY <- bi_lst$y %>% 
#    group_by(time) %>%
#    mutate(
#      q025 = quantile(value, 0.025),
#      q25 = quantile(value, 0.25),
#      q50 = quantile(value, 0.5),
#      q75 = quantile(value, 0.75),
#      q975 = quantile(value, 0.975)
#    ) %>% ungroup() %>%
#    left_join(y %>% rename(Y = value))
#  
#  g1 <- ggplot(data = fitY) +
#    geom_ribbon(aes(x = time, ymin = q25, ymax = q75), alpha = 0.3) +
#    geom_ribbon(aes(x = time, ymin = q025, ymax = q975), alpha = 0.3) +
#    geom_line(aes(x = time, y = q50)) +
#    geom_point(aes(x = time, y = Y), colour = "Red") +
#    ylab("Daily new H1N1 clinical cases") +
#    xlab("Time-Day")
#  
#  plot_df <- bi_lst$x %>% mutate(value = exp(value)) %>%
#    group_by(time) %>%
#    mutate(
#      q025 = quantile(value, 0.025),
#      q25 = quantile(value, 0.25),
#      q50 = quantile(value, 0.5),
#      q75 = quantile(value, 0.75),
#      q975 = quantile(value, 0.975)
#    ) %>% ungroup()
#  
#  g2 <- ggplot(data = plot_df) +
#    geom_ribbon(aes(x = time, ymin = q25, ymax = q75), alpha = 0.3) +
#   geom_ribbon(aes(x = time, ymin = q025, ymax = q975), alpha = 0.3) +
#   geom_line(aes(x = time, y = q50)) +
#   ylab(TeX("Transmissibility ($\\beta(t)$)")) +
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
#   ylab(TeX("Relative trans. ($\\beta(t)-\\beta(0)$)")) +
#   xlab("Time-Day")
# 
# 
# ggarrange(g1, g2, g3, ncol = 1, nrow = 3, align = "v")