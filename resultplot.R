

fitY <- read.csv("5winfection.csv", header=TRUE, stringsAsFactors=FALSE)
 g1 <- ggplot(data = fitY) +
   geom_ribbon(aes(x = time, ymin = q25, ymax = q75), alpha = 0.3) +
   geom_ribbon(aes(x = time, ymin = q025, ymax = q975), alpha = 0.3) +
   geom_line(aes(x = time, y = q50)) +
   geom_point(aes(x = time, y = Y), colour = "Red") +
   ylab("Daily new H1N1 clinical cases") +
   xlab("Time-Day")

 plot_df <- read.csv("5wbeta.csv", header=TRUE, stringsAsFactors=FALSE)

 g2 <- ggplot(data = plot_df) +
   geom_ribbon(aes(x = time, ymin = q25, ymax = q75), alpha = 0.3) +
  geom_ribbon(aes(x = time, ymin = q025, ymax = q975), alpha = 0.3) +
  geom_line(aes(x = time, y = q50)) +
  ylab(TeX("Transmissibility ($\\beta(t)$)")) +
  xlab("Time-Day")

plot_df1 <- read.csv("5wbeta1.csv", header=TRUE, stringsAsFactors=FALSE)

g3 <- ggplot(data = plot_df1) +
  geom_ribbon(aes(x = time, ymin = q25, ymax = q75), alpha = 0.3) +
  geom_ribbon(aes(x = time, ymin = q025, ymax = q975), alpha = 0.3) +
  geom_line(aes(x = time, y = q50)) +
  ylab(TeX("Relative trans. ($\\beta(t)-\\beta(0)$)")) +
  xlab("Time-Day")



ggarrange(g1, g2, g3, ncol = 1, nrow = 3, align = "v")

alpha<-read.csv("5walpha.csv", header=TRUE, stringsAsFactors=FALSE)
gamma<-read.csv("5wgamma.csv", header=TRUE, stringsAsFactors=FALSE)
par(mfrow=c(2,1))
plot(alpha,type='l',main=TeX("Sampled onset London COVID-19 symptoms rate ($\\alpha$)"),xlab="PMCMC iterations after first 5000 burn-in and thinning by 5",ylab="alpha")
plot(gamma,type='l',main=TeX("Sampled London COVID-19 Recovery rate ($\\gamma$)"),xlab="PMCMC iterations after first 5000 burn-in and thinning by 5",ylab="gamma")