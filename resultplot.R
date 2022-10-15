

fitY <- read.csv("cv37h1n1y.csv", header=TRUE, stringsAsFactors=FALSE)
 g1 <- ggplot(data = fitY) +
   geom_ribbon(aes(x = time, ymin = q25, ymax = q75), alpha = 0.3) +
   geom_ribbon(aes(x = time, ymin = q025, ymax = q975), alpha = 0.3) +
   geom_line(aes(x = time, y = q50)) +
   geom_point(aes(x = time, y = Y), colour = "Red") +
   ylab("Daily new H1N1 clinical cases") +
   xlab("Time-Day")

 plot_df <- read.csv("cv37h1n1beta.csv", header=TRUE, stringsAsFactors=FALSE)

 g2 <- ggplot(data = plot_df) +
   geom_ribbon(aes(x = time, ymin = q25, ymax = q75), alpha = 0.3) +
  geom_ribbon(aes(x = time, ymin = q025, ymax = q975), alpha = 0.3) +
  geom_line(aes(x = time, y = q50)) +
  ylab(TeX("Transmissibility ($\\beta(t)$)")) +
  xlab("Time-Day")

plot_df1 <- read.csv("cv37h1n1beta1.csv", header=TRUE, stringsAsFactors=FALSE)

g3 <- ggplot(data = plot_df1) +
  geom_ribbon(aes(x = time, ymin = q25, ymax = q75), alpha = 0.3) +
  geom_ribbon(aes(x = time, ymin = q025, ymax = q975), alpha = 0.3) +
  geom_line(aes(x = time, y = q50)) +
  ylab(TeX("Relative trans. ($\\beta(t)-\\beta(0)$)")) +
  xlab("Time-Day")



ggarrange(g1, g2, g3, ncol = 1, nrow = 3, align = "v")

alpha<-read.csv("cv37h1n1alpha.csv", header=TRUE, stringsAsFactors=FALSE)
gamma<-read.csv("cv37h1n1gamma.csv", header=TRUE, stringsAsFactors=FALSE)
par(mfrow=c(2,1))
plot(alpha,type='l',main=TeX("Sampled onset H1N1-2009 England symptoms rate ($\\alpha$)"),xlab="PMCMC iterations after first 5000 burn-in and thinning by 5",ylab="alpha")
plot(gamma,type='l',main=TeX("Sampled H1N1-2009 England Recovery rate ($\\gamma$)"),xlab="PMCMC iterations after first 5000 burn-in and thinning by 5",ylab="gamma")
