

fitY <- read.csv("covid365_y1.csv", header=TRUE, stringsAsFactors=FALSE)
g1 <- ggplot(data = fitY) +
  geom_ribbon(aes(x = time, ymin = q25, ymax = q75), alpha = 0.3) +
  geom_ribbon(aes(x = time, ymin = q025, ymax = q975), alpha = 0.3) +
  geom_line(aes(x = time, y = q50)) +
  geom_point(aes(x = time, y = Y), colour = "Red") +
  ylab("Daily new confirmed cases") +
  xlab("Time-Day")
fitS <- read.csv("covid365_S1.csv", header=TRUE, stringsAsFactors=FALSE)
g2 <- ggplot(data = fitS) +
  geom_ribbon(aes(x = time, ymin = q25, ymax = q75), alpha = 0.3) +
  geom_ribbon(aes(x = time, ymin = q025, ymax = q975), alpha = 0.3) +
  geom_line(aes(x = time, y = q50)) +
  geom_point(aes(x = time, y = S), colour = "Red")  +
  ylab("S")
fitE <- read.csv("covid259_E2.csv", header=TRUE, stringsAsFactors=FALSE)
g3 <- ggplot(data = fitE) +
  geom_ribbon(aes(x = time, ymin = q25, ymax = q75), alpha = 0.3) +
  geom_ribbon(aes(x = time, ymin = q025, ymax = q975), alpha = 0.3) +
  geom_line(aes(x = time, y = q50)) +
  geom_point(aes(x = time, y = E), colour = "Red") +
  ylab("E")
fitI <- read.csv("covid259_I2.csv", header=TRUE, stringsAsFactors=FALSE)
g4 <- ggplot(data = fitI) +
  geom_ribbon(aes(x = time, ymin = q25, ymax = q75), alpha = 0.3) +
  geom_ribbon(aes(x = time, ymin = q025, ymax = q975), alpha = 0.3) +
  geom_line(aes(x = time, y = q50)) +
  geom_point(aes(x = time, y = I), colour = "Red")  +
  ylab("I")
fitR <- read.csv("covid259_R2.csv", header=TRUE, stringsAsFactors=FALSE)
g5 <- ggplot(data = fitR) +
  geom_ribbon(aes(x = time, ymin = q25, ymax = q75), alpha = 0.3) +
  geom_ribbon(aes(x = time, ymin = q025, ymax = q975), alpha = 0.3) +
  geom_line(aes(x = time, y = q50)) +
  geom_point(aes(x = time, y = R), colour = "Red")  +
  ylab("R")
fitM <- read.csv("MM.csv", header=TRUE, stringsAsFactors=FALSE)
g6 <- ggplot(data = fitM) +
  geom_ribbon(aes(x = time, ymin = q25, ymax = q75), alpha = 0.3) +
  geom_ribbon(aes(x = time, ymin = q025, ymax = q975), alpha = 0.3) +
  geom_line(aes(x = time, y = q50)) +
  geom_point(aes(x = time, y = M), colour = "Red")




 plot_df <- read.csv("covid259_beta2.csv", header=TRUE, stringsAsFactors=FALSE)

 g6 <- ggplot(data = plot_df) +
   geom_ribbon(aes(x = time, ymin = q25, ymax = q75), alpha = 0.3) +
  geom_ribbon(aes(x = time, ymin = q025, ymax = q975), alpha = 0.3) +
  geom_line(aes(x = time, y = q50)) +
  geom_point(aes(x = time, y = B), colour = "Red")  +
  ylab(TeX("Transmissibility ($\\beta(t)$)")) +
  xlab("Time-Day")

plot_df1 <- read.csv(".csv", header=TRUE, stringsAsFactors=FALSE)

g3 <- ggplot(data = plot_df1) +
  geom_ribbon(aes(x = time, ymin = q25, ymax = q75), alpha = 0.3) +
  geom_ribbon(aes(x = time, ymin = q025, ymax = q975), alpha = 0.3) +
  geom_line(aes(x = time, y = q50)) +
  ylab(TeX("Relative trans. ($\\beta(t)-\\beta(0)$)")) +
  xlab("Time-Day")

alpha<-read.csv("covid259_alpha2.csv", header=TRUE, stringsAsFactors=FALSE)
gamma<-read.csv("covid259_gamma2.csv", header=TRUE, stringsAsFactors=FALSE)
par(mfrow=c(2,1))
plot(alpha,type='l',main=TeX("Trace plot of ($\\alpha$)"),xlab="PMCMC iterations after first 1000 burn-in and following 10000 thinning by 5",ylab="sigma")
abline(h=0.3, col="red")
abline(h=mean(alpha[,2]),col="blue")
plot(gamma,type='l',main=TeX("Trace plot of ($\\gamma$)"),xlab="PMCMC iterations after first 1000 burn-in and following 10000 thinning by 5",ylab="gamma")
abline(h=7, col="red")
abline(h=mean(gamma[,2]),col="blue")


ggarrange(g1, g2, g3, ncol = 1, nrow = 3, align = "v")
abline(h = x)
sigma<-read.csv("longdays_sigma1.csv", header=TRUE, stringsAsFactors=FALSE)
gamma<-read.csv("longdays_gamma1.csv", header=TRUE, stringsAsFactors=FALSE)
beta<-read.csv("longdays_beta1.csv", header=TRUE, stringsAsFactors=FALSE)
mu<- read.csv("longdays_mu1.csv", header=TRUE, stringsAsFactors=FALSE)
par(mfrow=c(4,1))
plot(sigma,type='l',main=TeX("Trace plot of ($\\sigma$)"),xlab="PMCMC iterations after first 5000 burn-in and following 10000 thinning by 5",ylab="sigma")
abline(h=0.25, col="red")
abline(h=mean(sigma[,2]),col="blue")
plot(gamma,type='l',main=TeX("Trace plot of ($\\gamma$)"),xlab="PMCMC iterations after first 5000 burn-in and following 10000 thinning by 5",ylab="gamma")
abline(h=0.2, col="red")
abline(h=mean(gamma[,2]),col="blue")
plot(beta,type='l',main=TeX("Trace plot of ($\\beta$)"),xlab="PMCMC iterations after first 5000 burn-in and following 10000 thinning by 5",ylab="beta")
abline(h=0.5, col="red")
abline(h=mean(beta[,2]),col="blue")
plot(mu,type='l',main=TeX("Trace plot of ($\\mu$)"),xlab="PMCMC iterations after first 5000 burn-in and following 10000 thinning by 5",ylab="mu")
abline(h=0.001, col="red")
abline(h=mean(mu[,2]),col="blue")

gamma<-read.csv("h1n1_gamma.csv", header=TRUE, stringsAsFactors=FALSE)
plot(gamma,type='l')
abline(h=1/1.08, col="red")
abline(h=mean(gamma[,2]),col="blue")











g2=rtruncnorm(1000,a=0, b=Inf, mean = 0.5, sd = 0.3)
G2<-density(g2)
g3<- rtruncnorm(1000,a=0, b=Inf, mean = 0.001, sd = 0.25)
G3<-density(g3)

plot(G,type='l',ylim=c(0,8),xlim=c(0,2))
lines(G1,col='red') #
lines(G2,col='blue') #beta
lines(G3)
