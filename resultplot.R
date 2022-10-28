

fitY <- read.csv("SEIRy.csv", header=TRUE, stringsAsFactors=FALSE)
g1 <- ggplot(data = fitY) +
  geom_ribbon(aes(x = time, ymin = q25, ymax = q75), alpha = 0.3) +
  geom_ribbon(aes(x = time, ymin = q025, ymax = q975), alpha = 0.3) +
  geom_line(aes(x = time, y = q50)) +
  geom_point(aes(x = time, y = Y), colour = "Red") +
  ylab("Daily new confirmed cases") +
  xlab("Time-Day")
#plot(g1)
#lines(Z[,1]/5,type='l',col='red')

 plot_df <- read.csv("60wbeta3.csv", header=TRUE, stringsAsFactors=FALSE)

 g2 <- ggplot(data = plot_df) +
   geom_ribbon(aes(x = time, ymin = q25, ymax = q75), alpha = 0.3) +
  geom_ribbon(aes(x = time, ymin = q025, ymax = q975), alpha = 0.3) +
  geom_line(aes(x = time, y = q50)) +
  ylab(TeX("Transmissibility ($\\beta(t)$)")) +
  xlab("Time-Day")

plot_df1 <- read.csv("60wbeta13.csv", header=TRUE, stringsAsFactors=FALSE)

g3 <- ggplot(data = plot_df1) +
  geom_ribbon(aes(x = time, ymin = q25, ymax = q75), alpha = 0.3) +
  geom_ribbon(aes(x = time, ymin = q025, ymax = q975), alpha = 0.3) +
  geom_line(aes(x = time, y = q50)) +
  ylab(TeX("Relative trans. ($\\beta(t)-\\beta(0)$)")) +
  xlab("Time-Day")



ggarrange(g1, g2, g3, ncol = 1, nrow = 3, align = "v")
abline(h = x)
sigma<-read.csv("SEIRsigma.csv", header=TRUE, stringsAsFactors=FALSE)
gamma<-read.csv("SEIRgamma.csv", header=TRUE, stringsAsFactors=FALSE)
beta<-read.csv("SEIRbeta.csv", header=TRUE, stringsAsFactors=FALSE)
par(mfrow=c(3,1))
plot(sigma,type='l',main=TeX("Sampled onset COVID-19 England symptoms rate ($\\sigma$)"),xlab="PMCMC iterations after first 5000 burn-in and thinning by 5",ylab="sigma")
abline(h=0.25, col="red")
abline(h=mean(sigma[,2]),col="blue")
plot(gamma,type='l',main=TeX("Sampled COVID-19 England Recovery rate ($\\gamma$)"),xlab="PMCMC iterations after first 5000 burn-in and thinning by 5",ylab="gamma",ylim=c(0.18,max(gamma[,2])))
abline(h=0.2, col="red")
abline(h=mean(gamma[,2]),col="blue")
plot(beta,type='l',main=TeX("Sampled COVID-19 England Contact rate ($\\beta$)"),xlab="PMCMC iterations after first 5000 burn-in and thinning by 5",ylab="beta",ylim=c(0.4,0.7))
abline(h=0.5, col="red")
abline(h=mean(beta[,2]),col="blue")











g2=rtruncnorm(1000,a=0, b=Inf, mean = 0.5, sd = 0.3)
G2<-density(g2)
g3<- rtruncnorm(1000,a=0, b=Inf, mean = 0.001, sd = 0.25)
G3<-density(g3)

plot(G,type='l',ylim=c(0,8),xlim=c(0,2))
lines(G1,col='red') #
lines(G2,col='blue') #beta
lines(G3)
