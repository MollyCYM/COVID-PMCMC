---
  title: "Random-walk model"
author: "Akira Endo"
date: "12 December 2018"
output:
html_document: default
pdf_document: default
---
  
  ```{r setup, include=FALSE}
knitr::opts_chunk$set(echo = TRUE)
```

## Random-walk model
This is a simple toy model (random-walk model) used as an example to in the main text. We assume that \(x_t\) evolves following a simple Gaussian kernel \(x_t\sim \mathcal N(x_{t-1},\sigma)\) with an initial state \(x_1=0\). The data is assumed to be obtained as a rounded integer with a Gaussian measurement error \(y_t=round(\mathcal N(x_{t-1},0.1))\). Consider that a five-step observation (\(T=5\)) with \(\sigma=0.5\) yielded \(y_{1:5}=(0,1,1,1,2)\). Without any observations, the possible trajectories of \(x_{1:5}\) are very diverse (Figure 1, Left panel). Trajectories \(x_{1:5}\) consistent with the observation are only a fractions of those (Figure 1, Right panel), and they need to be efficiently sampled to compute the marginal likelihood \(p(y_{1:T}|\sigma)\).

```{r prep}
## Preparation
yt<-c(0,1,1,1,2) #observed data
set.seed(122018)
#SMC by the bootstrap filter
BF_sim<-function(yt,sig,merror,size){
  tlen<-length(yt)
  xt<-matrix(0,tlen,size)
  xt_prop<-matrix(0,tlen,size)
  avg_wt<-numeric(tlen)
  avg_wt[1]=1
  for(t in 2:tlen){
    xt_prop[t,]<-rnorm(size,xt[t-1,],sig) #proposal particles
    weights=pnorm(yt[t]+0.5,xt_prop[t,],merror)-pnorm(yt[t]-0.5,xt_prop[t,],merror) #weighted by p(yt|xt_prop)
    avg_wt[t]=sum(weights)/size
    if(sum(weights)==0)weights=weights+1 #marginal likelihood 0: avoid error in sample() below due to weights=0
    ind<-sample(size,size,replace=T,prob=weights) #resample
    xt[1:t,]=rbind(xt[1:(t-1),ind],xt_prop[t,ind])
    
  }
  return(list(xt=xt,xt_prop=xt_prop,avg_wt=avg_wt))
}
sig=0.5
merror=0.1
xt_prior<-apply(rbind(0,matrix(rnorm(10000,0,sig),4)),2,cumsum) #prior samples of xt
xt_BF<-BF_sim(yt,sig,merror,10000)
xt_posterior<-xt_BF$xt#posterior samples of x_t
xt_prop<-xt_BF$xt_prop
```

```{r plot-traj, fig.width=8,fig.height=4}
## Figure 1. Trajectories of x_t
par(mfrow=c(1,2))
matplot(xt_prior,type="l",lty=1,ylim=c(-3,3),ylab="hidden state",xlab="time")
matplot(xt_posterior[,1:500],type="l",lty=1,ylim=c(-3,3),main="",ylab="hidden state",xlab="time")
```
```{r plot-pos, fig.width=12,fig.height=3}
## Figure 2. Posterior distribution p(x|y)
par(mfrow=c(1,4),cex=1)
for(i in 2:5){
  hist(xt_posterior[i,],xlab=paste0("x",i),main="",breaks=15,col="gray",freq=F,ylim=c(0,2))
  abline(v=c(0,1,1,1,2)[i],col="red")
}
```

## SMC
The bootstrap filter is one of the most popular algorithms for sequantial Monte Carlo (SMC). At each time step, candidate particles are generated from the time-evolution process \(f_\theta(x_t|x_{t-1})\). Particles are then filetered by the importance sampling resampling (ISR) based on \(g_\theta(y_t|x_t)\) to yield samples from the target distribution \(p_\theta(x_t|y_{1:t})\propto g_\theta(y_t|x_t)f_\theta(x_t|x_{t-1},y_{1:t-1})\).

```{r smc}
set.seed(122018)
xt_BF<-BF_sim(yt,sig,merror=0.1,200) #SMC
xt_posterior<-xt_BF$xt #posterior samples of xt
xt_prop<-xt_BF$xt_prop #proposed particles at each time step
```
```{r plot-smc, fig.width=12, fig.height=6}
## Figure 3. Sequential sampling of trajectories x
par(mfrow=c(2,4))
for(t in 2:5){ #plot trajectories
  matplot(rbind(xt_posterior[1:(t-1),],xt_prop[t,],matrix(NA,5-t,ncol(xt_posterior))),type="l",lty=1,ylim=c(-1,3),ylab="hidden state",xlab="time")
  matplot(rbind(xt_posterior[1:t,],matrix(NA,5-t,ncol(xt_posterior))),type="l",lty=1,ylim=c(-1,3),ylab="hidden state",xlab="time")
}
```

##PMCMC
Particle Markov-chain Monte Carlo is a method to estimate parameter \(\theta\) of a hidden-Markov process by efficiently sampling the hidden variable \(x\) with SMC. SMC is used as a sub-algorithm that returns the approximated marginal likelihood \(\hat p(y_{1:T}|\theta)\). Here, the parameter in the Gaussian kernel \(\sigma\) is estimated with an improper prior over the positive real line \(p(\sigma)=\mathcal U(0,\infty)\).

```{r pmcmc}
set.seed(122018)
tlen<-5 # total time steps
mcmclen<-5000 # Iterations in MCMC
nparticle<-200 # Particles used in SMC 
merror=0.1 #measurement error
# Initialise MCMC
sig_mcmc<-numeric(mcmclen) #sigma 
sig_mcmc[1]=1
ll<-numeric(mcmclen) #log-likelihood
ll[1]=-Inf
x_mcmc<-matrix(0,mcmclen,tlen) #x_t
prop_sd<-0.5 #sd of proposal distribution
# PMCMC
for(n in 2:mcmclen){
  sig_prop<-sig_mcmc[n-1]+rnorm(1,0,prop_sd)
  sig_prop=abs(sig_prop)
  BF<-BF_sim(yt,sig_prop,merror,nparticle)
  ll_prop<-sum(log(BF$avg_wt))
  prob_update<-min(1,exp(ll_prop-ll[n-1]))
  if(runif(1)<prob_update){
    sig_mcmc[n]=sig_prop
    x_mcmc[n,]=BF$xt[,sample(nparticle,1)]
    ll[n]=ll_prop
  }else{
    sig_mcmc[n]=sig_mcmc[n-1]
    x_mcmc[n,]=x_mcmc[n-1,]
    ll[n]=ll[n-1]
  }
}
```
```{r plot-pmcmc, fig.width=12, fig.height=8}
## Figure 4. PMCMC results
par(mfrow=c(2,2))
# MCMC chain for sigma
plot(sig_mcmc,type="l",xlab="iteration",ylab="sigma")
# Posterior distribution for sigma
hist(sig_mcmc[-(1:1000)],breaks=100,freq=F,main="",xlab="sigma",border="dimgray")
# Median and 95% Credible intervals
abline(v=quantile(sig_mcmc[-(1:1000)],c(0.025,0.5,0.975)),lty=c(2,1,2),lwd=2,col=1)
# Posterior samples of X_1:5
matplot(t(x_mcmc[-(1:1000),]),type="l",lty=1,ylim=c(-1,3),ylab="hidden state",xlab="time")
```