
rm(list=ls())
require(deSolve)


#Simulate an O-U Process
ornstein_uhlenbeck <- function(n,theta,sigma,x0,lens,mu,e){
  dt <- 1
  x<-vector(length=lens)
  for (i in 1:n) {
    if (i==1){x[i]=x0}
    else{
      x[i]  <-  x[i-1] + theta*(mu[i]-x[i])*dt + sigma*e[i]}
  }
  
  return(x);
}



gen_ode_model <- function(rep_time = 200, ori_seed = 2023, tims = 365, 
                          n = 52196381, pars = c(k=2, gamma=9),
                          ini_state = c(t=1,S=N-1, E=1, I=0, R=0),
                          tau = 0.1){
  
  times = 1:tims
  N = n
  
  mu = vector(length=tims) 
  for (t in 1:tims) {
    if (t <= 120){
      mu[t]=-0.02      }                    #No lock-down policy
    else {mu[t]=-0.2-0.02 }                 #Lock-down policy
  }
  mu = ts(mu)
  
  Zs = setNames(data.frame(matrix(NA, nrow = tims, ncol = rep_time)), paste0("Z", 1:rep_time))
  Ys = setNames(Zs, paste0("Y", 1:rep_time))
  dSs = setNames(Zs, paste0("dS", 1:rep_time))
  dEs = setNames(Zs, paste0("dE", 1:rep_time))
  dIs = setNames(Zs, paste0("dI", 1:rep_time))
  dRs = setNames(Zs, paste0("dR", 1:rep_time)) 
  
  params = pars   #c(k=2, gamma=9)
  initial_state = ini_state #ini_state = c(t=1,S=N-1, E=1, I=0, R=0)
  tau1 = tau #tau = 0.1
  
  #Main ODE Model
  Covid_OU <- function(time, current_state, params){
    
    with(as.list(c(current_state, params)),{
      
      # N <- S+E+I+R
      dt <- 1
      dS <- -exp(x[t])*S*(E+0.1*I)/N
      dE <- exp(x[t])*S*(E+0.1*I)/N - E*(1/k+1/gamma)
      dI <- E/k-I*(1/gamma+0.0087)
      dR <- (I+E)/gamma+0.0087*I
      
      return(list(c( dt, dS, dE, dI, dR)))
    })
  }
  
  ## rep 200 times 
  for (rep_i in 1:rep_time) {
    set.seed(ori_seed + 1001*rep_i) ###
    
    e = rnorm(n = length(times) , sd = 1)
    
    x = ornstein_uhlenbeck(tims, 0.05, sqrt(0.004), 0, tims,mu,e) ###
    
    
    model4 = ode(initial_state, times, Covid_OU, params)
    
    Z1 = model4[,4]/2
    
    Y1 = vector(length = tims)
    
    for (i in 1:tims){
      Y1[i] = rlnorm(1,log(Z1[i]/5),tau1)
    }
    
    dSi = model4[,3]
    dEi = model4[,4]
    dIi = model4[,5]
    dRi = model4[,6]
    # print(Z1)
    Zs[, rep_i] = Z1
    Ys[, rep_i] = Y1
    dSs[, rep_i] = dSi
    dEs[, rep_i] = dEi
    dIs[, rep_i] = dIi
    dRs[, rep_i] = dRi
  }
  
  # output
  return(list(    
    Z = Zs,
    Y = Ys,
    dS = dSs,
    dE = dEs,
    dI = dIs,
    dR = dRs))
}


ode_example <- gen_ode_model(200)

## for z

zs = ncol(ode_example$Z)

for (i in 1:zs) {
  if (i == 1){
    plot(ode_example$Z[,1], col="darkblue", type = "l", 
         ylim = c(0, max(ode_example$Z)),
         xlab = 'Times', ylab = 'Z')
    }
  else{
    lines(ode_example$Z[,i], col="darkblue")
  }
}

## for y

ys = ncol(ode_example$Y)

for (i in 1:ys) {
  if (i == 1){
    plot(ode_example$Y[,1], col="darkred", type = "l", 
         ylim = c(0, max(ode_example$Y)),
         xlab = 'Times', ylab = 'Y')
  }
  else{
    lines(ode_example$Y[,i], col="darkred")
  }
}

## for s

ss = ncol(ode_example$dS)

for (i in 1:ss) {
  if (i == 1){
    plot(ode_example$dS[,1], col="skyblue", type = "l", 
         ylim = c(0, max(ode_example$dS)),
         xlab = 'Times', ylab = 'S')
  }
  else{
    lines(ode_example$dS[,i], col="skyblue")
  }
}


## for e

es = ncol(ode_example$dE)

for (i in 1:es) {
  if (i == 1){
    plot(ode_example$dE[,1], col="yellow3", type = "l", 
         ylim = c(0, max(ode_example$dE)),
         xlab = 'Times', ylab = 'E')
  }
  else{
    lines(ode_example$dE[,i], col="yellow3")
  }
}


## for i

is = ncol(ode_example$dI)

for (i in 1:is) {
  if (i == 1){
    plot(ode_example$dI[,1], col="pink3", type = "l", 
         ylim = c(0, max(ode_example$dI)),
         xlab = 'Times', ylab = 'I')
  }
  else{
    lines(ode_example$dI[,i], col="pink3")
  }
}


## for R

rs = ncol(ode_example$dR)

for (i in 1:rs) {
  if (i == 1){
    plot(ode_example$dR[,1], col="darkgreen", type = "l", 
         ylim = c(0, max(ode_example$dR)),
         xlab = 'Times', ylab = 'R')
  }
  else{
    lines(ode_example$dR[,i], col="darkgreen")
  }
}



## for all

## for s


lims = c(0, max(max(ode_example$dS), max(ode_example$dE), 
           max(ode_example$dI), max(ode_example$dR)))

ss = ncol(ode_example$dS)

for (i in 1:ss) {
  if (i == 1){
    plot(ode_example$dS[,1], col="skyblue", type = "l", 
         ylim = lims,
         xlab = 'Times', ylab = 'S')
  }
  else{
    lines(ode_example$dS[,i], col="skyblue")
  }
  lines(ode_example$dE[,i], col="yellow")
  lines(ode_example$dI[,i], col="pink")
  lines(ode_example$dR[,i], col="green")
}



