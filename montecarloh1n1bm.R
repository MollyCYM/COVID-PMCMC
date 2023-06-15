

rm(list=ls())
require(deSolve)

gen_ode_model <- function(rep_time = 200, ori_seed = 2023, tims = 365, 
                          n = 52196381, pars = c(k=1.59, gamma=1.08),
                          ini_state = c(t=1,S=N-1, E=1, I=0, R=0),
                          tau = 0.1){
  
  times = 1:tims
  N = n
  sigma <- 0.07
  
  H1N1 <- function(time, current_state, params){
    
    with(as.list(c(current_state, params)),{
      
      dt <- 1
      dS <- -beta[t]*S*I/N
      dE <- beta[t]*S*I/N - E/k
      dI <- E/k - I/gamma
      dR <- I/gamma
      
      return(list(c(dt, dS, dE, dI, dR)))
    })
  }
  
  Zs = setNames(data.frame(matrix(NA, nrow = tims, ncol = rep_time)), paste0("Z", 1:rep_time))
  Ys = setNames(Zs, paste0("Y", 1:rep_time))
  dSs = setNames(Zs, paste0("dS", 1:rep_time))
  dEs = setNames(Zs, paste0("dE", 1:rep_time))
  dIs = setNames(Zs, paste0("dI", 1:rep_time))
  dRs = setNames(Zs, paste0("dR", 1:rep_time)) 
  
  params = pars   
  initial_state = ini_state
  tau1 = tau
  
  
  ## rep 200 times 
  for (rep_i in 1:rep_time) {
    set.seed(ori_seed + 1000*rep_i)
    
    e = rnorm(n = length(times)-1 , sd = 1)
    ## now compute their cumulative sum
    e = c(0, cumsum(e))
    beta = exp(sigma*e)
    model1 <- ode(initial_state, times, H1N1, params)
    
    Z1 <-model1[,4]/1.59
    
    Y1 = vector(length = tims)
    
    for (i in 1:tims){
      Y1[i]<- rlnorm(1,log(Z1[i]/5),tau1)
    }
    
    dSi = model1[,3]
    dEi = model1[,4]
    dIi = model1[,5]
    dRi = model1[,6]
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

## simulation

library('truncnorm')
R0 <- rtruncnorm(1, a=0, b=1, mean = 0.15, sd = 0.15)
E0 <-runif(1,-16, -9)
I0 <-runif(1,-16, -9)
S <- 52196381
R <- R0*S
S <- S - R

E <- exp(E0 + log(S))
S <- S - E
I <- exp(I0 + log(S))
S <- S - I

initial_state <- c(t=1, S=S, E=E, I=I, R=R)

ode_example <- gen_ode_model(rep_time = 200, ini_state = initial_state)






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







