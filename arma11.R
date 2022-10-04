#Simulate an Arma(1,1) time series data
Arma11<-arima.sim(n=1000,model=list(order=c(1,0,1),ar=0.7,ma=0.4))
Arma11
#Plot time series data
#ts.plot(Arma11)  
#Calculate the Sample Autocorrelation Function
#Arma11.acf<-acf(Arma11,type="correlation",plot=T)
#Arma11.acf

#Simulate an Arma(1,1) time series data
Arma11<-arima.sim(n=1000,model=list(order=c(1,0,1),ar=0.7,ma=-0.4))
Arma11
#Plot time series data
#ts.plot(Arma11)  
#Calculate the Sample Autocorrelation Function
#Arma11.acf<-acf(Arma11,type="correlation",plot=T)
#Arma11.acf


MA.1<-arima.sim(model=list(ma=c(0.5,0.4)), n=100)
#acf(MA.1)
#acf(MA.1, plot=F)

MA.2<-arima.sim(model=list(ma=c(1.2,-0.7)), n=100)
#acf(MA.2)
#acf(MA.2, plot=F)

MA.3<-arima.sim(model=list(ma=c(-1,-0.6)), n=100)
#acf(MA.3)
#acf(MA.3, plot=F)

AR.1<-arima.sim(model=list(ar=c(0.6,0.3)), n=100)
#acf(AR.1)

AR.2<-arima.sim(model=list(ar=c(-0.4,0.5)), n=100)
#acf(AR.2)
