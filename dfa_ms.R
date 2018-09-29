# This function computes mean-square DFA-solutions
# L is the length of the MA filter,
# periodogram is the frequency weighting function in the DFA
# Gamma is the transferfunction of the symmetric filter (target) and
# Lag is the lag-parameter: Lag=0 implies real-time filtering, Lag=L/2
#     implies symmetric filter
# The function returns optimal coefficients as well as the transfer function of the
#     optimized real-time filter
dfa_ms<-function(L,periodogram,Lag,Gamma)
{
  
  K<-length(periodogram)-1
  X<-exp(-1.i*Lag*pi*(0:(K))/(K))*rep(1,K+1)*sqrt(periodogram)
  X_y<-exp(-1.i*Lag*pi*(0:(K))/(K))*rep(1,K+1)
  for (l in 2:L)          #l<-L<-21
  {
    X<-cbind(X,(cos((l-1-Lag)*pi*(0:(K))/(K))+
                  1.i*sin((l-1-Lag)*pi*(0:(K))/(K)))*sqrt(periodogram))
    X_y<-cbind(X_y,(cos((l-1-Lag)*pi*(0:(K))/(K))+
                      1.i*sin((l-1-Lag)*pi*(0:(K))/(K))))
  }
  xtx<-t(Re(X))%*%Re(X)+t(Im(X))%*%Im(X)
  # MA-Filtercoefficients
  b<-as.vector(solve(xtx)%*%(t(Re(X_y))%*%(Gamma*periodogram)))
  # Transferfunction
  trffkt<-1:(K+1)
  trffkt[1]<-sum(b)
  for (k in 1:(K))#k<-1
  {
    trffkt[k+1]<-(b%*%exp(1.i*k*(0:(length(b)-1))*pi/(K)))
  }
  return(list(b=b,trffkt=trffkt))
}


