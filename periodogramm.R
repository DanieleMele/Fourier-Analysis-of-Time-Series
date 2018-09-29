per<-function(x,plot_T)
{
  len<-length(x)
  per<-0:(len/2)
  DFT<-per
  
  for (k in 0:(len/2))
  {
    cexp <- complex(arg=-(1:len)*2*pi*k/len)
    DFT[k+1]<-sum(cexp*x*sqrt(1/(2*pi*len)))
  }
  per<-abs(DFT)^2
  if (plot_T)
  {
    par(mfrow=c(2,1))
    plot(per,type="l",axes=F,xlab="Frequency",ylab="Periodogram",
         main="Periodogram")
    axis(1,at=1+0:6*len/12,labels=c("0","pi/6","2pi/6","3pi/6",
                                    "4pi/6","5pi/6","pi"))
    axis(2)
    box()
    plot(log(per),type="l",axes=F,xlab="Frequency",ylab="Log-periodogram",
         main="Log-periodogram")
    axis(1,at=1+0:6*len/12,labels=c("0","pi/6","2pi/6","3pi/6",
                                    "4pi/6","5pi/6","pi"))
    axis(2)
    box()
  }
  return(list(DFT=DFT,per=per))
}

