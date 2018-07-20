library(splines);

#make matrix for basis spline with 5 knots
Fknots <- 5
bs2 <- t(matrix(bs(1:numYr,df=Fknots,intercept=T),ncol=Fknots)) 

#make matrices for storing rsults
LEN <- scwts <- sswts <- matrix(NA, nrow=numAges,ncol=numYr)

####################
# PHASE 1
####################

#make nll for smoother, takes parameters and stock object
nll_phase1 <- function(par,stock){
  log_Linf  <- exp(par[1:5]) %*% bs2;

  for (t in 1:numYr){
    for (a in 1:numAges){
      LEN[a,t] <<-  exp(log_Linf[t])*(1-exp(-exp(par[6])*(a + par[7])))
    }
  }  
  scwts <<-  0.00001508 * LEN ^ 3.090 #1986 cefas report on literature   
  sswts <<-  scwts * exp(par[8])

  nllCWT  <- -sum(dnorm(landings.wt(stock), scwts,exp(par[9]),T),na.rm=T)
  nllSWT  <- -sum(dnorm(stock.wt(stock), sswts,exp(par[10]),T),na.rm=T)
                                  
  return(nllCWT+nllSWT)  
}


####################
# PHASE 2
####################

#make nll for smoother, takes parameters and stock object
nll_phase2 <- function(par,stock){
  log_Linf  <- exp(par[1:5]) %*% bs2;
  
  for (t in 1:numYr){
    for (a in 1:numAges){
      LEN[a,t] <<-  exp(log_Linf[t])*(1-exp(-exp(par[6])*(a + par[7])))
    }
  }  
  scwts <<-  0.00001508 * LEN ^ 3.090 #1986 cefas report on literature   
  sswts <<-  scwts * exp(par[8])
  
  nllCWT  <- -sum(dnorm(landings.wt(stock), scwts,exp(par[9:(8+numAges)]),T),na.rm=T)
  nllSWT  <- -sum(dnorm(stock.wt(stock), sswts,exp(par[(9+numAges):(8+2*numAges)]),T),na.rm=T)
  
  return(nllCWT+nllSWT)  
}

