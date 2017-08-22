# --------------------------------------------------------------------------------------
# Functions for the diagnostics of the North Sea plaice assessment
#
# Author  : Jan Jaap Poos
#           (+David Miller)
#
# Last edited: APRIL 2011
# --------------------------------------------------------------------------------------
########---------------------------+INDEX+------------------------------########
## DATA
#  data_list: Creates data frame with list of values (good for xyplot)
#  totalStk: Adds up totals for stock objects - computes landings, discards, catch, stock totals and sets units

## Assessment (XSA)
#  diagnostics: Runs XSA diagnostics
#  diagnosticsCEFAS: Runs (CEFAS) XSA diagnostics             #diagnosticsCEFAS(xsa.stock, xsa.indices, xsa.control)
#  survDiagCEFAS: Runs (CEFAS) XSA survivors diagnostics      #survDiagCEFAS(xsa.stock, xsa.res, xsa.control)
#  stock_Sum: Creates a stock summary tables for ICES WG report
#  write_Stock: Outputs a .txt file of slots from an FLR stock object
#  write_Indices: Outputs a .txt file of slots from an FLR indices object

## STF
#  scanSTF: Runs STFs for a range of Fmult values
#  summTableSTF: Creates a summary from 'scanSTF' outputs
#  inputTableSTF: Creates a table of inputs to the STF
#  ageTableSTF: Creates a table by year AND AGE for the projection period 

# FLSAM with pin turned on
FLSAM_pin <- function (stck, tun, ctrl, run.dir = tempdir(), batch.mode = FALSE) 
{
  inputSAM <- new("FLSAMinput")
  inputSAM@stck <- stck
  inputSAM@tun <- tun
  inputSAM@ctrl <- ctrl
  if (any(c(validObject(stck), validObject(tun), validObject(ctrl), 
            validObject(inputSAM)) == F)) 
    stop("Validity check failed")
  FLR2SAM(stck, tun, ctrl, run.dir)
  rtn <- runSAM(ctrl, run.dir, use.pin=T)
  if (rtn != 0) {
    if (batch.mode) {
      return(NULL)
    }
    else {
      stop(sprintf("An error occurred while running ADMB. Return code %s.", 
                   rtn))
    }
  }
  res <- SAM2FLR(ctrl, run.dir)
  return(res)
}


########--------------------------+++++++++-----------------------------########
####                                 INDICES                                   ####
########--------------------------+++++++++-----------------------------########
cor.tun.save <- function(stk.tun, pltInd=c(1), saveDir = getwd()){ 
  #for(i in names(stk.tun)) {
  for(i in pltInd) {
      if(dim(stk.tun[[i]]@index)[1]>1) {
      x11()
      plot(stk.tun[[i]],type="internal",main=name(stk.tun[[i]]))
      savePlot(filename=paste(saveDir,"\\Indices_intCor_",i, sep=""), type=c("wmf"), device=dev.cur())
      dev.off()
     }    
    }
  }              

########--------------------------+++++++++-----------------------------########
####                                 DATA                                   ####
########--------------------------+++++++++-----------------------------########
# NAME: totalStk
# DOES: Adds up totals for stock objects - computes landings, discards, catch, stock totals and sets units
totalStk <- function(stk, Units){
    landings(stk) <- computeLandings(stk)
    discards(stk) <- computeDiscards(stk)
    catch.n(stk)  <- landings.n(stk)+discards.n(stk)
    catch.wt(stk) <-(landings.n(stk)*landings.wt(stk)+discards.n(stk)*discards.wt(stk))/catch.n(stk)
    catch(stk)    <- computeCatch(stk)
    stock(stk)    <- computeStock(stk)
    units(stk)[1:17] <- as.list(c(rep(Units,4), "NA", "NA", "f", "NA", "NA"))
    return(stk)
    }

########--------------------------+++++++++-----------------------------########
# NAME: data_list
# DOES: Creates data frame with list of values (good for xyplot)
setGeneric("data_list", function(object, ...){
	standardGeneric("data_list")
	}
)
setMethod("data_list", signature("list"), function(object, ...){
	# names
	if(!is.null(names(object))){
		flqnames <- names(object)
	} else {
		flqnames <- paste("v", 1:length(object), sep="")
	}

	# data.frames
	flqs.lst <- lapply(object, as.data.frame)
	flqs.lst <- lapply(flqs.lst, function(x){x[,1] <- as.character(x[,1]); x})
	flqs.nlst <- lapply(flqs.lst, nrow)
	flqs.df <- do.call("rbind", flqs.lst)
	flqs.df$qname <- rep(flqnames, unlist(flqs.nlst))
	flqs.df

})  # }}}

########--------------------------+++++++++-----------------------------########
####                                 XSA                                    ####
########--------------------------+++++++++-----------------------------########
# NAME: diagnostics
# DOES: Runs XSA diagnostics
# Needs to be checked against original XSA program outputs (do a direct comparison!!!)
# INCORRECTLY CALCULATES WEIGHT OF INDICES ON RECENT YEAR CLASSES
diagnostics <- function(object){

    indices<-new("FLIndices")
    for (i in 1:length(object@index))
        {
        indices[[i]] <- FLIndex(index=object@index[[i]])
        indices[[i]]@name  <-object@index.name[i]
        indices[[i]]@range <-object@index.range[[i]]
        }
    control <-object@control #<- eval(parse(text=xsa@call[4]))

    cat("FLR XSA Diagnostics ",  as.character(Sys.time()),"\n\nCPUE data from ", object@call[3], "\n\n",
        "Catch data for ", dims(object@stock.n)$year," years. ", dims(object@stock.n)$minyear," to ",
        dims(object@stock.n)$maxyear, ". Ages ",dims(object@stock.n)$min," to ",dims(object@stock.n)$max,".\n\n", sep="")

    # print general tuning series information
    idx.info  <- NULL
    for (i in 1:length(object@index)) {
       idx.info  <-  rbind(idx.info,c(indices[[i]]@name, (dims(object@index[[i]]))$min,
          (dims(object@index[[i]]))$max,(dims(object@index[[i]]))$minyear,(dims(object@index[[i]]))$maxyear,
          indices[[i]]@range["startf"], indices[[i]]@range["endf"]))
    }

    dimnames(idx.info) <- list(NULL,c("fleet","first age","last age","first year","last year","alpha","beta"))
    print(as.data.frame(idx.info))
    cat("\n\n","Time series weights :\n\n")
    cat(ifelse(control@tsrange==0|control@tspower==0, "   Tapered time weighting not applied\n\n",
        paste("   Tapered time weighting applied\n", "  Power =  ",control@tspower,"over ",control@tsrange,
        "years\n\n", sep=" ")))

    cat("Catchability analysis :\n\n")
    cat(ifelse(as.numeric(control@rage) < dims(object@stock.n)$min, "    Catchability independent of size for all ages\n\n",
        paste("    Catchability independent of size for ages >  ",control@rage,"\n\n",sep=" ")))
    cat(ifelse(as.numeric(control@qage) < dims(object@stock.n)$min, "    Catchability independent of age for all ages\n\n",
        paste("    Catchability independent of age for ages >=  ",control@qage,"\n\n",sep=" ")))

    cat("Terminal population estimation :\n\n")
    cat(ifelse(control@shk.f, paste("    Survivor estimates shrunk towards the mean F\n",
        "   of the final  ",control@shk.yrs,"years or the ",control@shk.ages,"oldest ages.\n\n",
        "   S.E. of the mean to which the estimates are shrunk =  ", control@fse,"\n",sep=" "),
        "    Final estimates not shrunk towards mean F\n"))
    cat(ifelse(as.numeric(control@min.nse)==0, "\n", paste("\n", "   Minimum standard error for population\n",
        "   estimates derived from each fleet = ",control@min.nse,"\n\n", sep=" ")))
    cat("   prior weighting not applied\n\n")

    cat("Regression weights\n")
    ### Calculation of time series weighting
    yr.range <- (dims(object@harvest)$maxyear-9):dims(object@harvest)$maxyear
    regWt <- FLQuant(dimnames=list(age = 'all', year = yr.range))
    for(y in yr.range) regWt[,as.character(y)] <- (1-((max(yr.range)-y)/control@tsrange)^control@tspower)^control@tspower
    print(matrix(round(regWt,3),dims(regWt)$age,dimnames=list(age="all",year=yr.range)))

    cat("\n\n Fishing mortalities\n")
    print(matrix(round(trim(object@harvest,year=yr.range),3), dims(object@harvest)$age,
        dimnames=list(age=dims(object@harvest)$min:dims(object@harvest)$max, year=yr.range)))

    cat("\n\n XSA population number (Thousand)\n")
    print(t(matrix(round(trim(object@stock.n,year=yr.range),0), dims(object@stock.n)$age,
        dimnames=list(age=dims(object@stock.n)$min:dims(object@stock.n)$max, year=yr.range))))

    nextyear  <- dims(object@survivors)$maxyear
    cat("\n\n Estimated population abundance at 1st Jan ",nextyear,"\n")
    print(t(matrix(round(object@survivors[,as.character(nextyear)]),
        dimnames=list(age=dims(object@survivors)$min:dims(object@survivors)$max, year=nextyear))))

    ## tuning info
    for (f in 1:length(object@index)) {
        cat("\n\n Fleet: ",indices[[f]]@name,"\n\n","Log catchability residuals.\n\n")

        print(matrix(round(object@index.res[[f]],3), nrow=dims(object@index.res[[f]])$age,
            dimnames=list(age=dimnames(object@index.res[[f]])$age, year=dimnames(object@index.res[[f]])$year)))

        if (control@rage < dims(object@index[[f]])$max){
          cat("\n\n Mean log catchability and standard error of ages with catchability \n",
              "independent of year class strength and constant w.r.t. time \n\n")

          q.tab <- rbind(Mean_Logq=round(log(object@q.hat[[f]]),4), S.E_Logq=round(sd(matrix(object@index.res[[f]],
              dim(object@index.res[[f]])[2],dim(object@index.res[[f]])[1],byrow=T),na.rm=T),4))
          colnames(q.tab) <- dimnames(object@q.hat[[f]])$age

          if (dims(object@index[[f]])$min <= control@rage ) {
              print(q.tab[,as.character((control@rage+1):max(as.numeric(dimnames(object@index[[f]])$age)))])
          } else {print(q.tab)}
        }
        # print reg stats powermodel, note that maximum printed age is min(rage, max(age in tun series))
        if (dims(object@index[[f]])$min <= control@rage ) {
            cat("\n Regression statistics \n", "Ages with q dependent on year class strength \n")
            print(cbind((matrix(object@q2.hat[[f]][as.character(dims(object@index[[f]])$min:(min(control@rage,dims(object@index[[f]])$max)))],dimnames=list(age=paste("Age ",dims(object@index[[f]])$min:(min(control@rage,dims(object@index[[f]])$max)),sep=""),"slope"))),
            (matrix(object@q.hat[[f]][as.character(dims(object@index[[f]])$min:(min(control@rage,dims(object@index[[f]])$max)))],dimnames=list(age=dims(object@index[[f]])$min:(min(control@rage,dims(object@index[[f]])$max)),"intercept")))))
        }
    }

    cat("\n\n Terminal year survivor and F summaries: \n ")
    for ( age in sort(unique(object@diagnostics$age))){
        cat("\n Age ",age, " Year class =",  max(object@diagnostics$year) - age ," \n\n","source \n", sep="")
        weights <- object@diagnostics[(object@diagnostics$age==age) & (object@diagnostics$year== max(object@diagnostics$year)),]
        # calc surivors and scaled wts
        weights$survivors <- round(exp(weights$nhat))
        weights$scaledWts <- round(weights$w / sum(weights$w) ,3)
        row.names(weights) <- weights$source
        print(weights[ ,c("scaledWts","survivors","yrcls") ])
    }
    invisible()
}

########--------------------------+++++++++-----------------------------########
# NAME: diagnosticsCEFAS
# DOES: Runs (CEFAS) XSA diagnostics
# Needs to be checked against original XSA program outputs (do a direct comparison!!!)
diagnosticsCEFAS <- function(stock.obj, indices.obj,xsa.cntrl, ...){
     if(!validObject(stock.obj)) {stop("Invalid FLStock object")}
     if(!validObject(indices.obj)) {stop("Invalid FLIndices object")}
     pwer <- xsa.cntrl@rage
     q.plat <- xsa.cntrl@qage
     plot.count <- 0
     windows()      
     par(mfrow=c(2,2))
     
     for (f in c(1:length(indices.obj))){
         minYr <- max(indices.obj[[f]]@range[4],as.numeric(min(attributes(stock.obj@stock.n)$dimnames$year)))
         maxYr <- min(indices.obj[[f]]@range[5],as.numeric(max(attributes(stock.obj@stock.n)$dimnames$year)))
         minAge <- max(indices.obj[[f]]@range[1],as.numeric(min(attributes(stock.obj@stock.n)$dimnames$age)))
         maxAge <-  min(indices.obj[[f]]@range[2],as.numeric(max(attributes(stock.obj@stock.n)$dimnames$age)))

         fleet <- trim(indices.obj[[f]]@index,age=minAge:maxAge,year=minYr:maxYr)
         vpa <- trim(stock.obj@stock.n,age=minAge:maxAge,year=minYr:maxYr)
         q. <-  log(fleet)/log(vpa)
         
         ###### First Transform fleet data to begining of year #######
         alpha. <- FLQuant(indices.obj[[f]]@range[6],dim=c(length(minAge:maxAge), length(minYr:maxYr)))
         beta.  <- FLQuant(indices.obj[[f]]@range[7],dim=c(length(minAge:maxAge), length(minYr:maxYr)))
         Z.     <- trim(stock.obj@harvest,age=minAge:maxAge, year=minYr:maxYr) + trim(stock.obj@m,age=minAge:maxAge,year=minYr:maxYr)
         Raise. <- (exp(-alpha.*Z.) - exp(-beta.*Z.))/((beta.-alpha.)*Z.)
         fleet.adjusted <- fleet/Raise.         
       
         #windows()
         #par(mfrow=c(4,3))
         count<-0
         reg.sts.output <- NULL
         q.resids       <- NULL
         mean.qs        <- NULL
         SE.logqs       <- NULL 
         for (a in c(minAge:maxAge)){
              count <- count+1
     
              ###### Get required data for loop
              vpa.a1 <- as.vector(trim(log(vpa),age=a))              
              fleet.a1 <- as.vector(trim(log(fleet.adjusted),age=a))
              model.dat <- rbind(vpa.a1,fleet.a1) 
              model.dat2 <- model.dat[,which(is.na(model.dat[2,])==F)]
              vpa.a   <-  model.dat2[1,]
              fleet.a <-  model.dat2[2,]
              
              ##### perform calibration regression , NOTE - only implemented in model for ages with power model, otherwise uses mean.
              model <- lm(fleet.a~vpa.a )
              transpose.model <- as.data.frame(cbind("fleet" = fleet.a,                                                                
                                       "pred.VPA" = (fleet.a/summary(model)$coefficients[2,1])-(summary(model)$coefficients[1,1]/summary(model)$coefficients[2,1])) )
                  
              model2 <- lm(transpose.model$pred.VPA~ fleet.a)                         
         
              fleet.mean <-   mean(fleet.a)
              vpa.mean <- mean(vpa.a)                                                                    
              Sxx = sum((fleet.a-fleet.mean)^2)
              Sxy = sum((fleet.a-fleet.mean)*(vpa.a-vpa.mean))
              Syy = sum((vpa.a-vpa.mean)^2)
              Gradient =  summary(model2)$coefficients[2,1]
              Intercept = summary(model2)$coefficients[1,1]    
              
              ye <- vpa.a - transpose.model$pred.VPA
              ye2 <- ye ^ 2
              Se <- sqrt(sum(ye2)/(length(transpose.model$pred.VPA) - 2))        
        
              reg.sts <- cbind("Age"       = a,
                              "Model used?" = ifelse(a <= pwer,"Yes","No"),
                              "slope"     = round(Gradient,2),
                              #"t-value"   = round(Gradient/Se,2),
                              "Intercept" = round(Intercept,2),
                              "RSquare"   = round(Sxy*Sxy/(Sxx*Syy),2),
                              "Num Pts"   = length(model$model[,1]),
                              "Reg s.e"   = round(Se,2),
                              "Mean Q"    = round(mean(log(exp(fleet.a)/exp(vpa.a))),2))
              reg.sts.output <- rbind(reg.sts.output,reg.sts)
                 
              if (a <= pwer){
                  q.pwr   <-  log(exp(fleet.a)/exp(vpa.a)) 
                  q.pwr2  <-  log(exp(fleet.a)/exp(transpose.model$pred.VPA))
                  q.resid <-  q.pwr-q.pwr2
                  missing. <-   which(is.na(model.dat[2,]))
                  for (its in missing.){
                       q.resid <- round(c(q.resid[1:(its-1)],NA,q.resid[its:length(q.resid)]),2)
                       }                       
                  mean.q  <- "NA"
                  SE.logq <- "NA"          
                  }
              if (a > pwer & a <= q.plat){
                  q.   <-  log(exp(fleet.a1)/exp(vpa.a1))      
                  mean.q  <-  round(mean(q.[is.finite(q.)]),4)                          
                  q.resid <-  round(q. - mean.q,2)  
                  SE.logq <-  round(sqrt(sum(q.resid[is.finite(q.resid)]^2)/(length(q.resid[is.finite(q.resid)])-1)),4)                                                                   
                  }
              if (a > q.plat){
                  vpa.a.qhat <- as.vector(trim(log(vpa),age=q.plat))              
                  fleet.a.qhat <- as.vector(trim(log(fleet.adjusted),age=q.plat))
                  
                  q.qhat   <-  log(exp(fleet.a.qhat)/exp(vpa.a.qhat))      
                  mean.q.qhat  <-  round(mean(q.qhat[is.finite(q.qhat)]),4)
                                          
                  vpa.a. <- as.vector(trim(log(vpa),age=a))              
                  fleet.a. <- as.vector(trim(log(fleet.adjusted),age=a))
                  q.    <-  log(exp(fleet.a.)/exp(vpa.a.)) 
                                              
                  q.resid <-  round(q. - mean.q.qhat,2)  
                  SE.logq <-  round(sqrt(sum(q.resid[is.finite(q.resid)]^2)/(length(q.resid[is.finite(q.resid)])-1)),4)                                               
                  }                                      
                  q.resids <- rbind(q.resids,q.resid)                                 
                  mean.qs  <- cbind(mean.qs,mean.q)
                  SE.logqs  <- cbind(SE.logqs,SE.logq)
            }
         ###making tables etc pretty (ish)    
         rownames(mean.qs) <- "Mean log q"
         rownames(SE.logqs) <- "S.E. log q"
         log.qs <- rbind(mean.qs,SE.logqs)
         colnames(log.qs) <- c(minAge:maxAge)
     
         rownames(q.resids) <- c(minAge:maxAge)
         colnames(q.resids) <- c(minYr:maxYr)
     
         rownames(reg.sts.output) <- reg.sts.output[,1]
         reg.sts.output2 <- reg.sts.output[,-1]
     
         ###printing to screen to mimic XSA output
         cat("Fleet = ", indices.obj@names[f],"\n","\n","Catchability residuals:","\n","\n")
         print(round(q.resids,2))     
         cat("\n","\n","Mean log catchability and standard error of ages with", "\n", "independant of year class strength and constant w.r.t time:","\n","\n")
         print(log.qs)
         cat("\n","\n","Regression Statistics:","\n","\n")
         print(reg.sts.output2)
         cat("\n","\n","\n","\n") 
     
         ############# produce mean log q plot #################
         plot.count <- plot.count+1
         if(plot.count ==5){
                            windows()
                            par(mfrow=c(2,2))                            
                            plot.count <- 1
                            }
         plot.dat <- log.qs[,colnames(log.qs) > pwer]
         if(length(plot.dat) > 2){
         plot.dat2 <- exp(as.numeric(plot.dat[1,])+(2*as.numeric(plot.dat[2,])))
         plot.dat3 <- exp(as.numeric(plot.dat[1,])-(2*as.numeric(plot.dat[2,])))
         plot.dat4 <- exp(as.numeric(reg.sts.output2[as.numeric(rownames(reg.sts.output2)) > pwer,7]))
         plot(colnames(plot.dat), plot.dat4, type="l", lwd=2, xlab="Age", ylab="Catchability (+/- one s.d.)", 
              ylim=c(min(plot.dat3 ),max(plot.dat2)),xlim=c(0,maxAge),col="red",
              main = paste("Mean catchability of ",indices.obj@names[f],"\nblack line includes q plateau"))
         lines(colnames(plot.dat),plot.dat2, lty=2)
         lines(colnames(plot.dat),plot.dat3, lty=2)
         lines(colnames(plot.dat),exp(as.numeric(plot.dat[1,])), lty=1, col="black", lwd=2)
         }else{
              plot(c(0:1), c(0:1), type="n", lwd=2, xlab="Age", ylab="Catchability (+/- one s.d.)", 
              main = paste("Mean catchability of ",indices.obj@names[f],"\nblack line includes q plateau"))
              text(0.5,0.5, "NOT ENOUGH DATA TO PLOT")
              }
     }
     
}

########--------------------------+++++++++-----------------------------########
# NAME: survDiagCEFAS
# DOES: Runs (CEFAS) XSA survivors diagnostics
# Needs to be checked against original XSA program outputs (do a direct comparison!!!)
survDiagCEFAS <- function(stock.obj, xsa.object,xsa.cntrl,surv.plots=T,...){
    ow <- options("warn")
    options(warn=-1)
    
    xsa.output <- xsa.object@diagnostics
    minYr <- min(unlist(lapply(xsa.object@index.range, function(x) return(x[4]))))
    maxYr <- max(unlist(lapply(xsa.object@index.range, function(x) return(x[5]))))
    maxYrcl <- max(unlist(lapply(xsa.object@index.range, function(x) return(x[5]))))- min(xsa.object@diagnostics$age)
    minYrcl <- max(unlist(lapply(xsa.object@index.range, function(x) return(x[5]))))- max(xsa.object@diagnostics$age)    
    minAge <- min(xsa.object@diagnostics$age)
    maxAge <- max(xsa.object@diagnostics$age) 
    fl.surv.matrix <- fl.wgts.matrix <- matrix(NA,ncol=length(unique(xsa.output$source)),nrow=length(minAge:maxAge),dimnames=list(c(minAge:maxAge),c(sort(unique(xsa.output$source))))) 

    
    fleets.n <-   sort(unique(xsa.output$source)) 
    Surv. <- cbind("Fleet"=rep(fleets.n,2),"Stat"=c(rep(colnames(xsa.output)[1],length(fleets.n)),rep(colnames(xsa.output)[2],length(fleets.n))))
    Surv.4 <- NULL   
    final.surv.tot <- NULL
    
    for (yrcls. in c(maxYrcl:minYrcl)){         
         req.yrcls <-  xsa.output[xsa.output$yrcls == yrcls.,]
         age.for.output <-  max(req.yrcls$age)
         pwer <- xsa.cntrl@rage
         title. <- ifelse(age.for.output > pwer, "Catchability constand w.r.t. time and dependant on age","Catchability dependant on age and yearclass strength")                                                  
         cat("\n\nAge = ",age.for.output,". ",title.)
         cat("\n","Year class = ", yrcls.,"\n\n")
         
         fleet.out2 <- NULL
         ##############                                      ##############
         ############## To get Estimated survivors per fleet ##############
         ##############                                      ##############        
         for (fl. in sort(unique(req.yrcls$source))){
              fleet.dat <- req.yrcls[req.yrcls$source==fl.,]
              est.surv <- round(exp(sum(fleet.dat$nhat*fleet.dat$w)/sum(fleet.dat$w)))
              ext.se.1   <- (fleet.dat$nhat-log(est.surv))^2
              ext.se.2   <- sqrt((sum(ext.se.1*fleet.dat$w)/sum(fleet.dat$w)))
              ext.se     <- round(sqrt((1/(dim(fleet.dat)[1] - 1)))*ext.se.2 , 3)
              n          <- dim(fleet.dat)[1]
        
              ###### int.err is a bit more complicated as need to calculated ECF ########
              ECF.1 <- FLQuant(NA, dim=dim(xsa.object@harvest), dimnames=dimnames(xsa.object@harvest))
              ECF.1[dim(ECF.1)[1],] <- xsa.object@harvest[dim(ECF.1)[1],]
              ECF.1[,dim(ECF.1)[2]] <- xsa.object@harvest[,dim(ECF.1)[2]]
              for (cl. in c((dim(ECF.1)[2]-1):1)){
                   for (rw. in c((dim(ECF.1)[1]-1):1)){
                        ECF.1[rw.,cl.] <- xsa.object@harvest[rw.,cl.] + ECF.1[rw.+1,cl.+1]                        
                        }
                   }
              ECF.2 <- exp(ECF.1)
              w2.b <- NULL
              for (it.1 in c(1:dim(fleet.dat)[1])){
                        w2.a<- ECF.2[which(rownames(ECF.2)==fleet.dat$age[it.1]),which(colnames(ECF.2)==fleet.dat$year[it.1])]                                               
                        w2.b <- cbind(w2.b,w2.a) 
                        }
                   
              int.se <- round(sqrt(sum(fleet.dat$w * (1/w2.b))/sum(fleet.dat$w)^2),3)
              var.ratio <- round(ext.se/int.se,3)
              scaled.wgts <- round(sum(fleet.dat$w)/sum(req.yrcls$w),3)
              
              ####### Use Popes approx to get N, then use that to calc Z... hence F
              m. <- stock.obj@m[ac(age.for.output),ac(maxYr),,,,] 
              c. <- stock.obj@catch.n[ac(age.for.output),ac(maxYr),,,,] 
              N. <- (est.surv*exp(m.))+(c.*exp(m./2))
              est.f <- round(log(N./est.surv)-m.,3)
              
              fleet.out <- cbind("Fleet"       = fl.,
                                 "Est.Suvivors"= est.surv,
                                 "Int. s.e."   = int.se,
                                 "Ext. s.e."   = ext.se,
                                 "Var Ratio"   = var.ratio,
                                 "N"           = n,
                                 "Scaled Wgts" = scaled.wgts,
                                 "Estimated F" = est.f)
              fleet.out2 <- rbind(fleet.out2,fleet.out)
              
              
              ######## Tuck away some results to use later for Fleet surv ratios and wgts plot
              fl.surv.matrix[maxYr-yrcls. ,fl.] <-  est.surv
              fl.wgts.matrix[maxYr-yrcls. ,fl.] <-  scaled.wgts          
           
              ##############                                                ##############
              ##############    To get Survivors down the cohort per fleet  ##############
              ##############                                                ##############        
           
             if(dim(fleet.dat)[1] > 0){
                Surv.2  <- rbind("Survivors" =round(exp(fleet.dat$nhat)) ,
                                 "Raw weights" = c(round(fleet.dat$w,3)))
                colnames(Surv.2) <- fleet.dat$age
                #Surv.3 <- cbind("YearClass"=rep(yrcls.,dim(Surv.2)[2]), "fleet"=rep(fl.,dim(Surv.2)[2]),"age"=colnames(Surv.2), "Survivors"=Surv.2[1,],"Raw weights"=Surv.2[2,])               
                }else{
                      Surv.2  <- rbind("Survivors" = rep(0,dim(fleet.dat)[1]),
                                       "Raw weights" = rep(0,dim(fleet.dat)[1]))
                      colnames(Surv.2) <- fleet.dat$age
                      #Surv.3 <- cbind("YearClass"=rep(yrcls.,dim(Surv.2)[2]), "fleet"=rep(fl.,dim(Surv.2)[2]),"age"=colnames(Surv.2), "Survivors"=Surv.2[1,],"Raw weights"=Surv.2[2,])  
                      }
               #Surv.4 <- rbind(Surv.4,Surv.3)            ####can create and export his if require to get dataframe of surv est down the cohort... currently not set up to work.
               cat("Fleet = ", fl., "\n")
               print(Surv.2) 
               cat("\n\n")
             }                                 
      print(fleet.out2)
      ##############                                                ##############
      ##############    To get final weighted Survivors est.        ##############
      ##############                                                ##############     
      final.surv <- xsa.object@survivors[as.character(age.for.output+1),as.character(maxYr+1),,,,]
      final.int.se <- round(xsa.object@se.int[age.for.output],2)
      final.ext.se <- round(xsa.object@se.ext[age.for.output],2)
      final.n <- dim(Surv.4)[1]
      final.var.ratio <- round(final.ext.se/final.int.se,3)
      final.F <-  round(xsa.object@harvest[as.character(age.for.output),as.character(maxYr),,,,],3)
      fleet.out <- cbind("Suvivors"    = round(final.surv),
                         "Int.s.e."    = "",     #final.int.se,
                         "Ext.s.e."    = "",     #final.ext.se,
                         "Var.Ratio"   = "",     #final.var.ratio,
                         "N"           = final.n,
                         "F"           = final.F)
                         
      ##### store values for plots
      final.surv.tot <- cbind(final.surv.tot,final.surv)
      cat("\n\n","Weighted prediction:","\n\n")
      print(fleet.out)
     }
     
     
##############                                                ##############
##############    To get surv and wgt plots.                  ##############
##############                                                ##############               
windows(12,10)
#par(mfcol=c(1,2))
layout(cbind(rbind(1,2),rbind(3,4)), height=c(3,1))
fl.surv.matrix2 <- sweep(log(fl.surv.matrix),1,log(final.surv.tot),"/")
plot(rownames(fl.surv.matrix2), fl.surv.matrix2[,1],pch=1,ylim=c(0,2),main="Fleet Survivors ratio",xlab="Age",ylab="ln(fleet est)/ln(weighted surv est)")
for(it in c(2:dim(fl.surv.matrix2)[2])){
    lines(rownames(fl.surv.matrix2), fl.surv.matrix2[,it], pch=it,type="p")
    }
plot.new()
par(mar=c(0,0,0,0))
legend("center", legend=colnames(fl.surv.matrix2), pch= c(1:dim(fl.surv.matrix2)[2]),bty="n",ncol=2,title = "Key",cex=1.5)

#fl.wgts.matrix <- as.data.frame(fl.wgts.matrix)
if(sum(is.na(fl.wgts.matrix)) >1){
   fl.wgts.matrix2 <- replace(fl.wgts.matrix,is.na(fl.wgts.matrix),0)
   }else{fl.wgts.matrix2 <- fl.wgts.matrix}
cols. <- colorRampPalette(c("grey5", "grey99"),space="Lab")
cols.2 <- cols.(dim(fl.wgts.matrix2)[1])    
par(mar=c( 5.1, 4.1, 4.1, 2.1))
barplot(t(fl.wgts.matrix2),col=cols.2,xlab="Age",ylab="Proportion of weight",main="Fleet weights")
plot.new()
par(mar=c(0,0,0,0))
legend("center", legend=colnames(fl.surv.matrix2), pch= 15,col=cols.2, bty="n",ncol=2,title = "Key",pt.cex=3,cex=1.5)
      
options(ow)
}  

########--------------------------+++++++++-----------------------------########
# NAME: stock_Sum
# DOES: Creates a stock summary tables for ICES WG report
stock_Sum <- function(object, rage=1, fbar.min=2, fbar.max=6, fbar.dis.max=3, type="missing"){
  if (!inherits(object, "FLStock")) stop("Object must be an 'FLStock' object!")
   d.f <- NULL
      year.  <- dimnames(object@stock.n)[[2]]
      r.     <- as.data.frame(round(object@stock.n[as.character(rage),,,,],0))["data"]
      ssb.   <- as.data.frame(round(ssb(object),0))["data"]
      f.     <- as.data.frame(round(apply(slot(object, "harvest")[fbar.min:fbar.max,], 2:5, mean), 2))["data"]
      f.hc   <- as.data.frame(round(quantMeans((slot(object, "harvest") * ( object@landings.n /object@catch.n) )[as.character(fbar.min:fbar.max),]), 2))["data"]
      f.d    <- as.data.frame(round(quantMeans((slot(object, "harvest") * ( object@discards.n /object@catch.n) )[as.character(fbar.min:fbar.dis.max),]), 2))["data"]
      c.     <- as.data.frame(round(object@catch,0))["data"]
      l.     <- as.data.frame(round(object@landings,0))["data"]
      d.     <- as.data.frame(round(object@discards,0))["data"]
      ydssb. <- as.data.frame(round(object@landings/ssb(object),2))["data"]
      res  <- rbind(d.f, cbind(r., ssb., c., l., d., f., f.hc,f.d, ydssb.))
      colnames(res) <- c("recruits","ssb","catch","landings","discards",paste("fbar",fbar.min,"-",fbar.max,sep=""),paste("fbar hc",fbar.min,"-",fbar.max,sep=""),paste("fbar dis",fbar.min,"-",fbar.dis.max,sep=""),"Y/ssb")
      rownames(res) <- year.
     return(res)

 }

########--------------------------+++++++++-----------------------------########
# NAME: write_Stock
# DOES: Outputs a .txt file of slots from an FLR stock object
write_Stock <- function(number, x, slot, path=getwd(), round=0){
  capture.output(cat("table", number,".", x@name,".", slot, "\n", as.character(Sys.time()), " units=",
    units(slot(x,slot)), "\n"), file=paste(path,"/",number,"_",slot,".txt", sep=""))
  capture.output(matrix(round(slot(x, slot),round), dims(slot(x,slot))$year, byrow=T,dimnames=list(year =
    dims(slot(x, slot))$minyear:dims(slot(x, slot))$maxyear, age=dims(slot(x, slot))$min:dims(slot(x,
    slot))$max)), file=paste(path,"/",number,"_",slot,".txt", sep=""), append=T)
}

########--------------------------+++++++++-----------------------------########
# NAME: write_Indices
# DOES: Outputs a .txt file of slots from an FLR indices object
write_Indices <- function(number, x, name, filename, path=getwd(), round=0){
capture.output(cat("table", number, ".", name, "\n", as.character(Sys.time()), file=paste(path,"/",filename, sep="")))
for (i in 1:length(x)){
    capture.output(print(paste(x[[i]]@name," units=", units(x[[i]]@index)),quote=F), file=paste(path,"/",filename, sep=""),append=T)
    capture.output(print(cbind(x[[i]]@effort,matrix(x[[i]]@index,dims(x[[i]]@index)$year,byrow=T, dimnames=list(year=dims(x[[i]]@index)$minyear:dims(x[[i]]@index)$maxyear, age=dims(x[[i]]@index)$min:dims(x[[i]]@index)$max))), width=5, digits=3), file=paste(path,"/",filename, sep=""), append=T)
    }
}

########--------------------------+++++++++-----------------------------########
####                                 STF                                    ####
########--------------------------+++++++++-----------------------------########
# NAME: scanSTF        
# DOES: Runs STFs for a range of Fmult or F (mean) values  (if fmult=T, fscan used us fmult, else fscan used as absolute levels of F0
# Recoded for R 2.12.2 (D.Miller)
#stfObj <- stf_Object; projControl<- stf_Control; fscan=c(0.8, 1.0, 1.2); SRR=srr
scanSTF <- function(stfObj, projControl="missing", SRR=NA, fscan=c(0.8, 1.0, 1.2), fmult=T){
    if (missing(projControl)) stop("no projControl")
    if (!inherits(stfObj, "FLStock") ) stop("stfObj must be 'FLStock' or 'FLSTF' object")
    
    FSQ <- as.numeric(fbar(stfObj)[,ac(range(stfObj)["maxyear"])])
    
    res <- list()
    rIt <- 0
    for (it in fscan) {
       rIt <- rIt+1
       cat("Fmult: ", it, "\n")
       if (fmult) itVal <- it*FSQ else itVal <- it
       if (length(projControl@target) == 1) {
         projControl@target[1,"val"] <- itVal  
         } else {
         projControl@target[2:length(projControl@target),"val"] <-  itVal      
         }
       tmp <- project(stfObj, projControl, srr)
       res[[paste(rIt,it,sep="_")]] <- window(tmp, start=stfObj@range["maxyear"]-2, end=stfObj@range["maxyear"])
    }
    return(res)
}
########--------------------------+++++++++-----------------------------########
# NAME: summTableSTF
# DOES: Creates a summary from 'scanSTF' outputs
# object <- fwd_Fsq_Options
summTableSTF  <- function(object){
    d.f <- NULL
    for (i in 1:length(object) ) {
        fmult <- as.numeric(strsplit(names(object)[i],"_")[[1]][2])
        fbar.min <- as.numeric(range(object[[i]])["minfbar"])
        fbar.max <- as.numeric(range(object[[i]])["maxfbar"])
        ssb. <- as.data.frame(round(ssb(object[[i]]),0))
        f.   <- as.matrix(c(round(apply(slot(object[[i]], "harvest")[as.character(fbar.min:fbar.max),,,,], 2:5, mean),3)))
        r.   <- as.matrix(c(round(object[[i]]@stock.n[1,],0)))
        c.   <- as.matrix(c(round(object[[i]]@catch,0)))    
        f.d  <- as.matrix(c(round(apply((slot(object[[i]], "harvest") * object[[i]]@discards.n/object[[i]]@catch.n)[as.character(2:3),], 2:5, mean),2)))
        f.l  <- as.matrix(c(round(apply((slot(object[[i]], "harvest") * object[[i]]@landings.n/object[[i]]@catch.n)[as.character(fbar.min:fbar.max),,,,], 2:5, mean),2)))
        l.   <- as.matrix(c(round(apply(object[[i]]@landings.n*object[[i]]@landings.wt,2,sum, na.rm=TRUE),0)))
        d.   <- as.matrix(c(round(apply(object[[i]]@discards.n*object[[i]]@discards.wt,2,sum, na.rm=TRUE),0)))
        d.f  <- rbind(d.f, cbind(fmult, ssb., f., f.d, f.l, r., c., l., d.))
    }
    colnames(d.f) <- c("fmult", "age","year","unit","season","area","iter","ssb",
    paste("f",fbar.min,"-",fbar.max,sep=""),"f_dis2-3", "f_hc2-6","recruit", "catch", "landings", "discards")
    return(d.f)
}

########--------------------------+++++++++-----------------------------########
# NAME: inputTableSTF
# DOES: Creates a table of inputs to the STF
# E.G.# object <- window(ple4.stf,start=finalyear+1, end=finalyear+3)
inputTableSTF <- function(object){
  #note: for reasons of convenience I have removed calculation and output of fdisc and fland
  f.  <- as.data.frame(round(object@harvest,3))
  fd.   <- as.data.frame(round(object@harvest*(object@discards.n/object@catch.n),2))[,"data"]
  fl.   <- as.data.frame(round(object@harvest*(object@landings.n/object@catch.n),2))[,"data"]
  n.  <- as.data.frame(round(object@stock.n,0))[,"data"]
  cw. <- as.data.frame(round(object@catch.wt,2))[,"data"]
  lw. <- as.data.frame(round(object@landings.wt,2))[,"data"]
  dw. <- as.data.frame(round(object@discards.wt,2))[,"data"]
  sw. <- as.data.frame(round(object@stock.wt,2))[,"data"]
  mat.<- as.data.frame(object@mat)[,"data"]
  m.  <- as.data.frame(object@m)[,"data"]
  d.f <- cbind(f., fd., fl., n., cw., lw., dw.,sw., mat., m.)[,c(-3,-4,-5)]
  colnames(d.f) <- c("age",	"year",	"NA","f","f.disc","f.land",  "stock.n"	, "catch.wt", "landings.wt", "discards.wt","stock.wt", "mat", "M")
  return(d.f)
}
########--------------------------+++++++++-----------------------------########
# NAME: ageTableSTF
# DOES: Creates a table by year AND AGE for the projection period
ageTableSTF <- function(object){
  catch.n    <- as.data.frame(round(object@catch.n))[,"data"]
  catch      <- as.data.frame(round(object@catch.n*object@catch.wt,0))[,"data"]
  landings.n <- as.data.frame(round(object@landings.n,0))[,"data"]
  landings   <- as.data.frame(round(object@landings.n*object@landings.wt,0))[,"data"]
  discards.n <- as.data.frame(round(object@discards.n,0))[,"data"]
  discards   <- as.data.frame(round(object@discards.n*object@discards.wt,0))[,"data"]
  SSB        <- as.data.frame(round(object@stock.n*object@stock.wt*object@mat,0))[,"data"]
  TSB        <- as.data.frame(round(object@stock.n*object@stock.wt,0))[,"data"]
  inputtable <- inputTableSTF(object)
  d.f <- cbind(inputtable, catch.n, catch, landings.n, landings, discards.n, discards, SSB, TSB)
  return(d.f)
}

########--------------------------+++++++++-----------------------------########
# NAME: plotinternal
# DOES: have our own plotinternal function that shows more info

panel.ci1 <- function(x, y, interval='prediction',...) {	

        fit     <- lm(y ~ x, data=data.frame(x=x, y=y))
        rsq <-  sprintf("%4.3f",round(summary(fit)$r.squared,3))
	newX    <- data.frame(x=seq(min(x)-5,max(x)+5,0.01))
 	fitPred <- predict.lm(fit, newdata=newX, interval=interval, ...)
	panel.lmline(x, y, ..., identifier = "lmline")
	panel.xyplot(newX$x, fitPred[,2], type="a", lty=2, ...)
	panel.xyplot(newX$x, fitPred[,3], type="a", lty=2, ...)
        grid::grid.text(label = rsq,x = unit(0, "npc") + unit(0.25,"lines"), y = unit(1, "npc") - unit(0.25,"lines"),just="left",gp=gpar(cex=0.8)) 

}

plotinternal <- function(x, marklast=T ,... ) 
 { 

   skip <- matrix(1,nrow=dims(x@index)$age-1, ncol=dims(x@index)$age-1) 
   xydf <- NULL 
   for (age in dims(x@index)$min:(dims(x@index)$max-1) ) 
   { 
     for (inc in 1:(dims(x@index)$max-age)) 
       { 
         years <- dims(x@index)$minyear:(dims(x@index)$maxyear-inc) 
         xd <- as.numeric(x@index[as.character(age),as.character(years),]) 
         yd <- as.numeric(x@index[as.character(age+inc),as.character(years+inc),]) 
         d <- paste("age",age,"vs",age+inc) 
         xydf <- rbind(xydf,cbind(as.data.frame(cbind(xd,yd)),d)) 
         skip[dims(x@index)$max-dims(x@index)$min+1-inc, age-dims(x@index)$min + 1] <-0 
       } 
   } 
   xydf <- xydf[xydf$xd != 0 & xydf$yd != 0 & !is.na(xydf$xd) & !is.na(xydf$yd),] 

   if (marklast){
   	print(xyplot(yd~xd|d,outer=F, col="black", data=xydf, panel=
   		function(x,y,marklast,...) {
      			panel.xyplot(x,y, ...)
      			panel.xyplot(x[length(x)], y[length(y)], pch=16, ...)
      			panel.ci1(x,y,...)}, 
     scales=list(log="e",relation="sliced",draw=FALSE), layout=c(dims(x@index)$age-1, 
     dims(x@index)$age-1), xlab="log index", ylab="log index", skip=as.numeric(skip),  
     main=x@name))
   } else {
   	print(xyplot(yd~xd|d,outer=F, col="black", data=xydf, panel=
   		function(x,y,marklast,...) {
      			panel.xyplot(x,y, ...)
      			panel.ci(x,y,...)
                        }, 
     scales=list(log="e",relation="sliced",draw=FALSE), layout=c(dims(x@index)$age-1, 
     dims(x@index)$age-1), xlab="log index", ylab="log index", skip=as.numeric(skip),  
     main=x@name))
   }
 
 } # }}} 

# NAME: plotinternal
# DOES: have our own plot internal consistency plot function that shows more info

panel.ci2 <- function(x, y, interval='prediction',...) {	

        fit     <- lm(y ~ x, data=data.frame(x=x, y=y))
        rsq <-  sprintf("%4.3f",round(summary(fit)$r.squared,3))
	newX    <- data.frame(x=seq(min(x)-5,max(x)+5,0.01))
 	fitPred <- predict.lm(fit, newdata=newX, interval=interval, ...)
	panel.lmline(x, y, ..., identifier = "lmline")
	panel.xyplot(newX$x, fitPred[,2], type="a", col="black", lty=2, ...)
	panel.xyplot(newX$x, fitPred[,3], type="a", col="black", lty=2, ...)

}


# internal consistency  {{{ 
 plotInternalConsistency <-  function(idx,log.scales=TRUE, 
   cols=c("white", "yellow", "red"),use.rsq=TRUE,mark.significant=FALSE,mark.last=F,...) 
   { 
 
 
   # Define colour function 
   if(!length(cols)>0) stop("Colour definitions do not contain sufficient number of colours (at least one required)") 
   if(length(cols)==1) cols <- rep(cols,2) 
   colFn <- colorRamp(colors=cols) 
 

   #Number of ages 
   ages <- dimnames(idx@index)[[1]] 
 
 
   #Convert to Cohorts, reshape into appropriate format for splom 
   flc <-  if(log.scales) {log10(idx@index)}  
     else {idx@index} 
   flc  <- as.data.frame(FLCohort(flc)) 
   flc.wide <-  reshape(flc,direction="wide",timevar=names(flc)[1],idvar=names(flc)[2:6]) 
   names(flc.wide) <-  gsub("data.","",names(flc.wide)) 
    
   #Default plot settings 
   plot.args <- list(~flc.wide[ages],data=flc.wide, pscales=0,varname.font=2,varname.cex=1.5, 
       xlab = if(log.scales) {expression(paste(Log[10]," (Index Value)"))}  
         else {"Index Value"}, 
       ylab = if(log.scales) {expression(paste(Log[10]," (Index Value)"))} 
         else { "Index Value"}, 
       sub=list(if(use.rsq) {expression(paste("Lower right panels show the Coefficient of Determination (",italic(r^2),")"))} 
         else { expression(paste("Lower right panels show the Coefficient of Correlation (",italic(r),")"))},cex=0.7), 
       upper.panel=function(x,y,...) 
       { 
         # Filter out NAs 
         both.points  <-  is.finite(x) & is.finite(y) 
         x.filtered <-  x[both.points] 
         y.filtered <-  y[both.points] 
         # Only plot lmline if there is more than one point - colour panel according to rsq. 
         if(length(x.filtered)>2) 
         { 
           r <-  cor(y.filtered,x.filtered)     
           if(use.rsq) { 
             panel.colour <- r^2           #Colour & number panel based on the coefficient of determination (r^2) 
           } else { 
             panel.colour <- 0.5*r+0.5     #Colour & number panel based on the correlation coefficient (r) 
           } 
           if(is.numeric(mark.significant) | identical(TRUE,mark.significant) ) { 
               lm.model <- lm(y.filtered ~ x.filtered) 
               p.value  <- summary(lm.model)$coefficients["x.filtered",4]/2    #Halve the p-value, as we are doing a one sided test, not a two 
               slope    <- summary(lm.model)$coefficients["x.filtered",1] 
               signif.level <- 0.05 
               if(is.numeric(mark.significant)) signif.level <- mark.significant  
               if(p.value < signif.level & slope >0) {  #If marking significance, only fill panel and draw line when its significant 
                 number.format <- "%4.3f*"      #If its a significant correlation, mark with a * 
                 panel.fill(col = rgb(colFn(panel.colour),maxColorValue=255))   #Colour panel based on the coefficient of determination (r^2) 
                 panel.lmline(x.filtered,y.filtered,lwd=2) 
               }
               if (mark.last) panel.xyplot(x=x.filtered[length(x.filtered)], y=y.filtered[length(y.filtered)], pch=15) 
                 
           } else {  #If not marking significance, always fill panel and draw best fit line 

               panel.fill(col = rgb(colFn(panel.colour),maxColorValue=255))   #Colour panel based on the coefficient of determination (r^2) 
               panel.ci2(x.filtered,y.filtered,lwd=2) 
               if (mark.last) panel.xyplot(x=x.filtered[length(x.filtered)], y=y.filtered[length(y.filtered)], pch=15,...) 
           } 
         } 
         panel.splom(x.filtered,y.filtered,col="black",...) 
       }, 
       lower.panel=function(x, y, ...) 
       { 
         #Filter out NAs 
         both.points  <-  is.finite(x) & is.finite(y) 
         x.filtered <-  x[both.points] 
         y.filtered <-  y[both.points] 
          
         #Calculate r squared - but only if there is enough data to do so 
         if(length(x.filtered)>2) 
         { 
           r <-  cor(y.filtered,x.filtered) 
           if(use.rsq) { 
             panel.colour <- r^2           #Colour & number panel based on the coefficient of determination (r^2) 
             panel.number <- round(r^2,3) 
           } else { 
             panel.colour <- 0.5*r+0.5  #Colour & number panel based on the correlation coefficient (r) 
             panel.number <- round(r,3) 
           } 
           number.format <- "%4.3f" 
           if(is.numeric(mark.significant) | identical(TRUE,mark.significant) ) { 
               lm.model <- lm(y.filtered ~ x.filtered) 
               p.value  <- summary(lm.model)$coefficients["x.filtered",4]/2 
               slope    <- summary(lm.model)$coefficients["x.filtered",1] 
               signif.level <- 0.05 
               if(is.numeric(mark.significant)) signif.level <- mark.significant  
              if(p.value < signif.level & slope > 0) {  #If marking significance, only fill panel when its significant & positive 
                 number.format <- "%4.3f*"      #If its a significant correlation, mark with a * 
                 panel.fill(col = rgb(colFn(panel.colour),maxColorValue=255))   #Colour panel based on the coefficient of determination (r^2) 
               }                 
           } else {  #If not marking significance, always fill panel  
               panel.fill(col = rgb(colFn(panel.colour),maxColorValue=255))   #Colour panel based on the coefficient of determination (r^2) 
           } 
           grid::grid.text(label =sprintf(number.format,panel.number),x = unit(0.5, "npc"), 
             y = unit(0.5,"npc"),just="center")}}) 
 
 
   #Passed settings 
   passed.args   <- list(...) 
   plot.args[names(passed.args)] <- passed.args 
    
   #Do plot 
   p <- do.call(splom,plot.args) 
   print(p) 
   return(p) 
}   # }}} 








