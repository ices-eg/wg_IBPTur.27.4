rm(list=ls())
library(msy)
source("D:/ICES_WG/wg_IBPTur.27.4/refpoints/eqsim/eqsim_functions.R")
library(FLSAM)
load("D:/ICES_WG/wg_IBPTur.27.4/refpoints/final_finalassessmentOut.rdata")
library(FLCore)
name(TUR) <- "turbot"

TUR <- TUR + TUR.sam
stock.n(TUR) <- stock.n(TUR.sam)
harvest(TUR) <- harvest(TUR.sam)

catch(TUR) <- landings(TUR)
discards(TUR)[] <- 0


#step 1 do BLim and Flim
FIT101 <- eqsr_fit(TUR, nsamp = 1000, models = c( "Segreg"))

eqsr_plot(FIT,n=1000)

print(Blim <-  FIT101[["sr.det"]][,"b"])

print(Bpa <- 1.4 * Blim)


#segreg3 <- function(ab,sbb) log(ifelse(ssb>=207287, ab$a*207287, ab$a *ssb ))
FIT <- eqsr_fit(TUR,
                nsamp = 3000,   ###CC changed this number from 2000 to 3000, if this number changes, the result of SIM101 also changes
                models = c("Segreg"),  remove.years=c(ac(2013,2014,2015)))
print(Blim <-  FIT[["sr.det"]][,"b"][1])

SIM101 <- eqsim_run(FIT,  bio.years = c(2006, 2015), bio.const = FALSE,
                    sel.years = c(2006, 2015), sel.const = FALSE,
                    Fcv=0, Fphi=0,
                    Btrigger = 0,Blim=Blim,Bpa=NA,
                    Fscan = seq(0,1.2,len=61),verbose=FALSE)

eqsim_plot(SIM101,catch="FALSE")
Coby.fit(SIM101,outfile='aa')
# from this table get F50, catF
print(Flim <- SIM101$Refs2[1,3])

###### NEXT get BPA: go for 1.4 *Blim and Fpa (go for Flim /1.4)

print(Fpa <- Flim/1.4)

##### NOW estimate FMSY initially: calculated as the F that maximises median long-term yield in stochastic simulation, under constant F exploitation (i.e. without MSY Btrigger);   if this FMSY value is > Fpa, reduce it to Fpa (i.e. FMSY can not exceed Fpa).
##### usE fit from three S-Rs

FIT1 <- eqsr_fit(TUR,
                 nsamp = 2000, 
                 models = c("Ricker", "Segreg", "Bevholt"),  remove.years=c(ac(2013,2014,2015)))


eqsr_plot(FIT1,n=2000)

SIM1 <- eqsim_run(FIT1,  bio.years = c(2011,2015), bio.const = FALSE,
                  sel.years = c(2011,2015), sel.const = FALSE,
                  Fcv=0.212, Fphi=0.423,  # these are defauts, taken from WKMSYREF4, as used in Saithe assessments
                  Btrigger = 0,Blim=Blim,Bpa=Bpa,Fscan = seq(0,1.0,len=51),verbose=FALSE)

###########Chun's code########################### tried several options, 
#FIT1 <- eqsr_fit(TUR,
#                 nsamp = 1e4, 
#                 models = c("Ricker", "Segreg", "Bevholt"),remove.years=c(ac(2013,2014,2015)))
#eqsr_plot(FIT1,n=1e4)
FIT1 <- eqsr_fit(TUR,
                 nsamp = 3000, 
                 models = c("Ricker", "Segreg", "Bevholt"),remove.years=c(ac(2013,2014,2015)))
eqsr_plot(FIT1,n=3000)

SIM1 <- eqsim_run(FIT1,  bio.years = c(2006,2015), bio.const = FALSE,
                  sel.years = c(2006,2015), sel.const = FALSE,
                  Fcv=0.212, Fphi=0.423,  # these are defauts, taken from WKMSYREF4, as used in Saithe assessments
                  Btrigger = 0,Blim=Blim,Bpa=Bpa,Fscan = seq(0,1.0,len=51),verbose=FALSE)
#################################################

#eqsim_plot(SIM1,catch="FALSE")

Coby.fit(SIM1,outfile='ple sim1')
#get median MSY from lanF
print(Fmsy <- SIM1$Refs2[2,4])
#also get F05 from catF
print(F05 <- SIM1$Refs2[1,1])


# #########################################################
# # get Btrigger
# ##########################################################
#rerun with Fcv and Fphi are zero, needed to do the Btrigger calc
SIM2 <- eqsim_run(FIT,  bio.years = c(2006, 2015), bio.const = FALSE,
                  sel.years = c(2006, 2015), sel.const = FALSE,
                  Fcv=0, Fphi=0,  # these are defauts, taken from WKMSYREF4, as used in Saithe assessments
                  Btrigger = 0,Blim=Blim,Bpa=Bpa,Fscan = seq(0.1,0.3,len=21),verbose=FALSE)
# 
# Coby.fit(SIM2,outfile='saithe full time series no Btrigger but Fcv and Fphi')
#below is proxy for BMSY trigger, but note that this has to be manually set to Fmsy values because of rounding 
SIM2$rbp[,4][SIM2$rbp$variable=='Spawning stock biomass'& SIM2$rbp$Ftarget==0.21]
#
#
# Btrigger 564599 (current SSB/1.4) ??

# #if 5 or more years of fishing at Fmsy then next, otherwise Btrig= BPa 
# #if 5% of BFMSY > Bpa  then next otherwise Btrigger = Bpa
# #if 5% if BMSY > current trigger then next
# #is greater than 5% of SSB
# 
# #if indeed new Btrigger
# 
SIM3 <- eqsim_run(FIT,  bio.years = c(2006, 2015), bio.const = FALSE,
                  sel.years = c(2006, 2015), sel.const = FALSE,
                  Fcv=0.212, Fphi=0.423,  # these are defauts, taken from WKMSYREF4, as used in Saithe assessments
                  Btrigger = 564599 ,Blim=Blim,Bpa=Bpa,Fscan = seq(0,1.2,len=61),verbose=FALSE)
# 
# Coby.fit(SIM3,outfile='saithe full time series no Btrigger but Fcv and Fphi')
print(Fmsy <- SIM3$Refs2[2,4])

print(F05 <- SIM3$Refs2[1,1])

