rm(list=ls())
#library(devtools)
#install_git("https://github.com/ices-tools-prod/msy")
#install.packages(repos=NULL, pkgs= "D://wg_IBPTur.27.4/software/FLSAM_2.01.zip")

library(msy)
source("D:/wg_IBPTur.27.4/refpoints/eqsim/eqsim_functions.R")
library(FLSAM)
load("D:/wg_IBPTur.27.4/assessment runs/final_runs_2017/final_finalassessmentOut.rdata")
library(FLCore)
name(TUR) <- "turbot"

stock.n(TUR) <- stock.n(TUR.sam)
harvest(TUR) <- harvest(TUR.sam)

catch(TUR) <- landings(TUR)
discards(TUR)[] <- 0

###BLim and Bpa
FITlim <- eqsr_fit(TUR, nsamp = 5000, models = c( "Segreg"),
                   remove.years=c(2014, 2015, 2016))

eqsr_plot(FITlim,n=5000)

print(Blim <-  FITlim[["sr.det"]][,"b"])

print(Bpa <- 1.4 * Blim)

SIM101 <- eqsim_run(FITlim,  bio.years = c(2007, 2016), bio.const = FALSE,
                    sel.years = c(2007, 2016), sel.const = FALSE,
                    Fcv=0, Fphi=0,
                    Btrigger = 0,Blim=Blim,Bpa=NA,
                    Fscan = seq(0,1.2,len=61),verbose=FALSE)

print(Flim <- SIM101$Refs2[1,3])
print(Fpa <- Flim/1.4)

eqsim_plot(SIM101,catch="FALSE")
Coby.fit(SIM101,outfile='aa')


### explore fits
FITs <- eqsr_fit(TUR,
                 nsamp = 5000,
                 models = c("Segreg"),
                 remove.years=c(2014, 2015, 2016))

FITr <- eqsr_fit(TUR,
                 nsamp = 5000,   
                 models = c("Ricker"),
                 remove.years=c(2014, 2015, 2016))

FITrs <- eqsr_fit(TUR,
                   nsamp = 5000,
                   models = c("Ricker", "Segreg"),
                   remove.years=c(2014, 2015, 2016))

eqsr_plot(FITrs,n=5000)

##### NOW estimate FMSY initially: calculated as the F that maximises median long-term yield in stochastic simulation, under constant F exploitation (i.e. without MSY Btrigger);   if this FMSY value is > Fpa, reduce it to Fpa (i.e. FMSY can not exceed Fpa).
##### usE fit from ricker and segreg  S-Rs

SIM1rs <- eqsim_run(FITrs,  bio.years = c(2012,2016), bio.const = FALSE,
                  sel.years = c(2012,2016), sel.const = FALSE,
                  Fcv=0.212, Fphi=0.423,  # these are defauts, taken from WKMSYREF4, as used in Saithe assessments
                  Btrigger = 0,Blim=Blim, Bpa=Bpa, Fscan = seq(0,1.0,len=51),verbose=FALSE)

SIM1s <- eqsim_run(FITs,  bio.years = c(2012,2016), bio.const = FALSE,
                    sel.years = c(2012,2016), sel.const = FALSE,
                    Fcv=0.212, Fphi=0.423,  # these are defauts, taken from WKMSYREF4, as used in Saithe assessments
                    Btrigger = 0,Blim=Blim, Bpa=Bpa, Fscan = seq(0,1.0,len=51),verbose=FALSE)


########################### tried several options, 

eqsim_plot(SIM1rs,catch="FALSE")
Coby.fit(SIM1rs,outfile='tur sim1rs')

#get median MSY from lanF
print(FmsyRS <- SIM1rs$Refs2[2,4])
#also get F05 from catF
print(F05RS <- SIM1rs$Refs2[1,1])

#also the one based on segreg alon
eqsim_plot(SIM1s,catch="FALSE")
Coby.fit(SIM1s,outfile='tur sim1s')

#get median MSY from lanF
print(FmsyS <- SIM1s$Refs2[2,4])
#also get F05 from catF
print(F05S <- SIM1s$Refs2[1,1])

# #########################################################
# # get Btrigger
# ##########################################################
#rerun with Fcv and Fphi are zero, needed to do the Btrigger calc
SIM2rs <- eqsim_run(FITrs,  bio.years = c(2012, 2016), bio.const = FALSE,
                  sel.years = c(2012, 2016), sel.const = FALSE,
                  Fcv=0, Fphi=0,  # these are defauts, taken from WKMSYREF4, as used in Saithe assessments
                  Btrigger = 0,Blim=Blim,Bpa=Bpa,Fscan = seq(0.1,0.5,len=41),verbose=FALSE)

#rerun with Fcv and Fphi are zero, needed to do the Btrigger calc
SIM2s <- eqsim_run(FITs,  bio.years = c(2012, 2016), bio.const = FALSE,
                  sel.years = c(2012, 2016), sel.const = FALSE,
                  Fcv=0, Fphi=0,  # these are defauts, taken from WKMSYREF4, as used in Saithe assessments
                  Btrigger = 0,Blim=Blim,Bpa=Bpa,Fscan = seq(0.1,0.5,len=41),verbose=FALSE)

# 
# Coby.fit(SIM2,outfile='saithe full time series no Btrigger but Fcv and Fphi')
#below is proxy for BMSY trigger, but note that this has to be manually set to Fmsy values because of rounding 
print(Btrigrs <- SIM2rs$rbp[,4][SIM2rs$rbp$variable=='Spawning stock biomass'& SIM2rs$rbp$Ftarget==0.40])

print(Btrigs <- SIM2s$rbp[,4][SIM2s$rbp$variable=='Spawning stock biomass'& SIM2s$rbp$Ftarget==0.40])

# #if 5 or more years of fishing at Fmsy then next, otherwise Btrig= BPa 
# #if 5% of BFMSY > Bpa  then next otherwise Btrigger = Bpa
# #if 5% if BMSY > current trigger then next
# #is greater than 5% of SSB
# 
# #if indeed new Btrigger
# 
SIM3rs <- eqsim_run(FITrs,  bio.years = c(2012, 2016), bio.const = FALSE,
                  sel.years = c(2012, 2016), sel.const = FALSE,
                  Fcv=0.212, Fphi=0.423,  # these are defauts, taken from WKMSYREF4, as used in Saithe assessments
                  Btrigger = Btrigrs ,Blim=Blim,Bpa=Bpa,Fscan = seq(0,1.2,len=61),verbose=FALSE)

SIM3s <- eqsim_run(FITs,  bio.years = c(2012, 2016), bio.const = FALSE,
                    sel.years = c(2012, 2016), sel.const = FALSE,
                    Fcv=0.212, Fphi=0.423,  # these are defauts, taken from WKMSYREF4, as used in Saithe assessments
                    Btrigger = Btrigs ,Blim=Blim,Bpa=Bpa,Fscan = seq(0,1.2,len=61),verbose=FALSE)


# 
Coby.fit(SIM3rs,outfile='turbot RS  with Btrigger and Fcv and Fphi')
print(Fmsyrs <- SIM3rs$Refs2[2,4])
print(F05rs <- SIM3rs$Refs2[1,1])

Coby.fit(SIM3rs,outfile='turbot S  with Btrigger and Fcv and Fphi')
print(Fmsys <- SIM3s$Refs2[2,4])
print(F05s <- SIM3s$Refs2[1,1])

