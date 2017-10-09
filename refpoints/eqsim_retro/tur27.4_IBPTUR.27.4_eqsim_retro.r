rm(list=ls())
#library(devtools)
#install_git("https://github.com/ices-tools-prod/msy")
#install.packages(repos=NULL, pkgs= "D://wg_IBPTur.27.4/software/FLSAM_2.01.zip")

library(msy)
source("D:/wg_IBPTur.27.4/refpoints/eqsim/eqsim_functions.R")
library(FLSAM)
load("D:/wg_IBPTur.27.4/assessment runs/final_runs_2017/final_finalassessmentOut.rdata")
library(FLCore)

years <- 2016:2011
SRFits <- list()
SIMlim_pa <- list()
SIMfmsy <- list()
nits <- 5000
setwd("D:/wg_IBPTur.27.4/refpoints/eqsim_retro/")
iRetro <- 3
for(iRetro in 1:6 ){
  name(TUR) <- "turbot"
  TUR <- window(TUR, end = years[iRetro])
  stock.n(TUR) <- stock.n(TUR.retro[[iRetro]])
  harvest(TUR) <- harvest(TUR.retro[[iRetro]])
  
  catch(TUR) <- landings(TUR)
  discards(TUR)[] <- 0
  name(TUR) <- "turbot"
  
  ### segreg fit
  FITs <- eqsr_fit(TUR, nsamp = nits, models = c( "Segreg"),
                     remove.years=c((dims(TUR)$maxyear - 2):dims(TUR)$maxyear))
  
  eqsr_plot(FITs,n=nits)
  
  print(Blim <-  FITs[["sr.det"]][,"b"])
  
  print(Bpa <- 1.4 * Blim)
  
  ### Flim/Blim - Fpa/Bpa
  SIMlim <- eqsim_run(FITs,  bio.years = c((dims(TUR)$maxyear - 10), dims(TUR)$maxyear), bio.const = FALSE,
                      sel.years = c((dims(TUR)$maxyear - 10), dims(TUR)$maxyear), sel.const = FALSE,
                      Fcv=0, Fphi=0,
                      Btrigger = 0,Blim=Blim,Bpa=NA,
                      Fscan = seq(0,1.2,len=61),verbose=FALSE)
  
  print(Flim <- SIMlim$Refs2[1,3])
  print(Fpa <- Flim/1.4)
  
  eqsim_plot(SIMlim,catch = "FALSE")
  Coby.fit(SIMlim,outfile = 'aa')

  SIMlim_pa[[years[iRetro]]] <- SIMlim
  
  print(Flim <- SIMlim$Refs2[1,3])
  print(Fpa <- Flim/1.4)
  
  ### Fmsy - proxy
  SIM_fmsy <- eqsim_run(FITs,  bio.years = c((dims(TUR)$maxyear - 5), dims(TUR)$maxyear), bio.const = FALSE,
                     sel.years = c((dims(TUR)$maxyear - 5), dims(TUR)$maxyear), sel.const = FALSE,
                     Fcv=0.212, Fphi=0.423,  # these are defauts, taken from WKMSYREF4, as used in Saithe assessments
                     Btrigger = 0,Blim=Blim, Bpa=Bpa, Fscan = seq(0,1.0,len=51),verbose=FALSE)
  
  SIMfmsy[[years[iRetro]]] <- SIM_fmsy
  
  eqsim_plot(SIM_fmsy,catch="FALSE")
  Coby.fit(SIM_fmsy,outfile='tur sim1s')
  
  #get median MSY from lanF
  print(FmsyS <- SIM_fmsy$Refs2[2,4])
  #also get F05 from catF
  print(F05S <- SIM_fmsy$Refs2[1,1])
  
  ### Btrigger
  SIM2s <- eqsim_run(FITs,  bio.years = c((dims(TUR)$maxyear - 5), dims(TUR)$maxyear), bio.const = FALSE,
                     sel.years = c((dims(TUR)$maxyear - 5), dims(TUR)$maxyear), sel.const = FALSE,
                     Fcv=0, Fphi=0,  # these are defauts, taken from WKMSYREF4, as used in Saithe assessments
                     Btrigger = 0,Blim=Blim,Bpa=Bpa,Fscan = seq(0.1,0.5,len=41),verbose=FALSE)
  
  print(Btrigs <- SIM2s$rbp[,4][SIM2s$rbp$variable=='Spawning stock biomass'& SIM2s$rbp$Ftarget==0.40])
  
  ### Fmsy
  SIM3s <- eqsim_run(FITs,  bio.years = c((dims(TUR)$maxyear - 5), dims(TUR)$maxyear), bio.const = FALSE,
                     sel.years = c((dims(TUR)$maxyear - 5), dims(TUR)$maxyear), sel.const = FALSE,
                     Fcv=0.212, Fphi=0.423,  # these are defauts, taken from WKMSYREF4, as used in Saithe assessments
                     Btrigger = Btrigs ,Blim=Blim,Bpa=Bpa,Fscan = seq(0,1.2,len=61),verbose=FALSE)
  
  Coby.fit(SIM3s,outfile='turbot S  with Btrigger and Fcv and Fphi')
  print(Fmsys <- SIM3s$Refs2[2,4])
  print(F05s <- SIM3s$Refs2[1,1])
  
}
