#-------------------------------------------------------------------------------
# IBPNEW 2012
# Code to generate the .dat file for input into the ADMB model
# Uses Lowestoft format input files
# David Miller, some code adapted from FLR-project
# Date: 3 Oct 2012
#-------------------------------------------------------------------------------
rm(list=ls())

# FLR
# install.packages(repos="http://flr-project.org/R")
library(FLCore);library(mgcv)
library(FLAssess);
library(stockassessment)
library(FLEDA); library(splines);
library(scales); library(gplots);library(grid); library(gridExtra); library(latticeExtra)
library(sas7bdat)
library(TMB); library(FLSAM)

# Set paths to folders
Path      <- "D:/Repository/Turbot/assessment runs/"
dataPath  <- paste(Path,"Lowestoft files/",sep="")
outPath   <- paste(Path,"trial_runs_2017/Output/",sep="")
codePath  <- paste(Path,"Trial_runs_2017/",sep="")

#- Load runs and merge them into FLSAMs objects
runs      <- data.frame(name=c("base","trim89","trim98","addUKLPUE","lowPlusgroup6","lowPlusgroup7","lowPlusgroup8","lowPlusgroup9","lowPlusgroup10",
                               "addIBTS_CPUE_EB","addIBTS_CPUE","addIBTS_NpH","addDGBTS","noCat6LpUE","noDNKCaA","addAllSurveys",
                               "lowAgeSurvey3","lowAgeSurvey4","lowAgeSurvey5","lowAgeSurvey6","lowAgeSurvey7","noUKCaA","addNewNLLPUEbase","addNewNLLPUEModelA"
                               ,"addNewNLLPUEModelB","addNewNLLPUEModelC","addNewNLLPUEModelD"),
                        runID=1:27,
                        description=c("WGNSSK base run","Start CaA in 1989","Start CaA in 1998 + other surveys start in 1998",
                                      "Add the UK LPUE time-series","Reduce plusgroup to age 6","Reduce plusgroup to age 7","Reduce plusgroup to age 8",
                                      "Reduce plusgroup to age 9","Reduce plusgroup to age 10","Add IBTS survey CPUE_EB",
                                      "Add IBTS survey CPUE","Add IBTS survey Numbers per hour",
                                      "Replace BTS-ISIS with Delta-Gam BTS","Exclude age 1-2 from NL LPUE","Remove DNK CaA in 2014-106",
                                      "Add all surveys (but keeping unique data per survey)","Lower pg of the surveys to 3","Lower pg of the surveys to 4",
                                      "Lower pg of the surveys to 5","Lower pg of the surveys to 6","Lower pg of the surveys to 7",
                                      "Start CaA in 1981 and exclude 2000-2002","Add new standardized NL LPUE"),
                        datafile=c("base_assessmentOut","trim89_assessmentOut","trim98_assessmentOut","addUKLPUE_assessmentOut","pgroup_pg6assessmentOut",
                                   "pgroup_pg7assessmentOut","pgroup_pg8assessmentOut","pgroup_pg9assessmentOut","pgroup_pg10assessmentOut",
                                   "addIBTS_IBTS_CPUE_EBassessmentOut","addIBTS_IBTS_CPUEassessmentOut","addIBTS_IBTS_NpHassessmentOut","addBTSDG_assessmentOut",
                                   "noCat6LpUE_assessmentOut","noDNKCaA_assessmentOut","allSurveys_assessmentOut","pgroupSurveys_pgSurvey3assessmentOut",
                                   "pgroupSurveys_pgSurvey4assessmentOut","pgroupSurveys_pgSurvey5assessmentOut","pgroupSurveys_pgSurvey6assessmentOut",
                                   "pgroupSurveys_pgSurvey7assessmentOut","noUKCaA_assessmentOut","addNewNLLPUE_baseassessmentOut","addNewNLLPUE_ModelAassessmentOut",
                                   "addNewNLLPUE_ModelBassessmentOut","addNewNLLPUE_ModelCassessmentOut","addNewNLLPUE_ModelDassessmentOut"))

#- Function to load runs and give them the correct names
loadRuns <- function(runID,path,runs){

  sams  <- new("FLSAMs")
  sams.retro <- list()
  for(iRun in runID){
    load(file.path(path,paste0(runs[which(runs$runID==iRun),"datafile"],".RData")))
    sams[[iRun]] <- TUR.sam
    sams.retro[[iRun]] <- TUR.retro
  }
  names(sams) <- runs$name[match(runs$runID,runID)]
  names(sams.retro) <- runs$name[match(runs$runID,runID)]
  
  return(list(sams=sams,retros=sams.retro))}
  
### ------------------------------------------------------------------------------------------------------
###   Examples
### ------------------------------------------------------------------------------------------------------

#- change in CaA data
res <- loadRuns(runID=c(1,2,3,15,22),path=outPath,runs=runs)
plot(res[["sams"]])

#- change in pg of stock or survey
res <- loadRuns(runID=c(1,5:9,17:21),path=outPath,runs=runs)
plot(res[["sams"]])

#- Change the surveys that are included
res <- loadRuns(runID=c(1,4,10,11,12,13,16,23:27,path=outPath,runs=runs)

