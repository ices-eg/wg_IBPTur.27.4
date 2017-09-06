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
library(knitr)

# Set paths to folders
Path      <- "D:/Repository/Turbot/assessment runs/"
dataPath  <- paste(Path,"Lowestoft files/",sep="")
outPath   <- paste(Path,"trial_runs_2017/Output/",sep="")
codePath  <- paste(Path,"Trial_runs_2017/",sep="")

#- Source the run overview
source(file.path(codePath,"runOverview.r"))
kable(runs)

|name               | runID|description                                          |datafile                             |
|:------------------|-----:|:----------------------------------------------------|:------------------------------------|
|base               |     1|WGNSSK base run                                      |base_assessmentOut                   |
|trim89             |     2|Start CaA in 1989                                    |trim89_assessmentOut                 |
|trim98             |     3|Start CaA in 1998 + other surveys start in 1998      |trim98_assessmentOut                 |
|addUKLPUE          |     4|Add the UK LPUE time-series                          |addUKLPUE_assessmentOut              |
|lowPlusgroup6      |     5|Reduce plusgroup to age 6                            |pgroup_pg6assessmentOut              |
|lowPlusgroup7      |     6|Reduce plusgroup to age 7                            |pgroup_pg7assessmentOut              |
|lowPlusgroup8      |     7|Reduce plusgroup to age 8                            |pgroup_pg8assessmentOut              |
|lowPlusgroup9      |     8|Reduce plusgroup to age 9                            |pgroup_pg9assessmentOut              |
|lowPlusgroup10     |     9|Reduce plusgroup to age 10                           |pgroup_pg10assessmentOut             |
|addIBTS_CPUE_EB    |    10|Add IBTS survey CPUE_EB                              |addIBTS_IBTS_CPUE_EBassessmentOut    |
|addIBTS_CPUE       |    11|Add IBTS survey CPUE                                 |addIBTS_IBTS_CPUEassessmentOut       |
|addIBTS_NpH        |    12|Add IBTS survey Numbers per hour                     |addIBTS_IBTS_NpHassessmentOut        |
|addDGBTS           |    13|Replace BTS-ISIS with Delta-Gam BTS                  |addBTSDG_assessmentOut               |
|noCat6LpUE         |    14|Exclude age 1-2 from NL LPUE                         |noCat6LpUE_assessmentOut             |
|noDNKCaA           |    15|Remove DNK CaA in 2014-106                           |noDNKCaA_assessmentOut               |
|addAllSurveys      |    16|Add all surveys (but keeping unique data per survey) |allSurveys_assessmentOut             |
|lowAgeSurvey3      |    17|Lower pg of the surveys to 3                         |pgroupSurveys_pgSurvey3assessmentOut |
|lowAgeSurvey4      |    18|Lower pg of the surveys to 4                         |pgroupSurveys_pgSurvey4assessmentOut |
|lowAgeSurvey5      |    19|Lower pg of the surveys to 5                         |pgroupSurveys_pgSurvey5assessmentOut |
|lowAgeSurvey6      |    20|Lower pg of the surveys to 6                         |pgroupSurveys_pgSurvey6assessmentOut |
|lowAgeSurvey7      |    21|Lower pg of the surveys to 7                         |pgroupSurveys_pgSurvey7assessmentOut |
|noUKCaA            |    22|Start CaA in 1981 and exclude 2000-2002              |noUKCaA_assessmentOut                |
|addNewNLLPUEbase   |    23|Add new standardized NL LPUE base case               |addNewNLLPUE_baseassessmentOut       |
|addNewNLLPUEModelA |    24|Add new standardized NL LPUE ModelA                  |addNewNLLPUE_ModelAassessmentOut     |
|addNewNLLPUEModelB |    25|Add new standardized NL LPUE ModelB                  |addNewNLLPUE_ModelBassessmentOut     |
|addNewNLLPUEModelC |    26|Add new standardized NL LPUE ModelC                  |addNewNLLPUE_ModelCassessmentOut     |
|addNewNLLPUEModelD |    27|Add new standardized NL LPUE ModelD                  |addNewNLLPUE_ModelDassessmentOut     |
|noGerCaA           |    28|Remove German Weber data                             |noGerCaA_assessmentOut               |


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
res <- loadRuns(runID=c(1,4,10,11,12,13,16,23:27),path=outPath,runs=runs)
plot(res[["sams"]])

res <- loadRuns(runID=c(1,15,22,28),path=outPath,runs=runs)
plot(res[["sams"]])