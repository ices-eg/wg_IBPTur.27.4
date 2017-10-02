### spict run for WG_IBPTur.27.4
### code adapted from WGNSSK 2017
rm(list=ls())

library(FLCore)
### devtools::install_github("mawp/spict/spict")
library(spict)
library(scales)
library(ellipse)
library(ggplot2)

dataPath  <- "D:/wg_IBPTur.27.4/assessment runs/Lowestoft files/"
outPath   <- "D:/wg_IBPTur.27.4/refpoints/spict/"

load("D:/wg_IBPTur.27.4/assessment runs/Final_runs_2017/final_finalassessmentOut.RData")

### make indices age-aggregated - make biomass using stock weights
obsI <- list()
obsI$index1 <- c(apply(index(TUR.tun[[2]])[2:7,] * stock.wt(TUR)[ac(2:7),ac(1991:2016)],
                     FUN = sum, 2)) ##BTS-ISIS
obsI$index2 <- c(TUR.tun[[3]]@index)#NL_LPUE
obsI$index3 <- c(apply(index(TUR.tun[[1]])[2:6,] * stock.wt(TUR)[ac(2:6),ac(2004:2016)],
                       FUN = sum, 2))##SNS

timeI <- list()
timeI$index1 <- seq(1991.75,2016.75,1)
timeI$index2 <- seq(1995,2016,1)
timeI$index3 <- seq(2004.75,2016.75,1)

inp <- list(obsC=c(landings(TUR)),
         obsI=obsI,
         timeC=seq(1981,2016,1),
         timeI=timeI)

### Robust likelihood
#inp$robflagc=1 ## 0=not robust catch likelihood, 1=robust
#inp$robflagi=1 ## 0=not robust index likelihood, 1=robust
#inp$phases$logitpp=1; ## -1: don't estimate robust mixture proportion, >=1 estimate in phase xx
#inp$phases$logp1robfac=1; ## -1: don't estimate robust mixture sd.dev scaling , >=1 estimate in phase xx

## Fixed parameter values inis
#inp$ini$logbeta <- log(0.29)
#inp$ini$logr  <- log(0.3954)
#inp$ini$logm  <- log(4.41037e+3)
i#np$ini$logK  <- log(3.67625e+4)

#inp$ini$logalpha <- log(1)

## Phases
# -1: not estimated
#  1: estimated
#inp$phases$logbeta <- -1
#inp$phases$logalpha <- 1

## Predict catch or indices until this time
#inp$timepredc <- 1990
#inp$timepredi <- 1990

## Priors
## (value, std.dev, 0/1=Off/On)
inp$priors$logn <- c(log(2), 1, 0)
inp$priors$logalpha <- c(log(2), 3, 0)
inp$priors$logbeta <- c(log(2), 1, 0)

##no priors
inp_no <- inp
inp_no$priors$logn <- c(1, 1, 0)
inp_no$priors$logalpha <- c(1, 1, 0)
inp_no$priors$logbeta <- c(1, 1, 0)


inp <- check.inp(inp)
fit <- fit.spict(inp)
if(fit$opt$convergence!=0) stop("Error: model did not converge.");
fit <- calc.osa.resid(fit)
windows(10,10)
plot(fit)
summary(fit)

windows(10,10)
spict::plotspict.ffmsy(fit, ylim=c(0,3))
windows(10,10)
spict::plotspict.biomass(fit)

inp_no <- check.inp(inp_no)
fit_no <- fit.spict(inp_no)
if(fit_no$opt$convergence!=0) stop("Error: model did not converge.");
fit_no <- calc.osa.resid(fit_no)
windows(10,10)
plot(fit_no)
summary(fit_no)

windows(10,10)
spict::plotspict.ffmsy(fit_no, ylim=c(0,3))
windows(10,10)
spict::plotspict.biomass(fit_no)

#####################
# perform retro
######################

plot_tmsrs <- function(input, label = "model"){
  ### plot
  p <- ggplot(data = input, aes(x = as.numeric(as.character(year)), y = as.numeric(as.character(est)))) +
    scale_colour_discrete(label) + facet_wrap(~ quant, scale = "free", labeller = label_parsed) +
    theme_bw(base_size = 9) +
    #theme(panel.grid = element_line(linetype = 0)) +
    labs(x = "year", y = "estimate") +
    expand_limits(y = 0)
  
  ### plot colours if required
  if(length(unique(input$model)) > 1) {
    p <- p + geom_line(aes(colour = model))
  } else {
    p <- p + geom_line()
  }
  p
}
### ------------------------------------------------------------------------ ###
### extract estimated parameters, format, plot
### ------------------------------------------------------------------------ ###

plot_list_params <- function(model_list, ### list with model fits
                             correct_names = TRUE,
                             list_label = "", ### label of elements for plotting
                             use_last_year_unconditionally = TRUE,
                             use_label_as_axis = FALSE, ### x-axis or none
                             use_label_numeric = FALSE, ### continuous axis
                             plot_error = FALSE ### error bar
){
  
  ### extract parameters from model list
  pars <- get_pars(model_list, correct_names = correct_names) 
  
  ### set variable
  pars <- lapply(seq_along(pars), FUN = function(x){
    res <- as.data.frame(pars[[x]])
    res$variable <- names(pars)[x]
    res$par <- row.names(res)
    res
  })
  ### concatenate
  pars <- do.call(rbind, pars)
  
  #remove some 
  pars <- pars[!rownames(pars) %in% c("pp","robfac"),]
  
  ### use last year for B/F estimates
  if(isTRUE(use_last_year_unconditionally)){
    pars$par <- sub(x = pars$par, pattern = "\\d{4}.{0,1}\\d{0,2}", "last_year")
  }
  
  ### create missing levels for plotting
  pars <- create_missing_levels(data_frame = pars, col_split = "variable", col_par = "par")
  
  ### plot parameters
  if(isTRUE(use_label_numeric)){
    p <- ggplot(data = pars, aes(fill = variable, group = variable, y = estimate, x = as.numeric(as.character(variable))))
  } else {
    p <- ggplot(data = pars, aes(fill = variable, group = variable, y = estimate, x = variable))   
  }
  p <- p + geom_bar(stat = "identity", position = "dodge", width = 1) +
    scale_fill_discrete(list_label) +
    facet_wrap(~ par, scale = "free") +
    labs(x = "", y = "estimate") +
    theme_bw(base_size = 9)
  
  ### label x-axis or not
  if(!isTRUE(use_label_as_axis)){
    p <- p + theme(axis.text.x = element_blank(), axis.ticks.x = element_blank())  
  }     
  
  ### plot errorbar if desired
  if(isTRUE(plot_error)){
    p <- p + geom_errorbar(aes(ymin = cilow, ymax = ciupp))
  }
  p
}


extrct_tmsrs <- function(model, ### SPiCT object
                         get_B = TRUE,     ### biomass
                         get_F = TRUE,     ### fishing mortality
                         get_BBmsy = TRUE,  ### B / Bmsy
                         get_FFmsty = TRUE ### F / Fmsy
){
  
  ### extract values
  res <- as.data.frame(rbind(
    cbind(get.par("logB", model, exp = TRUE), quant = "B"),
    cbind(get.par("logBBmsy", model, exp = TRUE), quant = "B/B[MSY]"),
    cbind(get.par("logF", model, exp = TRUE), quant = "F"),
    cbind(get.par("logFFmsy", model, exp = TRUE), quant = "F/F[MSY]")
  ))
  
  ### set year
  res$year <- row.names(res)
  row.names(res) <- NULL
  
  ### return
  res
  
}

### wrapper for list
extrct_tmsrs_lst <- function(models, ### list of SPiCT objects
                             get_B = TRUE,     ### biomass
                             get_F = TRUE,     ### fishing mortality
                             get_BBmsy = TRUE,  ### B / Bmsy
                             get_FFmsty = TRUE ### F / Fmsy
){
  
  ### check if list is result from retro analysis
  if(is(models, "spictcls") & is(models$retro, "list")){
    
    ### get final input year
    last_year <- models$inp$timeC[length(models$inp$timeC)]
    
    ### create list with last years assessment first
    res <- list(models)
    ### add retro runs
    res[2:c(length(models$retro)+1)] <- models$retro
    
    ### set names
    names(res) <- rev(as.character(seq(to = last_year, length.out = length(res))))
    models <- res
  }
  
  ### extract timeseries
  res <- lapply(seq_along(models), function(x){
    ### extract values
    res_temp <- extrct_tmsrs(models[[x]])
    ### set names
    res_temp$model <- names(models)[x]
    ### return
    res_temp
  })
  
  ### combine
  res <- do.call(rbind, res)
  
  ### return
  return(res)
  
}
### function for extracting parameters from list of model outputs
get_pars <- function(models, correct_names = FALSE){
  
  ### extract model parameter estimates
  pars <- lapply(models, FUN = function(x){
    rbind(sumspict.parest(x), 
          sumspict.states(x), 
          sumspict.srefpoints(x)[, -5],
          sumspict.drefpoints(x))
  })
  
  survey_names <- names(models)
  
  ### modify names so that they match
  if(isTRUE(correct_names)){
    
    pars <- lapply(seq_along(survey_names), function(x){
      
      ### current survey name(s)
      names_temp <- survey_names[x]
      names_temp <- unlist(strsplit(x = names_temp, split = "\\."))
      
      ### temp pars object
      pars_temp <- pars[[x]]
      
      ### add survey names to param names
      if(length(names_temp) == 1){
        
        ### correct parameter names
        row.names(pars_temp)[row.names(pars_temp) %in% c("alpha", "q", "sdi")] <-
          paste0(row.names(pars_temp)[row.names(pars_temp) %in% c("alpha", "q", "sdi")], "_", names_temp)
      } else {
        ### get parameter name positions to be replaced
        pos <- grepl(x = row.names(pars_temp), pattern = "*[a-zA-Z][0-9]{1,3}$")
        
        ### "loop" through surveys
        for(y in seq_along(names_temp)){
          
          ### get current survey names
          name_temp <- names_temp[y]
          
          ### parameter name positions with this survey
          pos_temp <- grepl(x = row.names(pars_temp)[pos], pattern = paste0("[a-zA-Z]", y, "$"))
          ### replace those
          row.names(pars_temp)[pos][pos_temp] <-                
            gsub(x = row.names(pars_temp)[pos][pos_temp], 
                 pattern = paste0(y, "$"),
                 replacement = paste0("_", name_temp))
        }
      }
      ### return
      pars_temp
    })}
  
  ### set names
  names(pars) <- survey_names
  ### return
  pars
}

### function for creating empty fields for barplot
create_missing_levels <- function(data_frame, 
                                  col_split, ### column name defining the split 
                                  col_par ### column name with parameter names
){
  
  splits <- unique(data_frame[, col_split])
  checks <- unique(data_frame[, col_par])
  for(x in seq_along(splits)){
    ### available parameters
    available <- data_frame[data_frame[col_split] == splits[x], col_par]
    ### missing parameters
    missing_levels <- setdiff(checks, available)
    ### if missing levels exist, add them
    if(length(missing_levels) > 0){
      ### create template
      data_add <- data_frame[seq_along(missing_levels),]
      ### fill with NAs
      data_add[] <- NA
      ### add splitting variable
      data_add[col_split] <- splits[x]
      ### add parameter names
      data_add[col_par] <- missing_levels
      
      ### add to data frame
      data_frame <- rbind(data_frame, data_add)
    }
  }
  ### return data frame
  return(data_frame)
}

#################################################
# main part of retro for model without prior r
#################################################
model_list_retro <- spict::retro(fit, nretroyear = 4)

### create list, including base model
model_list <- model_list_retro$retro
model_list[[length(model_list)+1]] <- model_list_retro
### names (last input year)
names(model_list) <- lapply(model_list, function(x){
  max(x$inp$timeC)
})

### check convergence
lapply(model_list, function(x){ x$opt$convergence})
lapply(model_list, function(x){ x$opt$message})

### plot time series
df_retro <- extrct_tmsrs_lst(model_list_retro)

dev.new(width = 10, height = 10)
plot_tmsrs(df_retro, label = "retro year")  

### plot parameter estimates       
dev.new(width = 10, height = 10)
plot_list_params(model_list = model_list, list_label = "retro year",
                 use_label_as_axis = TRUE, use_label_numeric = TRUE,
                 plot_error = FALSE)
### with errorbar      
dev.new(width = 10, height = 10)
plot_list_params(model_list = model_list, list_label = "retro year",
                 use_label_as_axis = TRUE, use_label_numeric = TRUE,
                 plot_error = TRUE)
