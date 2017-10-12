
library(spict)
load("D:/wg_IBPTur.27.4/refpoints/spict/inp_spict.rdata")

inp <- check.inp(inp)

## Three indices so logsdi prior must be given as a list with three vectors.
## LPUE index is the second vector, only enable this one.
inp$priors$logsdi = list(  c(1,1,0), c(0.1,0.5,1), c(1,1,0))

fit <- fit.spict(inp)

plot(fit)

summary(fit)
