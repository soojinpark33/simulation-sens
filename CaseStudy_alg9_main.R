## Clear workspace
rm(list = ls())

## Load packages
library(geepack)
library(DynTxRegime)
library(boot)
library(rpart)
library(rpart.plot)
library(parallel)
library(dplyr)
library(tidyverse)
library(mice)
library(readxl)
library(distr) # AbscontDistribution()
library(rlang) #exprs()


source("CaseStudy_alg9_Functions.R")

#::::::::::::::::::::::::::
#:::: Data processing :::::
#::::::::::::::::::::::::::
load("data_imputed.RData")

# Select only Black and White students
data <- data_impute %>%
  filter(race %in% c("black", "white")) %>%
  mutate(female = scale(as.numeric(female), center = TRUE),
         S1M8GRADE = as.factor(S1M8GRADE)) #center baseline covariates

# determine variables for C,R,X,H1, M, and Y
cvars <- exprs(female)
rvar <- "race" # it should be factor, @Karen, if it is not factor, switch it to factor in the function.
xvars <- exprs(S1M8GRADE,ses,X1PAR1OCC_STEM1,X1PAR2OCC_STEM1,P1EDUEXPECT,P1EDUASPIRE,M1SEX,M1INTEREST,locale,
               X1SCHOOLCLI,X2PROBLEM,M1UNPREPPCT,requirements,schextracur,friend,
                 X1MTHINT,X1MTHUTI,X1MTHEFF,X1SCHOOLENG,X1SCHOOLBEL)
h1vars <- exprs(S1M8GRADE,X1MTHEFF,X1MTHINT)
mvar <- "algebra1" # it should be factor
yvar <- "X2TXMTH" # it should be continuous

## Setting (Change to numbers after #, we used the samllest number to reduce the run time)
num_cores = 16
nsim = 20
B = 500
Iternum = 15
k.m = k.y = 1 
# k.m = 1; k.y = -1 

## Unadjusted Results
results_unadj <- ind.decomp(B = B, data = data, cvars = cvars, rvar = rvar, xvars = xvars,
                          mvar = mvar, yvar = yvar, h1vars = h1vars, cluster = "SCH_ID")

## Adjusted Results
results_adj <- ind.decompsense(k.y, k.m, Iternum, data, cvars, rvar, xvars, mvar, yvar,
                        nsim, num_cores, B, h1vars, benchmark = "ses", cluster = "SCH_ID") 
  
  
  


