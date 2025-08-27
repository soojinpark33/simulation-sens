## Clear workspace
rm(list = ls())

## Load packages
library(causal.decomp)
library(tidyverse)

## Load data
load("data_imputed.RData")

## Select only Black and White students
data <- data_impute %>%
  filter(race %in% c("black", "white")) %>%
  droplevels() %>%
  mutate(female = scale(as.numeric(female), center = TRUE),
         S1M8GRADE = as.factor(S1M8GRADE)) # center baseline covariates

## Unadjusted Results
results_unadj <- ind.decomp(outcome = "X2TXMTH",      # continuous
                            group = "race",           # factor
                            group.ref = "white",
                            risk.factor = "algebra1", # factor
                            intermediates = c(
                              "S1M8GRADE", "ses", "X1PAR1OCC_STEM1", "X1PAR2OCC_STEM1", "P1EDUEXPECT",
                              "P1EDUASPIRE", "M1SEX", "M1INTEREST", "locale", "X1SCHOOLCLI",
                              "X2PROBLEM", "M1UNPREPPCT", "requirements", "schextracur", "friend",
                              "X1MTHINT", "X1MTHUTI", "X1MTHEFF", "X1SCHOOLENG", "X1SCHOOLBEL"
                            ),
                            moderators = c("S1M8GRADE", "X1MTHEFF", "X1MTHINT"),
                            covariates = "female",
                            data = data, B = 500, cluster = "SCH_ID")

## Adjusted Results
results_adj <- ind.sens(k.y = 1, k.m = 1, Iternum = 15,
                        outcome = "X2TXMTH",      # continuous
                        group = "race",           # factor
                        group.ref = "white",
                        risk.factor = "algebra1", # factor
                        intermediates = c(
                          "S1M8GRADE", "ses", "X1PAR1OCC_STEM1", "X1PAR2OCC_STEM1", "P1EDUEXPECT",
                          "P1EDUASPIRE", "M1SEX", "M1INTEREST", "locale", "X1SCHOOLCLI",
                          "X2PROBLEM", "M1UNPREPPCT", "requirements", "schextracur", "friend",
                          "X1MTHINT", "X1MTHUTI", "X1MTHEFF", "X1SCHOOLENG", "X1SCHOOLBEL"
                        ),
                        moderators = c("S1M8GRADE", "X1MTHEFF", "X1MTHINT"),
                        benchmark = "ses",
                        covariates = "female",
                        data = data, B = 500, cluster = "SCH_ID", nsim = 20, mc.cores = 8) 
