# Generation of U based on the stochastic EM algorithm
genU = function(k.m, k.y, benchmark, Iternum = 15, data, cvars, rvar, xvars, mvar, yvar, nsim = 20, num_cores = 1) {
  
  ## Generate b.m and b.y
  fit.y.rxc = lm(paste0(yvar, " ~ ", paste(c(rvar, xvars, cvars), collapse = " + ")), data = data)
  # Y ~ R + X + M + C
  fit.y.rxmc = lm(paste0(yvar, " ~ ", paste(c(rvar, xvars, mvar, cvars), collapse = " + ")), data = data)
  # M ~ R + X + C
  fit.m.rxc = glm(as.formula(paste0(mvar, " ~ ", paste(c(rvar, xvars, cvars), collapse = " + "))), family = binomial(logit), data = data)
  b.m.xj = coef(fit.m.rxc)[benchmark]
  b.y.xj = coef(fit.y.rxc)[benchmark]
  b.y.m = coef(fit.y.rxmc)[paste0(mvar, levels(data[[mvar]])[2])]
  sigma.m = sd(residuals(fit.m.rxc))
  sigma.xj = sd(residuals(lm(paste0(benchmark, " ~ ", paste(c(rvar, xvars, cvars), collapse = " + ")), data = data)))
  
  DATA.xj1 = DATA.xj0 = data
  DATA.xj1[, benchmark] = 1
  DATA.xj0[, benchmark] = 0
  
  b.m = log(k.m) + b.m.xj
  coef.m = c(coef(fit.m.rxc), U = b.m)
  u_fake = rnorm(length(data[[yvar]]), mean = 0, sd = sigma.xj)
  pre.m1 = 1 / (1 + exp(- (cbind(model.matrix(fit.m.rxc), U = u_fake + 1) %*% coef.m)))
  pre.m0 = 1 / (1 + exp(- (cbind(model.matrix(fit.m.rxc), U = u_fake) %*% coef.m)))
  diff.pi = mean(pre.m1 - pre.m0)
  rsq.u.m = (log(k.m) + b.m.xj)/sqrt((log(k.m) + b.m.xj)^2 + pi^2/(3 * sigma.m^2))
  b.y = (k.y * b.y.xj - b.y.m * diff.pi) *
    sigma.m / (sigma.m - abs(rsq.u.m/sqrt(1 - rsq.u.m)) * sigma.xj * diff.pi)
  
  ## Generate continuous U from its prior distribution first
  U = rnorm(nrow(data), 0, 0.5) #sigma.u
  #sigma.urxc: sens para
  
  coef.y.updated = NULL
  coef.m.updated = NULL
  
  ## EM steps
  for(iter in 1:Iternum) {
    cat(iter, "\r")
    
    ## scale.y = "continuous"
    l.y = lm(as.formula(paste0(yvar, " ~ ", paste(c(cvars,rvar, xvars, mvar, paste0(mvar, ":", h1vars)), collapse = " + "))),
             offset = b.y * U, data = data)
    coef.y = c(l.y$coef, U = b.y)
    sd.y = sigma(l.y) # This is the sd of Y conditional on t, m, C, X, and U.
    df.y = summary(l.y)$df[2]
    
    ## scale.m = "binary"
    l.m = glm(as.formula(paste0(mvar, " ~ ", paste(c(rvar, xvars, cvars), collapse = " + "))),
              offset = b.m * U, data = data, family = binomial(link = "logit"))
    coef.m = c(l.m$coef, U = b.m)
    
    ## For continuous U
    # Calculate the integral in the denominator
    integral_res <- mclapply(1:nrow(data), FUN = function(i) {
      integrand = function(u){
        mvar_numeric <- as.numeric(data[i, mvar]) - 1
        dnorm(data[i, yvar], mean = c(model.matrix(l.y)[i, ] %*% l.y$coef) + u * b.y, sd = sd.y) *
          (1/(1 + exp(-(c(model.matrix(l.m)[i, ] %*% l.m$coef) + u * b.m))))^mvar_numeric *
          (1 - 1/(1 + exp(-(c(model.matrix(l.m)[i, ] %*% l.m$coef) + u * b.m))))^(1 - mvar_numeric) *
          dnorm(u, mean = 0, sd = 1)
      }
      res = integrate(integrand, lower = -10, upper = 10)$value
    }, mc.cores = num_cores)
    integral = unlist(integral_res)
    
    # Generate random values of U
    U_res <- mclapply(1:nrow(data), FUN = function(i) {
      # = NULL
      #U.final.nsim = matrix(NA, nrow(data), nsim)
      # Obtain the condition probability of U
      conditional.u = function(u){
        mvar_numeric <- as.numeric(data[i, mvar]) - 1
        dnorm(data[i, yvar], mean = c(model.matrix(l.y)[i, ] %*% l.y$coef) + u * b.y, sd = sd.y) *
          (1/(1 + exp(-(c(model.matrix(l.m)[i, ] %*% l.m$coef) + u * b.m))))^mvar_numeric *
          (1 - 1/(1 + exp(-(c(model.matrix(l.m)[i, ] %*% l.m$coef) + u * b.m))))^(1 - mvar_numeric) *
          dnorm(u, mean = 0, sd = 1)/integral[i] #sigma.u
      }
      dist = AbscontDistribution(d = conditional.u)  # signature for a dist with pdf ~ conditional.u
      rdist = r(dist)# function to create random variates from conditional.u
      if(iter < Iternum){
        U = rdist(1)
        U.final.nsim = NA
      } else if (iter == Iternum){
        U = rdist(1)
        U.final.nsim = rdist(nsim)
      }
      return(list(U = U, U.final.nsim = U.final.nsim))
    }, mc.cores = num_cores)
    
    U = sapply(U_res, "[[", "U")
    U.final.nsim = sapply(U_res, "[[", "U.final.nsim")
    data$U = U #apply(U_res, 1, mean)
    
  }
  coef.y.updated = rbind(coef.y.updated, coef.y)
  coef.m.updated = rbind(coef.m.updated, coef.m)
  
  return(list(U = U, b.m= b.m, b.y=b.y,
              coef.y.updated = coef.y.updated,
              coef.m.updated = coef.m.updated,
              U.final.nsim = U.final.nsim))
}

reg1_no_int <- function(b.m = NULL, b.y = NULL, B, data,cvars, rvar, xvars, mvar, yvar, h1vars, cluster){
  
  est1 <- function(data){  
    data_R0 <- subset(data, data[[rvar]] == "white") # change this to reference group, @Karen
    data_R1 <- subset(data, data[[rvar]] == "black") # change this to the remaining groups, @Karen
    # weighting
    moPropen1 <- buildModelObj(model = as.formula(paste0(" ~ ", paste(c(cvars,xvars, "U"), collapse = " + "))),
                               solver.method = "glm",
                               solver.args = list("family" = "binomial"),
                               predict.method = 'predict.glm',
                               predict.args = list(type = "response"))

    moClass1 <- buildModelObj(model = as.formula(paste0(" ~ ", paste(c(h1vars), collapse = " + "))),
                              solver.method = "rpart",
                              solver.args = list(method = "class"),
                              predict.args = list(type = 'class'))

    fitFS1 <- optimalClass(moPropen = moPropen1,
                           moClass = moClass1,
                           data = data.frame(data_R1), response = data_R1[[yvar]],
                           txName = mvar, verbose = FALSE)
    fitFS0 <- optimalClass(moPropen = moPropen1,
                           moClass = moClass1,
                           data = data.frame(data_R0), response = data_R0[[yvar]],
                           txName = mvar, verbose = FALSE)
    
    
    data$opt.M <- ifelse(data[[rvar]] == "white", (optTx(fitFS0)$optimalTx),(optTx(fitFS1)$optimalTx)) #@Karen, change "white" to the reference group

    ta <- table(data$opt.M)
    p.opt <- ta[2] / (ta[2] + ta[1]) * 100
    # Construct weight for mediator
    DATA1 <- data   
    DATA1 <- DATA1 %>%
      mutate(Ind = (DATA1[[mvar]] == opt.M))
    
    fit.m <- glm(as.formula(paste0(mvar, " ~ ", paste(c(rvar,cvars,xvars), collapse = " + "))), offset = b.m * DATA1$U,
                 family = binomial(logit), data = DATA1)
    coef.m <- c(fit.m$coef, U = b.m)
    pre.m <- 1 / (1 + exp(-(cbind(model.matrix(fit.m), U = DATA1$U) %*% coef.m)))
    p.med <- ifelse(DATA1[[mvar]] == 0, 1 - pre.m, pre.m)
    
    # Weights
    DATA1$w <- 1 / p.med
    
    # Calculate zeta_ICDE
    fit <- lm(paste0(yvar, " ~ ", paste(c(cvars,rvar), collapse = " + ")), data=DATA1, weights = Ind * w)
    estimated_zeta_ICDE <- coef(fit)[3] #@Karen, change this to capture the coefficient of race
    
    # IIE: Method 2
    fit.lm1 <- lm(paste0(yvar, " ~ ", paste(c(cvars), collapse = " + ")), data=data_R1)
    W_Y1 <- coef(fit.lm1)[1]
    
    fit.lm0 <- lm(paste0(yvar, " ~ ", paste(c(cvars), collapse = " + ")), data=data_R0)
    W_Y0 <- coef(fit.lm0)[1]
    
    fit.lm <- lm(paste0(yvar, " ~ ", paste(c(cvars,rvar,"Ind", paste0(rvar, ":Ind")), collapse = " + ")), data = DATA1, weights = w) 
    fit.Ind <- glm(paste0("Ind", " ~ ", paste(c(rvar), collapse = " + ")), data = DATA1, family = binomial(logit))
    
    alpha1 <- plogis(sum(coef(fit.Ind)[1:2])) - plogis(coef(fit.Ind)[1])
    
    coef_names <- names(coef(fit.lm))
    ind_coef <- coef_names[grepl(paste0("^", "Ind"), coef_names)]
    interaction_coef <- coef_names[grepl(paste0(":", "Ind"), coef_names)]
    
    reg_delta_IIE <- alpha1 * sum(coef(fit.lm)[c(ind_coef, interaction_coef)])
    reg_zeta_IIE <- (W_Y1 - W_Y0) - reg_delta_IIE
    
    return(list(fitFS1 = fitFS1, fitFS0 = fitFS0, p.opt = p.opt,
                estimated_zeta_ICDE = estimated_zeta_ICDE,
                reg_zeta_IIE = reg_zeta_IIE,
                reg_delta_IIE = reg_delta_IIE))
  }
  
  est.ori <- est1(data)
  
  reg_delta_IIE <- est.ori$reg_delta_IIE
  reg_zeta_IIE <- est.ori$reg_zeta_IIE
  estimated_zeta_ICDE <- est.ori$estimated_zeta_ICDE
  p.opt <- est.ori$p.opt
  fitFS1 <- est.ori$fitFS1
  fitFS0 <- est.ori$fitFS0
  
  regzeta_IIE.boot <- regdelta_IIE.boot <- zeta_ICDE.boot <- rep(NA, B)
  
  for(i in 1:B){
    clusters <- unique(data[,cluster]) # @Karen, if they did not specify anything for cluster, it should operate standard bootsrapping
    m <- length(clusters)
    units <- sample(clusters, size = length(clusters), replace = TRUE)
    df.bs <- sapply(units, function(x) which(data[,cluster] == x))
    sb <- unlist(df.bs)
    data.boot <- data[sb, ]
    
    est.boot <- est1(data.boot)
    
    zeta_ICDE.boot[i] <- est.boot$estimated_zeta_ICDE
    regzeta_IIE.boot[i] <- est.boot$reg_zeta_IIE
    regdelta_IIE.boot[i] <- est.boot$reg_delta_IIE
  }
  
  SE_zeta_ICDE <- sd(zeta_ICDE.boot)
  SE_regdelta_IIE <- sd(regdelta_IIE.boot)
  SE_regzeta_IIE <- sd(regzeta_IIE.boot)
  
  results <- list(
    p.opt = p.opt, fitFS0 = fitFS0, fitFS1 = fitFS1,
    estimated_zeta_ICDE = as.numeric(estimated_zeta_ICDE), SE_zeta_ICDE = SE_zeta_ICDE,
    reg_delta_IIE = as.numeric(reg_delta_IIE), SE_regdelta_IIE = SE_regdelta_IIE, 
    reg_zeta_IIE = as.numeric(reg_zeta_IIE), SE_regzeta_IIE = SE_regzeta_IIE
  )
  
  return(results)
}

ind.decomp <- function(B, data,cvars, rvar, xvars, mvar, yvar, h1vars, cluster){
  
  est1 <- function(data){  
    
    data_R0 <- subset(data, data[[rvar]] == "white") # change this to reference group, @KAren
    data_R1 <- subset(data, data[[rvar]] == "black") # change this to the remaining groups, @KAren
    
    # weighting
    moPropen1 <- buildModelObj(model = as.formula(paste0(" ~ ", paste(c(cvars,xvars), collapse = " + "))),
                               solver.method = "glm",
                               solver.args = list("family" = "binomial", control = glm.control(maxit = 50)),
                               predict.method = 'predict.glm',
                               predict.args = list(type = "response"))

    moClass1 <- buildModelObj(model = as.formula(paste0(" ~ ", paste(c(h1vars), collapse = " + "))),
                              solver.method = "rpart",
                              solver.args = list(method = "class",  control = rpart.control(minsplit = 200, cp = 0.01)),
                              predict.args = list(type = 'class'))

    fitFS1 <- optimalClass(moPropen = moPropen1,
                           moClass = moClass1,
                           data = data.frame(data_R1), response = data_R1[[yvar]],
                           txName = mvar, verbose = FALSE)
    
    fitFS0 <- optimalClass(moPropen = moPropen1,
                          moClass = moClass1,
                          data = data.frame(data_R0), response = data_R0[[yvar]],
                          txName = mvar, verbose = FALSE)
    

    data$opt.M <- ifelse(data[[rvar]] == "white", optTx(fitFS0)$optimalTx,optTx(fitFS1)$optimalTx) #@Karen, change "white" to the reference group
    ta <- table(data$opt.M)
    p.opt <- ta[2]/(ta[2] + ta[1]) * 100
    # rpart.plot(classif(object=fitFS0))
    DATA1 <- data   
    DATA1 <- DATA1 %>%
      mutate(Ind = (DATA1[[mvar]] == DATA1$opt.M))
    
    fit.m <- glm(as.formula(paste0(mvar, " ~ ", paste(c(rvar,cvars,xvars), collapse = " + "))),
                 family = binomial(logit), data = DATA1)
    coef.m <- c(fit.m$coef)
    pre.m <- 1 / (1 + exp(- (cbind(model.matrix(fit.m)) %*% coef.m)))
    p.med <- ifelse(DATA1[[mvar]] == 0, 1 - pre.m, pre.m)
    
    DATA1$w <- 1 / p.med
    
    #zeta_ICDE
    fit <- lm(paste0(yvar, " ~ ", paste(c(cvars,rvar), collapse = " + ")), data=DATA1, weights = Ind * w)
    estimated_zeta_ICDE <- coef(fit)[3] #@KAren, change this to capture the coefficient of race
    
    #zeta_IIE
    fit.lm1 <- lm(paste0(yvar, " ~ ", paste(c(cvars), collapse = " + ")), data=data_R1)
    W_Y1 <- coef(fit.lm1)[1]
    
    fit.lm0 <- lm(paste0(yvar, " ~ ", paste(c(cvars), collapse = " + ")), data=data_R0)
    W_Y0 <- coef(fit.lm0)[1]
    
    fit.lm <- lm(paste0(yvar, " ~ ", paste(c(cvars,rvar,"Ind", paste0(rvar, ":Ind")), collapse = " + ")), data = DATA1, weights = w) 
    fit.Ind <- glm(paste0("Ind", " ~ ", paste(c(rvar), collapse = " + ")), data = DATA1, family = binomial(logit))
    
    alpha1 <- plogis(sum(coef(fit.Ind)[1:2])) - plogis(coef(fit.Ind)[1])
    
    coef_names <- names(coef(fit.lm))
    ind_coef <- coef_names[grepl(paste0("^", "Ind"), coef_names)]
    interaction_coef <- coef_names[grepl(paste0(":", "Ind"), coef_names)]
    
    reg_delta_IIE <- alpha1 * sum(coef(fit.lm)[c(ind_coef, interaction_coef)])
    reg_zeta_IIE <- (W_Y1 - W_Y0) - reg_delta_IIE
    
    return(list(fitFS0 = fitFS0,fitFS1 = fitFS1, p.opt = p.opt,
                estimated_zeta_ICDE = estimated_zeta_ICDE,
                reg_zeta_IIE = reg_zeta_IIE,
                reg_delta_IIE = reg_delta_IIE))
  }
  
  est.ori <- est1(data)
  
  reg_delta_IIE <- est.ori$reg_delta_IIE
  reg_zeta_IIE <- est.ori$reg_zeta_IIE
  estimated_zeta_ICDE <- est.ori$estimated_zeta_ICDE
  p.opt <- est.ori$p.opt
  fitFS1 <- est.ori$fitFS1
  fitFS0 <- est.ori$fitFS0
  
  regzeta_IIE.boot <- regdelta_IIE.boot <- zeta_ICDE.boot <- rep(NA, B)
  
  for(i in 1:B){
    clusters <- unique(data[, cluster]) # @ Karen, please make an option to do clustered boostrap or standard boostrap (null for cluster)
    units <- sample(clusters, size = length(clusters), replace = TRUE)
    df.bs <- sapply(units, function(x) which(data[, cluster] == x))
    sb <- unlist(df.bs)
    data.boot <- data[sb, ]
    
    est.boot <- est1(data.boot)
    
    zeta_ICDE.boot[i] <- est.boot$estimated_zeta_ICDE
    regzeta_IIE.boot[i] <- est.boot$reg_zeta_IIE
    regdelta_IIE.boot[i] <- est.boot$reg_delta_IIE
  }
  
  SE_zeta_ICDE <- sd(zeta_ICDE.boot)
  SE_regdelta_IIE <- sd(regdelta_IIE.boot)
  SE_regzeta_IIE <- sd(regzeta_IIE.boot)
  
  #Initial disparity
  results_ini <- lm(as.formula(paste0(yvar, " ~ ", paste(c(rvar,cvars), collapse = " + "))), data=data)
  estimated_tau <- coef(results_ini)[2]
  SE_tau <- summary(results_ini)$coefficients[2, 2]
  
  results <- list(# result table
    p.opt = p.opt, fitFS1 = fitFS1, fitFS0 = fitFS0,
    estimated_tau = estimated_tau, SE_tau = SE_tau, #Initial disparity
    estimated_zeta_ICDE = as.numeric(estimated_zeta_ICDE), SE_zeta_ICDE = SE_zeta_ICDE, #disparity remaining_ICDE
    reg_delta_IIE = as.numeric(reg_delta_IIE), SE_regdelta_IIE = SE_regdelta_IIE, #disparity reduction_IIE
    reg_zeta_IIE = as.numeric(reg_zeta_IIE), SE_regzeta_IIE = SE_regzeta_IIE #disparity remaining_IIE
  )
  # @ Karen, Create a table similar to Table 2, add % of recommended as Table 4
  return(results)
}

ind.decompsense <- function(k.y, k.m, Iternum, data, cvars, rvar, xvars, mvar, yvar,
                        nsim, num_cores, B, h1vars, benchmark = "ses", cluster = "SCH_ID") {
  # 1. Generate U
  genU.res <- genU(
    k.y = k.y, k.m = k.m, Iternum = Iternum,
    benchmark = benchmark,
    data = data, cvars = cvars, rvar = rvar, xvars = xvars,
    mvar = mvar, yvar = yvar, nsim = nsim, num_cores = num_cores
  )
  
  # 2. Define a helper function to run the adjusted causal decomposition on each bootstrap sample
  simOne <- function(ss) {
    DATA.newU <- data
    # Insert the generated U for this iteration
    DATA.newU$U <- genU.res$U.final.nsim[ss, ]
    b.m <- genU.res$b.m
    b.y <- genU.res$b.y
    
    # Run your adjusted causal decomposition
    res_reg <- reg1_no_int(
      b.y = b.y, b.m = b.m, B = B,
      data = DATA.newU, cvars = cvars, rvar = rvar, xvars = xvars,
      mvar = mvar, yvar = yvar, h1vars = h1vars, cluster = cluster
    )
    return(res_reg)
  }
  
  # 3. Run simOne() 'nsim' times in parallel
  #    This captures the uncertainty of generating U
  results_ind <- parallel::mclapply(
    X = 1:nsim,
    FUN = function(ss) {
      res <- try(simOne(ss), silent = TRUE)
      if (inherits(res, "try-error")) return(NULL) else return(res)
    },
    mc.cores = num_cores
  )
  
  # 4. Aggregate results into results_hom1
  #    We calculate mean estimates and standard errors
  #    using the standard formula with the extra (1 + 1/nsim) factor for variance
  results_hom1 <- list(
    p.opt = mean(
      sapply(results_ind, function(x) if (!is.null(x[["p.opt"]])) x[["p.opt"]] else NA),
      na.rm = TRUE
    ),
    estimated_zeta_ICDE = mean(
      sapply(results_ind, function(x) if (!is.null(x[["estimated_zeta_ICDE"]])) x[["estimated_zeta_ICDE"]] else NA),
      na.rm = TRUE
    ),
    estimated_delta_IIE = mean(
      sapply(results_ind, function(x) if (!is.null(x[["reg_delta_IIE"]])) x[["reg_delta_IIE"]] else NA),
      na.rm = TRUE
    ),
    estimated_zeta_IIE = mean(
      sapply(results_ind, function(x) if (!is.null(x[["reg_zeta_IIE"]])) x[["reg_zeta_IIE"]] else NA),
      na.rm = TRUE
    ),
    
    SE_zeta_ICDE = sqrt(
      mean(
        sapply(results_ind, function(x) if (!is.null(x[["SE_zeta_ICDE"]])) x[["SE_zeta_ICDE"]]^2 else NA),
        na.rm = TRUE
      ) +
        (1 + 1 / nsim) * var(
          sapply(results_ind, function(x) if (!is.null(x[["SE_zeta_ICDE"]])) x[["SE_zeta_ICDE"]] else NA),
          na.rm = TRUE
        )
    ),
    SE_delta_IIE = sqrt(
      mean(
        sapply(results_ind, function(x) if (!is.null(x[["SE_regdelta_IIE"]])) x[["SE_regdelta_IIE"]]^2 else NA),
        na.rm = TRUE
      ) +
        (1 + 1 / nsim) * var(
          sapply(results_ind, function(x) if (!is.null(x[["SE_regdelta_IIE"]])) x[["SE_regdelta_IIE"]] else NA),
          na.rm = TRUE
        )
    ),
    SE_zeta_IIE = sqrt(
      mean(
        sapply(results_ind, function(x) if (!is.null(x[["SE_regzeta_IIE"]])) x[["SE_regzeta_IIE"]]^2 else NA),
        na.rm = TRUE
      ) +
        (1 + 1 / nsim) * var(
          sapply(results_ind, function(x) if (!is.null(x[["SE_regzeta_IIE"]])) x[["SE_regzeta_IIE"]] else NA),
          na.rm = TRUE
        )
    )
  )
  
  # 5. Return the final list of results
  # @ Karen, Create a table similar to Table 4
  return(results_hom1)
}
