# Function to compute the coverage rate using pre-computed empirical confidence bounds
calculate_coverage_rate <- function(lower_bounds, upper_bounds, true_value, REP) {
  
  # Check if the true value is contained within each confidence interval
  contains_true_value <- (true_value >= lower_bounds) & (true_value <= upper_bounds)
  REP.nonNA <- sum(!is.na(contains_true_value))
  
  # Calculate the 95% coverage rate
  coverage_rate <- sum(contains_true_value, na.rm = TRUE) / REP.nonNA * 100

  return(coverage_rate)
}

reg1 = function(Y, R, M, X1, X2, X3, C, b.m = NULL, b.y = NULL, B, data){
  
  est1 = function(data){  
 # weighting
    moPropen1 <- buildModelObj(model = ~ R + C + X1 + X2 + X3 +U,
                               solver.method = "glm",
                               solver.args = list("family"="binomial"),
                               predict.method = 'predict.glm',
                               predict.args = list(type="response"))

    moClass1 <- buildModelObj(model = ~R + C + X1 + X2 + X3  ,
                              solver.method = "rpart",
                              solver.args = list(method="class"),
                              predict.args = list(type='class'))
    #
    fitFS <- optimalClass(moPropen = moPropen1,
                              #moMain = moMain, moCont = moCont,
                              moClass = moClass1,
                              data=data, response = data$Y, txName = "M", verbose = F)

    #:::::::::::::::::::::::::::::::#  
    #::::::: ICDE:Regression::::::::#
    #:::::::::::::::::::::::::::::::# 
    ### Assign optimal decision ###
    data$opt_M <- optTx(fitFS)$optimalTx
    a<- table(data$opt_M, data$Mopt) 
    accuracy <- (a[1,1]+a[2,2])/length(data$opt_M)
    DATA1 <- data   
    DATA1 <- DATA1 %>% mutate(Ind = (M == opt_M))

    ### construct weight for med ###
    fit.m <- glm(M ~ R + C + X1 +X2 + X3, offset = b.m * DATA1$U, family = binomial(logit), data = DATA1)
    coef.m = c(fit.m$coef, U = b.m)
    pre.m = 1/(1 + exp(-(cbind(model.matrix(fit.m), U = DATA1$U) %*% coef.m)))
    p.med <- ifelse(DATA1$M == 0,
                     1 - pre.m,
                     pre.m)

    ### Weights ###
    DATA1$w <- 1/p.med

    # Calculate zeta_icde #
    fit <- lm(Y ~ R * C, data=DATA1, weights = Ind*w)
    estimated_zeta_ICDE <- coef(fit)[2]

    #:::::::::::::::::::::::::::::::::::::::#
    #:::::::::::: IIE: Method 2 ::::::::::::#
    #:::::::::::::::::::::::::::::::::::::::#
    DATA_R0 <- subset(DATA1, R==0) # white
    DATA_R1 <- subset(DATA1, R==1) # black
    
    # E [Y | R=1, C=c] #
    W_Y1 <- DATA1 %>% filter(R == 1, C == 0) %>% summarise(W_Y1=mean(Y))
    W_Y0 <- DATA1 %>% filter(R == 0, C == 0) %>% summarise(W_Y0=mean(Y)) 
    
    # fit.lm1 <- lm(Y ~  C, data=DATA_R1)
    # W_Y1 <- coef(fit.lm1)[1]
    # # E [Y | R=0, C=c] #
    # fit.lm0 <- lm(Y ~ C, data=DATA_R0)
    # W_Y0 <- coef(fit.lm0)[1]
    
    ### Step 1 ###
    # Among white (R=0) #
    fit.I <- glm(Ind ~ C, data = DATA_R0, family = binomial(logit))

    # P [ I(M=d.opt)=1 | R=0, C=c ] # 
    pi_I1 <- plogis(coef(fit.I)[1])
    # P [ I(M=d.opt)=0 | R=0, C=c ] # 
    pi_I0 <- 1 - pi_I1

    ### Step 2 ###
    W_IIE.1 <- DATA1$Ind * DATA1$w
    W_IIE.0 <- (!DATA1$Ind) * DATA1$w
    
    # (theta1) = (1) #
    idx_R1 <- which(DATA1$R == 1)
    fit.lm1 <- lm(Y ~ C, data=DATA_R1, weights = W_IIE.1[idx_R1])
    term1 <- pi_I1 * coef(fit.lm1)[1]
    # (theta1) = (0) #
    fit.lm0 <- lm(Y ~ C, data=DATA_R1, weights = W_IIE.0[idx_R1])
    term0 <- pi_I0 * coef(fit.lm0)[1]
    
    W_iie.Y <- term1 + term0 
    
    estimated_delta_IIE <- W_Y1 - W_iie.Y 
    estimated_zeta_IIE <- W_iie.Y - W_Y0
    
    #regression
    fit.lm <- lm(Y ~ R + C + Ind , data=DATA1, weights=w)
    fit.Ind <- glm(Ind ~ R * C, data = DATA1, family = binomial(logit))

    alpha1 <- plogis(sum(coef(fit.Ind)[1:2])) - plogis(coef(fit.Ind)[1])


    reg_delta_IIE <- (alpha1*coef(fit.lm)["IndTRUE"])
    reg_zeta_IIE <- (W_Y1 - W_Y0) - reg_delta_IIE
  
    return(list(accuracy=accuracy, estimated_zeta_ICDE= estimated_zeta_ICDE, estimated_delta_IIE=estimated_delta_IIE,
                estimated_zeta_IIE=estimated_zeta_IIE, reg_zeta_IIE= reg_zeta_IIE, reg_delta_IIE=reg_delta_IIE   ))
}

est.ori = est1(data)

accuracy = est.ori$accuracy
reg_delta_IIE = est.ori$reg_delta_IIE
reg_zeta_IIE = est.ori$reg_zeta_IIE
estimated_delta_IIE = est.ori$estimated_delta_IIE
estimated_zeta_IIE = est.ori$estimated_zeta_IIE
estimated_zeta_ICDE = est.ori$estimated_zeta_ICDE


regzeta_IIE.boot <- regdelta_IIE.boot <- zeta_ICDE.boot <-zeta_IIE.boot <-delta_IIE.boot <- rep(NA, B)

for(i in 1:B){
  m <- nrow(data)
  data.boot = data[sample(1:nrow(data), m, replace = T), ]
  est.boot = est1(data.boot)

  delta_IIE.boot[i] = est.boot$estimated_delta_IIE
  zeta_IIE.boot[i] = est.boot$estimated_zeta_IIE  
  zeta_ICDE.boot[i] = est.boot$estimated_zeta_ICDE
  regzeta_IIE.boot[i] = est.boot$reg_zeta_IIE
  regdelta_IIE.boot[i] = est.boot$reg_delta_IIE
}

SE_delta_IIE = sd(unlist(delta_IIE.boot), na.rm=TRUE)
SE_zeta_IIE = sd(unlist(zeta_IIE.boot), na.rm=TRUE)
SE_zeta_ICDE = sd(unlist(zeta_ICDE.boot), na.rm=TRUE)
SE_regdelta_IIE = sd(unlist(regdelta_IIE.boot), na.rm=TRUE)
SE_regzeta_IIE = sd(unlist(regzeta_IIE.boot), na.rm=TRUE)


results = c(accuracy = as.numeric(accuracy), estimated_delta_IIE = as.numeric(estimated_delta_IIE), SE_delta_IIE = SE_delta_IIE, 
            estimated_zeta_IIE = as.numeric(estimated_zeta_IIE), SE_zeta_IIE = SE_zeta_IIE,
            estimated_zeta_ICDE= as.numeric(estimated_zeta_ICDE), SE_zeta_ICDE=SE_zeta_ICDE,
            reg_delta_IIE = as.numeric(reg_delta_IIE), SE_regdelta_IIE = SE_regdelta_IIE, 
            reg_zeta_IIE = as.numeric(reg_zeta_IIE), SE_regzeta_IIE = SE_regzeta_IIE
            )

return(results)
}

reg2 = function(Y, R, M, X1, X2, X3, C, b.m = NULL, b.y = NULL, B, data){
  
  est1 = function(data){  
    # weighting
    moPropen1 <- buildModelObj(model = ~ R + C + X1 + X2 + X3,
                               solver.method = "glm",
                               solver.args = list("family"="binomial"),
                               predict.method = 'predict.glm',
                               predict.args = list(type="response"))

    moClass1 <- buildModelObj(model = ~R + C + X1 + X2 + X3 ,
                              solver.method = "rpart",
                              solver.args = list(method="class"),
                              predict.args = list(type='class'))
    #
    fitFS <- optimalClass(moPropen = moPropen1,
                              #moMain = moMain, moCont = moCont,
                              moClass = moClass1,
                              data=data, response = data$Y, txName = "M", verbose = F)

    
    #:::::::::::::::::::::::::::::::#  
    #::::::: ICDE:Regression::::::::#
    #:::::::::::::::::::::::::::::::# 
    ### Assign optimal decision ###
    data$opt_M <- optTx(fitFS)$optimalTx
    a<- table(data$opt_M, data$Mopt) 
    accuracy <- (a[1,1]+a[2,2])/length(data$opt_M)
    DATA1 <- data   
    DATA1 <- DATA1 %>% mutate(Ind = (M == opt_M))
    
    ### construct weight for med ###
    fit.m <- glm(M ~ R + C + X1 +X2 + X3, family = binomial(logit), data = DATA1)
    coef.m = c(fit.m$coef)
    pre.m = 1/(1 + exp(-(cbind(model.matrix(fit.m)) %*% coef.m)))
    p.med <- ifelse(DATA1$M == 0,
                    1 - pre.m,
                    pre.m)
    
    ### Weights ###
    DATA1$w <- 1/p.med
    
    # Calculate zeta_icde #
    fit <- lm(Y ~ R * C, data=DATA1, weights = Ind*w)
    estimated_zeta_ICDE <- coef(fit)[2]
    
    #:::::::::::::::::::::::::::::::::::::::#
    #:::::::::::: IIE: Method 2 ::::::::::::#
    #:::::::::::::::::::::::::::::::::::::::#
    DATA_R0 <- subset(DATA1, R==0) # white
    DATA_R1 <- subset(DATA1, R==1) # black
    
    # E [Y | R=1, C=c] #
    W_Y1 <- DATA1 %>% filter(R == 1, C == 0) %>% summarise(W_Y1=mean(Y))
    W_Y0 <- DATA1 %>% filter(R == 0, C == 0) %>% summarise(W_Y0=mean(Y)) 
    # fit.lm1 <- lm(Y ~  C, data=DATA_R1)
    # W_Y1 <- coef(fit.lm1)[1]
    # # E [Y | R=0, C=c] #
    # fit.lm0 <- lm(Y ~ C, data=DATA_R0)
    # W_Y0 <- coef(fit.lm0)[1]
    # 
    ### Step 1 ###
    # Among white (R=0) #
    fit.I <- glm(Ind ~ C, data = DATA_R0, family = binomial(logit))
    
    # P [ I(M=d.opt)=1 | R=0, C=c ] # 
    pi_I1 <- plogis(coef(fit.I)[1])
    # P [ I(M=d.opt)=0 | R=0, C=c ] # 
    pi_I0 <- 1 - pi_I1
    
    ### Step 2 ###
    W_IIE.1 <- DATA1$Ind * DATA1$w
    W_IIE.0 <- (!DATA1$Ind) * DATA1$w
    
    # (theta1) = (1) #
    idx_R1 <- which(DATA1$R == 1)
    fit.lm1 <- lm(Y ~ C, data=DATA_R1, weights = W_IIE.1[idx_R1])
    term1 <- pi_I1 * coef(fit.lm1)[1]
    # (theta1) = (0) #
    fit.lm0 <- lm(Y ~ C, data=DATA_R1, weights = W_IIE.0[idx_R1])
    term0 <- pi_I0 * coef(fit.lm0)[1]
    
    W_iie.Y <- term1 + term0 
    
    estimated_delta_IIE <- W_Y1 - W_iie.Y 
    estimated_zeta_IIE <- W_iie.Y - W_Y0
    
    #regression
    fit.lm <- lm(Y ~ R + C + Ind , data=DATA1, weights=w)
    fit.Ind <- glm(Ind ~ R * C, data = DATA1, family = binomial(logit))
    
    alpha1 <- plogis(sum(coef(fit.Ind)[1:2])) - plogis(coef(fit.Ind)[1])
    
    
    reg_delta_IIE <- (alpha1*coef(fit.lm)["IndTRUE"])
    reg_zeta_IIE <- (W_Y1 - W_Y0) - reg_delta_IIE
    
    return(list(accuracy=accuracy, estimated_zeta_ICDE= estimated_zeta_ICDE, estimated_delta_IIE=estimated_delta_IIE,
                estimated_zeta_IIE=estimated_zeta_IIE, reg_zeta_IIE= reg_zeta_IIE, reg_delta_IIE=reg_delta_IIE   ))
  }
  
  est.ori = est1(data)
  
  accuracy = est.ori$accuracy
  reg_delta_IIE = est.ori$reg_delta_IIE
  reg_zeta_IIE = est.ori$reg_zeta_IIE
  estimated_delta_IIE = est.ori$estimated_delta_IIE
  estimated_zeta_IIE = est.ori$estimated_zeta_IIE
  estimated_zeta_ICDE = est.ori$estimated_zeta_ICDE
  
  
  regzeta_IIE.boot <- regdelta_IIE.boot <- zeta_ICDE.boot <-zeta_IIE.boot <-delta_IIE.boot <- rep(NA, B)
  
  for(i in 1:B){
    m <- nrow(data)
    data.boot = data[sample(1:nrow(data), m, replace = T), ]
    est.boot = est1(data.boot)
    
    delta_IIE.boot[i] = est.boot$estimated_delta_IIE
    zeta_IIE.boot[i] = est.boot$estimated_zeta_IIE  
    zeta_ICDE.boot[i] = est.boot$estimated_zeta_ICDE
    regzeta_IIE.boot[i] = est.boot$reg_zeta_IIE
    regdelta_IIE.boot[i] = est.boot$reg_delta_IIE
  }
  
  SE_delta_IIE = sd(unlist(delta_IIE.boot), na.rm=TRUE)
  SE_zeta_IIE = sd(unlist(zeta_IIE.boot), na.rm=TRUE)
  SE_zeta_ICDE = sd(unlist(zeta_ICDE.boot), na.rm=TRUE)
  SE_regdelta_IIE = sd(unlist(regdelta_IIE.boot), na.rm=TRUE)
  SE_regzeta_IIE = sd(unlist(regzeta_IIE.boot), na.rm=TRUE)
  
  
  results = c(accuracy=accuracy, estimated_delta_IIE = as.numeric(estimated_delta_IIE), SE_delta_IIE = SE_delta_IIE, 
              estimated_zeta_IIE = as.numeric(estimated_zeta_IIE), SE_zeta_IIE = SE_zeta_IIE,
              estimated_zeta_ICDE= as.numeric(estimated_zeta_ICDE), SE_zeta_ICDE=SE_zeta_ICDE,
              reg_delta_IIE = as.numeric(reg_delta_IIE), SE_regdelta_IIE = SE_regdelta_IIE, 
              reg_zeta_IIE = as.numeric(reg_zeta_IIE), SE_regzeta_IIE = SE_regzeta_IIE
  )
  
  return(results)
}
########################################################################################################################################################
##############################################                                                        ##################################################
##############################################            Define genU function to generate U          ##################################################
##############################################                     (with E-M steps)                   ##################################################
##############################################                                                        ##################################################
########################################################################################################################################################

# Generation of U based on the stochastic EM algorithm (When U is dependent of RXC)

genU1 = function(Y, R, M, X1, X2, X3, C, b.y, b.m, p.u, Iternum, data) {


  p.u=0.5 #Initial value of P(U)

  data$U = rbinom(nrow(data), 1, p.u) #start from U~binomial(0.5)
  coef.y.updated = NULL
  coef.m.updated = NULL

  # The EM steps
  for (iter in 1:Iternum) {
    # Outcome model
    l.y = lm(as.formula(paste("Y", "~", "R","+", "X1", "*", "M", "*","X2", "+", paste("C", "X3", collapse = "+", sep = "+"))),
             offset = b.y * data$U, data = data)
    coef.y = c(l.y$coef, U = b.y)
    sd.y = sigma(l.y) # This is the sd of Y conditional on t, m, C, X, and U.

    # Mediator model
    l.m = glm(as.formula(paste("M", "~", "R", "+", paste("C", "X1", "X2", "X3", collapse = "+", sep = "+"))),
              offset = b.m * data$U, data = data, family=binomial(link="logit"))
    coef.m = c(l.m$coef, U = b.m)

    # Conditional probability of Y
    mean.yu1 = cbind(model.matrix(l.y), U = 1) %*% coef.y
    sd.yu1 = sd.y
    mean.yu0 = cbind(model.matrix(l.y), U = 0) %*% coef.y
    sd.yu0 = sd.y
    pyu1 = dnorm(data[, "Y"], mean = mean.yu1, sd = sd.yu1)
    pyu0 = dnorm(data[, "Y"], mean = mean.yu0, sd = sd.yu0)

    # Conditional probability of M
    pm1u1 = 1/(1 + exp(-(cbind(model.matrix(l.m), U = 1) %*% coef.m)))
    pm1u0 = 1/(1 + exp(-(cbind(model.matrix(l.m), U = 0) %*% coef.m)))
    pmu1 = pm1u1^data[, "M"] * (1 - pm1u1)^(1 - data[, "M"])
    pmu0 = pm1u0^data[, "M"] * (1 - pm1u0)^(1 - data[, "M"])

    # New conditional probability of U
    p = pyu1 * pmu1 * p.u/(pyu1 * pmu1 * p.u + pyu0 * pmu0 * (1 - p.u))

    data$U <- rbinom(nrow(data), 1, p)

  }


  # Check the convergence of the coefficients
  coef.y.updated = rbind(coef.y.updated, coef.y)
  coef.m.updated = rbind(coef.m.updated, coef.m)

  return(list(new_U = data$U, coef.y.updated = coef.y.updated, coef.m.updated = coef.m.updated))

} #end of genU1 function


########################################################################################################################################################
########################################################################################################################################################




