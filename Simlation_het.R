rm(list = ls())
ns <- c(500, 1000, 2000)

for(nn in 1:3){
  
  rm(list = setdiff(ls(), c("ns", "nn")))
  
  nsample <- ns[nn]
  
  library(DynTxRegime)
  library(rpart)
  library(dplyr)
  library(parallel)
  
  source("functions_het.R")
  
  ################################################################################################
  ############################------ Pre-set these values ------##################################
  ################################################################################################
  
  Iternum <- 15 #100 #Number of iterations in E-M algorithm
  
  B<- 300 #100 #number of bootstrapping for each "nsim" 
  
  nsim<- 30 #50 # To account for the uncertainty of U and to obtain more accurate estimates(following Qin and Yang ; same as "K" in Qin's paper)
  
  REP<- 500    #1000 #Number of repetitions, at the end we get "REP" number of estimates.
  
  num_cores <- 48 #100  # Number of cores to use. Change if necessary.
  
  
  # Set important parameters
  b.m <- b.y <-1.5 #vary to 1 or 0.5
  beta_m <- 1
  # Generate population data
  set.seed(7124)
  C <-rbinom(1000000,1,0.4)
  U <- rbinom(1000000,1,0.5)
  pR1_C <-exp(1-0.5*C)/(1+exp(1-0.5*C)) #inverse logit function to get P(R=1|C).
  R <- rbinom(1000000,1,pR1_C) #Create population level of R using P(R=1|C) as defined above.
  X1 <- -.8+1*R+1.5*C+rnorm(1000000,0,1) 
  X2 <- .5+.5*R+.5*C+rnorm(1000000,0,1) 
  X3 <- -.8-1*R+.5*C+rnorm(1000000,0,1) 
  pM1_RXCU <- 1/(1+exp(-(.5-0.5*R+0.2*X1- .5*X2 -.2*X3 +.5*C+b.m*U))) #inverse logit function to get P(M=1|R,X,C,U). 
  M <- rbinom(1000000,1,pM1_RXCU) #Create population level of M using P(M=1|R,X,C,U) as defined above.
  Mopt <- (X1 > .1)*(U > 0.5)
  Y <- 0.5 + 0.25*X1 + 0.25*X2 -0.25*X3 -beta_m*(M-Mopt)^2 + 0.25*C -0.5*R +b.y*U + rnorm(1000000,0,1) 
  true.Ind =(M==Mopt)
  popdata <- data.frame(Y,R,M,X1, X2, X3,C,U, true.Ind, Mopt)
  
  a <- lm(Y~ R+ X1 +M + C+ X2 + X3 + M*(X1+ U), data=popdata )
  eb.y <- a$coef["U"]
  eb.yu <- a$coef["M:U"]
  
  #:::::::::::::::::::::::::#
  ### Calculate True ICDE ###
  #:::::::::::::::::::::::::#
  EY_R1 <- 0.5 + 0.25*mean((X1 + X2 -X3)[R==1 & C==0]) -0.5*1 +b.y*mean(U[R==1 & C==0])
  EY_R0 <- 0.5 + 0.25*mean((X1 + X2 -X3)[R==0 & C==0]) -0.5*0 +b.y*mean(U[R==0 & C==0])
  zeta_ICDE <- EY_R1 - EY_R0
  
  #:::::::::::::::::::::::::#
  ### Calculate True IIE ###
  #:::::::::::::::::::::::::#
  
  EY_R1C <- popdata %>% filter(R == 1, C == 0) %>% summarise(EY_R1c=mean(Y))
  EY_R0C <- popdata %>% filter(R == 0, C == 0) %>% summarise(EY_R0c=mean(Y)) 
  
  wm.Y1 <- popdata %>% filter(R == 1, C == 0) %>% summarise(wm.Y11=mean(0.5 + 0.25*(X1 + X2 - X3) - 0.5 + b.y*U))
  I1 <- popdata %>% filter(R == 0 , C == 0) %>% summarise(mean(true.Ind) )
  wm.Y0 <- popdata %>% filter(R == 1, C == 0) %>% summarise(wm.Y10=mean(0.5 + 0.25*(X1 + X2 - X3) - beta_m - 0.5 + b.y*U))
  I0 <- popdata %>% filter(R == 0 , C == 0) %>% summarise(mean(!true.Ind))
  
  EY_K1K2 <- wm.Y1*I1 + wm.Y0*I0 
  delta_IIE <- as.numeric(EY_R1C - EY_K1K2)
  zeta_IIE <- as.numeric(EY_K1K2 - EY_R0C)
  
  ################################################################################################
  ################################################################################################
  
  #Below is a modified code to use 'mclapply' function. Outer-loop is parallelized.
  
  
  ########################################################################################################################################################
  ##############################################################                                 #########################################################
  ##############################################################  Not addressing uncertainty   #########################################################
  ##############################################################                                 #########################################################
  ########################################################################################################################################################
  
  
  # Define the core functionality inside the for-loop as a function
  simulate_one_ind <- function(j, popdata, nsample, nsim, Y, R, M, X1, X2, X3, C, b.y, b.m, p.u, Iternum, B) {
    # Randomly select 'nsample' amount of data from population.
    random_rows <- sample(nrow(popdata), nsample, replace = FALSE)
    data <- popdata[random_rows, ] 
    results_ind  = NULL
    results_unadj  = NULL
    
    for(k in 1:nsim) {
      data$U = genU1(Y=Y, R=R, M=M, X1=X1, X2=X2, X3=X3, C=C, eb.y=eb.y, eb.yu=eb.yu, b.m=b.m, p.u=p.u, Iternum = Iternum, data= data)$new_U
      results_ind = rbind(results_ind, reg1(Y=Y, R=R, M=M, X1=X1, X2=X2, X3=X3, C=C, b.y=b.y, b.m=b.m, B=B, data = data))
      results_unadj = rbind(results_unadj, reg2(Y=Y, R=R, M=M, X1=X1, X2=X2, X3=X3, C=C, b.y=b.y, b.m=b.m, B=B, data = data))
    }
    
    list(accuracy = mean(results_ind[,"accuracy"]),
         estimated_zeta_ICDE = mean(results_ind[, "estimated_zeta_ICDE"]),
         estimated_delta_IIE = mean(results_ind[, "estimated_delta_IIE"]),
         estimated_zeta_IIE = mean(results_ind[, "estimated_zeta_IIE"]),
         reg_delta_IIE =mean(results_ind[, "reg_delta_IIE"]),
         reg_zeta_IIE=mean(results_ind[, "reg_zeta_IIE"]),
         SE_zeta_ICDE = sqrt(mean(results_ind[, "SE_zeta_ICDE"]^2)+(1 + 1/nsim) * var(results_ind[, "estimated_zeta_ICDE"])),
         SE_delta_IIE = sqrt(mean(results_ind[, "SE_delta_IIE"]^2)+(1 + 1/nsim) * var(results_ind[, "estimated_delta_IIE"])),
         SE_zeta_IIE = sqrt(mean(results_ind[, "SE_zeta_IIE"]^2)+(1 + 1/nsim) * var(results_ind[, "estimated_zeta_IIE"])),
         SE_regdelta_IIE = sqrt(mean(results_ind[, "SE_regdelta_IIE"]^2)+(1 + 1/nsim) * var(results_ind[, "reg_delta_IIE"])),
         SE_regzeta_IIE = sqrt(mean(results_ind[, "SE_regzeta_IIE"]^2)+(1 + 1/nsim) * var(results_ind[, "reg_zeta_IIE"])),
         unadj_accuracy =mean(results_unadj[,"accuracy"]),
         unadj_zeta_ICDE = mean(results_unadj[, "estimated_zeta_ICDE"]),
         unadj_delta_IIE = mean(results_unadj[, "estimated_delta_IIE"]),
         unadj_zeta_IIE = mean(results_unadj[, "estimated_zeta_IIE"]),
         unadj_reg_delta_IIE =mean(results_unadj[, "reg_delta_IIE"]),
         unadj_reg_zeta_IIE=mean(results_unadj[, "reg_zeta_IIE"]),
         unadj_SE_zeta_ICDE = sqrt(mean(results_unadj[, "SE_zeta_ICDE"]^2)+(1 + 1/nsim) * var(results_unadj[, "estimated_zeta_ICDE"])),
         unadj_SE_delta_IIE = sqrt(mean(results_unadj[, "SE_delta_IIE"]^2)+(1 + 1/nsim) * var(results_unadj[, "estimated_delta_IIE"])),
         unadj_SE_zeta_IIE = sqrt(mean(results_unadj[, "SE_zeta_IIE"]^2)+(1 + 1/nsim) * var(results_unadj[, "estimated_zeta_IIE"])),
         unadj_SE_regdelta_IIE = sqrt(mean(results_unadj[, "SE_regdelta_IIE"]^2)+(1 + 1/nsim) * var(results_unadj[, "reg_delta_IIE"])),
         unadj_SE_regzeta_IIE = sqrt(mean(results_unadj[, "SE_regzeta_IIE"]^2)+(1 + 1/nsim) * var(results_unadj[, "reg_zeta_IIE"]))
    )
  }
  
  # Now, use mclapply to run the above function in parallel
  
  system.time(results <- mclapply(1:REP, FUN = function(i) {
    res <- try(simulate_one_ind(i, popdata, nsample, nsim, Y, R, M, X1, X2, X3, C, b.y, b.m, p.u, Iternum, B), silent=TRUE)
    if (inherits(res, "try-error")) return(NULL) else return(res)
  }, mc.cores = num_cores))
  
  
  # Extract the results_ind(Need to check if this codes works)
  saved_accuracy <- sapply(results, "[[", "accuracy")
  saved_estimated_zeta_ICDE <- sapply(results, "[[", "estimated_zeta_ICDE")
  saved_estimated_delta_IIE <- sapply(results, "[[", "estimated_delta_IIE")
  saved_estimated_zeta_IIE <- sapply(results, "[[", "estimated_zeta_IIE")
  saved_reg_delta_IIE <- sapply(results, "[[", "reg_delta_IIE")
  saved_reg_zeta_IIE <- sapply(results, "[[", "reg_zeta_IIE")
  SE_zeta_ICDE.all <- sapply(results, "[[", "SE_zeta_ICDE")
  SE_delta_IIE.all <- sapply(results, "[[", "SE_delta_IIE")
  SE_zeta_IIE.all <- sapply(results, "[[", "SE_zeta_IIE")
  SE_regdelta_IIE.all <- sapply(results, "[[", "SE_regdelta_IIE")
  SE_regzeta_IIE.all <- sapply(results, "[[", "SE_regzeta_IIE")
  
  unadj_accuracy <- sapply(results, "[[", "unadj_accuracy")
  unadj_zeta_ICDE <- sapply(results, "[[", "unadj_zeta_ICDE")
  unadj_delta_IIE <- sapply(results, "[[", "unadj_delta_IIE")
  unadj_zeta_IIE <- sapply(results, "[[", "unadj_zeta_IIE")
  unadj_reg_delta_IIE <- sapply(results, "[[", "unadj_reg_delta_IIE")
  unadj_reg_zeta_IIE <- sapply(results, "[[", "unadj_reg_zeta_IIE")
  unadj_SE_zeta_ICDE.all <- sapply(results, "[[", "unadj_SE_zeta_ICDE")
  unadj_SE_delta_IIE.all <- sapply(results, "[[", "unadj_SE_delta_IIE")
  unadj_SE_zeta_IIE.all <- sapply(results, "[[", "unadj_SE_zeta_IIE")
  unadj_SE_regdelta_IIE.all <- sapply(results, "[[", "unadj_SE_regdelta_IIE")
  unadj_SE_regzeta_IIE.all <- sapply(results, "[[", "unadj_SE_regzeta_IIE")
  
  #end of simulation codes
  rm(popdata)
  save.image(paste0("indivdualized_het", nsample, ".1.5.RData"))
  
}





