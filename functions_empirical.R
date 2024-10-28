############################################################################################
###################### Empirical type I error or empirical power ###########################
############################################################################################

# Function to compute empirical power or empirical type I error under cross-sectional design, HTE test 
empirical_CS_HTE <- function(nullcase=F, parameter, nsims=2000){
  # Argument:
  # nullcase: calculate empirical type I error or power. The former if nullcase=TRUE
  # parameter: the combination of design parameters of interest 
  # nsims: number of data replications
  
  I <- as.numeric(parameter[9])
  J <- as.numeric(parameter[2])
  K <- as.numeric(parameter[1])
  CV <- as.numeric(parameter[3])
  alpha0 <- as.numeric(parameter[4])
  alpha1 <- as.numeric(parameter[5])
  rho0 <- as.numeric(parameter[6])
  rho1 <- as.numeric(parameter[7])
  
  beta2 <- 0.6
  beta4 <- delta_HTE
  
  if (nullcase==T){
    beta4 <- 0
  }
  
  b.est <- b.var <- pvalue <- NULL
  count <- NULL
  for (i in 1:nsims){
    set.seed(i)
    simdata <- gendata_CS(beta2=beta2, beta4=beta4, rho0=rho0, rho1=rho1, alpha0=alpha0, alpha1=alpha1, I=I, J=J, K=K, CV=CV)
    simdata$time <- factor(simdata$time)
    
    fit <- try(lme(Y ~ -1 + time + W + time:X + W:X, data=simdata, random=list(I_id = ~ 1, T_id = ~ 1)), silent=T)
    if(class(fit)=="try-error"){next}
    count[i] <- 1
    
    R <- c(rep(0,2*J+1),1)
    beta <- fit$coef$fixed
    b.est[i] <- as.numeric(t(R)%*%beta)
    b.var[i] <- as.numeric(t(R)%*%vcov(fit)%*%R)
    test.stat <- as.numeric((t(R)%*%beta)^2/(t(R)%*%vcov(fit)%*%R))
    pvalue[i] <- 2*(1-pt(sqrt(test.stat), df=I-2))
  }
  empirical <- mean(pvalue<0.05, na.rm=T)
  error.rate <- 1-sum(count, na.rm=T)/nsims
  
  data.results <- cbind(b.est, b.var, pvalue, empirical, error.rate)
  
  #return(c(empirical, error.rate))
  return(data.results)
}


# Function to compute empirical power or empirical type I error under closed-cohort design, HTE test 
empirical_CC_HTE <- function(nullcase=F, parameter, nsims=2000){
  # Argument:
  # nullcase: calculate empirical type I error or power. The former if nullcase=TRUE
  # parameter: the combination of design parameters of interest 
  # nsims: number of data replications
  
  I <- as.numeric(parameter[9])
  J <- as.numeric(parameter[2])
  K <- as.numeric(parameter[1])
  CV <- as.numeric(parameter[3])
  alpha0 <- as.numeric(parameter[4])
  alpha1 <- as.numeric(parameter[5])
  alpha2 <- as.numeric(parameter[6])
  rho0 <- as.numeric(parameter[7])
  
  beta2 <- 0.6
  beta4 <- delta_HTE
  
  if (nullcase==T){
    beta4 <- 0
  }
  
  b.est <- b.var <- pvalue <- NULL
  count <- NULL
  for (i in 1:nsims){
    set.seed(i)
    simdata <- gendata_CC(beta2=beta2, beta4=beta4, rho0=rho0, alpha0=alpha0, alpha1=alpha1, alpha2=alpha2, I=I, J=J, K=K, CV=CV)
    simdata$time <- factor(simdata$time)
    
    fit <- try(lmer(Y ~ -1 + time + W + time:X + W:X + (1|I_id) + (1|I_id:T_id) + (1|I_id:v_id), data=simdata), silent=T)
    if(class(fit)=="try-error"){next}
    count[i] <- 1
    
    R <- c(rep(0,2*J+1),1)
    beta <- fixef(fit)
    b.est[i] <- as.numeric(t(R)%*%beta)
    b.var[i] <- as.numeric(t(R)%*%vcov(fit)%*%R)
    test.stat <- as.numeric((t(R)%*%beta)^2/(t(R)%*%vcov(fit)%*%R))
    pvalue[i] <- 2*(1-pt(sqrt(test.stat), df=I-2))
  }
  empirical <- mean(pvalue<0.05, na.rm=T)
  error.rate <- 1-sum(count, na.rm=T)/nsims
  
  data.results <- cbind(b.est, b.var, pvalue, empirical, error.rate)
  
  #return(c(empirical, error.rate))
  return(data.results)
}


