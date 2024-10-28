############################################################################################
#############################      Power calculations      #################################
############################################################################################

# Function to calculate predicted number of clusters to get at least (1-beta) power, 
# the actual power that can be reached, and the analytical variance 
# under the cross-sectional design, HTE test
power_CS_HTE <- function(K=50, CV, p=0.5, I=50, sigma2y=1, sigma2x=1, eff=0.08, alpha0, alpha1, rho0, rho1, alpha=0.05){
  # Argument:
  # K: mean cluster-period size
  # CV: coefficient of variance of cluster-period size
  # p: proportion of clusters randomly assigned to the AB sequence
  # I: total number of clusters
  # sigma2y: total variance of outcome Y
  # sigma2x: total variance of the single covariate X
  # eff: effect size
  # alpha0: within-period outcome ICC
  # alpha1: between-period outcome ICC
  # rho0: within-period covariate ICC
  # rho1: between-period covariate ICC
  # alpha: type I error rate
  # 
  # Output:
  # power: predicted power
  
  lambda1 <- 1-alpha0
  lambda2 <- 1+(K-1)*alpha0-K*alpha1
  lambda3 <- 1+(K-1)*alpha0+K*alpha1
  zeta1 <- 1-rho0
  zeta2 <- 1+(K-1)*rho0-K*rho1
  zeta3 <- 1+(K-1)*rho0+K*rho1
  
  if (CV == 0){
    phi <- 0
    d <- 2*(K-1)*lambda1^{-1}*zeta1 + lambda2^{-1}*zeta3 + lambda3^{-1}*zeta2
  } else{
    phi <- (2*K-(2*K-3)*zeta1-zeta2)*lambda1*lambda2^{-2} - 
      2*(K-(K-1)*zeta1)*lambda1^2*lambda2^{-3} + (zeta1+zeta2)*lambda1*lambda3^{-2} - zeta1*(lambda2^{-1}+lambda3^{-1})
    d <- 2*K*lambda2^{-1} + 2*(K-1)*zeta1*(lambda1^{-1}-lambda2^{-1}) + zeta2*(lambda3^{-1}-lambda2^{-1})
  }
  
  Var4 <- sigma2y/sigma2x / (I*p*(1-p)) / (d-CV^2*phi)
  power <- pt( sqrt(eff^2/Var4)-qt(1-alpha/2,df=I-2) , df=I-2)
  
  return(power)
}


# Function to calculate predicted number of clusters to get at least (1-beta) power, 
# the actual power that can be reached, and the analytical variance 
# under the closed-cohort design, HTE test
power_CC_HTE <- function(K=50, CV, p=0.5, I=50, sigma2y=1, sigma2x=1, eff=0.08, alpha0, alpha1, alpha2, rho0, alpha=0.05){
  # Argument:
  # K: mean cluster-period size
  # CV: coefficient of variance of cluster-period size
  # p: proportion of clusters randomly assigned to the AB sequence
  # I: total number of clusters
  # sigma2y: total variance of outcome Y
  # sigma2x: total variance of the single covariate X
  # eff: effect size
  # alpha0: within-period outcome ICC
  # alpha1: between-period outcome ICC
  # alpha2: within-individual outcome ICC
  # rho0: within-period covariate ICC
  # alpha: type I error rate
  # 
  # Output:
  # power: actual predicted power
  
  tau1 <- 1-alpha0+alpha1-alpha2
  tau2 <- 1-alpha0-(alpha1-alpha2)
  tau3 <- 1+(K-1)*(alpha0-alpha1)-alpha2
  tau4 <- 1+(K-1)*alpha0+(K-1)*alpha1+alpha2
  eta1 <- 1-rho0
  eta2 <- 1+(K-1)*rho0
  
  if (CV == 0){
    phi <- 0
    d <- 2*(K-1)*tau1^{-1}*eta1 + 2*tau3^{-1}*eta2
  } else{
    phi <- 2*((K-(K-2)*eta1)*tau1*tau3^{-2} - (K-(K-1)*eta1)*tau1^2*tau3^{-3} - eta1*tau3^{-1})
    d <- 2*(K*tau3^{-1} + (K-1)*eta1*(tau1^{-1}-tau3^{-1}))
  }
  
  Var4 <- sigma2y/sigma2x / (I*p*(1-p)) / (d-CV^2*phi)
  power <- pt( sqrt(eff^2/Var4)-qt(1-alpha/2,df=I-2) , df=I-2)
  
  return(power)
}


