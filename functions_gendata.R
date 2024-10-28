############################################################################################
################################### Data generations #######################################
############################################################################################

# Function to generate the data with the cross-sectional design
gendata_CS <- function(beta2, beta4, rho0, rho1, alpha0, alpha1, sigma2x=1, sigma2y=1, I, J, K, CV){
  # Argument:
  # beta2: true parameter of the main effect of treatment
  # beta4: true parameter of the interaction effect of treatment and the univariate covariate
  # rho0: within-period covariate ICC
  # rho1: between-period covariate ICC
  # alpha0: within-period outcome ICC
  # alpha1: between-period outcome ICC
  # sigma2x: total variance of the single covariate X
  # sigma2y: total variance of outcome Y
  # I: number of clusters 
  # J: number of periods in each cluster
  # K: mean cluster-period size
  # CV: coefficient of variance of cluster-period size
  # 
  # Output:
  # data: simulated data frame
  
  nca <- I/2
  if (CV == 0){
    cpsize <- rep(K, I*J)
  } else{
    gamma_a <- CV^(-2)
    gamma_b <- 1/(K*CV^2)
    cpsize_h <- round(rgamma(I,shape=gamma_a,rate=gamma_b))
    cpsize_h[cpsize_h<=2] <- 2
    cpsize <- rep(cpsize_h, each=2)
  }
  
  Ksum <- sum(cpsize)
  
  K_id <- seq(1, Ksum, 1)
  T_id <- rep(1:(I*J), times=cpsize)
  I_id <- rep(rep(1:I, each=J), times=cpsize)
  time <- rep(rep(1:J, times=I), times=cpsize)
  data <- data.frame(cbind(K_id, T_id, I_id, time))
  
  # Generate beta1j and beta3j
  beta1j <- c(0.5, 0.6)
  beta3j <- c(0.2, 0.4)
  data$beta1 <- rep(rep(beta1j, times=I), times=cpsize)
  data$beta3 <- rep(rep(beta3j, times=I), times=cpsize)
  
  # Generate the treatment assignment indicators based on design parameters
  trtSeq <- matrix(c(1,0,0,1), ncol=2)
  W_assignment <- NULL
  for (i in 1:2){
    W_assignment <- c(W_assignment, rep(trtSeq[i,], nca))
  }
  data$W <- rep(W_assignment, times=cpsize)
  
  # Generate X (one continuous individual-level covariate)
  Xmean <- 1
  c <- rnorm(Ksum, 0, sqrt(sigma2x*(1-rho0)))
  b <- rep(rnorm(I*J, 0, sqrt(sigma2x*(rho0-rho1))), times=cpsize)
  a <- rep(rep(rnorm(I, 0, sqrt(sigma2x*rho1)), each=J), times=cpsize)
  data$X <- Xmean+a+b+c
  
  # Generate Y under null and alternative
  epsilon <- rnorm(Ksum, 0, sqrt(sigma2y*(1-alpha0)))
  u <- rep(rnorm(I*J, 0, sqrt(sigma2y*(alpha0-alpha1))), times=cpsize)
  gamma <- rep(rep(rnorm(I, 0, sqrt(sigma2y*alpha1)), each=J), times=cpsize)
  
  data$Y <- data$beta1 + beta2*data$W + data$beta3*data$X + beta4*data$W*data$X + gamma+u+epsilon
  
  return(data)
}

# Function to generate the data with the closed-cohort design
gendata_CC <- function(beta2, beta4, rho0, alpha0, alpha1, alpha2, sigma2x=1, sigma2y=1, I, J, K, CV){
  # Argument:
  # beta2: true parameter of the main effect of treatment
  # beta4: true parameter of the interaction effect of treatment and the univariate covariate
  # rho0: within-period covariate ICC
  # alpha0: within-period outcome ICC
  # alpha1: between-period outcome ICC
  # alpha2: within-individual outcome ICC
  # sigma2x: total variance of the single covariate X
  # sigma2y: total variance of outcome Y
  # I: number of clusters 
  # J: number of periods in each cluster
  # K: number of individuals in each period
  # 
  # Output:
  # data: simulated data frame
  
  nca <- I/2
  if (CV == 0){
    cpsize_h <- rep(K, I)
  } else{
    gamma_a <- CV^(-2)
    gamma_b <- 1/(K*CV^2)
    cpsize_h <- round(rgamma(I,shape=gamma_a,rate=gamma_b))
    cpsize_h[cpsize_h<=2] <- 2
  }
  cpsize <- rep(cpsize_h, each=2)
  
  Ksum <- sum(cpsize)
  Ksum_h <- sum(cpsize_h)
  
  K_id <- seq(1, Ksum, 1)
  T_id <- rep(1:(I*J), times=cpsize)
  I_id <- rep(rep(1:I, each=J), times=cpsize)
  time <- rep(rep(1:J, times=I), times=cpsize)

  # To fit random effect v, create v_id
  ik_num <- c(0, cumsum(cpsize_h))
  ik_vector <- seq(1, Ksum/2, 1)
  v_id <- NULL
  for (i in 1:length(cpsize_h)){
    v_id_i <- ik_vector[(ik_num[i]+1):ik_num[i+1]]
    v_id <- c(v_id, rep(v_id_i, time=2))
  }
  
  data <- data.frame(cbind(K_id, T_id, I_id, time, v_id))
  
  # Generate beta1j and beta3j
  beta1j <- c(0.5, 0.6)
  beta3j <- c(0.2, 0.4)
  data$beta1 <- rep(rep(beta1j, times=I), times=cpsize)
  data$beta3 <- rep(rep(beta3j, times=I), times=cpsize)
  
  # Generate the treatment assignment indicators based on design parameters
  trtSeq <- matrix(c(1,0,0,1), ncol=2)
  W_assignment <- NULL
  for (i in 1:2){
    W_assignment <- c(W_assignment, rep(trtSeq[i,], nca))
  }
  data$W <- rep(W_assignment, times=cpsize)
  
  # Generate X (one continuous individual-level covariate)
  Xmean <- 1
  
  # Create a_i, repeat for periods
  a <- rep(rnorm(I, 0, sqrt(sigma2x*rho0)), times=cpsize_h)
  c <- rnorm(Ksum/2, 0, sqrt(sigma2x*(1-rho0)))
  X_oneperiod <- Xmean+a+c

  X_vector <- NULL
  for (i in 1:length(cpsize_h)){
    X_vector_i <- X_oneperiod[(ik_num[i]+1):ik_num[i+1]]
    X_vector <- c(X_vector, rep(X_vector_i, time=2))
  }
  data$X <- X_vector
  
  # Generate Y under null and alternative
  # Create v_ik
  v_oneperiod <- rnorm(Ksum/2, 0, sqrt(sigma2y*(alpha2-alpha1)))
  v <- NULL
  for (i in 1:length(cpsize_h)){
    v_i <- v_oneperiod[(ik_num[i]+1):ik_num[i+1]]
    v <- c(v, rep(v_i, time=2))
  }
  
  epsilon <- rnorm(Ksum, 0, sqrt(sigma2y*(1-alpha0+alpha1-alpha2)))
  u <- rep(rnorm(I*J, 0, sqrt(sigma2y*(alpha0-alpha1))), times=cpsize)
  gamma <- rep(rep(rnorm(I, 0, sqrt(sigma2y*alpha1)), each=J), times=cpsize)
  
  data$Y <- data$beta1 + beta2*data$W + data$beta3*data$X + beta4*data$W*data$X + gamma+u+epsilon+v
  
  return(data)
}


