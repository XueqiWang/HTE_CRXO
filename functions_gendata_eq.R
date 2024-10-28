############################################################################################
################################### Data generations #######################################
############################################################################################

# Function to generate the data with the cross-sectional design
gendata_CS <- function(beta2, beta4, rho0, rho1, alpha0, alpha1, sigma2x=1, sigma2y=1, I, J, K){
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
  # K: number of individuals in each period
  # 
  # Output:
  # data: simulated data frame
  
  nca <- I/2
  K_id <- seq(1, (I*J*K), 1)
  T_id <- rep(1:(I*J), each=K)
  I_id <- rep(1:I, each=K*J)
  time <- rep(rep(1:J, each=K), I)
  data <- data.frame(cbind(K_id, T_id, I_id, time))
  
  # Generate beta1j and beta3j
  beta1j <- c(0.5, 0.6)
  beta3j <- c(0.2, 0.4)
  data$beta1 <- rep(rep(beta1j, each=K), I)
  data$beta3 <- rep(rep(beta3j, each=K), I)
  
  # Generate the treatment assignment indicators based on design parameters
  trtSeq <- matrix(c(1,0,0,1), ncol=2)
  W_assignment <- NULL
  for (i in 1:2){
    W_assignment <- c(W_assignment, rep(trtSeq[i,], nca))
  }
  data$W <- rep(W_assignment, each=K)
  
  # Generate X (one continuous individual-level covariate)
  Xmean <- 1
  c <- rnorm(I*J*K, 0, sqrt(sigma2x*(1-rho0)))
  b <- rep(rnorm(I*J, 0, sqrt(sigma2x*(rho0-rho1))), each=K)
  a <- rep(rnorm(I, 0, sqrt(sigma2x*rho1)), each=K*J)
  data$X <- Xmean+a+b+c
  
  # Generate Y under null and alternative
  epsilon <- rnorm(I*J*K, 0, sqrt(sigma2y*(1-alpha0)))
  u <- rep(rnorm(I*J, 0, sqrt(sigma2y*(alpha0-alpha1))), each=K)
  gamma <- rep(rnorm(I, 0, sqrt(sigma2y*alpha1)), each=K*J)
  
  data$Y <- data$beta1 + beta2*data$W + data$beta3*data$X + beta4*data$W*data$X + gamma+u+epsilon
  
  return(data)
}

# Function to generate the data with the closed-cohort design
gendata_CC <- function(beta2, beta4, rho0, alpha0, alpha1, alpha2, sigma2x=1, sigma2y=1, I, J, K){
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
  K_id <- seq(1, (I*J*K), 1)
  T_id <- rep(1:(I*J), each=K)
  I_id <- rep(1:I, each=K*J)
  time <- rep(rep(1:J, each=K), I)
  
  # To fit random effect v, create v_id
  ik_matrix <- matrix(1:(K*I), nrow=I, ncol=K, byrow=T)
  v_id_matrix <- NULL
  for (i in 1:J){
    v_id_matrix <- cbind(v_id_matrix, ik_matrix)
  }
  v_id <- NULL
  for (i in 1:I){
    v_id <- c(v_id, v_id_matrix[i,])
  }
  
  data <- data.frame(cbind(K_id, T_id, I_id, time, v_id))
  
  # Generate beta1j and beta3j
  beta1j <- c(0.5, 0.6)
  beta3j <- c(0.2, 0.4)
  data$beta1 <- rep(rep(beta1j, each=K), I)
  data$beta3 <- rep(rep(beta3j, each=K), I)
  
  # Generate the treatment assignment indicators based on design parameters
  trtSeq <- matrix(c(1,0,0,1), ncol=2)
  W_assignment <- NULL
  for (i in 1:2){
    W_assignment <- c(W_assignment, rep(trtSeq[i,], nca))
  }
  data$W <- rep(W_assignment, each=K)
  
  # Generate X (one continuous individual-level covariate)
  Xmean <- 1
  
  # Create a_i, repeat for periods
  a <- rep(rnorm(I, 0, sqrt(sigma2x*rho0)), each=K)
  c <- rnorm(I*K, 0, sqrt(sigma2x*(1-rho0)))
  X_oneperiod <- Xmean+a+c
  X_matrix <- matrix(X_oneperiod, ncol=K, nrow=I, byrow=T)
  
  X_vector <- NULL
  for (i in 1:I){
    X_vector <- c(X_vector, rep(X_matrix[i,], J))
  }
  data$X <- X_vector
  
  # Generate Y under null and alternative
  # Create v_ik
  v_matrix <- matrix(rnorm(I*K, 0, sqrt(sigma2y*(alpha2-alpha1))),
                     nrow=I, ncol=K, byrow=T)
  v_big_matrix <- NULL
  for (i in 1:J){
    v_big_matrix <- cbind(v_big_matrix, v_matrix)
  }
  v <- NULL
  for (i in 1:I){
    v <- c(v, v_big_matrix[i,])
  }
  
  epsilon <- rnorm(I*J*K, 0, sqrt(sigma2y*(1-alpha0+alpha1-alpha2)))
  u <- rep(rnorm(I*J, 0, sqrt(sigma2y*(alpha0-alpha1))), each=K)
  gamma <- rep(rnorm(I, 0, sqrt(sigma2y*alpha1)), each=K*J)
  
  data$Y <- data$beta1 + beta2*data$W + data$beta3*data$X + beta4*data$W*data$X + gamma+u+epsilon+v
  
  return(data)
}


