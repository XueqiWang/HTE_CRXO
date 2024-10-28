############################################################################################
###################### Simulations for the closed-cohort design ############################
############################################################################################

library(lme4)

source("functions_calc_ss.R")
source("functions_gendata.R")
source("functions_empirical.R")

delta_HTE <- 0.1

K <- c(rep(50, 16), rep(100, 16))
J <- 2
CV <- rep(c(rep(0, 4), rep(0.25, 4), rep(0.5, 4), rep(0.75, 4)), 2)
alpha0 <- rep(c(rep(0.015, 2), rep(0.1, 2)), 8)
alpha1 <- rep(c(rep(0.01, 2), rep(0.05, 2)), 8)
alpha2 <- rep(c(rep(0.2, 2), rep(0.5, 2)), 8)
rho0 <- rep(c(0.2, 0.5), 16)
table <- cbind(K, J, CV, alpha0, alpha1, alpha2, rho0)

I <- numeric(32)
nclusters_per_arm <- numeric(32)
pred.power <- numeric(32)
analytical_var4 <- numeric(32)

for (i in 1:nrow(table)){
  K.input <- as.numeric(table[i,][1])
  J.input <- as.numeric(table[i,][2])
  CV.input <- as.numeric(table[i,][3])
  alpha0.input <- as.numeric(table[i,][4])
  alpha1.input <- as.numeric(table[i,][5])
  alpha2.input <- as.numeric(table[i,][6])
  rho0.input <- as.numeric(table[i,][7])
  
  res_i <- calc_CC_HTE(eff=delta_HTE, K=K.input, CV=CV.input, alpha0=alpha0.input, alpha1=alpha1.input, alpha2=alpha2.input, rho0=rho0.input)
  nclusters_per_arm[i] <- res_i[1]
  I[i] <- res_i[2]
  pred.power[i] <- res_i[3]
  analytical_var4[i] <- res_i[4]
}

table <- cbind(table, nclusters_per_arm, I, analytical_var4)
setwd('conResults/S2')

empirical.tIe <- NULL
for (i in 1:nrow(table)){
  empirical.tIe.data <- empirical_CC_HTE(nullcase=T, parameter=table[i,])
  empirical.tIe <- rbind(empirical.tIe, empirical.tIe.data[1, 4:5])
  save(empirical.tIe.data,file=paste0("CC_",i,"_size_beta.RData"))
}

empirical.power <- NULL
for (i in 1:nrow(table)){
  empirical.power.data <- empirical_CC_HTE(parameter=table[i,])
  empirical.power <- rbind(empirical.power, empirical.power.data[1, 4:5])
  save(empirical.power.data,file=paste0("CC_",i,"_power_beta.RData"))
}

result <- data.frame(cbind(table, empirical.tIe, empirical.power, pred.power))
save(result, file="CC_HTE.RData")


