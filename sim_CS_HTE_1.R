############################################################################################
##################### Simulations for the cross-sectional design ###########################
############################################################################################

library(nlme)

source("functions_calc_ss.R")
source("functions_gendata.R")
source("functions_empirical.R")

delta_HTE <- 0.08

K <- c(rep(50, 24), rep(100, 24))
J <- 2
CV <- rep(c(rep(0, 6), rep(0.25, 6), rep(0.5, 6), rep(0.75, 6)), 2)
alpha0 <- rep(c(rep(0.015, 3), rep(0.1, 3)), 8)
alpha1 <- rep(c(rep(0.01, 3), rep(0.05, 3)), 8)
rho0 <- rep(c(0.15, 0.3, 0.5), 16)
rho1 <- rep(c(0.1, 0.15, 0.3), 16)
table <- cbind(K, J, CV, alpha0, alpha1, rho0, rho1)

I <- numeric(48)
nclusters_per_arm <- numeric(48)
pred.power <- numeric(48)
analytical_var4 <- numeric(48)

for (i in 1:nrow(table)){
  K.input <- as.numeric(table[i,][1])
  J.input <- as.numeric(table[i,][2])
  CV.input <- as.numeric(table[i,][3])
  alpha0.input <- as.numeric(table[i,][4])
  alpha1.input <- as.numeric(table[i,][5])
  rho0.input <- as.numeric(table[i,][6])
  rho1.input <- as.numeric(table[i,][7])
  
  res_i <- calc_CS_HTE(eff=delta_HTE, K=K.input, CV=CV.input, alpha0=alpha0.input, alpha1=alpha1.input, rho0=rho0.input, rho1=rho1.input)
  nclusters_per_arm[i] <- res_i[1]
  I[i] <- res_i[2]
  pred.power[i] <- res_i[3]
  analytical_var4[i] <- res_i[4]
}

table <- cbind(table, nclusters_per_arm, I, analytical_var4)
setwd('conResults/S1')

empirical.tIe <- NULL
for (i in 1:nrow(table)){
  empirical.tIe.data <- empirical_CS_HTE(nullcase=T, parameter=table[i,])
  empirical.tIe <- rbind(empirical.tIe, empirical.tIe.data[1, 4:5])
  save(empirical.tIe.data,file=paste0("CS_",i,"_size_beta.RData"))
}

empirical.power <- NULL
for (i in 1:nrow(table)){
  empirical.power.data <- empirical_CS_HTE(parameter=table[i,])
  empirical.power <- rbind(empirical.power, empirical.power.data[1, 4:5])
  save(empirical.power.data,file=paste0("CS_",i,"_power_beta.RData"))
}

result <- data.frame(cbind(table, empirical.tIe, empirical.power, pred.power))
save(result, file="CS_HTE.RData")


