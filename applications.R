############################################################################################
######################     Application to the ROC-PART study     ###########################
############################################################################################

source("functions_calc_ss.R")
source("functions_calc_power.R")

library(dplyr)
library(lme4)
library(performance)
library(ggplot2)
library(metR)
library(directlabels)
library(cowplot)

K <- 60; CV <- 0; sigma2x <- 0.24; eff <- 0.14
alpha0 <- 0.01
alpha1 <- 0.005
rho0 <- 0.2
rho1 <- 0.1
ss1 <- calc_CS_HTE(K=K, CV=CV, p=0.5, sigma2y=1, sigma2x=sigma2x, eff=eff, alpha0=alpha0, alpha1=alpha1, rho0=rho0, rho1=rho1, alpha=0.05, beta=0.2)
ss1
I <- ss1[2]
power_CS_HTE(K=K, CV=CV, p=0.5, I=I, sigma2y=1, sigma2x=sigma2x, eff=eff, alpha0=alpha0, alpha1=alpha1, rho0=rho0, rho1=rho1, alpha=0.05)
  
conPLOT <- function(alpha0_range, alpha1t0_range, rho_type){
  con_plot <- NULL
  rocpart <- NULL
  rho0 <- rho0_range[rho_type]
  rho1 <- rho1_range[rho_type]
  for (alpha0 in seq(alpha0_range[1], alpha0_range[2], 0.001)){
    for (alpha1t0 in seq(alpha1t0_range[1], alpha1t0_range[2], 0.05)){
      alpha1 <- alpha0*alpha1t0
      rocpart<-rbind(rocpart, c(alpha0, alpha1t0, rho0, rho1, power_CS_HTE(K=K, CV=CV, p=0.5, I=I, sigma2y=1, sigma2x=sigma2x, eff=eff, 
                                           alpha0=alpha0, alpha1=alpha1, rho0=rho0, rho1=rho1, alpha=0.05)))
    }
  }
  rocpart <- as.data.frame(rocpart)
  colnames(rocpart) <- c("alpha0", "alpha1t0", "rho0", "rho1", "power")
  fig <- ggplot()  +
    theme_bw() +
    ggtitle(bquote(paste(rho[0] == ~.(rho0), ", ", rho[1] == ~.(rho1)))) +
    xlab(expression(alpha[0])) +
    ylab(expression(alpha[1]/alpha[0])) +
    xlim(c(alpha0_range[1], alpha0_range[2])) +
    ylim(c(alpha1t0_range[1], alpha1t0_range[2])) +
    geom_contour(data = rocpart, aes(x = alpha0, y = alpha1t0, z = power, colour = ..level..), 
                 breaks = round(quantile(rocpart$power, seq(0, 1, 0.1)), 2), size = 1) +
    scale_color_continuous(name = "Power") +
    theme(plot.title = element_text(hjust = 0.5), 
          legend.justification=c(1, 0), legend.position=c(1, 0),
          text = element_text(size = 15))
  con_plot <- direct.label(fig, "top.pieces")
  return(con_plot)
}

rho0_range <- c(0.15, 0.2, 0.3)
rho1_range <- c(0.1, 0.1, 0.15)
alpha0_range <- c(0, 0.2)
alpha1t0_range <- c(0, 1)

rocpart_1<-conPLOT(alpha0_range, alpha1t0_range, 1)
rocpart_2<-conPLOT(alpha0_range, alpha1t0_range, 2)
rocpart_3<-conPLOT(alpha0_range, alpha1t0_range, 3)

plot_grid(rocpart_1, rocpart_2, rocpart_3, ncol=3)





############################################################################################
#######################     Application to the ROSTERS study     ###########################
############################################################################################

### Cross-sectional design

K <- 30; CV <- 0; sigma2x <- 0.23; eff <- 0.5
alpha0 <- 0.13
alpha1 <- 0.05
rho0 <- 0.026
rho1 <- 0.026
ss2 <- calc_CS_HTE(K=K, CV=CV, p=0.5, sigma2y=1, sigma2x=sigma2x, eff=eff, alpha0=alpha0, alpha1=alpha1, rho0=rho0, rho1=rho1, alpha=0.05, beta=0.2)
ss2
I <- ss2[2]
power_CS_HTE(K=K, CV=CV, p=0.5, I=I, sigma2y=1, sigma2x=sigma2x, eff=eff, alpha0=alpha0, alpha1=alpha1, rho0=rho0, rho1=rho1, alpha=0.05)

conPLOT <- function(alpha0_range, alpha1t0_range, rho_type){
  con_plot <- NULL
  rosters <- NULL
  rho0 <- rho0_range[rho_type]
  rho1 <- rho1_range[rho_type]
  for (alpha0 in seq(alpha0_range[1], alpha0_range[2], 0.001)){
    for (alpha1t0 in seq(alpha1t0_range[1], alpha1t0_range[2], 0.05)){
      alpha1 <- alpha0*alpha1t0
      rosters<-rbind(rosters, c(alpha0, alpha1t0, rho0, rho1, power_CS_HTE(K=K, CV=CV, p=0.5, I=I, sigma2y=1, sigma2x=sigma2x, eff=eff, 
                                                                           alpha0=alpha0, alpha1=alpha1, rho0=rho0, rho1=rho1, alpha=0.05)))
    }
  }
  rosters <- as.data.frame(rosters)
  colnames(rosters) <- c("alpha0", "alpha1t0", "rho0", "rho1", "power")
  fig <- ggplot()  +
    theme_bw() +
    ggtitle(bquote(paste(rho[0] == ~.(rho0), ", ", rho[1] == ~.(rho1)))) +
    xlab(expression(alpha[0])) +
    ylab(expression(alpha[1]/alpha[0])) +
    xlim(c(alpha0_range[1], alpha0_range[2])) +
    ylim(c(alpha1t0_range[1], alpha1t0_range[2])) +
    geom_contour(data = rosters, aes(x = alpha0, y = alpha1t0, z = power, colour = ..level..), 
                 breaks = round(quantile(rosters$power, seq(0, 1, 0.1)), 2), size = 1) +
    scale_color_continuous(name = "Power") +
    theme(plot.title = element_text(hjust = 0.5),
          legend.justification=c(1, 0), legend.position=c(1, 0),
          text = element_text(size = 15))
  con_plot <- direct.label(fig, "top.pieces")
  return(con_plot)
}

rho0_range <- c(0.026, 0.1, 0.15)
rho1_range <- c(0.026, 0.026, 0.1)
alpha0_range <- c(0, 0.2)
alpha1t0_range <- c(0, 1)

rosters_1<-conPLOT(alpha0_range, alpha1t0_range, 1)
rosters_2<-conPLOT(alpha0_range, alpha1t0_range, 2)
rosters_3<-conPLOT(alpha0_range, alpha1t0_range, 3)

plot_grid(rosters_1, rosters_2, rosters_3, ncol=3)


### Closed-cohort design

K <- 30; CV <- 0; sigma2x <- 0.23; eff <- 0.5
alpha0 <- 0.13
alpha1 <- 0.05
alpha2 <- 0.2
rho0 <- 0.026
ss2 <- calc_CC_HTE(K=K, CV=CV, p=0.5, sigma2y=1, sigma2x=sigma2x, eff=eff, alpha0=alpha0, alpha1=alpha1, alpha2=alpha2, rho0=rho0, alpha=0.05, beta=0.2)
ss2
I <- ss2[2]
power_CC_HTE(K=K, CV=CV, p=0.5, I=I, sigma2y=1, sigma2x=sigma2x, eff=eff, alpha0=alpha0, alpha1=alpha1, alpha2=alpha2, rho0=rho0, alpha=0.05)

conPLOT <- function(alpha0_range, alpha1t0_range, rho_type){
  con_plot <- NULL
  rosters <- NULL
  alpha2 <- alpha2_range[rho_type]
  rho0 <- rho0_range[rho_type]
  for (alpha0 in seq(alpha0_range[1], alpha0_range[2], 0.001)){
    for (alpha1t0 in seq(alpha1t0_range[1], alpha1t0_range[2], 0.05)){
      alpha1 <- alpha0*alpha1t0
      rosters<-rbind(rosters, c(alpha0, alpha1t0, alpha2, rho0, power_CC_HTE(K=K, CV=CV, p=0.5, I=I, sigma2y=1, sigma2x=sigma2x, eff=eff, 
                                                                           alpha0=alpha0, alpha1=alpha1, alpha2=alpha2, rho0=rho0, alpha=0.05)))
    }
  }
  rosters <- as.data.frame(rosters)
  colnames(rosters) <- c("alpha0", "alpha1t0", "alpha2", "rho0", "power")
  fig <- ggplot()  +
    theme_bw() +
    ggtitle(bquote(paste(alpha[2] == ~.(alpha2), ", ", rho[0] == ~.(rho0)))) +
    xlab(expression(alpha[0])) +
    ylab(expression(alpha[1]/alpha[0])) +
    xlim(c(alpha0_range[1], alpha0_range[2])) +
    ylim(c(alpha1t0_range[1], alpha1t0_range[2])) +
    geom_contour(data = rosters, aes(x = alpha0, y = alpha1t0, z = power, colour = ..level..), 
                 breaks = round(quantile(rosters$power, seq(0, 1, 0.1)), 2), size = 1) +
    scale_color_continuous(name = "Power") +
    theme(plot.title = element_text(hjust = 0.5),
          legend.justification=c(1, 0), legend.position=c(1, 0),
          text = element_text(size = 15))
  con_plot <- direct.label(fig, "top.pieces")
  return(con_plot)
}

alpha2_range <- c(0.15, 0.2, 0.25)
rho0_range <- c(0.026, 0.1, 0.15)
alpha0_range <- c(0, 0.2)
alpha1t0_range <- c(0, 1)

rosters_1<-conPLOT(alpha0_range, alpha1t0_range, 1)
rosters_2<-conPLOT(alpha0_range, alpha1t0_range, 2)
rosters_3<-conPLOT(alpha0_range, alpha1t0_range, 3)

plot_grid(rosters_1, rosters_2, rosters_3, ncol=3)


