################################
#        Numerical study
################################

source("functions_calc_power.R")

library(mvtnorm)
library(ggplot2)
library(gridExtra)

cv <- 0


############################
# CS
############################

# (A) Vary rho0 = 2*rho1 for different (alpha0, alpha1)
alpha0_range <- c(0.02, 0.05, 0.05, 0.1, 0.1)
alpha1_range <- c(0.01, 0.01, 0.05, 0.05, 0.1)
rho0_range <- seq(0, 1, 0.01) # x-axis
data_all <- NULL
for (i in 1:length(alpha0_range)){
  alpha0 <- alpha0_range[i]
  alpha1 <- alpha1_range[i]
  for (rho0 in rho0_range){
    rho1 <- rho0/2
    pred.power <- power_CS_HTE(CV=cv, alpha0=alpha0, alpha1=alpha1, rho0=rho0, rho1=rho1)
    data_all <- rbind(data_all, c(i=i, cv=cv, alpha0=alpha0, alpha1=alpha1, rho0=rho0, rho1=rho1, pred.power=pred.power))
  }
}
pdata_cs_1 <- as.data.frame(data_all)
min(pdata_cs_1[, 7]) # 0.2771762
max(pdata_cs_1[, 7]) # 0.828339

pplot_cs_1 <- ggplot(data=pdata_cs_1, aes(rho0, pred.power, colour = factor(i), shape = factor(i))) +
  geom_line(cex = 1.1, aes(linetype=factor(i))) + 
  scale_color_manual(values = c("darkorange2", "limegreen", "mediumpurple2", "lightskyblue", "palevioletred1"),
                     breaks = c("1", "2", "3", "4", "5"),
                     labels = c("(0.02, 0.01)", "(0.05, 0.01)", "(0.05, 0.05)", "(0.10, 0.05)", "(0.10, 0.10)")) +
  scale_linetype_manual(values = c(1, 5, 6, 2, 4),
                     breaks = c("1", "2", "3", "4", "5"),
                     labels = c("(0.02, 0.01)", "(0.05, 0.01)", "(0.05, 0.05)", "(0.10, 0.05)", "(0.10, 0.10)")) +
  labs(title = expression(paste("(A) Vary ", rho[0], "=2", rho[1])),
       y = "Power", x = expression(paste(rho[0],"=2",rho[1])),
       colour = " ", linetype = " ") +
  scale_x_continuous(breaks=seq(0, 1, 0.2), limits=c(0, 1)) +
  scale_y_continuous(breaks=seq(0.3, 0.9, 0.1), limits=c(0.27, 0.95)) +
  theme_bw() +
  theme(legend.position = "bottom", legend.box = "vertical", legend.margin = margin(), 
        text = element_text(size = 15), strip.text = element_text(size = 15), plot.title = element_text(size = 15)) +
  labs(color = expression(paste("(", alpha[0], ", ", alpha[1], ")")),
       linetype = expression(paste("(", alpha[0], ", ", alpha[1], ")"))) +
  guides(color = guide_legend(nrow=5, byrow=TRUE))
pplot_cs_1


#(B) Vary alpha0 = 2*alpha1 for different (rho0, rho1)
alpha0_range <- seq(0, 0.5, 0.01)  # x-axis
rho0_range <- c(0.15, 0.3, 0.3, 0.5, 0.5)
rho1_range <- c(0.1, 0.1, 0.3, 0.3, 0.5)
data_all <- NULL
for (i in 1:length(rho0_range)){
  rho0 <- rho0_range[i]
  rho1 <- rho1_range[i]
  for (alpha0 in alpha0_range){
    alpha1 <- alpha0/2
    pred.power <- power_CS_HTE(CV=cv, alpha0=alpha0, alpha1=alpha1, rho0=rho0, rho1=rho1)
    data_all <- rbind(data_all, c(i=i, cv=cv, alpha0=alpha0, alpha1=alpha1, rho0=rho0, rho1=rho1, pred.power=pred.power))
  }
}
pdata_cs_2 <- as.data.frame(data_all)
min(pdata_cs_2[, 7]) # 0.6198997
max(pdata_cs_2[, 7]) # 0.9475011

pplot_cs_2 <- ggplot(data=pdata_cs_2, aes(alpha0, pred.power, colour = factor(i), shape = factor(i))) +
  geom_line(cex = 1.1, aes(linetype=factor(i))) + 
  scale_color_manual(values = c("darkorange2", "limegreen", "mediumpurple2", "lightskyblue", "palevioletred1"),
                     breaks = c("1", "2", "3", "4", "5"),
                     labels = c("(0.15, 0.10)", "(0.30, 0.10)", "(0.30, 0.30)", "(0.50, 0.30)", "(0.50, 0.50)")) +
  scale_linetype_manual(values = c(1, 5, 6, 2, 4),
                        breaks = c("1", "2", "3", "4", "5"),
                        labels = c("(0.15, 0.10)", "(0.30, 0.10)", "(0.30, 0.30)", "(0.50, 0.30)", "(0.50, 0.50)")) +
  labs(title = expression(paste("(B) Vary ", alpha[0], "=2", alpha[1])),
       y = "Power", x = expression(paste(alpha[0],"=2",alpha[1])),
       colour = " ", linetype = " ") +
  scale_x_continuous(breaks=seq(0, 0.5, 0.1), limits=c(0, 0.5)) +
  scale_y_continuous(breaks=seq(0.3, 0.9, 0.1), limits=c(0.27, 0.95)) +
  theme_bw() +
  theme(legend.position = "bottom", legend.box = "vertical", legend.margin = margin(), 
        text = element_text(size = 15), strip.text = element_text(size = 15), plot.title = element_text(size = 15)) +
  labs(color = expression(paste("(", rho[0], ", ", rho[1], ")")),
       linetype = expression(paste("(", rho[0], ", ", rho[1], ")"))) +
  guides(color = guide_legend(nrow=5, byrow=TRUE))
pplot_cs_2


### final plot
pplot_cs <- grid.arrange(arrangeGrob(pplot_cs_1, pplot_cs_2, nrow=1))



############################
# CC
############################

# (A) Vary rho0 for (alpha0, alpha1, alpha2)
alpha0_range <- c(0.015, 0.1, 0.1, 0.1, 0.15)
alpha1_range <- c(0.01, 0.01, 0.05, 0.05, 0.05)
alpha2_range <- c(0.2, 0.2, 0.2, 0.5, 0.5)
rho0_range <- seq(0, 1, 0.01) # x-axis
data_all <- NULL
for (i in 1:length(alpha0_range)){
  alpha0 <- alpha0_range[i]
  alpha1 <- alpha1_range[i]
  alpha2 <- alpha2_range[i]
  for (rho0 in rho0_range){
    pred.power <- power_CC_HTE(CV=cv, alpha0=alpha0, alpha1=alpha1, alpha2=alpha2, rho0=rho0)
    data_all <- rbind(data_all, c(i=i, cv=cv, alpha0=alpha0, alpha1=alpha1, alpha2=alpha2, rho0=rho0, pred.power=pred.power))
  }
}
pdata_cc_1 <- as.data.frame(data_all)
min(pdata_cc_1[, 7]) # 0.2157035
max(pdata_cc_1[, 7]) # 0.9903196

pplot_cc_1 <- ggplot(data=pdata_cc_1, aes(rho0, pred.power, colour = factor(i), shape = factor(i))) +
  geom_line(cex = 1.1, aes(linetype=factor(i))) + 
  scale_color_manual(values = c("darkorange2", "limegreen", "mediumpurple2", "lightskyblue", "palevioletred1"),
                     breaks = c("1", "2", "3", "4", "5"),
                     labels = c("(0.015, 0.010, 0.200)", "(0.100, 0.010, 0.200)", "(0.100, 0.050, 0.200)", "(0.100, 0.050, 0.500)", "(0.150, 0.050, 0.500)")) +
  scale_linetype_manual(values = c(1, 5, 6, 2, 4),
                        breaks = c("1", "2", "3", "4", "5"),
                        labels = c("(0.015, 0.010, 0.200)", "(0.100, 0.010, 0.200)", "(0.100, 0.050, 0.200)", "(0.100, 0.050, 0.500)", "(0.150, 0.050, 0.500)")) +
  labs(title = expression(paste("(A) Vary ", rho[0])),
       y = "Power", x = expression(paste(rho[0])),
       colour = " ", linetype = " ") +
  scale_x_continuous(breaks=seq(0, 1, 0.2), limits=c(0, 1)) +
  scale_y_continuous(breaks=seq(0.2, 1.0, 0.1), limits=c(0.21, 1)) +
  theme_bw() +
  theme(legend.position = "bottom", legend.box = "vertical", legend.margin = margin(), 
        text = element_text(size = 15), strip.text = element_text(size = 15), plot.title = element_text(size = 15)) +
  labs(color = expression(paste("(", alpha[0], ", ", alpha[1], ", ", alpha[2], ")")),
       linetype = expression(paste("(", alpha[0], ", ", alpha[1], ", ", alpha[2], ")"))) +
  guides(color = guide_legend(nrow=5, byrow=TRUE))
pplot_cc_1


#(B) Vary alpha0 = 2*alpha1 for (alpha2, rho0)
alpha0_range <- seq(0, 0.5, 0.01) # x-axis
alpha2_range <- c(0.1, 0.2, 0.2, 0.5, 0.5)
rho0_range <- c(0.1, 0.1, 0.2, 0.2, 0.5)
data_all <- NULL
for (i in 1:length(alpha2_range)){
  alpha2 <- alpha2_range[i]
  rho0 <- rho0_range[i]
  for (alpha0 in alpha0_range){
    alpha1 <- alpha0/2
    pred.power <- power_CC_HTE(CV=cv, alpha0=alpha0, alpha1=alpha1, alpha2=alpha2, rho0=rho0)
    data_all <- rbind(data_all, c(i=i, cv=cv, alpha0=alpha0, alpha1=alpha1, alpha2=alpha2, rho0=rho0, pred.power=pred.power))
  }
}
pdata_cc_2 <- as.data.frame(data_all)
min(pdata_cc_2[, 7]) # 0.8153356
max(pdata_cc_2[, 7]) # 0.9979319

pplot_cc_2 <- ggplot(data=pdata_cc_2, aes(alpha0, pred.power, colour = factor(i), shape = factor(i))) +
  geom_line(cex = 1.1, aes(linetype=factor(i))) + 
  scale_color_manual(values = c("darkorange2", "limegreen", "mediumpurple2", "lightskyblue", "palevioletred1"),
                     breaks = c("1", "2", "3", "4", "5"),
                     labels = c("(0.1, 0.1)", "(0.2, 0.1)", "(0.2, 0.2)", "(0.5, 0.2)", "(0.5, 0.5)")) +
  scale_linetype_manual(values = c(1, 5, 6, 2, 4),
                        breaks = c("1", "2", "3", "4", "5"),
                        labels = c("(0.1, 0.1)", "(0.2, 0.1)", "(0.2, 0.2)", "(0.5, 0.2)", "(0.5, 0.5)")) +
  labs(title = expression(paste("(B) Vary ", alpha[0], "=2", alpha[1])),
       y = "Power", x = expression(paste(alpha[0],"=2",alpha[1])),
       colour = " ", linetype = " ") +
  scale_x_continuous(breaks=seq(0, 0.5, 0.1)) +
  scale_y_continuous(breaks=seq(0.2, 1.0, 0.1), limits=c(0.21, 1)) +
  theme_bw() +
  theme(legend.position = "bottom", legend.box = "vertical", legend.margin = margin(), 
        text = element_text(size = 15), strip.text = element_text(size = 15), plot.title = element_text(size = 15)) +
  labs(color = expression(paste("(", alpha[2], ", ", rho[0], ")")),
       linetype = expression(paste("(", alpha[2], ", ", rho[0], ")"))) +
  guides(color = guide_legend(nrow=5, byrow=TRUE))
pplot_cc_2


#(C) Vary alpha2 for (alpha0, alpha1, rho2)
alpha0_range <- c(0.015, 0.1, 0.1, 0.1, 0.15)
alpha1_range <- c(0.01, 0.01, 0.05, 0.05, 0.05)
alpha2_range <- seq(0, 0.5, 0.01)  # x-axis
rho0_range <- c(0.2, 0.2, 0.2, 0.5, 0.5)
data_all <- NULL
for (i in 1:length(alpha0_range)){
  alpha0 <- alpha0_range[i]
  alpha1 <- alpha1_range[i]
  rho0 <- rho0_range[i]
  for (alpha2 in alpha2_range){
    pred.power <- power_CC_HTE(CV=cv, alpha0=alpha0, alpha1=alpha1, alpha2=alpha2, rho0=rho0)
    data_all <- rbind(data_all, c(i=i, cv=cv, alpha0=alpha0, alpha1=alpha1, alpha2=alpha2, rho0=rho0, pred.power=pred.power))
  }
}
pdata_cc_3 <- as.data.frame(data_all)
min(pdata_cc_3[, 7]) # 0.592757
max(pdata_cc_3[, 7]) # 0.9711963

pplot_cc_3 <- ggplot(data=pdata_cc_3, aes(alpha2, pred.power, colour = factor(i), shape = factor(i))) +
  geom_line(cex = 1.1, aes(linetype=factor(i))) + 
  scale_color_manual(values = c("darkorange2", "limegreen", "mediumpurple2", "lightskyblue", "palevioletred1"),
                     breaks = c("1", "2", "3", "4", "5"),
                     labels = c("(0.015, 0.010, 0.200)", "(0.100, 0.010, 0.200)", "(0.100, 0.050, 0.200)", "(0.100, 0.050, 0.500)", "(0.150, 0.050, 0.500)")) +
  scale_linetype_manual(values = c(1, 5, 6, 2, 4),
                        breaks = c("1", "2", "3", "4", "5"),
                        labels = c("(0.015, 0.010, 0.200)", "(0.100, 0.010, 0.200)", "(0.100, 0.050, 0.200)", "(0.100, 0.050, 0.500)", "(0.150, 0.050, 0.500)")) +
  labs(title = expression(paste("(C) Vary ", alpha[2])),
       y = "Power", x = expression(paste(alpha[2])),
       colour = " ", linetype = " ") +
  scale_x_continuous(breaks=seq(0, 0.5, 0.1)) +
  scale_y_continuous(breaks=seq(0.2, 1.0, 0.1), limits=c(0.21, 1)) +
  theme_bw() +
  theme(legend.position = "bottom", legend.box = "vertical", legend.margin = margin(), 
        text = element_text(size = 15), strip.text = element_text(size = 15), plot.title = element_text(size = 15)) +
  labs(color = expression(paste("(", alpha[0], ", ", alpha[1], ", ", rho[0], ")")),
       linetype = expression(paste("(", alpha[0], ", ", alpha[1], ", ", rho[0], ")"))) +
  guides(color = guide_legend(nrow=5, byrow=TRUE))
pplot_cc_3


### final plot
pplot_cc <- grid.arrange(arrangeGrob(pplot_cc_1, pplot_cc_2, pplot_cc_3, nrow=1))


