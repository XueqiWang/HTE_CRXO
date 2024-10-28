############################################################################################
############## Figures of VIF for unequal to equal cluster-period sizes ####################
##############              for the cross-sectional design              ####################
############################################################################################

library(plotly)

VIF_CS_HTE <- function(K, CV, alpha0, alpha1=alpha0/2, rho0, rho1=rho0/2){
  # Argument:
  # K: mean cluster-period size
  # CV: coefficient of variation of cluster-period size
  # alpha0: within-period outcome ICC
  # alpha1: between-period outcome ICC
  # rho0: within-period covariate ICC
  # rho1: between-period covariate ICC
  
  lambda1 <- 1-alpha0
  lambda2 <- 1+(K-1)*alpha0-K*alpha1
  lambda3 <- 1+(K-1)*alpha0+K*alpha1
  zeta1 <- 1-rho0
  zeta2 <- 1+(K-1)*rho0-K*rho1
  zeta3 <- 1+(K-1)*rho0+K*rho1
  
  d0 <- 2*(K-1)*lambda1^{-1}*zeta1 + lambda2^{-1}*zeta3 + lambda3^{-1}*zeta2
  phi <- (2*K-(2*K-3)*zeta1-zeta2)*lambda1*lambda2^{-2} - 
    2*(K-(K-1)*zeta1)*lambda1^2*lambda2^{-3} + (zeta1+zeta2)*lambda1*lambda3^{-2} - zeta1*(lambda2^{-1}+lambda3^{-1})
  d <- 2*K*lambda2^{-1} + 2*(K-1)*zeta1*(lambda1^{-1}-lambda2^{-1}) + zeta2*(lambda3^{-1}-lambda2^{-1})
  
  return(d0 / (d-CV^2*phi))
}

x <- seq(0, 0.999, 0.001)
y <- seq(0, 0.2, 0.001)
k20cv25 <- k20cv50 <- k20cv75 <- matrix(NA, length(x), length(y), dimnames = list(x,y))
k50cv25 <- k50cv50 <- k50cv75 <- matrix(NA, length(x), length(y), dimnames = list(x,y))

for(i in 1:length(x)){
  for(j in 1:length(y)){
    k20cv25[i,j] <- VIF_CS_HTE(K=20, CV=0.25, alpha0=y[j], rho0=x[i])
    k20cv50[i,j] <- VIF_CS_HTE(K=20, CV=0.50, alpha0=y[j], rho0=x[i])
    k20cv75[i,j] <- VIF_CS_HTE(K=20, CV=0.75, alpha0=y[j], rho0=x[i])
    
    k50cv25[i,j] <- VIF_CS_HTE(K=50, CV=0.25, alpha0=y[j], rho0=x[i])
    k50cv50[i,j] <- VIF_CS_HTE(K=50, CV=0.50, alpha0=y[j], rho0=x[i])
    k50cv75[i,j] <- VIF_CS_HTE(K=50, CV=0.75, alpha0=y[j], rho0=x[i])
  }
}

common.range <- range(c(k20cv25, k20cv50, k20cv75, k50cv25, k50cv50, k50cv75))



################################################
############## Figures for K = 20 ##############
################################################

# p1: K = 20, CV = 0.25
fig <- plot_ly(x=y, y=x, z = k20cv25, type = "contour",colorscale = 'Viridis',
               line = list(width = 0.5, color = "black"),reversescale = T,
               contours = list(showlabels = TRUE,
                               labelfont = list(size = 20, color = 'white'),
                               coloring = 'heatmap'),
               colorbar = list(tickfont=list(size=20, color='black'),
                             tickformat = ".4n",title = "VIF",
                             titlefont=list(size=20))
)
fig <- fig %>% 
  #layout(title = list(text = "K = 20, CV = 0.25")) %>% 
  layout(xaxis = list(title = "\u03b1<sub><i>0</i></sub> = 2\u03b1<sub><i>1", titlefont = list(size=20), 
                      ticktext = c("0","0.05", "0.10", "0.15", "0.20"), 
                      tickvals = c(0, 0.05, 0.10, 0.15, 0.20), 
                      tickmode = "array", 
                      tickfont = list(size=20)),
         yaxis = list(title = "\u03c1<sub><i>0</i></sub> = 2\u03c1<sub><i>1", titlefont = list(size=20), 
                      ticktext = c("0", "0.10", "0.20", "0.30", "0.40", "0.50","0.60","0.70","0.80","0.90","1"),
                      tickvals = c(0, 0.10, 0.20, 0.30, 0.40, 0.50, 0.60, 0.70, 0.80, 0.90, 1), 
                      tickmode = "array",
                      tickfont = list(size=20)))
p_k20cv25 <- colorbar(fig, limits = common.range)
p_k20cv25


# p2: K = 20, CV = 0.5
fig <- plot_ly(x=y, y=x, z = k20cv50, type = "contour",colorscale = 'Viridis',
               line = list(width = 0.5, color = "black"),reversescale = T,
               contours = list(showlabels = TRUE,
                               labelfont = list(size = 20, color = 'white'),
                               coloring = 'heatmap'),
               colorbar = list(tickfont=list(size=20, color='black'),
                               tickformat = ".4n",title = "VIF",
                               titlefont=list(size=20))
)
fig <- fig %>% 
  #layout(title = list(text = "K = 20, CV = 0.50")) %>% 
  layout(xaxis = list(title = "\u03b1<sub><i>0</i></sub> = 2\u03b1<sub><i>1", titlefont = list(size=20), 
                      ticktext = c("0","0.05", "0.10", "0.15", "0.20"), 
                      tickvals = c(0, 0.05, 0.10, 0.15, 0.20), 
                      tickmode = "array", 
                      tickfont = list(size=20)),
         yaxis = list(title = "\u03c1<sub><i>0</i></sub> = 2\u03c1<sub><i>1", titlefont = list(size=20), 
                      ticktext = c("0", "0.10", "0.20", "0.30", "0.40", "0.50","0.60","0.70","0.80","0.90","1"),
                      tickvals = c(0, 0.10, 0.20, 0.30, 0.40, 0.50, 0.60, 0.70, 0.80, 0.90, 1), 
                      tickmode = "array",
                      tickfont = list(size=20)))
p_k20cv50 <- colorbar(fig, limits = common.range)
p_k20cv50


# p3: K = 20, CV = 0.75
fig <- plot_ly(x=y, y=x, z = k20cv75, type = "contour",colorscale = 'Viridis',
               line = list(width = 0.5, color = "black"),reversescale = T,
               contours = list(showlabels = TRUE,
                               labelfont = list(size = 20, color = 'white'),
                               coloring = 'heatmap'),
               colorbar = list(tickfont=list(size=20, color='black'),
                               tickformat = ".4n",title = "VIF",
                               titlefont=list(size=20))
)
fig <- fig %>% 
  #layout(title = list(text = "K = 20, CV = 0.75")) %>% 
  layout(xaxis = list(title = "\u03b1<sub><i>0</i></sub> = 2\u03b1<sub><i>1", titlefont = list(size=20), 
                      ticktext = c("0","0.05", "0.10", "0.15", "0.20"), 
                      tickvals = c(0, 0.05, 0.10, 0.15, 0.20), 
                      tickmode = "array", 
                      tickfont = list(size=20)),
         yaxis = list(title = "\u03c1<sub><i>0</i></sub> = 2\u03c1<sub><i>1", titlefont = list(size=20), 
                      ticktext = c("0", "0.10", "0.20", "0.30", "0.40", "0.50","0.60","0.70","0.80","0.90","1"),
                      tickvals = c(0, 0.10, 0.20, 0.30, 0.40, 0.50, 0.60, 0.70, 0.80, 0.90, 1), 
                      tickmode = "array",
                      tickfont = list(size=20)))
p_k20cv75 <- colorbar(fig, limits = common.range)
p_k20cv75



################################################
############## Figures for K = 50 ##############
################################################

# p1: K = 50, CV = 0.25
fig <- plot_ly(x=y, y=x, z = k50cv25, type = "contour",colorscale = 'Viridis',
               line = list(width = 0.5, color = "black"),reversescale = T,
               contours = list(showlabels = TRUE,
                               labelfont = list(size = 20, color = 'white'),
                               coloring = 'heatmap'),
               colorbar = list(tickfont=list(size=20, color='black'),
                               tickformat = ".4n",title = "VIF",
                               titlefont=list(size=20))
)
fig <- fig %>% 
  #layout(title = list(text = "K = 50, CV = 0.25")) %>% 
  layout(xaxis = list(title = "\u03b1<sub><i>0</i></sub> = 2\u03b1<sub><i>1", titlefont = list(size=20), 
                      ticktext = c("0","0.05", "0.10", "0.15", "0.20"), 
                      tickvals = c(0, 0.05, 0.10, 0.15, 0.20), 
                      tickmode = "array", 
                      tickfont = list(size=20)),
         yaxis = list(title = "\u03c1<sub><i>0</i></sub> = 2\u03c1<sub><i>1", titlefont = list(size=20), 
                      ticktext = c("0", "0.10", "0.20", "0.30", "0.40", "0.50","0.60","0.70","0.80","0.90","1"),
                      tickvals = c(0, 0.10, 0.20, 0.30, 0.40, 0.50, 0.60, 0.70, 0.80, 0.90, 1), 
                      tickmode = "array",
                      tickfont = list(size=20)))
p_k50cv25 <- colorbar(fig, limits = common.range)
p_k50cv25


# p2: K = 50, CV = 0.5
fig <- plot_ly(x=y, y=x, z = k50cv50, type = "contour",colorscale = 'Viridis',
               line = list(width = 0.5, color = "black"),reversescale = T,
               contours = list(showlabels = TRUE,
                               labelfont = list(size = 20, color = 'white'),
                               coloring = 'heatmap'),
               colorbar = list(tickfont=list(size=20, color='black'),
                               tickformat = ".4n",title = "VIF",
                               titlefont=list(size=20))
)
fig <- fig %>% 
  #layout(title = list(text = "K = 50, CV = 0.50")) %>% 
  layout(xaxis = list(title = "\u03b1<sub><i>0</i></sub> = 2\u03b1<sub><i>1", titlefont = list(size=20), 
                      ticktext = c("0","0.05", "0.10", "0.15", "0.20"), 
                      tickvals = c(0, 0.05, 0.10, 0.15, 0.20), 
                      tickmode = "array", 
                      tickfont = list(size=20)),
         yaxis = list(title = "\u03c1<sub><i>0</i></sub> = 2\u03c1<sub><i>1", titlefont = list(size=20), 
                      ticktext = c("0", "0.10", "0.20", "0.30", "0.40", "0.50","0.60","0.70","0.80","0.90","1"),
                      tickvals = c(0, 0.10, 0.20, 0.30, 0.40, 0.50, 0.60, 0.70, 0.80, 0.90, 1), 
                      tickmode = "array",
                      tickfont = list(size=20)))
p_k50cv50 <- colorbar(fig, limits = common.range)
p_k50cv50


# p3: K = 50, CV = 0.75
fig <- plot_ly(x=y, y=x, z = k50cv75, type = "contour",colorscale = 'Viridis',
               line = list(width = 0.5, color = "black"),reversescale = T,
               contours = list(showlabels = TRUE,
                               labelfont = list(size = 20, color = 'white'),
                               coloring = 'heatmap'),
               colorbar = list(tickfont=list(size=20, color='black'),
                               tickformat = ".4n",title = "VIF",
                               titlefont=list(size=20))
)
fig <- fig %>% 
  #layout(title = list(text = "K = 50, CV = 0.75")) %>% 
  layout(xaxis = list(title = "\u03b1<sub><i>0</i></sub> = 2\u03b1<sub><i>1", titlefont = list(size=20), 
                      ticktext = c("0","0.05", "0.10", "0.15", "0.20"), 
                      tickvals = c(0, 0.05, 0.10, 0.15, 0.20), 
                      tickmode = "array", 
                      tickfont = list(size=20)),
         yaxis = list(title = "\u03c1<sub><i>0</i></sub> = 2\u03c1<sub><i>1", titlefont = list(size=20), 
                      ticktext = c("0", "0.10", "0.20", "0.30", "0.40", "0.50","0.60","0.70","0.80","0.90","1"),
                      tickvals = c(0, 0.10, 0.20, 0.30, 0.40, 0.50, 0.60, 0.70, 0.80, 0.90, 1), 
                      tickmode = "array",
                      tickfont = list(size=20)))
p_k50cv75 <- colorbar(fig, limits = common.range)
p_k50cv75


# final plot
g_legend <- function(a.gplot){
  tmp <- ggplot_gtable(ggplot_build(a.gplot))
  leg <- which(sapply(tmp$grobs, function(x) x$name) == "guide-box")
  legend <- tmp$grobs[[leg]]
  return(legend)
}

mylegend <- g_legend(p_k20cv25)
pplot_om <- grid.arrange(arrangeGrob(p_k20cv25 + theme(legend.position="none"),
                                     p_k20cv50 + theme(legend.position="none"),
                                     p_k20cv75 + theme(legend.position="none"),
                                     p_k50cv25 + theme(legend.position="none"),
                                     p_k50cv50 + theme(legend.position="none"),
                                     p_k50cv75 + theme(legend.position="none"),
                                     nrow=2),
                         nrow=2, heights=c(12, 1))



fig <- subplot(p_k20cv25, p_k20cv50, p_k20cv75, p_k50cv25, p_k50cv50, p_k50cv75, nrows = 2) %>% 
  layout(title = 'Side By Side Subplots')
fig


