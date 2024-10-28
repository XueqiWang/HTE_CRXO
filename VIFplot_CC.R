############################################################################################
############## Figures of VIF for unequal to equal cluster-period sizes ####################
##############               for the closed-cohort design               ####################
############################################################################################

library(plotly)

VIF_CC_HTE <- function(K, CV, alpha0, alpha1=alpha0/2, alpha2=0.3, rho0){
  # Argument:
  # K: mean cluster-period size
  # CV: coefficient of variation of cluster-period size
  # alpha0: within-period outcome ICC
  # alpha1: between-period outcome ICC
  # alpha2: within-individual outcome ICC
  # rho0: within-period covariate ICC

  tau1 <- 1-alpha0+alpha1-alpha2
  tau2 <- 1-alpha0-(alpha1-alpha2)
  tau3 <- 1+(K-1)*(alpha0-alpha1)-alpha2
  tau4 <- 1+(K-1)*alpha0+(K-1)*alpha1+alpha2
  eta1 <- 1-rho0
  eta2 <- 1+(K-1)*rho0
  
  d0 <- 2*(K-1)*tau1^{-1}*eta1 + 2*tau3^{-1}*eta2
  phi <- 2*((K-(K-2)*eta1)*tau1*tau3^{-2} - (K-(K-1)*eta1)*tau1^2*tau3^{-3} - eta1*tau3^{-1})
  d <- 2*(K*tau3^{-1} + (K-1)*eta1*(tau1^{-1}-tau3^{-1}))
  
  return(d0 / (d-CV^2*phi))
}

x <- seq(0, 0.999, 0.001)
y <- seq(0, 0.2, 0.001)
k20cv25 <- k20cv50 <- k20cv75 <- matrix(NA, length(x), length(y), dimnames = list(x,y))
k50cv25 <- k50cv50 <- k50cv75 <- matrix(NA, length(x), length(y), dimnames = list(x,y))

for(i in 1:length(x)){
  for(j in 1:length(y)){
    k20cv25[i,j] <- VIF_CC_HTE(K=20, CV=0.25, alpha0=y[j], rho0=x[i])
    k20cv50[i,j] <- VIF_CC_HTE(K=20, CV=0.50, alpha0=y[j], rho0=x[i])
    k20cv75[i,j] <- VIF_CC_HTE(K=20, CV=0.75, alpha0=y[j], rho0=x[i])
    
    k50cv25[i,j] <- VIF_CC_HTE(K=50, CV=0.25, alpha0=y[j], rho0=x[i])
    k50cv50[i,j] <- VIF_CC_HTE(K=50, CV=0.50, alpha0=y[j], rho0=x[i])
    k50cv75[i,j] <- VIF_CC_HTE(K=50, CV=0.75, alpha0=y[j], rho0=x[i])
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
         yaxis = list(title = "\u03c1<sub><i>0", titlefont = list(size=20), 
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
         yaxis = list(title = "\u03c1<sub><i>0", titlefont = list(size=20), 
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
         yaxis = list(title = "\u03c1<sub><i>0", titlefont = list(size=20), 
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
         yaxis = list(title = "\u03c1<sub><i>0", titlefont = list(size=20), 
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
         yaxis = list(title = "\u03c1<sub><i>0", titlefont = list(size=20), 
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
         yaxis = list(title = "\u03c1<sub><i>0", titlefont = list(size=20), 
                      ticktext = c("0", "0.10", "0.20", "0.30", "0.40", "0.50","0.60","0.70","0.80","0.90","1"),
                      tickvals = c(0, 0.10, 0.20, 0.30, 0.40, 0.50, 0.60, 0.70, 0.80, 0.90, 1), 
                      tickmode = "array",
                      tickfont = list(size=20)))
p_k50cv75 <- colorbar(fig, limits = common.range)
p_k50cv75


