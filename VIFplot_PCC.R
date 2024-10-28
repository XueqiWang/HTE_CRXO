############################################################################################
################# Figures of RE for parallel CRT CRXO to closed-cohort #####################
############################################################################################

library(plotly)

EC_CC_P <- function(K=50, p=0.5, I=20, sigma2y=1, sigma2x=1, alpha0, alpha1=alpha0/2, alpha2=0.3, rho0){
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
  Var4_CC <- sigma2y/sigma2x / (I*p*(1-p)) / (d0)
  
  rhoy <- alpha0
  rhox <- rho0
  m <- 2*K
  Var4_P <- sigma2y*(1-rhoy)*(1+(m-1)*rhoy) / (0.5^2*sigma2x*m*I*(1+(m-2)*rhoy-(m-1)*rhox*rhoy))
  
  return(Var4_CC / Var4_P)
}

x <- seq(0, 0.999, 0.001)
y <- seq(0, 0.2, 0.001)
k20i20 <- k50i20 <- matrix(NA, length(x), length(y), dimnames = list(x,y))

for(i in 1:length(x)){
  for(j in 1:length(y)){
    k20i20[i,j] <- EC_CC_P(K=20, alpha0=y[j], rho0=x[i])

    k50i20[i,j] <- EC_CC_P(K=50, alpha0=y[j], rho0=x[i])
  }
}

common.range <- range(c(k20i20, k50i20))



################################################
############## Figures for K = 20 ##############
################################################

# p1: K = 20
fig <- plot_ly(x=y, y=x, z = k20i20, type = "contour",colorscale = 'Viridis',
               line = list(width = 0.5, color = "black"),reversescale = T,
               contours = list(showlabels = TRUE,
                               labelfont = list(size = 20, color = 'white'),
                               coloring = 'heatmap'),
               colorbar = list(tickfont=list(size=20, color='black'),
                               tickformat = ".4n",title = "RE",
                               titlefont=list(size=20))
)
fig <- fig %>% 
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
p_k20i20 <- colorbar(fig, limits = common.range)
p_k20i20



################################################
############## Figures for K = 50 ##############
################################################

# p1: K = 50
fig <- plot_ly(x=y, y=x, z = k50i20, type = "contour",colorscale = 'Viridis',
               line = list(width = 0.5, color = "black"),reversescale = T,
               contours = list(showlabels = TRUE,
                               labelfont = list(size = 20, color = 'white'),
                               coloring = 'heatmap'),
               colorbar = list(tickfont=list(size=20, color='black'),
                               tickformat = ".4n",title = "RE",
                               titlefont=list(size=20))
)
fig <- fig %>% 
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
p_k50i20 <- colorbar(fig, limits = common.range)
p_k50i20


