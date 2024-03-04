
install.packages("fields")
install.packages("desolve")

library(deSolve)
library(fields)


par(mfrow=c(1,1))
u <- 0.5 #m/day  # Vertical water velocity (advection velocity) settling velocity

d <- 100
n <- 100 # Number of grid cells along the water column
DeltaZ <- d/n # Vertical spacing between grid cells
kw <- 0.2
kp <- 0.05 #(umol N/m^3)^-1/m

gmax <- 1 #/day
I0 <- 350 #umol photons/m2/s(E)
kn <- 0.3
m <- 0.24 #/day
alpha <- 0.1

D <- 4.32 # m^2/day
z_values <- seq( 0.5*DeltaZ, d, by = DeltaZ)
grazing <- 0.1 #/day
remin <- 1.5 #(mmol N/m^3)^-1/day 
Av <- 4.32 #converted #When my vertical advection is 5 there is more phyto and when it is 1 there is less phyto

# Initial conditions 
Nb <- 5
P <- rep(5,n)
N <- rep(0.2,n)
De <- rep(0,n)


#--------------------------------------------------------
Light <- function(kw,kp,DeltaZ,P,I0){
  
  #Integral calculations - cumulative sums of the light in depth
  integral <- cumsum((kw+kp*P)*DeltaZ)+(-1/2*DeltaZ)*(kw+kp*P)
  
  #calculation the light values 
  I = I0*exp(-integral)
  
  return(I)
}

#------------------------------------------------------------

I = Light(kw,kp,DeltaZ,P,I0)
plot(z_values, I)

G.n <- N/(kn+N)
G.l <- (alpha*I)/sqrt(gmax^2+(alpha^2*I^2))
g <- gmax * pmin(G.n, G.l)- m

#------------------------------------------------
#Doing the diffusion and advection 
#------------------------------------------------

derivative <- function(times, y, params) {
  P <- y[1:n]
  N <- y[(n+1):(2*n)]
  De <- y[(2*n+1):(3*n)]
  
  J_a <- numeric(n + 1)
  J_d <- numeric(n + 1)
  J_n <- numeric(n + 1)
  J_D <- numeric(n + 1)
  
  for(i in 2:n){
    J_a[i] <-  u*P[i-1]
    J_d[i] <- -D*(P[i]-P[i-1])/DeltaZ
    J_n[i] <- -D*(N[i]-N[i-1])/DeltaZ 
    J_D[i] <- -D*(De[i]-De[i-1])/DeltaZ 
  }
  
  # Set boundary fluxes to 0
  J_a[1] <- 0 # Boundary condition
  J_d[1] <- 0
  J_n[1] <- 0 
  J_D[1] <- 0 
  J_a[n + 1] <- 0
  J_d[n + 1] <- 0
  J_D[n + 1] <- u*De[n]
  J_n[n + 1] <- -D * ((Nb-N[n])/DeltaZ)
  
  J = J_a + J_d
  
  jP= (-(J[2:(n+1)] - J[1:n])) / DeltaZ
  jN= -(J_n[2:(n+1)] - J_n[1:n]) / DeltaZ
  jD= -(J_D[2:(n+1)] - J_D[1:n]) / DeltaZ
 
  
  Light <- function(kw,kp,DeltaZ,P,I0){
    
    #Integral calculations - cumulative sums of the light in depth
    integral <- cumsum((kw+kp*P)*DeltaZ)+(-1/2*DeltaZ)*(kw+kp*P)
    
    #calculation the light values 
    I = I0*exp(-integral)
    
    return(I)
  }
  I = Light(kw,kp,DeltaZ,P,I0)
  
  G.n <- N/(kn+N)
  G.l <- (alpha*I)/sqrt(gmax^2+(alpha^2*I^2))
  g <- gmax * pmin(G.n, G.l)- m
  
  
  dP_dt <- numeric(n)
  dN_dt <- numeric(n)
  dD_dt <- numeric(n)
  
  for (i in 1:n){
    dP_dt[i] <- (g[i]*P[i])-(m*P[i]) - (grazing*P[i]^2) + jP[i]
    dN_dt[i] <- (-g[i]*P[i])+ (remin*De[i]) + jN[i]
    dD_dt[i] <- (m*P[i])+(grazing*P[i]^2) -(remin*De[i]) - (u*De[i]) + jD[i]
  }
  
  return(list(c(dP_dt, dN_dt,dD_dt)))
}


params <-
  list(
    n = n,
    DeltaZ = DeltaZ,
    u = u,
    D = D,
    m = m,
    g = g, 
    gmax = gmax, 
    I = I, 
    Nb = Nb, 
    grazing = grazing, 
    remin = remin,
    Av = Av,
    I0 = I0, 
    kw = kw,
    kp = kp, 
    alpha=alpha,
    kn = kn
    
  )
times <- seq(0, 1000, by = 1) # Time vector
result <- ode(y = c(P,N,De), times = times, func = derivative, parms = params)

P_out <- result[, 2:(n+1)]
N_out <- result[, (n+2):(2*n+1)]
D_out <- result[, (2*n+2):(3*n+1)]
Time <- result [,1]
z_values <- z_values


image.plot(Time, z_values, P_out[,ncol(P_out):1], col = hcl.colors(50, "viridis"), main = "Phytoplankton [mmol N/m^3]", xlab = "Days", ylab = "Distance from seabed [m]", cex.lab=1.25)
image.plot(Time, z_values, N_out[,ncol(P_out):1], col = hcl.colors(50, "viridis"), main = "Nutrients [mmol N/m^3]", xlab = "Days", ylab = "Distance from seabed [m]", cex.lab=1.25)
image.plot(Time, z_values, D_out[,ncol(P_out):1], col = hcl.colors(50, "viridis"), main = "Detritus [mmol N/m^3]", xlab = "Days", ylab = "Distance from seabed [m]", cex.lab=1.25)


plot(N_out[ dim(N_out)[2],], -z_values, type = "l", col = "peru", lwd=4, xlim = c(0,5), xlab = "mmol N/m^3", ylab = "Depth [m]", main = "NPD model - Steady state", cex.lab=1.25)
plot(P_out[ dim(P_out)[2],], -z_values, type = "l", col = "darkgreen", lwd=4)
lines(D_out[ dim(D_out)[2],], -z_values, type = "l", col = "brown4", lwd=4)
legend(
  "topright",
  #inset = c(-0.5, 0),
  legend = c("Phytoplankton", "Nutrients", "Detritus"),
  lty = c(1, 1,1),
  lwd = 3,
  col = c("darkgreen","peru", "brown4")
)



#Growth limitation plot nutrients and light 
g.l.L <- Light(kw, kp, DeltaZ, P_out[nrow(P_out), ], I0) / (alpha + Light(kw, kp, DeltaZ, P_out[nrow(P_out), ], I0))
g.l.N <- N_out[nrow(N_out), ] / (kn + N_out[nrow(N_out), ])

plot(
  x = g.l.L,
  y = n:1,
  type = "l",
  #xlim = c(0,1),
  main = "Growth limitation by light or nutrients",
  ylab = "Distance from seabed[m]",
  xlab = "Growth Limitation Factor",
  lwd = 3
)
lines(x = g.l.N, y = n:1, type = "l", lwd = 3, lty = 2)
lines(x=P_out[dim(P_out)[2],], y=n:1, type = "l", lty = 3, lwd=3, col = "darkgreen")
legend(
  "bottom",
  #inset = c(-0.5, 0),
  legend = c("Light", "Nutrients", "Phytoplankton"),
  lty = c(1, 2,3),
  lwd = 2
)


#---------with and without damping ---------------------------

Light <- function(kw,kp,DeltaZ,P,I0){
  
  #Integral calculations - cumulative sums of the light in depth
  integral <- cumsum((kw+kp*P)*DeltaZ)+(-1/2*DeltaZ)*(kw+kp*P)
  
  #calculation the light values 
  I = I0*exp(-integral)
  
  return(I)
}

#------------------------------------------------------------

I = Light(kw,kp,DeltaZ,P,I0)
plot(I, -z_values, type = "l", xlim = c(0, 300), lwd = 3, main = "Light attenuation with and without dampering", xlab = "Watts m^-2", ylab = "Depth [m]", col = "khaki3")

kp1=0

Light <- function(kw,kp1,DeltaZ,P,I0){
  
  #Integral calculations - cumulative sums of the light in depth
  integral <- cumsum((kw+kp1*P)*DeltaZ)+(-1/2*DeltaZ)*(kw+kp1*P)
  
  #calculation the light values 
  I2 = I0*exp(-integral)
  
  return(I2)
}


I2 = Light(kw,kp1,DeltaZ,P,I0)
lines(I2, -z_values, type = "l", lty=2, lwd = 3, col = "khaki3")
legend(
  "bottomright",
  legend = c("with damping", "without damping"),
  lty = c(1, 2),
  lwd = 3, 
  col = c("khaki3", "khaki3")
)

#Notes 
#When i make gmax smaller than the n and d become much bigger and raises longer
#Making Hi 30 will make it not as high but will come with a bump around 200
#The higher remineralisation give higher nutrients mmol N/m3, phytoplankton mmol N/m^3
#SInking velocity can change the phytoplankton. 
#Light absorption 


