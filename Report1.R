

library(deSolve)


par(mfrow=c(1,1))
u <- 0.5 #m/day  # Vertical water velocity (advection velocity) settling velocity

DeltaZ <- 1 # Vertical spacing between grid cells
kw <- 0.2
kp <- 0.05 #(umol N/m^3)^-1/m

gmax <- 1 #/day
I0 <- 350 #umol photons/m2/s(E)
yield <- 0.00005

Hi <- 20 #(20*(60*60*24))/1000
Hn <- 0.0425
kn <- 0.3
m <- 0.24 #/day
alpha <- 0.3

D <- 5 # m^2/day
d <- 100
n <- 100 # Number of grid cells along the water column
z_values <- seq( DeltaZ/2, d, by = DeltaZ)
e <- 0.03 # d^-1 
w <- 15
grazing <- 0.1 #/day
remin <- 1.5 #(mmol N/m^3)^-1/day 
Av <- 5
# Initial conditions 
Nb <- 50
P <- rep(30,n)
N <- rep(10,n)
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

g <- gmax * pmin(I/(I+Hi), N/(N+Hn)) - m 
#r <- g-m 

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
  J_D[n + 1] <- -w*De[n]
  J_n[n + 1] <- -D * ((Nb-N[n])/DeltaZ)
  
  J <- J_a + J_d
  

  Light <- function(kw,kp,DeltaZ,P,De,I0){
    
    #Integral calculations - cumulative sums of the light in depth
    integral <- cumsum((kw+kp*P)*DeltaZ)+(-1/2*DeltaZ)*(kw+kp*P)
    
    #calculation the light values 
    I = I0*exp(-integral)
    
    return(I)
  }
  I = Light(kw,kp,DeltaZ,P,I0,De)
  
  G.n <- N/(kn+N)
  G.l <- (alpha*I)/sqrt(gmax^2+(alpha^2*I^2))
  g <- gmax * pmin(G.n, G.l)- m
  
  dP_dt <- seq (1,n,1)
  dN_dt <- seq (1,n,1)
  dD_dt <- seq (1,n,1)
  
  for (i in 1:n){
    dP_dt[i] <- (gmax*g[i]*P[i])-(m*P[i]) - (grazing*P[i]^2) + Av * (-J[i+1]- J[i])/DeltaZ
    dN_dt[i] <- (-gmax*g[i]*P[i])+ (remin*De[i]) + Av * (-J_n[i+1]- J_n[i])/DeltaZ
    dD_dt[i] <- (m*P[i])+(grazing*P[i]^2) -(remin*De[i]) - (w*De[i]) + Av * (-J_D[i+1]+ J_D[i])/DeltaZ
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
    yield = yield, 
    Hi = Hi,
    Hn = Hn, 
    grazing = grazing, 
    remin = remin,
    Av = Av,
    I0 = I0, 
    kw = kw,
    kp = kp, 
    alpha=alpha,
    kn = kn
    
  )
times <- seq(0, 100, by = 1) # Time vector
result <- ode(y = c(P,N,De), times = times, func = derivative, parms = params)

P_out <- result[, 2:(n+1)]
N_out <- result[, (n+2):(2*n+1)]
D_out <- result[, (n+3):(3*n+1)]
Time <- result [,1]
Depth <- z_values


image(Time, Depth, P_out[,ncol(P_out):1], col = hcl.colors(30, "viridis"), main = "phytoplankton")
image(Time, Depth, N_out[,ncol(P_out):1], col = hcl.colors(50, "viridis"), main = "Nutrients")
image(Time, Depth, D_out[,ncol(P_out):1], col = hcl.colors(50, "viridis"), main = "Detritus")


plot(z_values, I, type = "l", col = "yellow", lwd=3) #ylim=c(0,50))
plot(Time, N_out[,2], type = "l", col = "darkorange", lwd=3) # ylim = c(0,50))
plot(Time, D_out[,2], type = "l", col = "brown", lwd=3)
plot(Time, P_out[,2], type = "l", col = "darkgreen", lwd=3)

#It seems like there is an issue with the code in the way it understands time and depth


#Growth limitation plot nutrients and light 
g.l.L <- Light(kw, kp, DeltaZ, P_out[nrow(P_out), ], I0) / (Hi + Light(kw, kp, DeltaZ, P_out[nrow(P_out), ], I0))
g.l.N <- N_out[nrow(N_out), ] / (Hn + N_out[nrow(N_out), ])

plot(
  x = g.l.L,
  y = n:1,
  type = "l",
  xlim = c(0, max(g.l.N)),
  main = "Limitation by light or nutrients",
  ylab = "Distance from seabed[m]",
  xlab = "growth limitiation factor",
  lwd = 3
)
lines(x = g.l.N, y = n:1, type = "l", lwd = 3, lty = 2)
#lines(x=P_out[,2], y=n:1)
legend(
  "topright",
  #inset = c(-0.5, 0),
  legend = c("Light", "nutrients"),
  lty = c(1, 2),
  lwd = 2
)

#Sensitivity analysis 

#Calculating for different values for D 
D_values <- c(5, 10, 30, 50)

# Create a new plot window
plot(times, N_out[,2], type = "l", col = 1, lwd = 2, 
     main = "Phytoplankton Concentration for D[i]", 
     xlab = "Time", ylab = "Phytoplankton Concentration")

# Iterate over each value of D
for (i in 1:length(D_values)) {
  D <- D_values[i]  # Update D value
  
  # Solve ODEs
  Result <- ode(y=c(P,N), times = times, func = derivative, 
                parms = list(u = u, D = D, kw = kw, kp = kp, Iin = Iin, 
                             HI = HI, gmax = gmax, m = m, HN = HN, Nb = Nb, yield = yield))
  
  N_out <- Result[, (n+1):(2*n+1)] # Extract phytoplankton output
  
  # Plot phytoplankton concentration for each D value
  lines(times, N_out[,2], type = "l", col = i+1, lwd = 2)
}



#Calculating the maximum concentration of phytoplankton 
find_phytoplankton_maximum <- function(P_out, z) {
  # Find the index of the maximum phytoplankton concentration
  max_index <- which.max(P_out)
  
  # Calculate the position and magnitude of the maximum
  max_depth <- z[max_index]
  max_concentration <- P_out[max_index]
  
  return(list(position = max_depth, magnitude = max_concentration))
}

max_phytoplankton <- find_phytoplankton_maximum(P_out, z)
print(max_phytoplankton)




