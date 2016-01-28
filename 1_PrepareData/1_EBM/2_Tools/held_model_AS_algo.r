#######################################################################################
# EBM
#######################################################################################
# Based on Geoffroy et al. (2013)
models_name = c("BCC","CCCMA","CNRM","CSIRO","GFDL","INM","IPSL","MIROC","MPIM","MRI","NCC")

F    = c(6.70, 7.60, 7.10, 5.10, 6.60, 6.20, 6.40, 8.50, 8.20, 6.60, 6.20)
lamb = c(1.21, 1.03, 1.08, 0.61, 1.34, 1.51, 0.78, 1.58, 1.14, 1.26, 1.11)
gamm = c(0.71, 0.69, 0.47, 0.92, 0.87, 0.66, 0.60, 0.64, 0.76, 0.73, 0.95)
c    = c(7.80, 7.30, 10.0, 6.20, 8.10, 8.90, 9.90, 8.20, 7.30, 8.90, 7.80)
c0   = c(52.0, 64.0, 130., 68.0, 110., 314., 98.0, 134., 70.0, 62.0, 100.)

parameters <- data.frame(models_name, F, lamb, gamm, c, c0)

# Two-box model for the climate system (ocean-atmosphere).
# c  dT/dt   = F - lamb T - gamm (T - T_0)
# c0 dT_0/dt = gamm (T - T_0)
# Based on Held et al. (2010)

#-- Arguments :
# forcing : step - linear - stabilization
# F_infty : amplitude of the forcing
# c       : atm/upper ocean heat capacity
# c0      : deep ocean heat capacity
# lamb    : feedback parameter
# gamm    : ocean heat uptake coefficient
# N       : number of points simulated   

# function [T, To, H] = held_model(forcing, F_infty, c, c0, lamb, gamm, N)
held_model  <- function (FF, model="CNRM"){
  model_param <- subset(parameters, models_name == model)
  list2env(model_param, environment())
  #-0- Forcing functions
  N <- length(FF)
  dt <- 1; #-- timestep (year)

  #-1- Numerical solutions (explicit scheme)
  T <- numeric(N+1);
  To <- numeric(N+1);
  H <- numeric(N+1);
  T[1] <- 0;
  To[1] <- 0;
  for(n in 2:(N+1)){
    T[n]  <-  (T[n-1] + dt/c*(FF[n-1] - lamb*T[n-1] - gamm*(T[n-1]-To[n-1])));
    To[n] <-  (To[n-1] + dt/c0*(gamm*(T[n-1]-To[n-1])));
    H[n]  <-  gamm*(T[n] - To[n]);
  }
  ans <- cbind(T, To, H)
  colnames(ans) <- c("T", "To", "H")
  ans
}

