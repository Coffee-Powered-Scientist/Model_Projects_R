
rm(list=ls())

library(deSolve)

## Units are as follows ##
  
#VM_C1: maximum uptake rate of nutrient 1 (g / g / day)
#VM_C2: maximum uptake rate of nutrient 2 (g / g / day)
#K1: half-saturation constant for nutrient 1 uptake (g / g)
#K2: half-saturation constant for nutrient 2 uptake (g / g)
#s1: stoichiometric constant for nutrient 1 accumulation (g nutrient 1 / g biomass)
#s2: stoichiometric constant for nutrient 2 accumulation (g nutrient 2 / g biomass)
#VI: inflow rate of water (L / day)
#VO: outflow rate of water (L / day)
#mu_w: specific growth rate of water (day^-1)
#CI: initial concentration of nutrients in water (g / L)

#X: biomass of the plant (g)
#C1: concentration of nutrient 1 in the water (g / L)
#C2: concentration of nutrient 2 in the water (g / L)
#V: volume of water in the system (L)

# Define parameter values
parameters <- c(VM_C1 = 0.05, VM_C2 = 0.05, K1 = 0.01, K2 = 0.01, s1 = 0.0012,
                s2 = 0.0012, VI = 1, VO = 2, mu_w = 2, CI1 = 0.01, CI2 = .01)

# Define initial state values
initial_state <- c(X = 0.01, XC1 = 0.001, XC2 = 0.001, V = 10, C1 = 0.05, C2 = 0.05)

# Define the ODE system
Micro_Greens <- function(time, state, parameters) {
  with(as.list(c(state, parameters)), {
    
    # Expressions to simplify DEs
    mu = min(VM_C1*C1/(K1+C1), VM_C2*C2/(K2+C2))
    
    # Define the system of ODE
    
    # Plant and Nutrients in Plant
    dXdt <- mu*X
    dXC1dt <- s1*dXdt
    dXC2dt <- s2*dXdt
    
    # Water
    #dVdt <- VI - VO- mu_w*dXdt
    
    # Water
    if (V > 5) {
      dVdt <- VI - VO - mu_w * dXdt
    } else {
      dVdt <- ifelse(V <= 5 & VI - VO - mu_w * dXdt < 0, 0, VI - VO - mu_w * dXdt)
    }
    
    # Nutrients
    dC1dt <- (CI1*VI - C1*VO - s1*dXdt)/V
    dC2dt <- (CI2*VI  - C2*VO - s2*dXdt)/V
    
    # Return the system of ODEs
    return(list(c(dXdt, dXC1dt, dXC2dt, dVdt, dC1dt, dC2dt)))
  })
}



# Define time points for simulation
times <- seq(0, 10, by = 0.1)

# Solve the ODE system using the deSolve package
output <- ode(y = initial_state, times = times, func = Micro_Greens, parms = parameters)

