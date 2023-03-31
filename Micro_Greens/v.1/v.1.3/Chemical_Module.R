# Load the deSolve package
library(deSolve)

# Define the function for the system of differential equations
carbonate_eqs <- function(t, y, parms) {
  
  # Unpack the state variables
  CO2 <- y[1]
  H2CO3 <- y[2]
  HCO3 <- y[3]
  CO3 <- y[4]
  H <- y[5]
  OH <- y[6]
  
  # Unpack the rate constants
  k1f <- parms[1]
  k1b <- parms[2]
  k11f <- parms[3]
  k11b <- parms[4]
  k12f <- parms[5]
  k12b <- parms[6]
  kwf <- parms[7]
  kwb <- parms[8]
  
  # Calculate the reaction rates
  r1f <- k1f * CO2
  r1b <- k1b * H2CO3
  r11f <- k11f * H2CO3
  r11b <- k11b * H * HCO3
  r12f <- k12f * HCO3 - k12b * H * CO3
  r12b <- k12b * H * CO3
  rwf <- kwf
  rwb <- kwb * H * OH
  
  # Calculate the derivatives
  dCO2 <- -r1f + r1b
  dH2CO3 <- r1f + r11b - r1b - r11f
  dHCO3 <- r11f + r12b - r11b - r12f
  dCO3 <- r12f - r12b
  dH <- rwf - rwb - r11b + r11f
  dOH <- rwf - rwb
  
  # Return the derivatives as a list
  list(c(dCO2, dH2CO3, dHCO3, dCO3, dH, dOH))
}

# Set the initial state variables and rate constants
y0 <- c(CO2 = 0.1, H2CO3 = 0.5, HCO3 = 0.0, CO3 = 0.00, H = 1e-6, OH = 1e-8)
parms <- c(k1f = 1e-3, k1b = 1e-2, k11f = 1e-4, k11b = 1e-5, k12f = 1e-7, k12b = 1e-8, kwf = 1e-14, kwb = 1e-14)

# Set the time interval to simulate
t <- seq(0, 100, by = 0.1)

# Solve the system of differential equations using the ode solver in deSolve
out <- ode(y = y0, times = t, func = carbonate_eqs, parms = parms)

# Plot the results
plot(out, xlab = "Time (s)", ylab = "Concentration (M)")
