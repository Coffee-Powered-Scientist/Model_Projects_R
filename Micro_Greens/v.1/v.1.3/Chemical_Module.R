# Load the deSolve package

setwd("F:/Model_Projects/Model_Projects_R/Micro_Greens/v.1/v.1.3")

library(deSolve)

k_data <- read.csv("Reaction_Parameters.csv")
species_data <- read.csv("Chemical_Nom.csv")

# Initialize an empty list to store the parameter values
parms <- list()

# Loop through the rows of the data frame
for(i in 1:nrow(k_data)) {
  
  # Extract the reaction code, kf, and kb values
  reaction_code <- k_data[i, "Reaction_Code"]
  kf <- k_data[i, "kf"]
  kb <- k_data[i, "kb"]
  
  # Assign the kf and kb values to dynamically created variables
  assign(paste0("k", reaction_code, "f"), kf)
  assign(paste0("k", reaction_code, "b"), kb)
  
  # Add the parameter values to the list
  parms[[paste0("k", reaction_code, "f")]] <- kf
  parms[[paste0("k", reaction_code, "b")]] <- kb
}

# initialize lists to hold state variable names and initial concentrations


# Define the function for the system of differential equations
carbonate_eqs <- function(t, y, parms) {
  
  # Unpack the state variables
  CO2 <- y[1]
  H2CO3 <- y[2]
  HCO3 <- y[3]
  CO3 <- y[4]
  H <- y[5]
  OH<- y[6]
  
  # Unpack the rate constants
  k1f <- parms$k1_0f
  k1b <- parms$k1_0b
  k11f <- parms$k1_1f
  k11b <- parms$k1_1b
  k12f <- parms$k1_2f
  k12b <- parms$k1_2b
  kwf <- parms$kwf
  kwb <- parms$kwb
  
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
#parms <- c(k1f = 1e-3, k1b = 1e-2, k11f = 1e-4, k11b = 1e-5, k12f = 1e-7, k12b = 1e-8, kwf = 1e-14, kwb = 1e-14)

# Set the time interval to simulate
t <- seq(0, 100, by = 0.1)

# Solve the system of differential equations using the ode solver in deSolve
out <- ode(y = y0, times = t, func = carbonate_eqs, parms = parms)

# Plot the results
plot(out, xlab = "Time (s)", ylab = "Concentration (M)")
