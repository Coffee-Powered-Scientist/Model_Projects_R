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
  k1f <- parms$k1_0f
  k1b <- parms$k1_0b
  k11f <- parms$k1_1f
  k11b <- parms$k1_1b
  k12f <- parms$k1_2f
  k12b <- parms$k1_2b
  kwf <- parms$kwf
  kwb <- parms$kwb
  kgf <- parms$kgf
  kgb <- parms$kgb
  
  # Calculate the reaction rates
  r1f <- kgf * CO2
  r1b <- k1b * H2CO3
  r11f <- k11f * H2CO3
  r11b <- k11b * H * HCO3
  r12f <- k12f * HCO3 - k12b * H * CO3
  r12b <- k12b * H * CO3
  rwf <- kwf
  rwb <- kwb * H * OH
  rgf <- kgf * CO2
  rgb <- kgb * H2CO3
  
  # Calculate the derivatives
  dCO2 <- rgf - r1f
  dH2CO3 <- r1f + r11b - r1b - r11f
  dHCO3 <- r11f + r12b - r11b - r12f - rgb
  dCO3 <- r12f - r12b + rgb
  dH <- rwf - rwb - r11b + r11f
  dOH <- rwf - rwb
  
  # Return the derivatives as a list
  list(c(dCO2, dH2CO3, dHCO3, dCO3, dH, dOH))
}