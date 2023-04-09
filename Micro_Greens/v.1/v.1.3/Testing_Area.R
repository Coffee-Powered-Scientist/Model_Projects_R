
rm(list = ls())

library(deSolve)

parms <- c(
  VI <- 0.1,
  VO <- 0.1,
  kag <- 1.22
  kgf <- 0.032, # forward rate constant for CO2 hydration (1/s)
  kgb <- 1.7e-4,# backward rate constant for CO2 hydration (1/s)
  k1f <- 4.45e-2, # forward rate constant for H2CO3 dissociation (1/s)
  k1b <- 5.6e-7, # backward rate constant for H2CO3 dissociation (1/s)
  k11f <- 5.6e-3, # forward rate constant for H2CO3 + H2O <=> H3O+ + HCO3- (1/s)
  k11b <- 5.6e-8, # backward rate constant for H2CO3 + H2O <=> H3O+ + HCO3- (1/s)
  k12f <- 5.01e-7, # forward rate constant for HCO3- + H2O <=> H3O+ + CO3^2- (1/s)
  k12b <- 5.01e-11, # backward rate constant for HCO3- + H2O <=> H3O+ + CO3^2- (1/s)
  PCO2 <- 1e-3, # partial pressure of CO2 in the system (atm)
  kwf <- 1e-7, # ion product constant for water (at 25°C)
  kwb <- 1e-7,
  CO2g_I = 0,
  CO2_I = 0,
  H2CO3=0,
  HCO3 = 0,
  CO3 = 0# ion product constant for water (at 25°C)
)

carbonate_eqs <- function(t, y, parms) {
CO2g<-y[8]
CO2 <- y[1]
H2CO3 <- y[2]
HCO3 <- y[3]
CO3 <- y[4]
H <- y[5]
OH <- y[6]
V<-y[7]

# Unpack Constants

VI <- parms["VI"]
VO <- parms["VO"]
kag <- parms["kag"]
kgf <- parms["kgf"]
kgb <- parms["kgb"]
k1f <- parms["k1f"]
k1b <- parms["k1b"]
k11f <- parms["k11f"]
k11b <- parms["k11b"]
k12f <- parms["k12f"]
k12b <- parms["k12b"]
PCO2 <- parms["PCO2"]
kwf <- parms["kwf"]
kwb <- parms["kwb"]

CO2g_I <- parms["CO2g_I"]
CO2_I <- parms["CO2_I"]
H2CO3_I <- parms["H2CO3_I"]
HCO3_I <- parms["HCO3_I"]
CO3_I <- parms["CO3_I"]
H_I <- parms["H_I"]
OH_I <- parms["OH_I"]


#dCO2g<- (VI*CO2g_I)/V + (kag*(PCO2-CO2g) - kgf*CO2g + kgb*CO2)/V - VO*CO2g/V
dCO2 <- VI*CO2_I/V+(- kgb*CO2 + k1b*H2CO3 - k1f*CO2)/V -VO*CO2/V
dH2CO3 <- VI*H2CO3_I/V + (k1f * CO2 - k1b * H2CO3 - k11f * H2CO3 + k11b * H * HCO3)/V - VO*H2CO3/V
dHCO3 <- VI*HCO3_I/V + (k11f * H2CO3 - k11b * H * HCO3 - k12f * HCO3 + k12b * H * CO3)/V - VO*HCO3/V
dCO3 <- VI*CO3_I/V + (k12f * HCO3 - k12b * H * CO3)/V - VI*CO3/V


H_Carb<- k12f*HCO3 - k12b*CO3*H + k11f*H2CO3 - k11b*H*HCO3

dH <- VI*H_I/V + (kwf - kwb * H * OH + H_Carb)/V  - VO*H/V
dOH <- VI*OH_I/V + (kwf - kwb * H * OH)/V - VO*OH/V

dV <- VI - VO


list(c(dCO2, dH2CO3, dHCO3, dCO3, dH, dOH, dV))
}
# Set the initial state variables and rate constants
y0 <- c(CO2 = 0, H2CO3 = 0, HCO3 = 0.0, CO3 = 0.00, H = 1e-7, OH = 1e-7, V=10)



# Set the time interval to simulate
t <- seq(0, 24, by = .1)

# Solve the system of differential equations using the ode solver in deSolve
out <- ode(y = y0, times = t, func = carbonate_eqs, parms = parms, method = "rk4")


plot(out, xlab = "Time (days)", ylab = "Concentration (M)")