
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



# Read in Water data

Water_Data<-read.csv("Water_Data.csv")
names(Water_Data) <- Water_Data[1,]
Water_Data <- Water_Data[-1,]

Water_Data <-Water_Data[ , colSums(is.na(Water_Data))==0]

Time_Conversion<-1/24 # seconds to days to hours
# Define the function for the system of differential equations
carbonate_eqs <- function(t, y, parms) {
  
  # Unpack the state variables
  
  # Water
  
  V<-y[7]
  
  # Chemistry
  CO2 <- y[1]
  H2CO3 <- y[2]
  HCO3 <- y[3]
  CO3 <- y[4]
  H <- y[5]
  OH <- y[6]
  
  # Unpack the rate constants
  k1f <- parms$k1_0f*Time_Conversion
  k1b <- parms$k1_0b*Time_Conversion
  k11f <- parms$k1_1f*Time_Conversion
  k11b <- parms$k1_1b*Time_Conversion
  k12f <- parms$k1_2f*Time_Conversion
  k12b <- parms$k1_2b*Time_Conversion
  kwf <- parms$kwf*Time_Conversion
  kwb <- parms$kwb*Time_Conversion
  kgf <- parms$kgf*Time_Conversion
  kgb <- parms$kgb*Time_Conversion
  
  # Manual Parameter Setting
  kH<-parms$kH<-1.7e-5#.034 # At 25°C and 1 atm pressure, the value of kH for CO2 is approximately 3.4 x 10^-2 mol/L/atm.
  PCO2<-parms$PCO2<-.0101 # PCO2 = (412 ppm) * (1/10^6) * (24.45 L/mol) = 0.0101 atm
  CO2_I<-0
  H2CO3_I<-0
  HCO3_I<-0
  CO3_I<-0
  H_I<-10^-7
  OH_I<-10^-7
  # Water
  
  Break_Condition_1<-parms$BC1<-as.numeric(Water_Data$Break_Condition1)
  VI<-parms$VI<-as.numeric(Water_Data$Inflow_Rate)
  VO<-parms$VO<-as.numeric(Water_Data$Outflow_Rate)
  
  
  # Calculate the reaction rates
  dCO2 <- VI*CO2_I/V+(kgf*(PCO2 - kH*CO2) - k1f*CO2 + k1b*H2CO3 - kgb*CO2)/V -VO*CO2/V
  dH2CO3 <- VI*H2CO3_I/V + (k1f * CO2 - k1b * H2CO3 - k11f * H2CO3 + k11b * H * HCO3)/V - VO*H2CO3/V
  dHCO3 <- VI*HCO3_I/V + (k11f * H2CO3 - k11b * H * HCO3 - k12f * HCO3 + k12b * H * CO3)/V - VO*HCO3/V
  dCO3 <- VI*CO3_I/V + (k12f * HCO3 - k12b * H * CO3)/V - VI*CO3/V
  dH <- VI*H_I/V + (kwf - kwb * H * OH - k11f*H2CO3 - k11b*H*HCO3 + k12f*HCO3 -k12b*CO3*H)/V - VO*H/V
  dOH <- VI*OH_I/V + (kwf - kwb * H * OH)/V - VO*OH/V
  
  # Water Equations
  if (V > 5) {
    dV <- VI - VO
  } else {
    dV <- ifelse(V <= 5 & VI - VO < 0, 0, VI - VO)
  }
  
  
  # Return the derivatives as a list
  list(c(dCO2, dH2CO3, dHCO3, dCO3, dH, dOH, dV))
}
# Set the initial state variables and rate constants
y0 <- c(CO2 = 0, H2CO3 = 0, HCO3 = 0.0, CO3 = 0.00, H = 1e-7, OH = 1e-7, V=10)
#parms <- c(k1f = 1e-3, k1b = 1e-2, k11f = 1e-4, k11b = 1e-5, k12f = 1e-7, k12b = 1e-8, kwf = 1e-14, kwb = 1e-14)

# Set the time interval to simulate
t <- seq(0, 240, by = 1)

# Solve the system of differential equations using the ode solver in deSolve
out <- ode(y = y0, times = t, func = carbonate_eqs, parms = parms, method = "lsoda")

# Plot the results
plot(out, xlab = "Time (days)", ylab = "Concentration (M)")

Model_Dataframe<-as.data.frame(out)
