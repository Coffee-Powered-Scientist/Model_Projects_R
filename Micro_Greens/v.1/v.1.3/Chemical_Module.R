
# Clear local environment

rm(list = ls())

# Load the deSolve package

setwd("F:/Model_Projects/Model_Projects_R/Micro_Greens/v.1/v.1.3")

library(deSolve)
library(ggplot2)

k_data <- read.csv("Reaction_Parameters.csv")
chem_data <- read.csv("Chemical_Nom.csv")

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

for(i in 1:nrow(chem_data)) {
  
  # Extract the species code and input concentration
  species_code <- chem_data[i, "Species_Code"]
  input_conc <- chem_data[i, "Input_Conc"]
  
  # Assign the input concentration to dynamically created variables
  assign(paste0(species_code, "_I"), input_conc)
  
  # Add the parameter value to the list
  parms[[paste0(species_code, "_I")]] <- input_conc
}

# Read in Water data

Water_Data<-read.csv("Water_Data.csv")
names(Water_Data) <- Water_Data[1,]
Water_Data <- Water_Data[-1,]

Water_Data <-Water_Data[ , colSums(is.na(Water_Data))==0]

Time_Conversion<-1 # seconds to hours 
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
  
  H2SO4<-y[8]
  HSO4<-y[9]
  SO4<-y[10]
  
  H3PO4<-y[11]
  H2PO4<-y[12]
  HPO4<-y[13]
  PO4<-y[14]
  
  NH4<-y[15]
  NH3<-y[16]
  HNO3<-y[17]
  NO3<-y[18]
  
  Ca<-y[19]
  CaOH<-y[20]
  CaOH2<-y[21]
  
  Mg<-y[22]
  MgOH<-y[23]
  MgOH2<-y[24]
  
  K<-y[25]
  KOH<-y[26]

  
  # Unpack the rate constants
  k1f <- parms$k1_0f*Time_Conversion
  k1b <- parms$k1_0b*Time_Conversion
  k11f <- parms$k1_1f*Time_Conversion
  k11b <- parms$k1_1b*Time_Conversion
  k12f <- parms$k1_2f*Time_Conversion
  k12b <- parms$k1_2b*Time_Conversion
  
  k21f <- parms$k2_1f*Time_Conversion
  k21b<-parms$k2_1b*Time_Conversion
  k22f <- parms$k2_2f*Time_Conversion
  k22b<-parms$k2_2b*Time_Conversion
  
  k31f<-parms$k3_1f*Time_Conversion
  k31b<-parms$k3_1b*Time_Conversion
  k32f<-parms$k3_2f*Time_Conversion
  k32b<-parms$k3_2b*Time_Conversion
  k33f<-parms$k3_3f*Time_Conversion
  k33b<-parms$k3_3b*Time_Conversion
  
  k4f<-parms$k4_0f*Time_Conversion
  k4b<-parms$k4_0b*Time_Conversion
  
  k5f<-parms$k5_0f*Time_Conversion
  k5b<-parms$k5_0b*Time_Conversion
  
  k61f<-parms$k5_0f*Time_Conversion
  k61b<-parms$k5_0b*Time_Conversion
  k62f<-parms$k5_0f*Time_Conversion
  k62b<-parms$k5_0b*Time_Conversion
  
  k7f<-parms$k5_0f*Time_Conversion
  k7b<-parms$k5_0b*Time_Conversion
  
  k81f<-parms$k5_0f*Time_Conversion
  k81b<-parms$k5_0b*Time_Conversion
  k82f<-parms$k5_0f*Time_Conversion
  k82b<-parms$k5_0b*Time_Conversion
  
  #CO2 and water
  kwf <- parms$kwf*Time_Conversion
  kwb <- parms$kwb*Time_Conversion
  kgf <- parms$kgf*Time_Conversion
  kgb <- parms$kgb*Time_Conversion
  
  # Manual Parameter Setting
  kH<-parms$kH<-1.7e-5#.034 # At 25°C and 1 atm pressure, the value of kH for CO2 is approximately 3.4 x 10^-2 mol/L/atm.
  PCO2<-parms$PCO2<-.0101 # PCO2 = (412 ppm) * (1/10^6) * (24.45 L/mol) = 0.0101 atm

  kag<-1.22*Time_Conversion
  # Water
  
  Break_Condition_1<-parms$BC1<-as.numeric(Water_Data$Break_Condition1)
  VI<-parms$VI<-as.numeric(Water_Data$Inflow_Rate)
  VO<-parms$VO<-as.numeric(Water_Data$Outflow_Rate)
  
  
  # Calculate the reaction rates
  ## Reaction 1: Bicarbonate
  dCO2 <- VI*CO2_I/V  + (k1b*H2CO3 - kgb*CO2 - k1f*CO2 +kag*.0101)/V -VO*CO2/V
  dH2CO3 <- VI*H2CO3_I/V + (k1f * CO2 - k1b * H2CO3 - k11f * H2CO3 + k11b * H * HCO3)/V - VO*H2CO3/V
  dHCO3 <- VI*HCO3_I/V + (k11f * H2CO3 - k11b * H * HCO3 - k12f * HCO3 + k12b * H * CO3)/V - VO*HCO3/V
  dCO3 <- VI*CO3_I/V + (k12f * HCO3 - k12b * H * CO3)/V - VI*CO3/V
  
  ## Reaction 2: Sulfuric Acid
  dH2SO4 <- VI*H2SO4_I/V + (-k21f*H2SO4 + k21b*HSO4*H)/V - VO*H2SO4/V
  dHSO4 <- VI*HSO4_I/V + (k21f*H2SO4 + k22b*SO4*H - k21b*HSO4*H - k22f*HSO4)/V - VO*HSO4/V
  dSO4 <- VI*SO4_I/V+ (k22f*HSO4 - k22b*SO4*H)/V - VO*SO4/V
  
  ## Reaction 3: Phosphate
  dH3PO4<- VI*H3PO4_I/V + (k31b*H2PO4*H - k31f*H3PO4)/V - VO*H3PO4
  dH2PO4<- VI*H2PO4_I/V + (k31f*H3PO4 + k32b*HPO4*H - k31b*H2PO4*H - k32f*H2PO4)/V - VO*H2PO4/V
  dHPO4<- VI*HPO4_I/V + (k32f*H2PO4 + k33b*PO4*H - k33f*HPO4 - k32b*HPO4*H)/V - VO*HPO4/V
  dPO4<- VI*PO4_I/V + (k33f*HPO4 - k33b*PO4*H)/V - VO*PO4/V
  
  ## Reaction 4: Ammonium/Ammonia
  #### considers NH4 as the more common ion, thus f-b relative to NH4
  
  dNH4<- VI*NH4_I/V + (k4f*NH3*H - k4b*NH4)/V - VO*NH4/V
  dNH3<- VI*NH3_I/V + (k4b*NH4 - k4f*NH3*H)/V - VO*NH3/V
  
  ## Reaction 5: Nitric acid/Nitrate
  
  dHNO3<- VI*HNO3/V + (k5b*NO3*H-k5f*HNO3)/V - VO*HNO3/V
  dNO3<- VI*NO3/V + (k5f*HNO3 - k5b*NO3*H)/V - VO*NO3/V
  
  ## Reaction 6: Calcium
  
  dCa <- VI*Ca_I/V + (k61b*CaOH - k61f*Ca*OH)/V - VO*Ca/V
  dCaOH <- VI*CaOH_I/V + (k61b*CaOH - k61f*Ca*OH)/V - VO*CaOH/V
  dCaOH2 <- VI*CaOH2_I/V + (k62f * CaOH * OH-k62b * CaOH2)/V - VO*CaOH2/V
  
  ## Reaction 7: Potassium
  
  dK<- VI*K/V + (k7b*KOH-k7f*K*OH)/V - VO*K/V
  dKOH<-VI*KOH/V + (k7f*K*OH - k7b*KOH)/V - VO*KOH/V
  
  ## Reaction 8: Magnesium
  dMg <- VI*Mg_I/V + (k81b*MgOH - k81f*Mg*OH)/V - VO*Mg/V
  dMgOH <- VI*MgOH_I/V + (k81b*MgOH - k81f*Mg*OH)/V - VO*MgOH/V
  dMgOH2 <- VI*MgOH2_I/V + (k82f * MgOH * OH-k82b * MgOH2)/V - VO*MgOH2/V
 
  ## Hydronium/Hydroxide Reactions
  
  ###Hydronium Equations per reaction
  H_Carb<- k12f*HCO3 - k12b*CO3*H + k11f*H2CO3 - k11b*H*HCO3
  H_Sulf<- k21f*H2SO4 - k22b*SO4*H + k22f*HSO4 -k21b*HSO4*H
  H_Phos<- k31f*H3PO4 + k32f*H2PO4 + k33f*HPO4 - k31b*H2PO4*H - k32b*HPO4*H - k33b*PO4*H
  H_Ammon<-k4b*NH4 - k4f*NH3*H
  H_Nitrate<-k5f*HNO3 - k5b*NO3*H
 
  ### Hydroxide Equaitons per reaction
  
  OH_Ca<-k61b*CaOH + k62b*CaOH2 - k61f*Ca*OH - k62f*CaOH*OH
  OH_Mg<-k81b*MgOH + k82b*MgOH2 - k81f*Mg*OH - k82f*MgOH*OH
  OH_K<-k7b*KOH - k7f*K*OH
  
  #### complete hydronium/hydroxide equations
    
  dH <- VI*H_I/V + (kwf - kwb * H * OH + H_Carb + H_Sulf + H_Phos + H_Ammon + H_Nitrate)/V  - VO*H/V
  dOH <- VI*OH_I/V + (kwf - kwb * H * OH + OH_Ca + OH_Mg + OH_K)/V - VO*OH/V
  
  
  # Water Equations
  if (V > 5) {
    dV <- VI - VO
  } else {
    dV <- ifelse(V <= 5 & VI - VO < 0, 0, VI - VO)
  }
  
  
  # Return the derivatives as a list
  list(c(dCO2, dH2CO3, dHCO3, dCO3, dH, dOH, dV, dH2SO4, dHSO4, dSO4, dH3PO4, dH2PO4, dHPO4, dPO4, dNH3, dNH4, dHNO3, dNO3,
         dCa, dCaOH, dCaOH2, dK, dKOH, dMg, dMgOH, dMgOH2))
}
# Set the initial state variables and rate constants
y0 <- c(CO2 = 0.0001, H2CO3 = 0, HCO3 = 0.0, CO3 = 0.00, H = 1e-7, OH = 1e-7, V=10,
        H2SO4 = 0, HSO4 = 0, SO4 = 0, H3PO4 = 0, H2PO4 = 0, HPO4 = 0, PO4 = 0,
        NH3 = 0, NH4 = 0, HNO3 = 0, NO3 = 0, Ca = 0, CaOH = 0, CaOH2 = 0, Mg = 0, MgOH=0, MgOH2 = 0, K = 0, KOH = 0)
#parms <- c(k1f = 1e-3, k1b = 1e-2, k11f = 1e-4, k11b = 1e-5, k12f = 1e-7, k12b = 1e-8, kwf = 1e-14, kwb = 1e-14)

# Set the time interval to simulate
t <- seq(0, 3600, by = .1)

# Solve the system of differential equations using the ode solver in deSolve
out <- ode(y = y0, times = t, func = carbonate_eqs, parms = parms, method = "rk4")

# Plot the results
plot(out, xlab = "Time (seconds)", ylab = "Concentration (M)")

Model_Dataframe<-as.data.frame(out)


