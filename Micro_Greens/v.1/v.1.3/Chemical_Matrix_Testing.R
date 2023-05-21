library(dplyr)
library(tidyr)

setwd("F:/Model_Projects/Model_Projects_R/Micro_Greens/v.1/v.1.3")


components <- read.csv("Component_Matrix.csv")


cations <- filter(components, Class == "Cation")
anions <- filter(components, Class == "Anion")

combinations <- expand.grid(cations$Name, anions$Name)

combinations <- expand.grid(cations$Name, anions$Name)

combinations <- rename(combinations, Cation = Var1, Anion = Var2)

combinations$Compound <- paste0(combinations$Cation, combinations$Anion)
combinations$Charge <- cations$Charge[cations$Name == combinations$Cation] + anions$Charge[anions$Name == combinations$Anion]







