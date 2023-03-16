Model Project: Micro_Greens

Objective: Develop a multi-system plant growth and nutrient model, for the purposes of developing a model capable to simulating nutrient limitation for a hydroponic crop system. 

Description of Model v.1:

Micro_Greens v.1 has two core submodels, the Plant Model and the Water Model. The Plant Model is a modified Monod equation, relating nutrient uptake from water to plant growth. The water model is a modified "tank mixing" problem, where water enters a limited chamber of some volume, adding nutrients with every time step to the water. 

Micro_Greens v.1.1:

Another submodel, the nutrient model, should be added, as the modified Monad equation requires some calcution of uptake amount.

Assumptions

### Plant Model

Here I define two ways of modelling plant growth. The first way is using a modified Monad/Michelis Mention Equation, which uses a differential equation to model plant growth. The second is a discrete version of the model, further modified to have more realistic dependencies with respect to growth.  We use Liebling's law of the minimum, which states that the nutrient of the lowest availability will govern growth rate. 

Our objectives for this model is to:

1). Trace the change in nutrient storage in the plant over time

2). model reduction and acceleration of plant growth due to nutrient availability



We say that generally, the change in biomass over time is:

dX/dt = μ*X

where X is the biomass of the plant and μ is the specific growth rate. However, we need to modify this equation to take into account the multiple nutrients that are available. One way to do this is to use the concept of Liebig's law of the minimum, which states that the growth rate of an organism is limited by the nutrient that is in the shortest supply relative to its demand.

dX/dt = μ*X, $\mu  =min(\frac{V_{max_{C_1}}C_1}{K_1+C_1}, \frac{V_{max_{C_2}}C_2}{K_2+C_2})$,

where $C_1$ and $C_2$ are the concentrations of nutrients 1 and 2.  This equation links biomass growth rate to nutrient concentration in the solution. This is merely a continuous model linking nutrient availability to biomass, it assumes that the plant is always able to extract (this is not the case in reality). 

We now would like to trace how much nutrient goes into the plant. We can do this by using a simple stoichiometric constant $s$ , which multiplies by the change in biomass to yield the change in nutrient accrual. 

$\frac{dX_{C_1}}{dt} = s_1*\frac{dX}{dt}$

$\frac{dX_{C_2}}{dt} = s_2*\frac{dX}{dt}$

Water Model

Now, we would like to include a model that allows for the concentrations of nutrients to change over time. We can consider that the plant's roots grow in a tank of a fixed volume of 10 liters. Let's assume the tank is a perfect cube, as this will later allow us to write coded dependencies related to water level in the tank. For simplicity, we say that the tank adds the same amount of water at every time point as it subtracts. In the discrete version of the model, an order to this process will have to be defined (does addition or subtraction occur first?).  We first assume no evapotranspiration, we add a water plant uptake term such that our nutrient and water model make more sense. 

The differential equations governing change in water volume in the tank are thus:

v.1.0:

$\frac{dV}{dt}= Inputs - Outputs = V_I - V_o -\mu_w*\frac{dX}{dt}$,

v.1.01: Added Dependency

The v.1.0 equation was edited such that if the volume of the tank is reduced below 5 liters, then the VO parameters is set to 0 (i.e., no outflow). This assumption was made assuming that water pours in from the top of the tank, and is pumped out from a outlet situated at the 5 liter mark on the 10 liter tank.  

where $\mu_w$ is the volume of water uptaken for every unit of dry biomass produced by the plant, $V_I, V_o$ are the input and output water terms (like pumping water in and out of the tank). 



Nutrient Model

We now extend the water model to a nutrient model. The nutrient model considers gain of nutrients from the pumping in of nutrient containing water, and loss of nutrients from the pumping out of nutrient containing water, along with the removal of nutrients due to plant growth. 

$\frac{dC_1}{dt} = (C_I*V_I)-\frac{dX_{C_1}}{dt} - (C_1*V_O)$, then we divide by the volume of the solution to get the change in the concentration over time:

$\frac{d[C_1]}{dt} =\frac{(C_I*V_I + C_1*V)-\frac{dX_{C_1}}{dt}}{V}$



Below is the copy pasted R code of the above model, with version history comments added throughout. 

```R

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
      ## Added for v.1.01
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

```

