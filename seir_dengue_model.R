# SEIR-SIR Model for Dengue Disease in R (Human + Vector)
# Required packages: deSolve, ggplot2, reshape2

# Install required packages if not already installed
if (!require(deSolve)) install.packages('deSolve', dependencies=TRUE)
if (!require(ggplot2)) install.packages('ggplot2', dependencies=TRUE)
if (!require(reshape2)) install.packages('reshape2', dependencies=TRUE)
library(deSolve)
library(ggplot2)
library(reshape2)

# SEIR-SIR model function (humans: SEIR, vectors: SIR)
dengue_seir_sir <- function(time, state, parameters) {
  with(as.list(c(state, parameters)), {
    # Human population
    N_h <- S_h + E_h + I_h + R_h
    # Vector population
    N_v <- S_v + I_v + R_v
    # ODEs
    dS_h <- -beta_hv * S_h * I_v / N_v
    dE_h <- beta_hv * S_h * I_v / N_v - sigma_h * E_h
    dI_h <- sigma_h * E_h - gamma_h * I_h
    dR_h <- gamma_h * I_h
    dS_v <- -beta_vh * S_v * I_h / N_h
    dI_v <- beta_vh * S_v * I_h / N_h - gamma_v * I_v
    dR_v <- gamma_v * I_v
    return(list(c(dS_h, dE_h, dI_h, dR_h, dS_v, dI_v, dR_v)))
  })
}

# Initial state values
init <- c(
  S_h = 990000,  # Susceptible humans
  E_h = 5000,    # Exposed humans
  I_h = 5000,    # Infectious humans
  R_h = 0,       # Recovered humans
  S_v = 100000,  # Susceptible vectors
  I_v = 1000,    # Infectious vectors
  R_v = 0        # Recovered vectors (often 0 for dengue)
)

# Parameters (example values)
params <- c(
  beta_hv = 0.5,    # Transmission rate: vector -> human
  beta_vh = 0.5,    # Transmission rate: human -> vector
  sigma_h = 1/5,    # Incubation rate in humans (1/incubation period)
  gamma_h = 1/7,    # Recovery rate in humans (1/infectious period)
  gamma_v = 0       # Recovery rate in vectors (0 for dengue, as mosquitoes do not recover)
)

# Time vector (days)
times <- seq(0, 180, by = 1)

# Solve ODE
out <- ode(y = init, times = times, func = dengue_seir_sir, parms = params)
out <- as.data.frame(out)

# Prepare data for plotting
out_long_h <- melt(out[, c('time', 'S_h', 'E_h', 'I_h', 'R_h')], id = 'time')
out_long_v <- melt(out[, c('time', 'S_v', 'I_v', 'R_v')], id = 'time')

# Plot human compartments
p1 <- ggplot(out_long_h, aes(x = time, y = value, color = variable)) +
  geom_line(size = 1) +
  labs(title = 'SEIR Model for Humans (Dengue)',
       x = 'Time (days)',
       y = 'Number of Humans',
       color = 'Compartment') +
  theme_minimal()

# Plot vector compartments
p2 <- ggplot(out_long_v, aes(x = time, y = value, color = variable)) +
  geom_line(size = 1) +
  labs(title = 'SIR Model for Vectors (Mosquitoes)',
       x = 'Time (days)',
       y = 'Number of Vectors',
       color = 'Compartment') +
  theme_minimal()

print(p1)
print(p2) 
