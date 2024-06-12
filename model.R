# Load necessary libraries
library(deSolve)
library(ggplot2)

data <- guecDat

# Data
date <- seq(as.Date('1976-08-25'), by = "day", length.out = 61)
cases <- c(
  1, 0, 0, 0, 0, 0, 0, 1, 1, 1, 5, 1, 1, 3, 2, 6, 10, 4, 9, 9, 8, 10, 6, 6, 5, 8,
  11, 6, 16, 12, 12, 11, 7, 7, 4, 5, 3, 8, 3, 3, 8, 5, 7, 5, 2, 4, 4, 3, 6, 2, 1,
  0, 0, 0, 1, 1, 2, 0, 0, 0, 1
)

data <- data.frame(date = date, cases = cases)


# Data
#date <- data['date']
#cases <- data['cases']

#data <- data.frame(date = date, cases = cases)







N <- 5000
I0 <- cases[1]
E0 <- 0  # Assuming initially there are no exposed individuals
R0 <- 0  # Initially no one has recovered
S0 <- N - I0 - E0 - R0

# SEIR model differential equations
seir_model <- function(t, y, parms) {
  S <- y[1]
  E <- y[2]
  I <- y[3]
  R <- y[4]
  
  beta <- parms["beta"]
  sigma <- parms["sigma"]
  gamma <- parms["gamma"]
  
  dSdt <- -beta * S * I / N
  dEdt <- beta * S * I / N - sigma * E
  dIdt <- sigma * E - gamma * I
  dRdt <- gamma * I
  
  list(c(dSdt, dEdt, dIdt, dRdt))
}

# Initial conditions vector
y0 <- c(S = S0, E = E0, I = I0, R = R0)

# Time grid
t <- seq(0, length(cases) - 1, by = 1)
# Assume some initial beta, sigma, and gamma

beta_calc <- function(
    R_0, mu = 0.02 / 365.25, sigma = 1 / 8, gamma = 1 / 5
) {
  R_0 / ((sigma / (mu + sigma)) * (1 / (mu + gamma)))
}

parms <- c(beta = beta_calc(1), sigma = 0.2, gamma = 0.15)

# Integrate the SEIR equations over the time grid, t
out <- ode(y = y0, times = t, func = seir_model, parms = parms)
out <- as.data.frame(out)

ggplot() +
  geom_line(data = out, aes(x = t, y = I), color = 'red', linetype = "dashed", size = 1, alpha = 0.7) +
  geom_point(data = data, aes(x = as.numeric(date - min(date)), y = cases), color = 'red') +
  labs(x = 'Date', y = 'Number of cases', title = 'SEIR Model Fit to Data') +
  scale_x_continuous(breaks = seq(0, length(cases) - 1, by = 5), labels = format(date, "%Y-%m-%d")[seq(1, length(date), by = 5)]) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

