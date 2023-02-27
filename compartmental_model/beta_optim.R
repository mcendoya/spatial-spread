library(deSolve)

##### Xf survey data #####

df <- read.table("xf_surveys.txt")

##### Compartmental Model #####

# state_values:
  # S: susceptible individuals
  # IA: asymptomatic infected
  # IS: symptomatic infected
# parameters:
  # beta: transmission rate
  # lmbda: IA transmission rate reduction
  # sigma: transition rate from IA to IS

si_model = function (current_timepoint, state_values, parameters)
{
  S = state_values [1]        
  IA = state_values [2]  
  IS = state_values [3]  
  with ( 
    as.list (parameters),     
    {
      dS = -((beta * lmbda * IA) + (beta * IS)) * S
      dIA = (((beta * lmbda * IA) + (beta * IS)) * S) - (sigma * IA)
      dIS = (sigma * IA)
      
      results = c (dS, dIA, dIS)
      list (results)
    }
  )
}

##### Function to solve ODE and run the model and obtain the sum of the squared error #####
  # beta: transmission rate (1/month)
  # lmbda: IA transmission rate reduction
  # se: asymptomatic period (time from IA to IS) sigma=1/se
  # IA: proportion of initial asymptomatic infected individuals 
  # yr: time for run the model (years)

ss <- function(beta, lmbda, se, IA, yr) {
  beta_value = beta
  lmbda_value = lmbda
  sigma_value = 1/se
  
  parameter_list = c (beta = beta_value, lmbda = lmbda_value, sigma = sigma_value)
  
  initial_values = c (S = 1-IA, IA = IA, IS = 0)
  timepoints = seq (1, yr*12, by=12)
  
  predictions <- ode (initial_values, timepoints, si_model, parameter_list)
  predictions <- as.data.frame(predictions)
  pred_I <- (predictions$IA+predictions$IS)[27:30]
  sum((pred_I - df$proportion[df$stat=="I"])^2)
}

##### "Optimal" beta based on sum of the squared error #####

# Total population
N = 282041
# Initial proportion of infected individuals
initial = 100/N
# beta values to run the model
beta_val <- seq(from = 0.001, to = 0.1, le = 100)
# run the model
ss_val <- sapply(beta_val, ss, lmbda=0.015, se=8, IA=initial, yr=30)
# minimum sum of the squared error
min_ss_val <- min(ss_val)
# Optimal beta
beta_hat <- beta_val[ss_val == min_ss_val]
beta_opt <- beta_hat

plot(beta_val, ss_val, type = "l", lwd = 2,
     xlab = expression(paste("transmission rate ", beta, " (",month^{-1}, ")")),
     ylab = "sum of squares")



