library(deSolve)
library(dplyr)
library(spaMM)

##### Xf survey data #####

data <- read.table("xf_surveys.txt")

df <- data %>% 
  mutate(stat = case_when(res==0 ~ "S",
                          .default = "I")) %>% 
  select(!res) %>% 
  group_by(year, stat) %>% 
  summarise(n = length(stat)) %>% 
  mutate(n_cum = n)

n_cum <- cumsum(df$n[df$stat=="I"])
df$n_cum[df$stat=="I"] <- n_cum

df <- df %>% 
  group_by(year) %>%
  mutate(N = sum(n_cum),
         proportion = n_cum/N)

###########

# Mean cumulative spatial correlation (Matérn correlation function) for susceptibles

df_Cij <- data.frame(yr = rep(c(2018,2019,2020,2021), each=4),
                     Range = rep(c(250,500,750,1000),4),
                     mean_Cor = NA)

r <- c(250,500,750,1000)
yr <- c(2018,2019,2020,2021)

for(y in yr){
  df_yr <- data[data$year==y,]
  locs_S <- df_yr[df_yr$res==0, 1:2]
  locs_I <- df_yr[df_yr$res==1, 1:2]
  for(rg in r){
    d_max <- rg +(1.5*rg)
    cor_mat <- matrix(NA, nrow = nrow(locs_S), ncol = nrow(locs_I))
    for(i in 1:nrow(locs_S)){
      for(j in 1:nrow(locs_I)){
        d <- sqrt((locs_S[i,1]-locs_I[j,1])^2 +(locs_S[i,2]-locs_I[j,2])^2)
        if(d>d_max){
          cor_mat[i,j] <- 0
        }else{

        }
        cor_mat[i,j] <- MaternCorr(d, rho = sqrt(8)/rg, nu = 1)
      }
    }
    df_Cij$mean_Cor[df_Cij$yr==y & df_Cij$Range==rg] <- mean(apply(cor_mat, 1, sum))
  }
}

dfC <- df_Cij %>% 
  group_by(Range) %>% 
  summarise(meanCor = mean(mean_Cor))

##### Compartmental Model #####

# state_values:
# S: susceptible individuals
# IA: asymptomatic infected
# IS: symptomatic infected
# parameters:
# beta: transmission rate
# lmbda: IA transmission rate reduction
# sigma: transition rate from IA to IS
# Cij: spatial correlation

siC_model = function (current_timepoint, state_values, parameters)
{
  S = state_values [1]        
  IA = state_values [2]  
  IS = state_values [3]  
  with ( 
    as.list (parameters),     
    {
      dS = -((beta * lmbda * Cij * IA) + (beta * Cij * IS)) * S
      dIA = (((beta * lmbda * Cij * IA) + (beta * Cij * IS)) * S) - (sigma * IA)
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
# Cij: spatial correlation
# yr: time for run the model (years)

ss_C <- function(beta, lmbda, se, Cij, IA, yr) {
  beta_value = beta
  lmbda_value = lmbda
  sigma_value = 1/se
  C_mat = Cij
  
  parameter_list = c (beta = beta_value, lmbda = lmbda_value, sigma = sigma_value, Cij = C_mat)
  
  initial_values = c (S = 1-IA, IA = IA, IS = 0)
  timepoints = seq (1, yr*12, by=12)
  
  predictions <- ode (initial_values, timepoints, siC_model, parameter_list)
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
# Values for the range parameter of the Matérn correlation for run the model
rg <- dfC$Range

# run the model
dfC$beta_opt <- NA
for(i in 1:length(dfC$meanCor)){
  ss_val <- sapply(beta_val, ss_C, lmbda=0.015, se=8, Cij = round(dfC$meanCor[i],2), IA=initial, yr=30)
  min_ss_val <- min(ss_val)
  beta_hat <- beta_val[ss_val == min_ss_val]
  dfC$beta_opt[dfC$Range==rg[i]] <- beta_hat
}

# Plot results
plot(beta_val, type="n", xlab = expression(paste("transmission rate ", beta, " (",month^{-1}, ")")),
     ylab = "sum of squares", ylim=c(0,1.5), xlim=c(0,0.1))
for(i in 1:length(ss_val)){
  lines(beta_val, ss_val[[i]], lwd = 2, col=i)
}
legend("topright", c(paste0("r = ", dfC$Range)), col=c(1:length(dfC$Range)), lwd=2,
       title = "Range parameter (m)")
