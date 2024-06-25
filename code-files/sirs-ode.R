



# Libs --------------------------------------------------------------------

library(pacman)
p_load(dplyr, ggplot2, tidyr, deSolve, stringr)


# Load data from ABM ------------------------------------------------------

abm_data = read.csv('data/p500/sirs-model.csv') |>
  select(time, S = Susceptible, I = Infected.infectious, R = Recovered) |>
  mutate(P = I/(S+I+R)) |>
  select(S, I, R, time, P)



# Functions --------------------------------------------------------------


#' Define the SIRS model with the beta parametrization
#'
#' @param time # a vector of time points at which to evaluate 
#' @param state a vector with numbr of S, I, R at a given point
#' @param parameters a vector of parameters: prob_infect, prob_recover, prob_reinfect, contact_rate
#'
#' @return a list with the dS, dI, dR
sirs_model <- function(time, state, parameters) {
  with(as.list(c(state, parameters)), {
    N <- S + I + R
    
    # Differential equations
    dS <- xi * R - beta * S * I / N
    dI <- beta * S * I / N - gamma * I
    dR <- gamma * I - xi * R
    
    list(c(dS, dI, dR))
  })
}

#' This is the SIRS model with c*P implementation
#'
#' @param time # a vector of time points at which to evaluate 
#' @param state a vector with numbr of S, I, R at a given point
#' @param parameters a vector of parameters: prob_infect, prob_recover, prob_reinfect, contact_rate
#'
#' @return a list with the dS, dI, dR
sirs_model_cp <- function(time, state, parameters) {
  
  with(as.list(c(state, parameters)), {
    N <- S + I + R
    beta = c*p
    
    # Differential equations
    dS <- xi * R - beta * S * I / N
    dI <- beta * S * I / N - gamma * I
    dR <- gamma * I - xi * R
    
    list(c(dS, dI, dR))
  })
}

#' Function for sampling from the ABM world dataset: returns either full or reduced sample
#'
#' @param simDat The data from the ABM 
#'
#' @return a dataset with time, positive cases, sample size, prevalence, lower and upper CI
sampleEpidemic <- function(simDat # Simulated "data" which we treat as real 
                           , sampleDates = seq(from = 1, to = 200, by = 1) # Sample every 3 years 
                           , numSamp = rep(1000, length(sampleDates)) # Number of individuals sampled at each time point
){
  prev_at_sample_times <- simDat |> pull(P)
  samp_size = lag(simDat$S) |> zoo::na.locf(fromLast = T)
  numSamp = samp_size
  
  numPos <- rbinom(length(numSamp), round(numSamp, 0), prev_at_sample_times)
  
  lci <- mapply(function(x,n) binom.test(x,n)$conf.int[1], x = numPos, n = round(numSamp, 0))
  uci <- mapply(function(x,n) binom.test(x,n)$conf.int[2], x = numPos, n = round(numSamp, 0)) 
  
  return(data.frame(time = sampleDates, numPos, numSamp, sampPrev =  numPos/numSamp,
                    lci = lci, uci = uci))
}

# myDat = sampleEpidemic(abm_data)

#' Function for computing the negative log likelihood
#'
#' @param parms the parameters of the SIRS model
#' @param obsDat the observation data i.e. a sample from the ABM data
#' @param abm_data # the actual ABM data
#'
#' @return numeric value: the sum of the negative log-likelihoods
nllikelihood <- function(parms = (parameters), obsDat=myDat, abm_data = abm_data) {
  
  initial_state = (abm_data)[1, 1:3] |> unlist()
  time = seq(from = 0, to = max(abm_data$time), by = 1)
  
  simDat <- ode(y = initial_state, times = time, func = sirs_model_cp, parms = parms) |> # sirs_model
    as.data.frame() |>
    mutate(P = I / (S + I + R))
  
  nlls <- -dbinom(obsDat$numPos, obsDat$numSamp, prob = simDat$P, log = T)
  return(sum(nlls))
}


nllikelihood(parms =  c(c = 1, p = .2, gamma = .4, xi = .5),
             obsDat=myDat, abm_data = abm_data) ## loglikelihood of the true parameters (which we usually never know)



# Fitting the model by MLE ------------------------------------------------

optim.vals <- optim(par = c(c = 3, p = .1, gamma = .4, xi = .5)
                    , nllikelihood
                    # , fixed.params = disease_params()
                    , obsDat = myDat
                    , abm_data = abm_data
                    , control = list(trace = 1, maxit = 1000)
                    , method = "Nelder-Mead")

# THE CHOSEN PARAMETERS
optim.vals$par

# Solving the system of equations: SIRS-ODE using the fitted parameters
out <- ode(y = (abm_data)[1, 1:3] |> unlist(),
           times = 1:200,
           func = sirs_model_cp,
           parms = optim.vals$par)

# Convert to data frame for easier handling
out <- as.data.frame(out)


# plottiong to confirm
ggplot(fsir |> mutate(simulation = factor(simulation)) |> filter(simulation == 1),
       aes(x = time)) + 
  geom_point(aes(y = I, col = simulation), cex = .3) + 
  geom_line(data = out, aes(x = time, y = I), col = 'red') + 
  theme_classic()
