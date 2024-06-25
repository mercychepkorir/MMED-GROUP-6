


# Libs --------------------------------------------------------------------

library(pacman)
p_load(dplyr, ggplot2, tidyr, deSolve, stringr)


# Load data from ABM ------------------------------------------------------

abm_data = read.csv('data/p500/sirs-model.csv') |>
  select(time, S = Susceptible, I = Infected.infectious, R = Recovered) |>
  mutate(P = I/(S+I+R)) |>
  select(S, I, R, time, P)



# Fitting an SIR model by MLE ---------------------------------------------

#' chain binomial SIR function
#'
#' @param time gives the time we want our simulation to run from
#' @param state gives the state values for S, I, R
#' @param parameters include: prob_infection, prob_recovery, prob_reinfection, contact_rate
#'
#' @return list with dS, dI, dR
sirs_chainbinomial <- function(time, state, parameters) {
  
  with(as.list(c(state, parameters)), {
    N <- S + I + R
    
    prob_infect_full = 1 - (1 - c*(I/N)*prob_infect)**c
    
    # Differential equations
    infected_individuals = rbinom(size = last(S), n = 1, prob = prob_infect_full)
    recovered_individuals = rbinom(size = last(I), n = 1, prob = prob_recover)
    reinfected_individuals = rbinom(size = last(R), n = 1, prob = prob_reinfect)
    
    dS = reinfected_individuals - infected_individuals
    dI = infected_individuals - recovered_individuals
    dR = recovered_individuals - reinfected_individuals
    
    list(c(dS, dI, dR))
  })
  
}



#' The function for generating data given time points
#'
#' @param time gives the time we want our simulation to run from
#' @param state gives the state values for S, I, R
#' @param parameters include: prob_infection, prob_recovery, prob_reinfection, contact_rate
#'
#' @return dataset with time series of S, I, R
solve_chainbinomial <- function(time, initial_state, parameters) {
  states = as.list(initial_state)
  N <- last(states$S) + last(states$I) + last(states$R)
  
  
  
  prob_infect = 1 - (1 - parameters[4]*(last(states$I)/N)*parameters[1])**parameters[4]
  prob_recover = parameters[2]
  prob_reinfect = parameters[3]
  
  for (i in time) {
    
    # Differential equations
    infected_individuals = rbinom(size = last(states$S), n = 1, prob = prob_infect)
    recovered_individuals = rbinom(size = last(states$I), n = 1, prob = prob_recover)
    reinfected_individuals = rbinom(size = last(states$R), n = 1, prob = prob_reinfect)
    
    states$S = c(states$S, last(states$S) - infected_individuals + reinfected_individuals)
    states$I = c(states$I, last(states$I) + infected_individuals - recovered_individuals)
    states$R = c(states$R, last(states$R) + recovered_individuals - reinfected_individuals)
  }
  simDat <- list2DF(states) |>
    as.data.frame() |>
    mutate(P = I / (S + I + R))
  return(simDat)
}

# trial
# sirs_chainbinomial(time = 1:200, state = c(S=1000, I = 90, R = 30), parameters = parms)


#' Function for sampling from the ABM world dataset: returns either full or reduced sample
#'
#' @param simDat The data from the ABM 
#'
#' @return a dataset with time, positive cases, sample size, prevalence, lower and upper CI
sampleEpidemic <- function(simDat) {
  
  prev_at_sample_times <- simDat |> pull(P)
  samp_size = lag(simDat$S) |> zoo::na.locf(fromLast = T)
  numSamp = samp_size
  
  numPos <- rbinom(length(numSamp), round(numSamp, 0), prev_at_sample_times)
  
  lci <- mapply(function(x,n) binom.test(x,n)$conf.int[1], x = numPos, n = round(numSamp, 0))
  uci <- mapply(function(x,n) binom.test(x,n)$conf.int[2], x = numPos, n = round(numSamp, 0)) 
  
  return(data.frame(time = simDat |> pull(time),
                    numPos, numSamp, sampPrev =  numPos/numSamp,
                    lci = lci, uci = uci))
}

# sample code
# myDat = sampleEpidemic(abm_data)



#' Function for computing the negative log likelihood
#'
#' @param parms the parameters of the SIRS model
#' @param obsDat the observation data i.e. a sample from the ABM data
#' @param abm_data # the actual ABM data
#'
#' @return numeric value: the sum of the negative log-likelihoods
nllikelihood <- function(parms, obsDat=myDat, abm_data = abm_data) {
  
  initial_state = (abm_data)[1, 1:3] |> unlist()
  time = seq(from = 1, to = max(abm_data$time), by = 1)
  
  simDat = solve_chainbinomial(time = time, 
                               initial_state = initial_state,
                               parameters = parms) |>
    mutate(P = I / (S + I + R))
  
  nlls <- -dbinom(obsDat$numPos, round(obsDat$numSamp, 0), prob = simDat$P, log = T)
  return(sum(nlls))
}

nllikelihood(parms = c(prob_infect = 0.3, prob_recover = 0.5, prob_reinfect = 0.4, c=2),
             obsDat=myDat,
             abm_data = abm_data) ## loglikelihood of the true parameters (which we usually never know)


# Fitting the model by MLE ------------------------------------------------

optim.vals <- optim(par = c(prob_infect = 0.2, prob_recover = .3, prob_reinfect = 0.1, c=2)
                    , nllikelihood
                    # , fixed.params = disease_params()
                    , obsDat = myDat
                    , abm_data = abm_data
                    , control = list(trace = 1, maxit = 1000)
                    , method = "BFGS")

# THE CHOSEN PARAMETERS
optim.vals$par


# plotting to confirm
out = solve_chainbinomial(time = time, initial_state = initial_state,
                          parameters = optim.vals$par) |> # optim.vals$par
  as.data.frame() |>
  mutate(time = row_number())

ggplot(abm_data, aes(x = time)) + 
  geom_point(aes(y = I), cex = .3, show.legend = F) + 
  geom_line(data = out, aes(x = time, y = I), col = 'red') + 
  theme_classic()





