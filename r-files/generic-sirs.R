

for (i in 1:9) {
  
  message('Running simulation: ', i)
  
  num_agents = 1000
  max_time = 200
  
  
  # generating initial pop dist
  popmat = matrix(NA, nrow = num_agents, ncol = max_time)
  initially_infected = which(rbinom(n = 1000, size = 1, prob = pars[1]) == 1)
  popmat[initially_infected, 1] = 'I'
  popmat[-initially_infected, 1] = 'S'
  initially_recovered = sample(which(popmat[, 1] == 'S'), 30)
  popmat[initially_recovered, 1] = 'R'
  
  for (t in 2:max_time) {
    
    num_susceptible = length(which(popmat[, t - 1] == 'S'))
    num_infected = length(which(popmat[, t - 1] == 'I'))
    num_recovered = length(which(popmat[, t - 1] == 'R'))
    
    # prob_infect = (rbinom( n = num_agents, size = num_susceptible, prob = rbeta(num_agents, 1, 8))) / num_susceptible
    # prob_recover = (rbinom(n = num_agents, size = num_infected, prob = rbeta(num_agents, 2, 6))) / num_infected
    # prob_susceptible = (rbinom(n = num_agents, size = num_recovered, prob = rbeta(num_agents, 3, 100))) / num_recovered

    prob_infect = (rbinom( n = num_agents, size = num_susceptible, prob = .3)) / num_susceptible
    prob_recover = (rbinom(n = num_agents, size = num_infected, prob = .5)) / num_infected
    prob_susceptible = (rbinom(n = num_agents, size = num_recovered, prob = .4)) / num_recovered

        
    popmat[, t] <- case_when(
      popmat[, t - 1] == 'S' ~ ifelse(prob_infect > runif(num_agents), 'I', 'S') ,
      popmat[, t - 1] == 'I' ~ ifelse(prob_recover > runif(num_agents), 'R', 'I'),
      popmat[, t - 1] == 'R' ~ ifelse(prob_susceptible > runif(num_agents), 'S', 'R'),
      TRUE ~  NA
    )
    
  }
  
  sir_data = popmat |>
    data.frame() |>
    apply(MARGIN = 2, FUN = \(x) paste((length(which(x == 'S'))), ' | ', (length(which(x == 'I'))), ' | ', (length(which(x == 'R'))))) |>
    data.frame() |>
    setNames('SIR') |>
    separate_wider_delim(col = 'SIR',
                         names = c('S', 'I', 'R'),
                         delim = ' | ') |>
    mutate(across(.cols = everything(), .fns = as.numeric),
           time = row_number())
  
  if (i == 1) fsir = sir_data |> mutate(simulation = 1)
  else fsir = rbind(fsir, sir_data |> mutate(simulation = i))
  
  
}

fsir = fsir |>
  mutate(P = I / (S + I + R))

# plottiong to confirm
ggplot(fsir |> mutate(simulation = factor(simulation)),
       aes(x = time)) + 
  geom_line(aes(y = S, col = simulation)) + 
  geom_line(aes(y = I, col = simulation)) + 
  geom_line(aes(y = R, col = simulation)) + 
  theme_classic()




# Fitting an SIR model by MLE ---------------------------------------------


# WE USE THE FSIR DATA

# Load necessary package
library(deSolve)

# Define the SIRS model
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


# once you get your data, do an I/N
abm_data = fsir |>
  filter(simulation == 1)

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

myDat = sampleEpidemic(abm_data)

## Return the negative log of likelihood by convention
nllikelihood <- function(parms = (parameters), obsDat=myDat, abm_data = abm_data) {
  
  initial_state = (abm_data)[1, 1:3] |> unlist()
  time = seq(from = 0, to = max(abm_data$time), by = 1)
  
  simDat <- ode(y = initial_state, times = time, func = sirs_model_cp, parms = parms) |> # sirs_model
    as.data.frame() |>
    mutate(P = I / (S + I + R))
  
  ## What are the rows from our simulation at which we have observed data?
  # matchedTimes <- simDat$time %in% obsDat$time
  nlls <- -dbinom(obsDat$numPos, obsDat$numSamp, prob = simDat$P, log = T)
  return(sum(nlls))
}


nllikelihood(parms =  c(c = 1, p = .2, gamma = .4, xi = .5),
             obsDat=myDat, abm_data = abm_data) ## loglikelihood of the true parameters (which we usually never know)



optim.vals <- optim(par = c(c = 3, p = .1, gamma = .4, xi = .5)
                    , nllikelihood
                    # , fixed.params = disease_params()
                    , obsDat = myDat
                    , abm_data = abm_data
                    , control = list(trace = 1, maxit = 1000)
                    , method = "Nelder-Mead")


optim.vals$par




# Solving the system of equations
out <- ode(y = (abm_data)[1, 1:3] |> unlist(),
           times = 1:200,
           func = sirs_model_cp,
           parms = optim.vals$par)

# Convert to data frame for easier handling
out <- as.data.frame(out)




# plottiong to confirm
ggplot(fsir |> mutate(simulation = factor(simulation)) |> filter(simulation == 1),
       aes(x = time)) + 
  # geom_point(aes(y = S, col = simulation), cex = .3) + 
  geom_point(aes(y = I, col = simulation), cex = .3) + 
  # geom_point(aes(y = R, col = simulation), cex = .3) + 
  
  # geom_line(data = out, aes(x = time, y = S), col = 'black') + 
  geom_line(data = out, aes(x = time, y = I), col = 'red') + 
  # geom_line(data = out, aes(x = time, y = R), col = 'blue') + 
  
  theme_classic()


abm_runs = 9
plts = list()
opt_vals = list()

for (i in 1:abm_runs) {
  
  
  optim.vals <- optim(par = c(beta = .3, gamma = .4, xi = .5)
                      , nllikelihood
                      # , fixed.params = disease_params()
                      , obsDat = sampleEpidemic(fsir |> filter(simulation == i))
                      , abm_data = (fsir |> filter(simulation == i))
                      , control = list(trace = 1, maxit = 1000)
                      , method = "Nelder-Mead")
  
  opt_vals[[i]] = optim.vals$par
  
  
  # Solving the system of equations
  out <- ode(y = (fsir |> filter(simulation == i))[1, 1:3] |> unlist(),
             times = 1:200,
             func = sirs_model,
             parms = optim.vals$par)
  
  # Convert to data frame for easier handling
  out <- as.data.frame(out)
  

  
  # plottiong to confirm
  plts[[i]] = ggplot(fsir |> mutate(simulation = factor(simulation)) |> filter(simulation == i),
         aes(x = time)) + 
    # geom_point(aes(y = S, col = simulation), cex = .3) + 
    geom_point(aes(y = I, col = simulation), cex = .3) + 
    # geom_point(aes(y = R, col = simulation), cex = .3) + 
    
    # geom_line(data = out, aes(x = time, y = S), col = 'black') + 
    geom_line(data = out, aes(x = time, y = I), col = 'red') + 
    # geom_line(data = out, aes(x = time, y = R), col = 'blue') + 
    
    labs(title = paste('Simulation: ', i)) + 
    
    theme_classic()
  
}

do.call(gridExtra::grid.arrange, plts)

as.data.frame(
  list2DF(opt_vals) |> t()
) |>
  setNames(c('beta', 'gamma', 'xi')) |>
  mutate(beta = beta / 10, gamma = 1 / gamma) |>
  pivot_longer(everything()) |>
  ggplot() + 
  geom_boxplot(aes(x = name, y = value)) + 
  geom_point(data = data.frame(name = c('beta', 'gamma', 'xi'),
                      value = c(.3, .5, .4)),
             aes(x = name, y = value, group = name),
             cex = 3, col = 'red')



