

# pars she used:
# beta=33%
# imunity: 11 days
# infectoius period = 6 days
# population = 2000
# initially infected = 1

# beta = 0.33; xi = 0.09090909; gamma = 0.1666667

abm_data = read.csv(file.choose()) |>
  select(time, S = Susceptible, I = Infected.infectious, R = Recovered) |>
  mutate(P = I/(S+I+R)) |>
  select(S, I, R, time, P)


myDat = sampleEpidemic(abm_data, sampleDates = seq(from=0, to=max(abm_data$time), by=1))

# nllikelihood(parms = c(beta = .4, gamma = 1/6, xi = 1/11),
#              obsDat = myDat, abm_data = abm_data)



optim.vals <- optim(par = c(beta = .4, gamma = .01, xi = .05)
                    , nllikelihood
                    # , fixed.params = disease_params()
                    , obsDat = myDat
                    , abm_data = abm_data
                    , control = list(trace = 1, maxit = 1000)
                    , method = "BFGS")

# USING CP
optim.vals <- optim(par = c(c = 1, p = .4, gamma = .01, xi = .05)
                    , nllikelihood
                    # , fixed.params = disease_params()
                    , obsDat = myDat
                    , abm_data = abm_data
                    , control = list(trace = 1, maxit = 1000)
                    , method = "BFGS")


(opar = abs(optim.vals$par))

# Solving the system of equations
out <- ode(y = (abm_data)[1, 1:3] |> unlist(),
           times = seq(from = 1, to = max(abm_data$time), by = 1),
           func = sirs_model_cp,
           parms = opar)

# Convert to data frame for easier handling
out <- as.data.frame(out)

# plotting the abm data
abm_data |>
  select(S, I, R, time) |>
  pivot_longer(-time) |>
  ggplot() + 
  geom_line(aes(x = time, y = value, col = name)) + 
  theme_classic()

# plottiong to confirm
ggplot(abm_data, aes(x = time)) + 
  # geom_point(aes(y = S, col = simulation), cex = .3) + 
  geom_point(aes(y = I), cex = .3) + 
  geom_line(aes(y = I)) + 
  # geom_point(aes(y = R, col = simulation), cex = .3) + 
  
  # geom_line(data = out, aes(x = time, y = S), col = 'black') + 
  geom_line(data = out, aes(x = time, y = I), col = 'red') + 
  # geom_line(data = out, aes(x = time, y = R), col = 'blue') + 
  
  theme_classic()



