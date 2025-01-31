prob_target <- function(t_target, t_step, k){
  return(1 - pnorm(t_target, t_step, k*t_step))
}

prob_action <- function(prob_target, c, g, a = 10){
  return(g + (1 - g)/(1 + exp(-a*(log(prob_target / (1 - prob_target)) + c))))
}

plot_action_prob <- function(t_target, nsteps, g, c, k, a){
  t_steps <- seq(0, t_target, t_target/nsteps)
  probability_action <- prob_action(prob_target(t_target, t_steps, k), c * t_target, g, a)
  # prob_no_action_yet <- prob_no_action_until_time(t_target, k, g, threshold, a, nsteps)
  plot(t_steps, probability_action)
}

hazard_function <- function(t, t_target, k, g, a, c){
  z = (1 - pnorm(t_target, t, k*t)) / pnorm(t_target, t, k*t)
  action_prob = g + (1 - g) / (1 + exp(-a*c*t_target) * z^-a)
  
  return(action_prob)
}

t_target = 1
t_steps = seq(0, 1, 0.01)
k = 1
g = 0.1
a = 1
c = 4

plot(t_steps, hazard_function(t_steps, t_target, k, g, a, c))

survival_function(t, t*1.2, 1, 0.01, 1, 1)

plot(t_steps, survival_function(t_steps, 1, 1, 0.1, 1, 1) * hazard_function(t_steps, t_target, 1, 0.01, 1, 1))

plot_action_prob(300, 100, 0.1, 0.001, 1, 2)

t_target <- 100
plot(seq(0, t_target, t_target/nsteps), prob_target(t_target, seq(0, t_target, t_target/nsteps), k))

t_target <- 200
t_steps <- seq(0, t_target, t_target/nsteps)
nsteps <- 100
g <- 0
c <- 1
k <- 1
a <- 1

plot(t_steps, prob_target(t_target, t_steps, k))

probability_action <- prob_action(prob_target(t_target, t_steps, k), c, g, a)
# prob_no_action_yet <- prob_no_action_until_time(t_target, k, g, threshold, a, nsteps)
plot(t_steps, probability_action)


prob_no_action_until_time <- function(t_target, k, g, c, a = 10, nsteps = 100){
  action_prob = 0
  
  steps = seq(0, t_target, t_target/nsteps)
  
  step_action_prob = cumsum(prob_action(prob_target(t_target, steps, k), g, c, a))
  
  prob_no_action = 1 - c(0, step_action_prob[-length(step_action_prob)])
  return(prob_no_action)
}

probability_no_action <- prob_no_action_until_time(t_target, k, -3, 0, a)
# prob_no_action_yet <- prob_no_action_until_time(t_target, k, g, threshold, a, nsteps)
plot(t_steps, probability_no_action)


prob_check_time <- function(t_target, k, g, threshold, a = 10, nsteps = 100){
  t_steps = seq(0, t_target, t_target/nsteps)
  
  prob_check_time = prob_no_action_until_time(t_target, k, g, threshold, a, nsteps) * prob_action(prob_target(t_target, t_steps, k), g, threshold, a)
  
  return(prob_check_time)
}

t_target <- 1
nsteps <- 100
t_steps <- seq(0, t_target, t_target/nsteps)
g <- 0.01
c <- 0
k <- 1
a <- 1

probability_action <- prob_action(prob_target(t_target, t_steps, k), c, g, a)
prob_no_action_yet <- prob_no_action_until_time(t_target, k, g, c, a, nsteps)
plot(t_steps, probability_action)
plot(t_steps, prob_no_action_yet)

resp <- prob_check_time(t_target,k, g, threshold, a, nsteps)

plot(t_steps[which(resp > 0)], resp[which(resp > 0)])
