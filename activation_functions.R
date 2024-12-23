prob_target <- function(t_target, t_step, k){
  return(1 - pnorm(t_target, t_step, k*t_step))
}

prob_action <- function(prob_target, g, threshold, a = 10){
  return(g + (1 - g)/(1 + exp(-a*(prob_target - threshold))))
}

prob_check_time <- function(t_target, k, g, threshold, a = 10, dstep = t_target/1000){
  t_steps = seq(0, t_target, dstep)
  
  prob_check_time = numeric(length(t_steps))
  
  for (istep in seq_along(t_steps)){
    prob_action_now = prob_action(prob_target(t_target, t_steps[istep], k), g, threshold, a)
    if (istep == 1){
      prob_check_time[istep] = prob_action_now
    } else {
      prob_no_previous_action = 1 - (sum(prob_action(prob_target(t_target, t_steps[1:(istep - 1)], k), g, threshold, a)) * t_steps[istep])
      prob_check_time[istep] = prob_action_now * prob_no_previous_action
    }
  }
  
  return(prob_check_time)
}

t_steps <- seq(0, 1, 0.01)
g <- 0.1
threshold <- 0.2
k <- 1
a <- 10
t_target <- 1

probability_action <- prob_action(prob_target(1, t_steps, k), g, threshold, a)
plot(t_steps, probability_action)

resp <- prob_check_time(t_target,k, g, threshold, a)

hist(resp)
