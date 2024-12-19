prob_target <- function(t_target, t_step, k){
  return(1 - pnorm(t_target, t_step, k*t_step))
}

prob_action <- function(prob_target, g, threshold, a = 1){
  return(g + (1-g)*(exp(a*(prob_target - threshold))/1+exp(a*(prob_target - threshold))))
}

prob_action(prob_target(1, seq(0, 1, 0.01), 1), 0.1, 0.3, 0.1)
