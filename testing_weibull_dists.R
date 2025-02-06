standard_weibull_pdf <- function(x, k, lambda){
  return((k/lambda)*(x/lambda)^(k-1)*exp(-(x/lambda))^k)
}

standard_weibull_hazard <- function(x, k, lambda){
  return((k/lambda)*(x/lambda)^(k-1))
}

prop_hazard_weibull_pdf <- function(x, k, b){
  lambda = b^(-1/k)
  
  standard_weibull_pdf(x, k, lambda)
}

prop_hazard_weibull_hazard <- function(x, k, b){
  lambda = b^(-1/k)
  
  standard_weibull_hazard(x, k, lambda)
}

three_pl_weibull_pdf <- function(x, a, b, c){
  term1 = log(a) + log(b + c * x) + (b - 1) * log(x) + c * x
  term2 = -a * x^b * exp(c * x)
  
  return(term1 * term2)
}


t_steps <- seq(0, 1, length.out = 100)
k = 4
lambda = 0.6
b = 200
# plot(t_steps, standard_weibull_pdf(t_steps, k,lambda))
# plot(t_steps, standard_weibull_hazard(t_steps, k, lambda))
# plot(t_steps, prop_hazard_weibull_hazard(t_steps, k, b))
plot(t_steps, prop_hazard_weibull_pdf(t_steps, k, b))

plot(t_steps, three_pl_weibull_pdf(t_steps, 2, 2, 0.1))

