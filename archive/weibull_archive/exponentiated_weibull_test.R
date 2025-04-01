surv_exp_survival <- function(t, k, rho, y){
  surv = 1 - (1 - exp(-(rho*t)^k))^y
}

t <- seq(0, 1, length.out = 100)

plot(t, surv_exp_survival(t, 1, 1, 1))
