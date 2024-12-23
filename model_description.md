# Clock Check Model

## Outline

There is a time estimation module inside the system with linearly scaling noise. This leads to uncertainty in time estimation. This time estimation and uncertainty is noted by the system, the probability of checking the clock is related to the probability of the target time being included in the estimated time plus uncertainty.

With $\theta$ being the parameters $t_{target}$ and $k$:

$$P(target | x, \theta) = 1 - \Phi(\frac{\t_{target} - x}{k x})$$

This probability of the target being passed already is mapped to a probability of acting via a sigmoid function with the parameters $g$, $\delta$, and $a$:

$$P(response | x, \theta) = g + \frac{1 - g}{1 + exp(-a*(1 - \Phi(\frac{\t_{target} - x}{k x}) - \delta))} $$

The likelihood function of a survival model with these probabilities is

$$P(x | response, \theta) = (1-\int_0^x P(response | t, \theta) dt) P(response | x, \theta)$$
