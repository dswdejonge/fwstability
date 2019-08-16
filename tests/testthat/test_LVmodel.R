# Lotka-Volterra model
# https://www.r-bloggers.com/lotka-volterra-model%C2%A0%C2%A0intro/
library(FME)

LotVmod <- function (Time, State, Pars) {
  with(as.list(c(State, Pars)), {
    dx = x*(alpha - beta*y)
    dy = -y*(gamma - delta*x)
    return(list(c(dx, dy)))
  })
}

Pars <- c(alpha = 2, beta = 0.5, gamma = 0.2, delta = 0.6)
State <- c(x = 10, y = 10)
Time <- seq(0, 100, by = 1)

out <- as.data.frame(ode(func = LotVmod, y = State, parms = Pars, times = Time))

library(ggplot2)
ggplot(out) +
  geom_point(aes(x = time, y = x, colour = "blue")) +
  geom_point(aes(x = time, y = y, colour = "red"))

jac <- jacobian.full(y=State, func=LotVmod, parms = Pars)
