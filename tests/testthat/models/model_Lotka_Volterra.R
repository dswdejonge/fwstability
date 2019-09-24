# Lotka-Volterra model
# https://www.r-bloggers.com/lotka-volterra-model%C2%A0%C2%A0intro/
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

model <- list(
  type = "ODE",
  func = LotVmod,
  y = State,
  parms = Pars
)
answer <- max(Re(eigen(rootSolve::jacobian.full(y=State, func=LotVmod, parms = Pars)
)$value))
