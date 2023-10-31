# Generated on 2023-10-27 11:01:50 by gEcon ver. 1.2.1 (2023-01-18)
# http://gecon.r-forge.r-project.org/

# Model name: RW

# info
info__ <- c("RW", "C:/Users/fm007/Documents/GitHub/OptimalRSW/RW.gcn", "2023-10-27 11:01:50", "false")

# index sets
index_sets__ <- list()

# variables
variables__ <- c("etapi",
                 "i",
                 "lambda__RW_1",
                 "lambda__RW_2",
                 "pi",
                 "y",
                 "U")

variables_tex__ <- c("{e\\!t\\!a\\!p\\!i}",
                     "i",
                     "\\lambda^{\\mathrm{RW}^{\\mathrm{1}}}",
                     "\\lambda^{\\mathrm{RW}^{\\mathrm{2}}}",
                     "\\pi",
                     "y",
                     "U")

# shocks
shocks__ <- c("epsilon_pi")

shocks_tex__ <- c("\\epsilon^{\\pi}")

# parameters
parameters__ <- c("beta",
                  "kappa",
                  "phi",
                  "pistar",
                  "sigma",
                  "theta")

parameters_tex__ <- c("\\beta",
                     "\\kappa",
                     "\\phi",
                     "{p\\!i\\!s\\!t\\!a\\!r}",
                     "\\sigma",
                     "\\theta")

# free parameters
parameters_free__ <- c("beta",
                       "kappa",
                       "phi",
                       "pistar",
                       "sigma",
                       "theta")

# free parameters' values
parameters_free_val__ <- c(0.99,
                           0.2465,
                           0.95,
                           0,
                           1,
                           6)

# equations
equations__ <- c("-etapi[] + exp(epsilon_pi[] + phi * log(etapi[-1])) = 0",
                 "-y[-1] + y[] - sigma * (i[-1] - pi[]) = 0",
                 "lambda__RW_2[] + beta * (kappa * E[][lambda__RW_1[1]] - E[][lambda__RW_2[1]]) - kappa * theta^-1 * y[] = 0",
                 "etapi[-1] - pi[-1] + beta * pi[] + kappa * y[-1] = 0",
                 "U[] + 0.5 * (-pistar + pi[])^2 - beta * E[][U[1]] + 0.5 * kappa * theta^-1 * y[]^2 = 0",
                 "pistar - pi[] + beta * lambda__RW_1[] - beta * E[][lambda__RW_1[1]] + sigma * lambda__RW_2[] = 0",
                 "-beta * sigma * E[][lambda__RW_2[1]] = 0")

# calibrating equations
calibr_equations__ <- character(0)

# variables / equations map
vareqmap__ <- sparseMatrix(i = c(1, 2, 2, 2, 3, 3, 3, 4, 4, 4,
                                 5, 5, 5, 6, 6, 6, 7),
                           j = c(1, 2, 5, 6, 3, 4, 6, 1, 5, 6,
                                 5, 6, 7, 3, 4, 5, 4),
                           x = c(3, 1, 2, 3, 4, 6, 2, 1, 3, 1,
                                 2, 2, 6, 6, 2, 2, 4),
                           dims = c(7, 7))

# variables / calibrating equations map
varcalibreqmap__ <- sparseMatrix(i = NULL, j = NULL, dims = c(0, 7))

# calibrated parameters / equations map
calibrpareqmap__ <- sparseMatrix(i = NULL, j = NULL, dims = c(7, 0))

# calibrated parameters / calibrating equations map
calibrparcalibreqmap__ <- sparseMatrix(i = NULL, j = NULL, dims = c(0, 0))

# free parameters / equations map
freepareqmap__ <- sparseMatrix(i = c(1, 2, 3, 3, 3, 4, 4, 5, 5, 5,
                                     5, 6, 6, 6, 7, 7),
                               j = c(3, 5, 1, 2, 6, 1, 2, 1, 2, 4,
                                     6, 1, 4, 5, 1, 5),
                               x = rep(1, 16), dims = c(7, 6))

# free parameters / calibrating equations map
freeparcalibreqmap__ <- sparseMatrix(i = NULL, j = NULL, dims = c(0, 6))

# shocks / equations map
shockeqmap__ <- sparseMatrix(i = c(1),
                             j = c(