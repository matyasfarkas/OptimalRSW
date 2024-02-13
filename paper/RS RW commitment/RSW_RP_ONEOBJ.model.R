# Generated on 2024-02-13 18:09:18 by gEcon ver. 1.2.1 (2023-01-18)
# http://gecon.r-forge.r-project.org/

# Model name: RSW_RP_ONEOBJ

# info
info__ <- c("RSW_RP_ONEOBJ", "C:/Users/fm007/Documents/GitHub/OptimalRSW/paper/RS RW commitment/RSW_RP_ONEOBJ.gcn", "2024-02-13 18:09:18", "false")

# index sets
index_sets__ <- list()

# variables
variables__ <- c("etapi",
                 "lambda__OPTIMALMP_1",
                 "lambda__OPTIMALMP_2",
                 "piH",
                 "piL",
                 "yH",
                 "yL",
                 "U")

variables_tex__ <- c("{e\\!t\\!a\\!p\\!i}",
                     "\\lambda^{\\mathrm{OPTIMALMP}^{\\mathrm{1}}}",
                     "\\lambda^{\\mathrm{OPTIMALMP}^{\\mathrm{2}}}",
                     "{p\\!i\\!H}",
                     "{p\\!i\\!L}",
                     "{y\\!H}",
                     "{y\\!L}",
                     "U")

# shocks
shocks__ <- c("epsilon_pi")

shocks_tex__ <- c("\\epsilon^{\\pi}")

# parameters
parameters__ <- c("beta",
                  "kappa",
                  "lambda",
                  "phi",
                  "pitH",
                  "pitCB",
                  "pitL",
                  "pH",
                  "pL",
                  "sigma",
                  "theta")

parameters_tex__ <- c("\\beta",
                     "\\kappa",
                     "\\lambda",
                     "\\phi",
                     "{p\\!i\\!t\\!H}",
                     "{p\\!i\\!t\\!C\\!B}",
                     "{p\\!i\\!t\\!L}",
                     "{p\\!H}",
                     "{p\\!L}",
                     "\\sigma",
                     "\\theta")

# free parameters
parameters_free__ <- c("beta",
                       "kappa",
                       "lambda",
                       "phi",
                       "pitH",
                       "pitCB",
                       "pitL",
                       "pH",
                       "pL",
                       "sigma",
                       "theta")

# free parameters' values
parameters_free_val__ <- c(0.99,
                           0.2465,
                           0.04106,
                           0.95,
                           2,
                           2,
                           4,
                           0.99,
                           0.99,
                           1,
                           6)

# equations
equations__ <- c("-etapi[] + exp(epsilon_pi[] + phi * log(etapi[-1])) = 0",
                 "-0.5 * lambda * yH[] + beta * kappa * E[][lambda__OPTIMALMP_1[1]] = 0",
                 "-0.5 * lambda * yL[] + beta * kappa * E[][lambda__OPTIMALMP_2[1]] = 0",
                 "-piH[-1] + log(etapi[-1]) + beta * piH[] + kappa * yH[-1] + beta * (1 - pH) * (-piH[] + piL[]) = 0",
                 "-piL[-1] + log(etapi[-1]) + kappa * yL[-1] + beta * pL * piL[] + beta * (1 - pL) * (piH[] - piL[]) = 0",
                 "-0.5 * pitH + 0.5 * pitCB - 0.5 * piH[] - beta * E[][lambda__OPTIMALMP_1[1]] + lambda__OPTIMALMP_1[] * (beta - beta * (1 - pH)) + beta * lambda__OPTIMALMP_2[] * (1 - pL) = 0",
                 "0.5 * pitCB - 0.5 * pitL - 0.5 * piL[] - beta * E[][lambda__OPTIMALMP_2[1]] + lambda__OPTIMALMP_2[] * (beta * pL - beta * (1 - pL)) + beta * lambda__OPTIMALMP_1[] * (1 - pH) = 0",
                 "U[] + 0.25 * (pitH - pitCB + piH[])^2 + 0.25 * (-pitCB + pitL + piL[])^2 - beta * E[][U[1]] + 0.25 * lambda * yH[]^2 + 0.25 * lambda * yL[]^2 = 0")

# calibrating equations
calibr_equations__ <- character(0)

# variables / equations map
vareqmap__ <- sparseMatrix(i = c(1, 2, 2, 3, 3, 4, 4, 4, 4, 5,
                                 5, 5, 5, 6, 6, 6, 7, 7, 7, 8,
                                 8, 8, 8, 8),
                           j = c(1, 2, 6, 3, 7, 1, 4, 5, 6, 1,
                                 4, 5, 7, 2, 3, 4, 2, 3, 5, 4,
                                 5, 6, 7, 8),
                           x = c(3, 4, 2, 4, 2, 1, 3, 2, 1, 1,
                                 2, 3, 1, 6, 2, 2, 2, 6, 2, 2,
                                 2, 2, 2, 6),
                           dims = c(8, 8))

# variables / calibrating equations map
varcalibreqmap__ <- sparseMatrix(i = NULL, j = NULL, dims = c(0, 8))

# calibrated parameters / equations map
calibrpareqmap__ <- sparseMatrix(i = NULL, j = NULL, dims = c(8, 0))

# calibrated parameters / calibrating equations map
calibrparcalibreqmap__ <- sparseMatrix(i = NULL, j = NULL, dims = c(0, 0))

# free parameters / equations map
freepareqmap__ <- sparseMatrix(i = c(1, 2, 2, 2, 3, 3, 3, 4, 4, 4,
                                     5, 5, 5, 6, 6, 6, 6, 6, 7, 7,
                                     7, 7, 7, 8, 8, 8, 8, 8),
                               j = c(