# Generated on 2025-04-15 00:18:35 by gEcon ver. 1.2.1 (2023-01-18)
# http://gecon.r-forge.r-project.org/

# Model name: NK_RS

# info
info__ <- c("NK_RS", "C:/Users/fm007/Documents/GitHub/OptimalRSW/paper/Endogenous_Switching/For Larry/NK_RS.gcn", "2025-04-15 00:18:35", "false")

# index sets
index_sets__ <- list()

# variables
variables__ <- c("epsilon_G",
                 "g_1",
                 "g_2",
                 "inflation_gap",
                 "lambda",
                 "mc",
                 "nu_p",
                 "percieved_pi_obj",
                 "pi",
                 "pi_star",
                 "pi_obj",
                 "pH",
                 "pL",
                 "q",
                 "r",
                 "B",
                 "C",
                 "Div",
                 "G",
                 "I",
                 "K_s",
                 "L_s",
                 "Q",
                 "R",
                 "T",
                 "U",
                 "W",
                 "Y",
                 "Y_j",
                 "Y_s",
                 "Z")

variables_tex__ <- c("\\epsilon^{\\mathrm{G}}",
                     "g^{\\mathrm{1}}",
                     "g^{\\mathrm{2}}",
                     "{i\\!n\\!f\\!l\\!a\\!t\\!i\\!o\\!n}^{\\mathrm{gap}}",
                     "\\lambda",
                     "{m\\!c}",
                     "\\nu^{\\mathrm{p}}",
                     "{p\\!e\\!r\\!c\\!i\\!e\\!v\\!e\\!d}^{\\pi^{\\mathrm{obj}}}",
                     "\\pi",
                     "\\pi^{\\star}",
                     "\\pi^{\\mathrm{obj}}",
                     "{p\\!H}",
                     "{p\\!L}",
                     "q",
                     "r",
                     "B",
                     "C",
                     "{D\\!i\\!v}",
                     "G",
                     "I",
                     "K^{\\mathrm{s}}",
                     "L^{\\mathrm{s}}",
                     "Q",
                     "R",
                     "T",
                     "U",
                     "W",
                     "Y",
                     "Y^{\\mathrm{j}}",
                     "Y^{\\mathrm{s}}",
                     "Z")

# shocks
shocks__ <- c("epsilon_Z",
              "eta_p",
              "eta_R",
              "eta_pi",
              "eta_G")

shocks_tex__ <- c("\\epsilon^{\\mathrm{Z}}",
                  "\\eta^{\\mathrm{p}}",
                  "\\eta^{\\mathrm{R}}",
                  "\\eta^{\\pi}",
                  "\\eta^{\\mathrm{G}}")

# parameters
parameters__ <- c("alpha",
                  "beta",
                  "calibr_pi",
                  "delta",
                  "eta",
                  "gamma_p",
                  "kappa",
                  "lambda_p",
                  "mu",
                  "pi_H",
                  "pLss",
                  "r_pi",
                  "r_Y",
                  "rho",
                  "rho_pi_bar",
                  "rho_G",
                  "rho_a",
                  "xi_p",
                  "G_bar")

parameters_tex__ <- c("\\alpha",
                     "\\beta",
                     "{c\\!a\\!l\\!i\\!b\\!r}^{\\pi}",
                     "\\delta",
                     "\\eta",
                     "\\gamma^{\\mathrm{p}}",
                     "\\kappa",
                     "\\lambda^{\\mathrm{p}}",
                     "\\mu",
                     "\\pi^{\\mathrm{H}}",
                     "{p\\!L\\!s\\!s}",
                     "r^{\\pi}",
                     "r^{\\mathrm{Y}}",
                     "\\rho",
                     "\\rho^{\\pi^{\\mathrm{bar}}}",
                     "\\rho^{\\mathrm{G}}",
                     "\\rho^{\\mathrm{a}}",
                     "\\xi^{\\mathrm{p}}",
                     "G^{\\mathrm{bar}}")

# free parameters
parameters_free__ <- c("alpha",
                       "beta",
                       "delta",
                       "eta",
                       "gamma_p",
                       "kappa",
                       "lambda_p",
                       "mu",
                       "pi_H",
                       "r_pi",
                       "r_Y",
                       "rho",
                       "rho_pi_bar",
                       "rho_G",
                       "rho_a",
                       "xi_p")

# free parameters' values
parameters_free_val__ <- c(0.3,
                           0.99,
                           0.025,
                           2,
                           0.469,
                           1,
                           0.5,
                           0.3,
                           1,
                           1.684,
                           0.099,
                           0.961,
                           0.924,
                           0.949,
                           0.823,
                           0.908)

# equations
equations__ <- c("-B[] = 0",
                 "-lambda[] + q[] = 0",
                 "-lambda[] + mu * C[]^(-1 + mu) * (1 - L_s[])^(1 - mu) * (C[]^mu * (1 - L_s[])^(1 - mu))^(-eta) = 0",
                 "-pL[] + (1 + exp(pLss - kappa * log(inflation_gap[])))^-1 = 0",
                 "-q[] + beta * ((1 - delta) * E[][q[1]] + E[][lambda[1] * r[1]]) = 0",
                 "-r[] + alpha * mc[] * Z[] * K_s[-1]^(-1 + alpha) * L_s[]^(1 - alpha) = 0",
                 "-G[] + G_bar * epsilon_G[] = 0",
                 "-Q[] + lambda[]^-1 * q[] = 0",
                 "-W[] + mc[] * Z[] * (1 - alpha) * K_s[-1]^alpha * L_s[]^(-alpha) = 0",
                 "-Y_j[] + Z[] * K_s[-1]^alpha * L_s[]^(1 - alpha) = 0",
                 "Y_j[] - Y_s[] = 0",
                 "Y_s[] - nu_p[] * Y[] = 0",
                 "-Z[] + exp(epsilon_Z[] + rho_a * log(Z[-1])) = 0",
                 "beta * E[][lambda[1] * pi[1]^-1] - lambda[] * R[]^-1 = 0",
                 "lambda[] * W[] + (-1 + mu) * C[]^mu * (1 - L_s[])^(-mu) * (C[]^mu * (1 - L_s[])^(1 - mu))^(-eta) = 0",
                 "-1 + xi_p * (pi[]^-1 * pi[-1]^gamma_p)^(-lambda_p^-1) + (1 - xi_p) * pi_star[]^(-lambda_p^-1) = 0",
                 "1 - pH[] - pL[] = 0",
                 "eta_p[] - g_1[] + g_2[] * (1 + lambda_p) = 0",
                 "eta_G[] - log(epsilon_G[]) + rho_G * log(epsilon_G[-1]) = 0",
                 "-g_1[] + lambda[] * pi_star[] * Y[] + beta * xi_p * pi_star[] * E[][g_1[1] * pi_star[1]^-1 * (pi[1]^-1 * pi[]^gamma_p)^(-lambda_p^-1)] = 0",
                 "-g_2[] + beta * xi_p * E[][g_2[1] * (pi[1]^-1 * pi[]^gamma_p)^(-lambda_p^-1 * (1 + lambda_p))] + lambda[] * mc[] * Y[] = 0",
                 "-nu_p[] + (1 - xi_p) * pi_star[]^(-lambda_p^-1 * (1 + lambda_p)) + xi_p * nu_p[-1] * (pi[]^-1 * pi[-1]^gamma_p)^(-lambda_p^-1 * (1 + lambda_p)) = 0",
                 "I[] - K_s[] + K_s[-1] * (1 - delta) = 0",
                 "U[] - beta * E[][U[1]] - (1 - eta)^-1 * (C[]^mu * (1 - L_s[])^(1 - mu))^(1 - eta) = 0",
                 "-log(inflation_gap[]) - log(percieved_pi_obj[]) + log(pi[]) = 0",
                 "-log(percieved_pi_obj[]) + pH[] * log(pi_H) + pL[] * log(pi[]) = 0",
                 "eta_pi[] - log(pi_obj[]) + rho_pi_bar * log(pi_obj[-1]) + log(percieved_pi_obj[]) * (1 - rho_pi_bar) = 0",
                 "-Div[] + Y[] - K_s[-1] * r[] - L_s[] * W[] = 0",
                 "-G[] + T[] - B[-1] * pi[]^-1 + B[] * R[]^-1 = 0",
                 "-calibr_pi + eta_R[] - log(R[ss]^-1 * R[]) + rho * log(R[ss]^-1 * R[-1]) + (1 - rho) * (log(pi_obj[]) + r_pi * (-log(pi_obj[]) + log(pi[ss]^-1 * pi[-1])) + r_Y * log(Y[ss]^-1 * Y[])) = 0",
                 "-C[] + Div[] - I[] - T[] + B[-1] * pi[]^-1 + K_s[-1] * r[] - B[] * R[]^-1 + L_s[] * W[] = 0")

# calibrating equations
calibr_equations__ <- c("-0.18 + G[ss] * Y[ss]^-1 = 0",
                        "-0.05 + pL[ss] = 0",
                        "pi[ss] - pi_obj[ss] = 0")

# variables / equations map
vareqmap__ <- sparseMatrix(i = c(1, 2, 2, 3, 3, 3, 4, 4, 5, 5,
                                 5, 6, 6, 6, 6, 6, 7, 7, 8, 8,
                                 8, 9, 9, 9, 9, 9, 10, 10, 10, 10,
                                 11, 11, 12, 12, 12, 13, 14, 14, 14, 15,
                                 15, 15, 15, 16, 16, 17, 17, 18, 18, 19,
                                 20, 20, 20, 20, 20, 21, 21, 21, 21, 21,
                                 22, 22, 22, 23, 23, 24, 24, 24, 25, 25,
                                 25, 26, 26, 26, 26, 27, 27, 28, 28, 28,
                                 28, 28, 28, 29, 29, 29, 29, 29, 30, 30,
                                 30, 30, 31, 31, 31, 31, 31, 31, 31, 31,
                                 31, 31, 31),
                           j = c(16, 5, 14, 5, 17, 22, 4, 13, 5, 14,
                                 15, 6, 15, 21, 22, 31, 1, 19, 5, 14,
                                 23, 6, 21, 22, 27, 31, 21, 22, 29, 31,
                                 29, 30, 7, 28, 30, 31, 5, 9, 24, 5,
                                 17, 22, 27, 9, 10, 12, 13, 2, 3, 1,
                                 2, 5, 9, 10, 28, 3, 5, 6, 9, 28,
                                 7, 9, 10, 20, 21, 17, 22, 26, 4, 8,
                                 9, 8, 9, 12, 13, 8, 11, 15, 18, 21,
                                 22, 27, 28, 9, 16, 19, 24, 25, 9, 11,
                                 24, 28, 9, 15, 16, 17, 18, 20, 21, 22,
                                 24, 25, 27),
                           x = c(2, 2, 2, 2, 2, 2, 2, 2, 4, 6,
                                 4, 2, 2, 1, 2, 2, 2, 2, 2, 2,
                                 2, 2, 1, 2, 2, 2, 1, 2, 2, 2,
                                 2, 2, 2, 2, 2, 3, 6, 4, 2, 2,
                                 2, 2, 2, 3, 2, 2, 2, 2, 2, 3,
                                 6, 2, 6, 6, 2, 6, 2, 2, 6, 2,
                                 3, 3, 2, 2, 3, 2, 2, 6, 2, 2,
                                 2, 2, 2, 2, 2, 2, 3, 2, 2, 1,
                                 2, 2, 2, 2, 3, 2, 2, 2, 9, 2,
                                 11, 10, 2, 2, 3, 2, 2, 2, 1, 2,
                                 2, 2, 2),
                           dims = c(31, 31))

# variables / calibrating equations map
varcalibreqmap__ <- sparseMatrix(i = c(1, 1, 2, 3, 3),
                                 j = c(19, 28, 13, 9, 11),
                                 x = rep(1, 5), dims = c(3, 31))

# calibrated parameters / equations map
calibrpareqmap__ <- sparseMatrix(i = c(4, 7, 30),
                                 j = c(2, 3, 1),
                                 x = rep(1, 3), dims = c(31, 3))

# calibrated parameters / calibrating equations map
calibrparcalibreqmap__ <- sparseMatrix(i = NULL, j = NULL, dims = c(3, 3))

# free parameters / equations map
freepareqmap__ <- sparseMatrix(i = c(3, 3, 4, 5, 5, 6, 9, 10, 13, 14,
                                     15, 15, 16, 16, 16, 18, 19, 20, 20, 20,
                                     20, 21, 21, 21, 21, 22, 22, 22, 23, 24,
                                     24, 24, 26, 27, 30, 30, 30),
                               j = c(4, 8, 6, 2, 3, 1, 1, 1, 15, 2,
                                     4, 8, 5, 7, 16, 7, 14, 2, 5, 7,
                                     16, 2, 5, 7, 16, 5, 7, 16, 3, 2,
                                     4, 8, 9, 13, 10, 11, 12),
                               x = rep(1, 37), dims = c(31, 16))

# free parameters / calibrating equations map
freeparcalibreqmap__ <- sparseMatrix(i = NULL, j = NULL, dims = c(3, 16))

# shocks / equations map
shockeqmap__ <- sparseMatrix(i = c(13, 18, 19, 27, 30),
                             j = c(