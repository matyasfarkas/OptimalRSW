# Generated on 2025-04-15 18:16:55 by gEcon ver. 1.2.1 (2023-01-18)
# http://gecon.r-forge.r-project.org/

# Model name: NK_RS

# info
info__ <- c("NK_RS", "C:/Users/fm007/Documents/GitHub/OptimalRSW/paper/Endogenous_Switching/For Larry/NK_RS.gcn", "2025-04-15 18:16:55", "false")

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
                  "tau",
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
                     "\\tau",
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
                       "tau",
                       "xi_p")

# free parameters' values
parameters_free_val__ <- c(0.3,
                           0.99,
                           0.025,
                           2,
                           0.469,
                           10,
                           0.5,
                           0.3,
                           1,
                           1.684,
                           0.099,
                           0.961,
                           0.9999,
                           0.949,
                           0.823,
                           0.085,
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
                 "-log(percieved_pi_obj[]) + pH[] * log(pi_H) + pL[] * (log(percieved_pi_obj[-1]) + tau * log(inflation_gap[])) = 0",
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
                                 9, 4, 8, 12, 13, 8, 11, 15, 18, 21,
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
                                 2, 2, 3, 2, 2, 2, 3, 2, 2, 1,
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
                                     24, 24, 26, 26, 27, 30, 30, 30),
                               j = c(4, 8, 6, 2, 3, 1, 1, 1, 15, 2,
                                     4, 8, 5, 7, 17, 7, 14, 2, 5, 7,
                                     17, 2, 5, 7, 17, 5, 7, 17, 3, 2,
                                     4, 8, 9, 16, 13, 10, 11, 12),
                               x = rep(1, 38), dims = c(31, 17))

# free parameters / calibrating equations map
freeparcalibreqmap__ <- sparseMatrix(i = NULL, j = NULL, dims = c(3, 17))

# shocks / equations map
shockeqmap__ <- sparseMatrix(i = c(13, 18, 19, 27, 30),
                             j = c(1, 2, 5, 4, 3),
                             x = rep(1, 5), dims = c(31, 5))

# steady state equations
ss_eq__ <- function(v, pc, pf)
{
    r <- numeric(31)
    r[1] = -v[16]
    r[2] = -v[5] + v[14]
    r[3] = -v[5] + pf[8] * v[17]^(-1 + pf[8]) * (1 - v[22])^(1 - pf[8]) * (v[17]^pf[8] * (1 - v[22])^(1 - pf[8]))^(-pf[4])
    r[4] = -v[13] + (1 + exp(pc[2] - pf[6] * log(v[4])))^-1
    r[5] = -v[14] + pf[2] * (v[5] * v[15] + v[14] * (1 - pf[3]))
    r[6] = -v[15] + pf[1] * v[6] * v[31] * v[21]^(-1 + pf[1]) * v[22]^(1 - pf[1])
    r[7] = -v[19] + pc[3] * v[1]
    r[8] = -v[23] + v[5]^-1 * v[14]
    r[9] = -v[27] + v[6] * v[31] * (1 - pf[1]) * v[21]^pf[1] * v[22]^(-pf[1])
    r[10] = -v[29] + v[31] * v[21]^pf[1] * v[22]^(1 - pf[1])
    r[11] = v[29] - v[30]
    r[12] = v[30] - v[7] * v[28]
    r[13] = -v[31] + exp(pf[15] * log(v[31]))
    r[14] = -v[5] * v[24]^-1 + pf[2] * v[5] * v[9]^-1
    r[15] = v[5] * v[27] + (-1 + pf[8]) * v[17]^pf[8] * (1 - v[22])^(-pf[8]) * (v[17]^pf[8] * (1 - v[22])^(1 - pf[8]))^(-pf[4])
    r[16] = -1 + pf[17] * (v[9]^-1 * v[9]^pf[5])^(-pf[7]^-1) + (1 - pf[17]) * v[10]^(-pf[7]^-1)
    r[17] = 1 - v[12] - v[13]
    r[18] = -v[2] + v[3] * (1 + pf[7])
    r[19] = -log(v[1]) + pf[14] * log(v[1])
    r[20] = -v[2] + v[5] * v[10] * v[28] + pf[2] * pf[17] * v[2] * (v[9]^-1 * v[9]^pf[5])^(-pf[7]^-1)
    r[21] = -v[3] + v[5] * v[6] * v[28] + pf[2] * pf[17] * v[3] * (v[9]^-1 * v[9]^pf[5])^(-pf[7]^-1 * (1 + pf[7]))
    r[22] = -v[7] + (1 - pf[17]) * v[10]^(-pf[7]^-1 * (1 + pf[7])) + pf[17] * v[7] * (v[9]^-1 * v[9]^pf[5])^(-pf[7]^-1 * (1 + pf[7]))
    r[23] = v[20] - v[21] + v[21] * (1 - pf[3])
    r[24] = v[26] - pf[2] * v[26] - (1 - pf[4])^-1 * (v[17]^pf[8] * (1 - v[22])^(1 - pf[8]))^(1 - pf[4])
    r[25] = -log(v[4]) - log(v[8]) + log(v[9])
    r[26] = -log(v[8]) + v[12] * log(pf[9]) + v[13] * (log(v[8]) + pf[16] * log(v[4]))
    r[27] = -log(v[11]) + pf[13] * log(v[11]) + log(v[8]) * (1 - pf[13])
    r[28] = -v[18] + v[28] - v[15] * v[21] - v[22] * v[27]
    r[29] = -v[19] + v[25] - v[9]^-1 * v[16] + v[16] * v[24]^-1
    r[30] = -pc[1] + (1 - pf[12]) * (log(v[11]) - pf[10] * log(v[11]))
    r[31] = -v[17] + v[18] - v[20] - v[25] + v[9]^-1 * v[16] + v[15] * v[21] - v[16] * v[24]^-1 + v[22] * v[27]

    return(r)
}

# calibrating equations
calibr_eq__ <- function(v, pc, pf)
{
    r <- numeric(3)
    r[1] = -0.18 + v[19] * v[28]^-1
    r[2] = -0.05 + v[13]
    r[3] = v[9] - v[11]

    return(r)
}

# steady state and calibrating equations Jacobian
ss_calibr_eq_jacob__ <- function(v, pc, pf)
{
    r <- numeric(3)
    jac <- numeric(108)
    jac[1] = -1
    jac[2] = -1
    jac[3] = 1
    jac[4] = -1
    jac[5] = pf[8] * (-1 + pf[8]) * v[17]^(-2 + pf[8]) * (1 - v[22])^(1 - pf[8]) * (v[17]^pf[8] * (1 - v[22])^(1 - pf[8]))^(-pf[4]) - pf[4] * pf[8]^2 * (v[17]^(-1 + pf[8]))^2 * ((1 - v[22])^(1 - pf[8]))^2 * (v[17]^pf[8] * (1 - v[22])^(1 - pf[8]))^(-1 - pf[4])
    jac[6] = pf[8] * (-1 + pf[8]) * v[17]^(-1 + pf[8]) * (1 - v[22])^(-pf[8]) * (v[17]^pf[8] * (1 - v[22])^(1 - pf[8]))^(-pf[4]) - pf[4] * pf[8] * (-1 + pf[8]) * v[17]^(-1 + 2 * pf[8]) * (1 - v[22])^(-pf[8]) * (1 - v[22])^(1 - pf[8]) * (v[17]^pf[8] * (1 - v[22])^(1 - pf[8]))^(-1 - pf[4])
    jac[7] = pf[6] * v[4]^-1 * exp(pc[2] - pf[6] * log(v[4])) * (1 + exp(pc[2] - pf[6] * log(v[4])))^-2
    jac[8] = -1
    jac[9] = -exp(pc[2] - pf[6] * log(v[4])) * (1 + exp(pc[2] - pf[6] * log(v[4])))^-2
    jac[10] = pf[2] * v[15]
    jac[11] = -1 + pf[2] * (1 - pf[3])
    jac[12] = pf[2] * v[5]
    jac[13] = pf[1] * v[31] * v[21]^(-1 + pf[1]) * v[22]^(1 - pf[1])
    jac[14] = -1
    jac[15] = pf[1] * v[6] * v[31] * (-1 + pf[1]) * v[21]^(-2 + pf[1]) * v[22]^(1 - pf[1])
    jac[16] = pf[1] * v[6] * v[31] * (1 - pf[1]) * v[21]^(-1 + pf[1]) * v[22]^(-pf[1])
    jac[17] = pf[1] * v[6] * v[21]^(-1 + pf[1]) * v[22]^(1 - pf[1])
    jac[18] = pc[3]
    jac[19] = -1
    jac[20] = v[1]
    jac[21] = -v[5]^-2 * v[14]
    jac[22] = v[5]^-1
    jac[23] = -1
    jac[24] = v[31] * (1 - pf[1]) * v[21]^pf[1] * v[22]^(-pf[1])
    jac[25] = pf[1] * v[6] * v[31] * (1 - pf[1]) * v[21]^(-1 + pf[1]) * v[22]^(-pf[1])
    jac[26] = -pf[1] * v[6] * v[31] * (1 - pf[1]) * v[21]^pf[1] * v[22]^(-1 - pf[1])
    jac[27] = -1
    jac[28] = v[6] * (1 - pf[1]) * v[21]^pf[1] * v[22]^(-pf[1])
    jac[29] = pf[1] * v[31] * v[21]^(-1 + pf[1]) * v[22]^(1 - pf[1])
    jac[30] = v[31] * (1 - pf[1]) * v[21]^pf[1] * v[22]^(-pf[1])
    jac[31] = -1
    jac[32] = v[21]^pf[1] * v[22]^(1 - pf[1])
    jac[33] = 1
    jac[34] = -1
    jac[35] = -v[28]
    jac[36] = -v[7]
    jac[37] = 1
    jac[38] = -1 + pf[15] * v[31]^-1 * exp(pf[15] * log(v[31]))
    jac[39] = -v[24]^-1 + pf[2] * v[9]^-1
    jac[40] = -pf[2] * v[5] * v[9]^-2
    jac[41] = v[5] * v[24]^-2
    jac[42] = v[27]
    jac[43] = pf[8] * (-1 + pf[8]) * v[17]^(-1 + pf[8]) * (1 - v[22])^(-pf[8]) * (v[17]^pf[8] * (1 - v[22])^(1 - pf[8]))^(-pf[4]) - pf[4] * pf[8] * (-1 + pf[8]) * v[17]^(-1 + 2 * pf[8]) * (1 - v[22])^(-pf[8]) * (1 - v[22])^(1 - pf[8]) * (v[17]^pf[8] * (1 - v[22])^(1 - pf[8]))^(-1 - pf[4])
    jac[44] = pf[8] * (-1 + pf[8]) * v[17]^pf[8] * (1 - v[22])^(-1 - pf[8]) * (v[17]^pf[8] * (1 - v[22])^(1 - pf[8]))^(-pf[4]) - pf[4] * (-1 + pf[8])^2 * (v[17]^pf[8])^2 * ((1 - v[22])^(-pf[8]))^2 * (v[17]^pf[8] * (1 - v[22])^(1 - pf[8]))^(-1 - pf[4])
    jac[45] = v[5]
    jac[46] = -pf[7]^-1 * pf[17] * (-v[9]^-2 * v[9]^pf[5] + pf[5] * v[9]^-1 * v[9]^(-1 + pf[5])) * (v[9]^-1 * v[9]^pf[5])^(-1 - pf[7]^-1)
    jac[47] = -pf[7]^-1 * (1 - pf[17]) * v[10]^(-1 - pf[7]^-1)
    jac[48] = -1
    jac[49] = -1
    jac[50] = -1
    jac[51] = 1 + pf[7]
    jac[52] = -v[1]^-1 + pf[14] * v[1]^-1
    jac[53] = -1 + pf[2] * pf[17] * (v[9]^-1 * v[9]^pf[5])^(-pf[7]^-1)
    jac[54] = v[10] * v[28]
    jac[55] = -pf[2] * pf[7]^-1 * pf[17] * v[2] * (-v[9]^-2 * v[9]^pf[5] + pf[5] * v[9]^-1 * v[9]^(-1 + pf[5])) * (v[9]^-1 * v[9]^pf[5])^(-1 - pf[7]^-1)
    jac[56] = v[5] * v[28]
    jac[57] = v[5] * v[10]
    jac[58] = -1 + pf[2] * pf[17] * (v[9]^-1 * v[9]^pf[5])^(-pf[7]^-1 * (1 + pf[7]))
    jac[59] = v[6] * v[28]
    jac[60] = v[5] * v[28]
    jac[61] = -pf[2] * pf[7]^-1 * pf[17] * v[3] * (1 + pf[7]) * (-v[9]^-2 * v[9]^pf[5] + pf[5] * v[9]^-1 * v[9]^(-1 + pf[5])) * (v[9]^-1 * v[9]^pf[5])^(-1 - pf[7]^-1 * (1 + pf[7]))
    jac[62] = v[5] * v[6]
    jac[63] = -1 + pf[17] * (v[9]^-1 * v[9]^pf[5])^(-pf[7]^-1 * (1 + pf[7]))
    jac[64] = -pf[7]^-1 * pf[17] * v[7] * (1 + pf[7]) * (-v[9]^-2 * v[9]^pf[5] + pf[5] * v[9]^-1 * v[9]^(-1 + pf[5])) * (v[9]^-1 * v[9]^pf[5])^(-1 - pf[7]^-1 * (1 + pf[7]))
    jac[65] = -pf[7]^-1 * (1 + pf[7]) * (1 - pf[17]) * v[10]^(-1 - pf[7]^-1 * (1 + pf[7]))
    jac[66] = 1
    jac[67] = -pf[3]
    jac[68] = -pf[8] * v[17]^(-1 + pf[8]) * (1 - v[22])^(1 - pf[8]) * (v[17]^pf[8] * (1 - v[22])^(1 - pf[8]))^(-pf[4])
    jac[69] = -(-1 + pf[8]) * v[17]^pf[8] * (1 - v[22])^(-pf[8]) * (v[17]^pf[8] * (1 - v[22])^(1 - pf[8]))^(-pf[4])
    jac[70] = 1 - pf[2]
    jac[71] = -v[4]^-1
    jac[72] = -v[8]^-1
    jac[73] = v[9]^-1
    jac[74] = pf[16] * v[4]^-1 * v[13]
    jac[75] = -v[8]^-1 + v[8]^-1 * v[13]
    jac[76] = log(pf[9])
    jac[77] = log(v[8]) + pf[16] * log(v[4])
    jac[78] = v[8]^-1 * (1 - pf[13])
    jac[79] = -v[11]^-1 + pf[13] * v[11]^-1
    jac[80] = -v[21]
    jac[81] = -1
    jac[82] = -v[15]
    jac[83] = -v[27]
    jac[84] = -v[22]
    jac[85] = 1
    jac[86] = v[9]^-2 * v[16]
    jac[87] = v[24]^-1 - v[9]^-1
    jac[88] = -1
    jac[89] = -v[16] * v[24]^-2
    jac[90] = 1
    jac[91] = (1 - pf[12]) * (v[11]^-1 - pf[10] * v[11]^-1)
    jac[92] = -1
    jac[93] = -v[9]^-2 * v[16]
    jac[94] = v[21]
    jac[95] = v[9]^-1 - v[24]^-1
    jac[96] = -1
    jac[97] = 1
    jac[98] = -1
    jac[99] = v[15]
    jac[100] = v[27]
    jac[101] = v[16] * v[24]^-2
    jac[102] = -1
    jac[103] = v[22]
    jac[104] = v[28]^-1
    jac[105] = -v[19] * v[28]^-2
    jac[106] = 1
    jac[107] = 1
    jac[108] = -1
    jacob <- sparseMatrix(i = c(1, 2, 2, 3, 3, 3, 4, 4, 4, 5,
                                5, 5, 6, 6, 6, 6, 6, 7, 7, 7,
                                8, 8, 8, 9, 9, 9, 9, 9, 10, 10,
                                10, 10, 11, 11, 12, 12, 12, 13, 14, 14,
                                14, 15, 15, 15, 15, 16, 16, 17, 17, 18,
                                18, 19, 20, 20, 20, 20, 20, 21, 21, 21,
                                21, 21, 22, 22, 22, 23, 23, 24, 24, 24,
                                25, 25, 25, 26, 26, 26, 26, 27, 27, 28,
                                28, 28, 28, 28, 28, 29, 29, 29, 29, 29,
                                30, 30, 31, 31, 31, 31, 31, 31, 31, 31,
                                31, 31, 31, 32, 32, 33, 34, 34),
                          j = c(16, 5, 14, 5, 17, 22, 4, 13, 33, 5,
                                14, 15, 6, 15, 21, 22, 31, 1, 19, 34,
                                5, 14, 23, 6, 21, 22, 27, 31, 21, 22,
                                29, 31, 29, 30, 7, 28, 30, 31, 5, 9,
                                24, 5, 17, 22, 27, 9, 10, 12, 13, 2,
                                3, 1, 2, 5, 9, 10, 28, 3, 5, 6,
                                9, 28, 7, 9, 10, 20, 21, 17, 22, 26,
                                4, 8, 9, 4, 8, 12, 13, 8, 11, 15,
                                18, 21, 22, 27, 28, 9, 16, 19, 24, 25,
                                11, 32, 9, 15, 16, 17, 18, 20, 21, 22,
                                24, 25, 27, 19, 28, 13, 9, 11),
                          x = jac, dims = c(34, 34))

    return(jacob)
}

# 1st order perturbation
pert1__ <- function(v, pc, pf)
{
    Atm1x <- numeric(17)
    Atm1x[1] = pf[1] * v[6] * v[31] * (-1 + pf[1]) * v[21]^(-2 + pf[1]) * v[22]^(1 - pf[1])
    Atm1x[2] = pf[1] * v[6] * v[31] * (1 - pf[1]) * v[21]^(-1 + pf[1]) * v[22]^(-pf[1])
    Atm1x[3] = pf[1] * v[31] * v[21]^(-1 + pf[1]) * v[22]^(1 - pf[1])
    Atm1x[4] = pf[15] * v[31]^-1 * exp(pf[15] * log(v[31]))
    Atm1x[5] = -pf[5] * pf[7]^-1 * pf[17] * v[9]^-1 * v[9]^(-1 + pf[5]) * (v[9]^-1 * v[9]^pf[5])^(-1 - pf[7]^-1)
    Atm1x[6] = pf[14] * v[1]^-1
    Atm1x[7] = pf[17] * (v[9]^-1 * v[9]^pf[5])^(-pf[7]^-1 * (1 + pf[7]))
    Atm1x[8] = -pf[5] * pf[7]^-1 * pf[17] * v[7] * v[9]^-1 * (1 + pf[7]) * v[9]^(-1 + pf[5]) * (v[9]^-1 * v[9]^pf[5])^(-1 - pf[7]^-1 * (1 + pf[7]))
    Atm1x[9] = 1 - pf[3]
    Atm1x[10] = v[8]^-1 * v[13]
    Atm1x[11] = pf[13] * v[11]^-1
    Atm1x[12] = -v[15]
    Atm1x[13] = -v[9]^-1
    Atm1x[14] = pf[10] * v[9]^-1 * (1 - pf[12])
    Atm1x[15] = pf[12] * v[24]^-1
    Atm1x[16] = v[9]^-1
    Atm1x[17] = v[15]
    Atm1 <- sparseMatrix(i = c(6, 9, 10, 13, 16, 19, 22, 22, 23, 26,
                               27, 28, 29, 30, 30, 31, 31),
                         j = c(21, 21, 21, 31, 9, 1, 7, 9, 21, 8,
                               11, 21, 16, 9, 24, 16, 21),
                         x = Atm1x, dims = c(31, 31))

    Atx <- numeric(94)
    Atx[1] = -1
    Atx[2] = -1
    Atx[3] = 1
    Atx[4] = -1
    Atx[5] = pf[8] * (-1 + pf[8]) * v[17]^(-2 + pf[8]) * (1 - v[22])^(1 - pf[8]) * (v[17]^pf[8] * (1 - v[22])^(1 - pf[8]))^(-pf[4]) - pf[4] * pf[8]^2 * v[17]^(-2 + 2 * pf[8]) * (1 - v[22])^(2 - 2 * pf[8]) * (v[17]^pf[8] * (1 - v[22])^(1 - pf[8]))^(-1 - pf[4])
    Atx[6] = pf[8] * (-1 + pf[8]) * v[17]^(-1 + pf[8]) * (1 - v[22])^(-pf[8]) * (v[17]^pf[8] * (1 - v[22])^(1 - pf[8]))^(-pf[4]) - pf[4] * pf[8] * (-1 + pf[8]) * v[17]^(-1 + 2 * pf[8]) * (1 - v[22])^(1 - 2 * pf[8]) * (v[17]^pf[8] * (1 - v[22])^(1 - pf[8]))^(-1 - pf[4])
    Atx[7] = pf[6] * v[4]^-1 * exp(pc[2] - pf[6] * log(v[4])) * (1 + exp(pc[2] - pf[6] * log(v[4])))^-2
    Atx[8] = -1
    Atx[9] = -1
    Atx[10] = pf[1] * v[31] * v[21]^(-1 + pf[1]) * v[22]^(1 - pf[1])
    Atx[11] = -1
    Atx[12] = pf[1] * v[6] * v[31] * (1 - pf[1]) * v[21]^(-1 + pf[1]) * v[22]^(-pf[1])
    Atx[13] = pf[1] * v[6] * v[21]^(-1 + pf[1]) * v[22]^(1 - pf[1])
    Atx[14] = pc[3]
    Atx[15] = -1
    Atx[16] = -v[5]^-2 * v[14]
    Atx[17] = v[5]^-1
    Atx[18] = -1
    Atx[19] = v[31] * (1 - pf[1]) * v[21]^pf[1] * v[22]^(-pf[1])
    Atx[20] = -pf[1] * v[6] * v[31] * (1 - pf[1]) * v[21]^pf[1] * v[22]^(-1 - pf[1])
    Atx[21] = -1
    Atx[22] = v[6] * (1 - pf[1]) * v[21]^pf[1] * v[22]^(-pf[1])
    Atx[23] = v[31] * (1 - pf[1]) * v[21]^pf[1] * v[22]^(-pf[1])
    Atx[24] = -1
    Atx[25] = v[21]^pf[1] * v[22]^(1 - pf[1])
    Atx[26] = 1
    Atx[27] = -1
    Atx[28] = -v[28]
    Atx[29] = -v[7]
    Atx[30] = 1
    Atx[31] = -1
    Atx[32] = -v[24]^-1
    Atx[33] = v[5] * v[24]^-2
    Atx[34] = v[27]
    Atx[35] = pf[8] * (-1 + pf[8]) * v[17]^(-1 + pf[8]) * (1 - v[22])^(-pf[8]) * (v[17]^pf[8] * (1 - v[22])^(1 - pf[8]))^(-pf[4]) - pf[4] * pf[8] * (-1 + pf[8]) * v[17]^(-1 + 2 * pf[8]) * (1 - v[22])^(1 - 2 * pf[8]) * (v[17]^pf[8] * (1 - v[22])^(1 - pf[8]))^(-1 - pf[4])
    Atx[36] = pf[8] * (-1 + pf[8]) * v[17]^pf[8] * (1 - v[22])^(-1 - pf[8]) * (v[17]^pf[8] * (1 - v[22])^(1 - pf[8]))^(-pf[4]) - pf[4] * (-1 + pf[8])^2 * v[17]^(2 * pf[8]) * (1 - v[22])^(-2 * pf[8]) * (v[17]^pf[8] * (1 - v[22])^(1 - pf[8]))^(-1 - pf[4])
    Atx[37] = v[5]
    Atx[38] = pf[7]^-1 * pf[17] * v[9]^-2 * v[9]^pf[5] * (v[9]^-1 * v[9]^pf[5])^(-1 - pf[7]^-1)
    Atx[39] = -pf[7]^-1 * (1 - pf[17]) * v[10]^(-1 - pf[7]^-1)
    Atx[40] = -1
    Atx[41] = -1
    Atx[42] = -1
    Atx[43] = 1 + pf[7]
    Atx[44] = -v[1]^-1
    Atx[45] = -1
    Atx[46] = v[10] * v[28]
    Atx[47] = -pf[2] * pf[5] * pf[7]^-1 * pf[17] * v[2] * v[9]^-1 * v[9]^(-1 + pf[5]) * (v[9]^-1 * v[9]^pf[5])^(-1 - pf[7]^-1)
    Atx[48] = v[5] * v[28] + pf[2] * pf[17] * v[2] * v[10]^-1 * (v[9]^-1 * v[9]^pf[5])^(-pf[7]^-1)
    Atx[49] = v[5] * v[10]
    Atx[50] = -1
    Atx[51] = v[6] * v[28]
    Atx[52] = v[5] * v[28]
    Atx[53] = -pf[2] * pf[5] * pf[7]^-1 * pf[17] * v[3] * v[9]^-1 * (1 + pf[7]) * v[9]^(-1 + pf[5]) * (v[9]^-1 * v[9]^pf[5])^(-1 - pf[7]^-1 * (1 + pf[7]))
    Atx[54] = v[5] * v[6]
    Atx[55] = -1
    Atx[56] = pf[7]^-1 * pf[17] * v[7] * v[9]^-2 * (1 + pf[7]) * v[9]^pf[5] * (v[9]^-1 * v[9]^pf[5])^(-1 - pf[7]^-1 * (1 + pf[7]))
    Atx[57] = -pf[7]^-1 * (1 + pf[7]) * (1 - pf[17]) * v[10]^(-1 - pf[7]^-1 * (1 + pf[7]))
    Atx[58] = 1
    Atx[59] = -1
    Atx[60] = -pf[8] * v[17]^(-1 + pf[8]) * (1 - v[22])^(1 - pf[8]) * (v[17]^pf[8] * (1 - v[22])^(1 - pf[8]))^(-pf[4])
    Atx[61] = (1 - pf[8]) * v[17]^pf[8] * (1 - v[22])^(-pf[8]) * (v[17]^pf[8] * (1 - v[22])^(1 - pf[8]))^(-pf[4])
    Atx[62] = 1
    Atx[63] = -v[4]^-1
    Atx[64] = -v[8]^-1
    Atx[65] = v[9]^-1
    Atx[66] = pf[16] * v[4]^-1 * v[13]
    Atx[67] = -v[8]^-1
    Atx[68] = log(pf[9])
    Atx[69] = log(v[8]) + pf[16] * log(v[4])
    Atx[70] = v[8]^-1 * (1 - pf[13])
    Atx[71] = -v[11]^-1
    Atx[72] = -v[21]
    Atx[73] = -1
    Atx[74] = -v[27]
    Atx[75] = -v[22]
    Atx[76] = 1
    Atx[77] = v[9]^-2 * v[16]
    Atx[78] = v[24]^-1
    Atx[79] = -1
    Atx[80] = -v[16] * v[24]^-2
    Atx[81] = 1
    Atx[82] = (1 - pf[12]) * (v[11]^-1 - pf[10] * v[11]^-1)
    Atx[83] = -v[24]^-1
    Atx[84] = pf[11] * v[28]^-1 * (1 - pf[12])
    Atx[85] = -v[9]^-2 * v[16]
    Atx[86] = v[21]
    Atx[87] = -v[24]^-1
    Atx[88] = -1
    Atx[89] = 1
    Atx[90] = -1
    Atx[91] = v[27]
    Atx[92] = v[16] * v[24]^-2
    Atx[93] = -1
    Atx[94] = v[22]
    At <- sparseMatrix(i = c(1, 2, 2, 3, 3, 3, 4, 4, 5, 6,
                             6, 6, 6, 7, 7, 8, 8, 8, 9, 9,
                             9, 9, 10, 10, 10, 11, 11, 12, 12, 12,
                             13, 14, 14, 15, 15, 15, 15, 16, 16, 17,
                             17, 18, 18, 19, 20, 20, 20, 20, 20, 21,
                             21, 21, 21, 21, 22, 22, 22, 23, 23, 24,
                             24, 24, 25, 25, 25, 26, 26, 26, 26, 27,
                             27, 28, 28, 28, 28, 28, 29, 29, 29, 29,
                             29, 30, 30, 30, 31, 31, 31, 31, 31, 31,
                             31, 31, 31, 31),
                       j = c(16, 5, 14, 5, 17, 22, 4, 13, 14, 6,
                             15, 22, 31, 1, 19, 5, 14, 23, 6, 22,
                             27, 31, 22, 29, 31, 29, 30, 7, 28, 30,
                             31, 5, 24, 5, 17, 22, 27, 9, 10, 12,
                             13, 2, 3, 1, 2, 5, 9, 10, 28, 3,
                             5, 6, 9, 28, 7, 9, 10, 20, 21, 17,
                             22, 26, 4, 8, 9, 4, 8, 12, 13, 8,
                             11, 15, 18, 22, 27, 28, 9, 16, 19, 24,
                             25, 11, 24, 28, 9, 15, 16, 17, 18, 20,
                             22, 24, 25, 27),
                         x = Atx, dims = c(31, 31))

    Atp1x <- numeric(11)
    Atp1x[1] = pf[2] * v[15]
    Atp1x[2] = pf[2] * (1 - pf[3])
    Atp1x[3] = pf[2] * v[5]
    Atp1x[4] = pf[2] * v[9]^-1
    Atp1x[5] = -pf[2] * v[5] * v[9]^-2
    Atp1x[6] = pf[2] * pf[17] * (v[9]^-1 * v[9]^pf[5])^(-pf[7]^-1)
    Atp1x[7] = pf[2] * pf[7]^-1 * pf[17] * v[2] * v[9]^-2 * v[9]^pf[5] * (v[9]^-1 * v[9]^pf[5])^(-1 - pf[7]^-1)
    Atp1x[8] = -pf[2] * pf[17] * v[2] * v[10]^-1 * (v[9]^-1 * v[9]^pf[5])^(-pf[7]^-1)
    Atp1x[9] = pf[2] * pf[17] * (v[9]^-1 * v[9]^pf[5])^(-pf[7]^-1 * (1 + pf[7]))
    Atp1x[10] = pf[2] * pf[7]^-1 * pf[17] * v[3] * v[9]^-2 * (1 + pf[7]) * v[9]^pf[5] * (v[9]^-1 * v[9]^pf[5])^(-1 - pf[7]^-1 * (1 + pf[7]))
    Atp1x[11] = -pf[2]
    Atp1 <- sparseMatrix(i = c(5, 5, 5, 14, 14, 20, 20, 20, 21, 21,
                               24),
                         j = c(5, 14, 15, 5, 9, 2, 9, 10, 3, 9,
                               26),
                         x = Atp1x, dims = c(31, 31))

    Aepsx <- numeric(5)
    Aepsx[1] = exp(pf[15] * log(v[31]))
    Aepsx[2] = 1
    Aepsx[3] = 1
    Aepsx[4] = 1
    Aepsx[5] = 1
    Aeps <- sparseMatrix(i = c(13, 18, 19, 27, 30),
                         j = c(1, 2, 5, 4, 3),
                         x = Aepsx, dims = c(31, 5))

    return(list(Atm1, At, Atp1, Aeps))
}

ext__ <- list()

# create model object
gecon_model(model_info = info__,
            index_sets = index_sets__,
            variables = variables__,
            variables_tex = variables_tex__,
            shocks = shocks__,
            shocks_tex = shocks_tex__,
            parameters = parameters__,
            parameters_tex = parameters_tex__,
            parameters_free = parameters_free__,
            parameters_free_val = parameters_free_val__,
            equations = equations__,
            calibr_equations = calibr_equations__,
            var_eq_map = vareqmap__,
            shock_eq_map = shockeqmap__,
            var_ceq_map = varcalibreqmap__,
            cpar_eq_map = calibrpareqmap__,
            cpar_ceq_map = calibrparcalibreqmap__,
            fpar_eq_map = freepareqmap__,
            fpar_ceq_map = freeparcalibreqmap__,
            ss_function = ss_eq__,
            calibr_function = calibr_eq__,
            ss_calibr_jac_function = ss_calibr_eq_jacob__,
            pert = pert1__,
            ext = ext__)
