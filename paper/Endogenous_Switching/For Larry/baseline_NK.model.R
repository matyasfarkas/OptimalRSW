# Generated on 2025-04-15 22:47:12 by gEcon ver. 1.2.1 (2023-01-18)
# http://gecon.r-forge.r-project.org/

# Model name: baseline_NK

# info
info__ <- c("baseline_NK", "C:/Users/fm007/Documents/GitHub/OptimalRSW/paper/Endogenous_Switching/For Larry/baseline_NK.gcn", "2025-04-15 22:47:12", "false")

# index sets
index_sets__ <- list()

# variables
variables__ <- c("epsilon_G",
                 "g_1",
                 "g_2",
                 "lambda",
                 "mc",
                 "nu_p",
                 "pi",
                 "pi_star",
                 "pi_obj",
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
                     "\\lambda",
                     "{m\\!c}",
                     "\\nu^{\\mathrm{p}}",
                     "\\pi",
                     "\\pi^{\\star}",
                     "\\pi^{\\mathrm{obj}}",
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
                  "calibr_pi_obj",
                  "delta",
                  "eta",
                  "gamma_p",
                  "lambda_p",
                  "mu",
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
                     "{c\\!a\\!l\\!i\\!b\\!r}^{\\pi^{\\mathrm{obj}}}",
                     "\\delta",
                     "\\eta",
                     "\\gamma^{\\mathrm{p}}",
                     "\\lambda^{\\mathrm{p}}",
                     "\\mu",
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
                       "lambda_p",
                       "mu",
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
                           0.5,
                           0.3,
                           1.684,
                           0.099,
                           0.961,
                           0.9999,
                           0.949,
                           0.823,
                           0.908)

# equations
equations__ <- c("-B[] = 0",
                 "-lambda[] + q[] = 0",
                 "-lambda[] + mu * C[]^(-1 + mu) * (1 - L_s[])^(1 - mu) * (C[]^mu * (1 - L_s[])^(1 - mu))^(-eta) = 0",
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
                 "eta_p[] - g_1[] + g_2[] * (1 + lambda_p) = 0",
                 "eta_G[] - log(epsilon_G[]) + rho_G * log(epsilon_G[-1]) = 0",
                 "-g_1[] + lambda[] * pi_star[] * Y[] + beta * xi_p * pi_star[] * E[][g_1[1] * pi_star[1]^-1 * (pi[1]^-1 * pi[]^gamma_p)^(-lambda_p^-1)] = 0",
                 "-g_2[] + beta * xi_p * E[][g_2[1] * (pi[1]^-1 * pi[]^gamma_p)^(-lambda_p^-1 * (1 + lambda_p))] + lambda[] * mc[] * Y[] = 0",
                 "-nu_p[] + (1 - xi_p) * pi_star[]^(-lambda_p^-1 * (1 + lambda_p)) + xi_p * nu_p[-1] * (pi[]^-1 * pi[-1]^gamma_p)^(-lambda_p^-1 * (1 + lambda_p)) = 0",
                 "I[] - K_s[] + K_s[-1] * (1 - delta) = 0",
                 "U[] - beta * E[][U[1]] - (1 - eta)^-1 * (C[]^mu * (1 - L_s[])^(1 - mu))^(1 - eta) = 0",
                 "eta_pi[] - log(pi_obj[]) + rho_pi_bar * log(pi_obj[-1]) + log(calibr_pi_obj) * (1 - rho_pi_bar) = 0",
                 "-Div[] + Y[] - K_s[-1] * r[] - L_s[] * W[] = 0",
                 "-G[] + T[] - B[-1] * pi[]^-1 + B[] * R[]^-1 = 0",
                 "-calibr_pi + eta_R[] - log(R[ss]^-1 * R[]) + rho * log(R[ss]^-1 * R[-1]) + (1 - rho) * (log(pi_obj[]) + r_pi * (-log(pi_obj[]) + log(pi[ss]^-1 * pi[-1])) + r_Y * log(Y[ss]^-1 * Y[])) = 0",
                 "-C[] + Div[] - I[] - T[] + B[-1] * pi[]^-1 + K_s[-1] * r[] - B[] * R[]^-1 + L_s[] * W[] = 0")

# calibrating equations
calibr_equations__ <- c("-1 + pi_obj[ss] = 0",
                        "-0.18 + G[ss] * Y[ss]^-1 = 0",
                        "pi[ss] - pi_obj[ss] = 0")

# variables / equations map
vareqmap__ <- sparseMatrix(i = c(1, 2, 2, 3, 3, 3, 4, 4, 4, 5,
                                 5, 5, 5, 5, 6, 6, 7, 7, 7, 8,
                                 8, 8, 8, 8, 9, 9, 9, 9, 10, 10,
                                 11, 11, 11, 12, 13, 13, 13, 14, 14, 14,
                                 14, 15, 15, 16, 16, 17, 18, 18, 18, 18,
                                 18, 19, 19, 19, 19, 19, 20, 20, 20, 21,
                                 21, 22, 22, 22, 23, 24, 24, 24, 24, 24,
                                 24, 25, 25, 25, 25, 25, 26, 26, 26, 26,
                                 27, 27, 27, 27, 27, 27, 27, 27, 27, 27,
                                 27),
                           j = c(12, 4, 10, 4, 13, 18, 4, 10, 11, 5,
                                 11, 17, 18, 27, 1, 15, 4, 10, 19, 5,
                                 17, 18, 23, 27, 17, 18, 25, 27, 25, 26,
                                 6, 24, 26, 27, 4, 7, 20, 4, 13, 18,
                                 23, 7, 8, 2, 3, 1, 2, 4, 7, 8,
                                 24, 3, 4, 5, 7, 24, 6, 7, 8, 16,
                                 17, 13, 18, 22, 9, 11, 14, 17, 18, 23,
                                 24, 7, 12, 15, 20, 21, 7, 9, 20, 24,
                                 7, 11, 12, 13, 14, 16, 17, 18, 20, 21,
                                 23),
                           x = c(2, 2, 2, 2, 2, 2, 4, 6, 4, 2,
                                 2, 1, 2, 2, 2, 2, 2, 2, 2, 2,
                                 1, 2, 2, 2, 1, 2, 2, 2, 2, 2,
                                 2, 2, 2, 3, 6, 4, 2, 2, 2, 2,
                                 2, 3, 2, 2, 2, 3, 6, 2, 6, 6,
                                 2, 6, 2, 2, 6, 2, 3, 3, 2, 2,
                                 3, 2, 2, 6, 3, 2, 2, 1, 2, 2,
                                 2, 2, 3, 2, 2, 2, 9, 2, 11, 10,
                                 2, 2, 3, 2, 2, 2, 1, 2, 2, 2,
                                 2),
                           dims = c(27, 27))

# variables / calibrating equations map
varcalibreqmap__ <- sparseMatrix(i = c(1, 2, 2, 3, 3),
                                 j = c(9, 15, 24, 7, 9),
                                 x = rep(1, 5), dims = c(3, 27))

# calibrated parameters / equations map
calibrpareqmap__ <- sparseMatrix(i = c(6, 23, 26),
                                 j = c(3, 2, 1),
                                 x = rep(1, 3), dims = c(27, 3))

# calibrated parameters / calibrating equations map
calibrparcalibreqmap__ <- sparseMatrix(i = NULL, j = NULL, dims = c(3, 3))

# free parameters / equations map
freepareqmap__ <- sparseMatrix(i = c(3, 3, 4, 4, 5, 8, 9, 12, 13, 14,
                                     14, 15, 15, 15, 16, 17, 18, 18, 18, 18,
                                     19, 19, 19, 19, 20, 20, 20, 21, 22, 22,
                                     22, 23, 26, 26, 26),
                               j = c(4, 7, 2, 3, 1, 1, 1, 13, 2, 4,
                                     7, 5, 6, 14, 6, 12, 2, 5, 6, 14,
                                     2, 5, 6, 14, 5, 6, 14, 3, 2, 4,
                                     7, 11, 8, 9, 10),
                               x = rep(1, 35), dims = c(27, 14))

# free parameters / calibrating equations map
freeparcalibreqmap__ <- sparseMatrix(i = NULL, j = NULL, dims = c(3, 14))

# shocks / equations map
shockeqmap__ <- sparseMatrix(i = c(12, 16, 17, 23, 26),
                             j = c(1, 2, 5, 4, 3),
                             x = rep(1, 5), dims = c(27, 5))

# steady state equations
ss_eq__ <- function(v, pc, pf)
{
    r <- numeric(27)
    r[1] = -v[12]
    r[2] = -v[4] + v[10]
    r[3] = -v[4] + pf[7] * v[13]^(-1 + pf[7]) * (1 - v[18])^(1 - pf[7]) * (v[13]^pf[7] * (1 - v[18])^(1 - pf[7]))^(-pf[4])
    r[4] = -v[10] + pf[2] * (v[4] * v[11] + v[10] * (1 - pf[3]))
    r[5] = -v[11] + pf[1] * v[5] * v[27] * v[17]^(-1 + pf[1]) * v[18]^(1 - pf[1])
    r[6] = -v[15] + pc[3] * v[1]
    r[7] = -v[19] + v[4]^-1 * v[10]
    r[8] = -v[23] + v[5] * v[27] * (1 - pf[1]) * v[17]^pf[1] * v[18]^(-pf[1])
    r[9] = -v[25] + v[27] * v[17]^pf[1] * v[18]^(1 - pf[1])
    r[10] = v[25] - v[26]
    r[11] = v[26] - v[6] * v[24]
    r[12] = -v[27] + exp(pf[13] * log(v[27]))
    r[13] = -v[4] * v[20]^-1 + pf[2] * v[4] * v[7]^-1
    r[14] = v[4] * v[23] + (-1 + pf[7]) * v[13]^pf[7] * (1 - v[18])^(-pf[7]) * (v[13]^pf[7] * (1 - v[18])^(1 - pf[7]))^(-pf[4])
    r[15] = -1 + pf[14] * (v[7]^-1 * v[7]^pf[5])^(-pf[6]^-1) + (1 - pf[14]) * v[8]^(-pf[6]^-1)
    r[16] = -v[2] + v[3] * (1 + pf[6])
    r[17] = -log(v[1]) + pf[12] * log(v[1])
    r[18] = -v[2] + v[4] * v[8] * v[24] + pf[2] * pf[14] * v[2] * (v[7]^-1 * v[7]^pf[5])^(-pf[6]^-1)
    r[19] = -v[3] + v[4] * v[5] * v[24] + pf[2] * pf[14] * v[3] * (v[7]^-1 * v[7]^pf[5])^(-pf[6]^-1 * (1 + pf[6]))
    r[20] = -v[6] + (1 - pf[14]) * v[8]^(-pf[6]^-1 * (1 + pf[6])) + pf[14] * v[6] * (v[7]^-1 * v[7]^pf[5])^(-pf[6]^-1 * (1 + pf[6]))
    r[21] = v[16] - v[17] + v[17] * (1 - pf[3])
    r[22] = v[22] - pf[2] * v[22] - (1 - pf[4])^-1 * (v[13]^pf[7] * (1 - v[18])^(1 - pf[7]))^(1 - pf[4])
    r[23] = -log(v[9]) + pf[11] * log(v[9]) + log(pc[2]) * (1 - pf[11])
    r[24] = -v[14] + v[24] - v[11] * v[17] - v[18] * v[23]
    r[25] = -v[15] + v[21] - v[7]^-1 * v[12] + v[12] * v[20]^-1
    r[26] = -pc[1] + (1 - pf[10]) * (log(v[9]) - pf[8] * log(v[9]))
    r[27] = -v[13] + v[14] - v[16] - v[21] + v[7]^-1 * v[12] + v[11] * v[17] - v[12] * v[20]^-1 + v[18] * v[23]

    return(r)
}

# calibrating equations
calibr_eq__ <- function(v, pc, pf)
{
    r <- numeric(3)
    r[1] = -1 + v[9]
    r[2] = -0.18 + v[15] * v[24]^-1
    r[3] = v[7] - v[9]

    return(r)
}

# steady state and calibrating equations Jacobian
ss_calibr_eq_jacob__ <- function(v, pc, pf)
{
    r <- numeric(3)
    jac <- numeric(96)
    jac[1] = -1
    jac[2] = -1
    jac[3] = 1
    jac[4] = -1
    jac[5] = pf[7] * (-1 + pf[7]) * v[13]^(-2 + pf[7]) * (1 - v[18])^(1 - pf[7]) * (v[13]^pf[7] * (1 - v[18])^(1 - pf[7]))^(-pf[4]) - pf[4] * pf[7]^2 * (v[13]^(-1 + pf[7]))^2 * ((1 - v[18])^(1 - pf[7]))^2 * (v[13]^pf[7] * (1 - v[18])^(1 - pf[7]))^(-1 - pf[4])
    jac[6] = pf[7] * (-1 + pf[7]) * v[13]^(-1 + pf[7]) * (1 - v[18])^(-pf[7]) * (v[13]^pf[7] * (1 - v[18])^(1 - pf[7]))^(-pf[4]) - pf[4] * pf[7] * (-1 + pf[7]) * v[13]^(-1 + 2 * pf[7]) * (1 - v[18])^(-pf[7]) * (1 - v[18])^(1 - pf[7]) * (v[13]^pf[7] * (1 - v[18])^(1 - pf[7]))^(-1 - pf[4])
    jac[7] = pf[2] * v[11]
    jac[8] = -1 + pf[2] * (1 - pf[3])
    jac[9] = pf[2] * v[4]
    jac[10] = pf[1] * v[27] * v[17]^(-1 + pf[1]) * v[18]^(1 - pf[1])
    jac[11] = -1
    jac[12] = pf[1] * v[5] * v[27] * (-1 + pf[1]) * v[17]^(-2 + pf[1]) * v[18]^(1 - pf[1])
    jac[13] = pf[1] * v[5] * v[27] * (1 - pf[1]) * v[17]^(-1 + pf[1]) * v[18]^(-pf[1])
    jac[14] = pf[1] * v[5] * v[17]^(-1 + pf[1]) * v[18]^(1 - pf[1])
    jac[15] = pc[3]
    jac[16] = -1
    jac[17] = v[1]
    jac[18] = -v[4]^-2 * v[10]
    jac[19] = v[4]^-1
    jac[20] = -1
    jac[21] = v[27] * (1 - pf[1]) * v[17]^pf[1] * v[18]^(-pf[1])
    jac[22] = pf[1] * v[5] * v[27] * (1 - pf[1]) * v[17]^(-1 + pf[1]) * v[18]^(-pf[1])
    jac[23] = -pf[1] * v[5] * v[27] * (1 - pf[1]) * v[17]^pf[1] * v[18]^(-1 - pf[1])
    jac[24] = -1
    jac[25] = v[5] * (1 - pf[1]) * v[17]^pf[1] * v[18]^(-pf[1])
    jac[26] = pf[1] * v[27] * v[17]^(-1 + pf[1]) * v[18]^(1 - pf[1])
    jac[27] = v[27] * (1 - pf[1]) * v[17]^pf[1] * v[18]^(-pf[1])
    jac[28] = -1
    jac[29] = v[17]^pf[1] * v[18]^(1 - pf[1])
    jac[30] = 1
    jac[31] = -1
    jac[32] = -v[24]
    jac[33] = -v[6]
    jac[34] = 1
    jac[35] = -1 + pf[13] * v[27]^-1 * exp(pf[13] * log(v[27]))
    jac[36] = -v[20]^-1 + pf[2] * v[7]^-1
    jac[37] = -pf[2] * v[4] * v[7]^-2
    jac[38] = v[4] * v[20]^-2
    jac[39] = v[23]
    jac[40] = pf[7] * (-1 + pf[7]) * v[13]^(-1 + pf[7]) * (1 - v[18])^(-pf[7]) * (v[13]^pf[7] * (1 - v[18])^(1 - pf[7]))^(-pf[4]) - pf[4] * pf[7] * (-1 + pf[7]) * v[13]^(-1 + 2 * pf[7]) * (1 - v[18])^(-pf[7]) * (1 - v[18])^(1 - pf[7]) * (v[13]^pf[7] * (1 - v[18])^(1 - pf[7]))^(-1 - pf[4])
    jac[41] = pf[7] * (-1 + pf[7]) * v[13]^pf[7] * (1 - v[18])^(-1 - pf[7]) * (v[13]^pf[7] * (1 - v[18])^(1 - pf[7]))^(-pf[4]) - pf[4] * (-1 + pf[7])^2 * (v[13]^pf[7])^2 * ((1 - v[18])^(-pf[7]))^2 * (v[13]^pf[7] * (1 - v[18])^(1 - pf[7]))^(-1 - pf[4])
    jac[42] = v[4]
    jac[43] = -pf[6]^-1 * pf[14] * (-v[7]^-2 * v[7]^pf[5] + pf[5] * v[7]^-1 * v[7]^(-1 + pf[5])) * (v[7]^-1 * v[7]^pf[5])^(-1 - pf[6]^-1)
    jac[44] = -pf[6]^-1 * (1 - pf[14]) * v[8]^(-1 - pf[6]^-1)
    jac[45] = -1
    jac[46] = 1 + pf[6]
    jac[47] = -v[1]^-1 + pf[12] * v[1]^-1
    jac[48] = -1 + pf[2] * pf[14] * (v[7]^-1 * v[7]^pf[5])^(-pf[6]^-1)
    jac[49] = v[8] * v[24]
    jac[50] = -pf[2] * pf[6]^-1 * pf[14] * v[2] * (-v[7]^-2 * v[7]^pf[5] + pf[5] * v[7]^-1 * v[7]^(-1 + pf[5])) * (v[7]^-1 * v[7]^pf[5])^(-1 - pf[6]^-1)
    jac[51] = v[4] * v[24]
    jac[52] = v[4] * v[8]
    jac[53] = -1 + pf[2] * pf[14] * (v[7]^-1 * v[7]^pf[5])^(-pf[6]^-1 * (1 + pf[6]))
    jac[54] = v[5] * v[24]
    jac[55] = v[4] * v[24]
    jac[56] = -pf[2] * pf[6]^-1 * pf[14] * v[3] * (1 + pf[6]) * (-v[7]^-2 * v[7]^pf[5] + pf[5] * v[7]^-1 * v[7]^(-1 + pf[5])) * (v[7]^-1 * v[7]^pf[5])^(-1 - pf[6]^-1 * (1 + pf[6]))
    jac[57] = v[4] * v[5]
    jac[58] = -1 + pf[14] * (v[7]^-1 * v[7]^pf[5])^(-pf[6]^-1 * (1 + pf[6]))
    jac[59] = -pf[6]^-1 * pf[14] * v[6] * (1 + pf[6]) * (-v[7]^-2 * v[7]^pf[5] + pf[5] * v[7]^-1 * v[7]^(-1 + pf[5])) * (v[7]^-1 * v[7]^pf[5])^(-1 - pf[6]^-1 * (1 + pf[6]))
    jac[60] = -pf[6]^-1 * (1 + pf[6]) * (1 - pf[14]) * v[8]^(-1 - pf[6]^-1 * (1 + pf[6]))
    jac[61] = 1
    jac[62] = -pf[3]
    jac[63] = -pf[7] * v[13]^(-1 + pf[7]) * (1 - v[18])^(1 - pf[7]) * (v[13]^pf[7] * (1 - v[18])^(1 - pf[7]))^(-pf[4])
    jac[64] = -(-1 + pf[7]) * v[13]^pf[7] * (1 - v[18])^(-pf[7]) * (v[13]^pf[7] * (1 - v[18])^(1 - pf[7]))^(-pf[4])
    jac[65] = 1 - pf[2]
    jac[66] = -v[9]^-1 + pf[11] * v[9]^-1
    jac[67] = pc[2]^-1 * (1 - pf[11])
    jac[68] = -v[17]
    jac[69] = -1
    jac[70] = -v[11]
    jac[71] = -v[23]
    jac[72] = -v[18]
    jac[73] = 1
    jac[74] = v[7]^-2 * v[12]
    jac[75] = v[20]^-1 - v[7]^-1
    jac[76] = -1
    jac[77] = -v[12] * v[20]^-2
    jac[78] = 1
    jac[79] = (1 - pf[10]) * (v[9]^-1 - pf[8] * v[9]^-1)
    jac[80] = -1
    jac[81] = -v[7]^-2 * v[12]
    jac[82] = v[17]
    jac[83] = v[7]^-1 - v[20]^-1
    jac[84] = -1
    jac[85] = 1
    jac[86] = -1
    jac[87] = v[11]
    jac[88] = v[23]
    jac[89] = v[12] * v[20]^-2
    jac[90] = -1
    jac[91] = v[18]
    jac[92] = 1
    jac[93] = v[24]^-1
    jac[94] = -v[15] * v[24]^-2
    jac[95] = 1
    jac[96] = -1
    jacob <- sparseMatrix(i = c(1, 2, 2, 3, 3, 3, 4, 4, 4, 5,
                                5, 5, 5, 5, 6, 6, 6, 7, 7, 7,
                                8, 8, 8, 8, 8, 9, 9, 9, 9, 10,
                                10, 11, 11, 11, 12, 13, 13, 13, 14, 14,
                                14, 14, 15, 15, 16, 16, 17, 18, 18, 18,
                                18, 18, 19, 19, 19, 19, 19, 20, 20, 20,
                                21, 21, 22, 22, 22, 23, 23, 24, 24, 24,
                                24, 24, 24, 25, 25, 25, 25, 25, 26, 26,
                                27, 27, 27, 27, 27, 27, 27, 27, 27, 27,
                                27, 28, 29, 29, 30, 30),
                          j = c(12, 4, 10, 4, 13, 18, 4, 10, 11, 5,
                                11, 17, 18, 27, 1, 15, 30, 4, 10, 19,
                                5, 17, 18, 23, 27, 17, 18, 25, 27, 25,
                                26, 6, 24, 26, 27, 4, 7, 20, 4, 13,
                                18, 23, 7, 8, 2, 3, 1, 2, 4, 7,
                                8, 24, 3, 4, 5, 7, 24, 6, 7, 8,
                                16, 17, 13, 18, 22, 9, 29, 11, 14, 17,
                                18, 23, 24, 7, 12, 15, 20, 21, 9, 28,
                                7, 11, 12, 13, 14, 16, 17, 18, 20, 21,
                                23, 9, 15, 24, 7, 9),
                          x = jac, dims = c(30, 30))

    return(jacob)
}

# 1st order perturbation
pert1__ <- function(v, pc, pf)
{
    Atm1x <- numeric(16)
    Atm1x[1] = pf[1] * v[5] * v[27] * (-1 + pf[1]) * v[17]^(-2 + pf[1]) * v[18]^(1 - pf[1])
    Atm1x[2] = pf[1] * v[5] * v[27] * (1 - pf[1]) * v[17]^(-1 + pf[1]) * v[18]^(-pf[1])
    Atm1x[3] = pf[1] * v[27] * v[17]^(-1 + pf[1]) * v[18]^(1 - pf[1])
    Atm1x[4] = pf[13] * v[27]^-1 * exp(pf[13] * log(v[27]))
    Atm1x[5] = -pf[5] * pf[6]^-1 * pf[14] * v[7]^-1 * v[7]^(-1 + pf[5]) * (v[7]^-1 * v[7]^pf[5])^(-1 - pf[6]^-1)
    Atm1x[6] = pf[12] * v[1]^-1
    Atm1x[7] = pf[14] * (v[7]^-1 * v[7]^pf[5])^(-pf[6]^-1 * (1 + pf[6]))
    Atm1x[8] = -pf[5] * pf[6]^-1 * pf[14] * v[6] * v[7]^-1 * (1 + pf[6]) * v[7]^(-1 + pf[5]) * (v[7]^-1 * v[7]^pf[5])^(-1 - pf[6]^-1 * (1 + pf[6]))
    Atm1x[9] = 1 - pf[3]
    Atm1x[10] = pf[11] * v[9]^-1
    Atm1x[11] = -v[11]
    Atm1x[12] = -v[7]^-1
    Atm1x[13] = pf[8] * v[7]^-1 * (1 - pf[10])
    Atm1x[14] = pf[10] * v[20]^-1
    Atm1x[15] = v[7]^-1
    Atm1x[16] = v[11]
    Atm1 <- sparseMatrix(i = c(5, 8, 9, 12, 15, 17, 20, 20, 21, 23,
                               24, 25, 26, 26, 27, 27),
                         j = c(17, 17, 17, 27, 7, 1, 6, 7, 17, 9,
                               17, 12, 7, 20, 12, 17),
                         x = Atm1x, dims = c(27, 27))

    Atx <- numeric(82)
    Atx[1] = -1
    Atx[2] = -1
    Atx[3] = 1
    Atx[4] = -1
    Atx[5] = pf[7] * (-1 + pf[7]) * v[13]^(-2 + pf[7]) * (1 - v[18])^(1 - pf[7]) * (v[13]^pf[7] * (1 - v[18])^(1 - pf[7]))^(-pf[4]) - pf[4] * pf[7]^2 * v[13]^(-2 + 2 * pf[7]) * (1 - v[18])^(2 - 2 * pf[7]) * (v[13]^pf[7] * (1 - v[18])^(1 - pf[7]))^(-1 - pf[4])
    Atx[6] = pf[7] * (-1 + pf[7]) * v[13]^(-1 + pf[7]) * (1 - v[18])^(-pf[7]) * (v[13]^pf[7] * (1 - v[18])^(1 - pf[7]))^(-pf[4]) - pf[4] * pf[7] * (-1 + pf[7]) * v[13]^(-1 + 2 * pf[7]) * (1 - v[18])^(1 - 2 * pf[7]) * (v[13]^pf[7] * (1 - v[18])^(1 - pf[7]))^(-1 - pf[4])
    Atx[7] = -1
    Atx[8] = pf[1] * v[27] * v[17]^(-1 + pf[1]) * v[18]^(1 - pf[1])
    Atx[9] = -1
    Atx[10] = pf[1] * v[5] * v[27] * (1 - pf[1]) * v[17]^(-1 + pf[1]) * v[18]^(-pf[1])
    Atx[11] = pf[1] * v[5] * v[17]^(-1 + pf[1]) * v[18]^(1 - pf[1])
    Atx[12] = pc[3]
    Atx[13] = -1
    Atx[14] = -v[4]^-2 * v[10]
    Atx[15] = v[4]^-1
    Atx[16] = -1
    Atx[17] = v[27] * (1 - pf[1]) * v[17]^pf[1] * v[18]^(-pf[1])
    Atx[18] = -pf[1] * v[5] * v[27] * (1 - pf[1]) * v[17]^pf[1] * v[18]^(-1 - pf[1])
    Atx[19] = -1
    Atx[20] = v[5] * (1 - pf[1]) * v[17]^pf[1] * v[18]^(-pf[1])
    Atx[21] = v[27] * (1 - pf[1]) * v[17]^pf[1] * v[18]^(-pf[1])
    Atx[22] = -1
    Atx[23] = v[17]^pf[1] * v[18]^(1 - pf[1])
    Atx[24] = 1
    Atx[25] = -1
    Atx[26] = -v[24]
    Atx[27] = -v[6]
    Atx[28] = 1
    Atx[29] = -1
    Atx[30] = -v[20]^-1
    Atx[31] = v[4] * v[20]^-2
    Atx[32] = v[23]
    Atx[33] = pf[7] * (-1 + pf[7]) * v[13]^(-1 + pf[7]) * (1 - v[18])^(-pf[7]) * (v[13]^pf[7] * (1 - v[18])^(1 - pf[7]))^(-pf[4]) - pf[4] * pf[7] * (-1 + pf[7]) * v[13]^(-1 + 2 * pf[7]) * (1 - v[18])^(1 - 2 * pf[7]) * (v[13]^pf[7] * (1 - v[18])^(1 - pf[7]))^(-1 - pf[4])
    Atx[34] = pf[7] * (-1 + pf[7]) * v[13]^pf[7] * (1 - v[18])^(-1 - pf[7]) * (v[13]^pf[7] * (1 - v[18])^(1 - pf[7]))^(-pf[4]) - pf[4] * (-1 + pf[7])^2 * v[13]^(2 * pf[7]) * (1 - v[18])^(-2 * pf[7]) * (v[13]^pf[7] * (1 - v[18])^(1 - pf[7]))^(-1 - pf[4])
    Atx[35] = v[4]
    Atx[36] = pf[6]^-1 * pf[14] * v[7]^-2 * v[7]^pf[5] * (v[7]^-1 * v[7]^pf[5])^(-1 - pf[6]^-1)
    Atx[37] = -pf[6]^-1 * (1 - pf[14]) * v[8]^(-1 - pf[6]^-1)
    Atx[38] = -1
    Atx[39] = 1 + pf[6]
    Atx[40] = -v[1]^-1
    Atx[41] = -1
    Atx[42] = v[8] * v[24]
    Atx[43] = -pf[2] * pf[5] * pf[6]^-1 * pf[14] * v[2] * v[7]^-1 * v[7]^(-1 + pf[5]) * (v[7]^-1 * v[7]^pf[5])^(-1 - pf[6]^-1)
    Atx[44] = v[4] * v[24] + pf[2] * pf[14] * v[2] * v[8]^-1 * (v[7]^-1 * v[7]^pf[5])^(-pf[6]^-1)
    Atx[45] = v[4] * v[8]
    Atx[46] = -1
    Atx[47] = v[5] * v[24]
    Atx[48] = v[4] * v[24]
    Atx[49] = -pf[2] * pf[5] * pf[6]^-1 * pf[14] * v[3] * v[7]^-1 * (1 + pf[6]) * v[7]^(-1 + pf[5]) * (v[7]^-1 * v[7]^pf[5])^(-1 - pf[6]^-1 * (1 + pf[6]))
    Atx[50] = v[4] * v[5]
    Atx[51] = -1
    Atx[52] = pf[6]^-1 * pf[14] * v[6] * v[7]^-2 * (1 + pf[6]) * v[7]^pf[5] * (v[7]^-1 * v[7]^pf[5])^(-1 - pf[6]^-1 * (1 + pf[6]))
    Atx[53] = -pf[6]^-1 * (1 + pf[6]) * (1 - pf[14]) * v[8]^(-1 - pf[6]^-1 * (1 + pf[6]))
    Atx[54] = 1
    Atx[55] = -1
    Atx[56] = -pf[7] * v[13]^(-1 + pf[7]) * (1 - v[18])^(1 - pf[7]) * (v[13]^pf[7] * (1 - v[18])^(1 - pf[7]))^(-pf[4])
    Atx[57] = (1 - pf[7]) * v[13]^pf[7] * (1 - v[18])^(-pf[7]) * (v[13]^pf[7] * (1 - v[18])^(1 - pf[7]))^(-pf[4])
    Atx[58] = 1
    Atx[59] = -v[9]^-1
    Atx[60] = -v[17]
    Atx[61] = -1
    Atx[62] = -v[23]
    Atx[63] = -v[18]
    Atx[64] = 1
    Atx[65] = v[7]^-2 * v[12]
    Atx[66] = v[20]^-1
    Atx[67] = -1
    Atx[68] = -v[12] * v[20]^-2
    Atx[69] = 1
    Atx[70] = (1 - pf[10]) * (v[9]^-1 - pf[8] * v[9]^-1)
    Atx[71] = -v[20]^-1
    Atx[72] = pf[9] * v[24]^-1 * (1 - pf[10])
    Atx[73] = -v[7]^-2 * v[12]
    Atx[74] = v[17]
    Atx[75] = -v[20]^-1
    Atx[76] = -1
    Atx[77] = 1
    Atx[78] = -1
    Atx[79] = v[23]
    Atx[80] = v[12] * v[20]^-2
    Atx[81] = -1
    Atx[82] = v[18]
    At <- sparseMatrix(i = c(1, 2, 2, 3, 3, 3, 4, 5, 5, 5,
                             5, 6, 6, 7, 7, 7, 8, 8, 8, 8,
                             9, 9, 9, 10, 10, 11, 11, 11, 12, 13,
                             13, 14, 14, 14, 14, 15, 15, 16, 16, 17,
                             18, 18, 18, 18, 18, 19, 19, 19, 19, 19,
                             20, 20, 20, 21, 21, 22, 22, 22, 23, 24,
                             24, 24, 24, 24, 25, 25, 25, 25, 25, 26,
                             26, 26, 27, 27, 27, 27, 27, 27, 27, 27,
                             27, 27),
                       j = c(12, 4, 10, 4, 13, 18, 10, 5, 11, 18,
                             27, 1, 15, 4, 10, 19, 5, 18, 23, 27,
                             18, 25, 27, 25, 26, 6, 24, 26, 27, 4,
                             20, 4, 13, 18, 23, 7, 8, 2, 3, 1,
                             2, 4, 7, 8, 24, 3, 4, 5, 7, 24,
                             6, 7, 8, 16, 17, 13, 18, 22, 9, 11,
                             14, 18, 23, 24, 7, 12, 15, 20, 21, 9,
                             20, 24, 7, 11, 12, 13, 14, 16, 18, 20,
                             21, 23),
                         x = Atx, dims = c(27, 27))

    Atp1x <- numeric(11)
    Atp1x[1] = pf[2] * v[11]
    Atp1x[2] = pf[2] * (1 - pf[3])
    Atp1x[3] = pf[2] * v[4]
    Atp1x[4] = pf[2] * v[7]^-1
    Atp1x[5] = -pf[2] * v[4] * v[7]^-2
    Atp1x[6] = pf[2] * pf[14] * (v[7]^-1 * v[7]^pf[5])^(-pf[6]^-1)
    Atp1x[7] = pf[2] * pf[6]^-1 * pf[14] * v[2] * v[7]^-2 * v[7]^pf[5] * (v[7]^-1 * v[7]^pf[5])^(-1 - pf[6]^-1)
    Atp1x[8] = -pf[2] * pf[14] * v[2] * v[8]^-1 * (v[7]^-1 * v[7]^pf[5])^(-pf[6]^-1)
    Atp1x[9] = pf[2] * pf[14] * (v[7]^-1 * v[7]^pf[5])^(-pf[6]^-1 * (1 + pf[6]))
    Atp1x[10] = pf[2] * pf[6]^-1 * pf[14] * v[3] * v[7]^-2 * (1 + pf[6]) * v[7]^pf[5] * (v[7]^-1 * v[7]^pf[5])^(-1 - pf[6]^-1 * (1 + pf[6]))
    Atp1x[11] = -pf[2]
    Atp1 <- sparseMatrix(i = c(4, 4, 4, 13, 13, 18, 18, 18, 19, 19,
                               22),
                         j = c(4, 10, 11, 4, 7, 2, 7, 8, 3, 7,
                               22),
                         x = Atp1x, dims = c(27, 27))

    Aepsx <- numeric(5)
    Aepsx[1] = exp(pf[13] * log(v[27]))
    Aepsx[2] = 1
    Aepsx[3] = 1
    Aepsx[4] = 1
    Aepsx[5] = 1
    Aeps <- sparseMatrix(i = c(12, 16, 17, 23, 26),
                         j = c(1, 2, 5, 4, 3),
                         x = Aepsx, dims = c(27, 5))

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
