# Generated on 2023-04-17 10:10:16 by gEcon ver. 1.2.1 (2023-01-18)
# http://gecon.r-forge.r-project.org/

# Model name: RSWTV

# info
info__ <- c("RSWTV", "C:/Users/fm007/Documents/GitHub/Optimal Gecon/RSWTV.gcn", "2023-04-17 10:10:16", "false")

# index sets
index_sets__ <- list()

# variables
variables__ <- c("iH",
                 "iL",
                 "lambda__HIGHREGIME_1",
                 "lambda__LOWREGIME_1",
                 "lambda__HIGHREGIME_2",
                 "lambda__LOWREGIME_2",
                 "piH",
                 "piL",
                 "rn",
                 "yH",
                 "yL",
                 "UH",
                 "UL")

variables_tex__ <- c("{i\\!H}",
                     "{i\\!L}",
                     "\\lambda^{\\mathrm{HIGHREGIME}^{\\mathrm{1}}}",
                     "\\lambda^{\\mathrm{LOWREGIME}^{\\mathrm{1}}}",
                     "\\lambda^{\\mathrm{HIGHREGIME}^{\\mathrm{2}}}",
                     "\\lambda^{\\mathrm{LOWREGIME}^{\\mathrm{2}}}",
                     "{p\\!i\\!H}",
                     "{p\\!i\\!L}",
                     "{r\\!n}",
                     "{y\\!H}",
                     "{y\\!L}",
                     "{U\\!H}",
                     "{U\\!L}")

# shocks
shocks__ <- c("epsilon_Z")

shocks_tex__ <- c("\\epsilon^{\\mathrm{Z}}")

# parameters
parameters__ <- c("beta",
                  "kappa",
                  "phi",
                  "pHc",
                  "pLc",
                  "sigma",
                  "theta")

parameters_tex__ <- c("\\beta",
                     "\\kappa",
                     "\\phi",
                     "{p\\!H\\!c}",
                     "{p\\!L\\!c}",
                     "\\sigma",
                     "\\theta")

# free parameters
parameters_free__ <- c("beta",
                       "kappa",
                       "phi",
                       "pHc",
                       "pLc",
                       "sigma",
                       "theta")

# free parameters' values
parameters_free_val__ <- c(0.99,
                           0.2465,
                           0.95,
                           0.99,
                           0.99,
                           1,
                           6)

# equations
equations__ <- c("-rn[] + exp(epsilon_Z[] + phi * log(rn[-1])) = 0",
                 "-piH[-1] + beta * (piH[] * (1 + exp(-pHc - piH[] + piL[]))^-1 + piL[] * (1 - (1 + exp(-pHc - piH[] + piL[]))^-1)) + kappa * yH[-1] = 0",
                 "-piL[-1] + beta * (piH[] * (1 - (1 + exp(-pLc + piH[] - piL[]))^-1) + piL[] * (1 + exp(-pLc + piH[] - piL[]))^-1) + kappa * yL[-1] = 0",
                 "lambda__HIGHREGIME_2[] * (1 + exp(-pHc - piH[] + piL[]))^-1 + beta * (1 + exp(-pHc - piH[] + piL[]))^-1 * E[][-lambda__HIGHREGIME_2[1] + kappa * lambda__HIGHREGIME_1[1]] - kappa * theta^-1 * yH[] = 0",
                 "lambda__LOWREGIME_2[] * (1 + exp(-pLc + piH[] - piL[]))^-1 + beta * (1 + exp(-pLc + piH[] - piL[]))^-1 * E[][-lambda__LOWREGIME_2[1] + kappa * lambda__LOWREGIME_1[1]] - kappa * theta^-1 * yL[] = 0",
                 "-yH[-1] - sigma * (iH[-1] - rn[-1] - piH[] * (1 + exp(-pHc - piH[] + piL[]))^-1 - piL[] * (1 - (1 + exp(-pHc - piH[] + piL[]))^-1)) + yH[] * (1 + exp(-pHc - piH[] + piL[]))^-1 + yL[] * (1 - (1 + exp(-pHc - piH[] + piL[]))^-1) = 0",
                 "-yL[-1] - sigma * (iL[-1] - rn[-1] - piH[] * (1 - (1 + exp(-pLc + piH[] - piL[]))^-1) - piL[] * (1 + exp(-pLc + piH[] - piL[]))^-1) + yH[] * (1 - (1 + exp(-pLc + piH[] - piL[]))^-1) + yL[] * (1 + exp(-pLc + piH[] - piL[]))^-1 = 0",
                 "UH[] + 0.5 * piH[]^2 - beta * ((1 + exp(-pHc - piH[] + piL[]))^-1 * E[][UH[1]] + (1 - (1 + exp(-pHc - piH[] + piL[]))^-1) * E[][UL[1]]) + 0.5 * kappa * theta^-1 * yH[]^2 = 0",
                 "UL[] + 0.5 * piL[]^2 - beta * ((1 + exp(-pLc + piH[] - piL[]))^-1 * E[][UL[1]] + (1 - (1 + exp(-pLc + piH[] - piL[]))^-1) * E[][UH[1]]) + 0.5 * kappa * theta^-1 * yL[]^2 = 0",
                 "-piH[] + beta * (exp(-pHc - piH[] + piL[]) * (1 + exp(-pHc - piH[] + piL[]))^-2 * E[][UH[1]] - exp(-pHc - piH[] + piL[]) * (1 + exp(-pHc - piH[] + piL[]))^-2 * E[][UL[1]]) + lambda__HIGHREGIME_2[] * (-sigma * (-(1 + exp(-pHc - piH[] + piL[]))^-1 - piH[] * exp(-pHc - piH[] + piL[]) * (1 + exp(-pHc - piH[] + piL[]))^-2 + piL[] * exp(-pHc - piH[] + piL[]) * (1 + exp(-pHc - piH[] + piL[]))^-2) + yH[] * exp(-pHc - piH[] + piL[]) * (1 + exp(-pHc - piH[] + piL[]))^-2 - yL[] * exp(-pHc - piH[] + piL[]) * (1 + exp(-pHc - piH[] + piL[]))^-2) + beta * lambda__HIGHREGIME_1[] * ((1 + exp(-pHc - piH[] + piL[]))^-1 + piH[] * exp(-pHc - piH[] + piL[]) * (1 + exp(-pHc - piH[] + piL[]))^-2 - piL[] * exp(-pHc - piH[] + piL[]) * (1 + exp(-pHc - piH[] + piL[]))^-2) - beta * (1 + exp(-pHc - piH[] + piL[]))^-1 * E[][lambda__HIGHREGIME_1[1]] = 0",
                 "-piL[] + beta * (exp(-pLc + piH[] - piL[]) * (1 + exp(-pLc + piH[] - piL[]))^-2 * E[][UL[1]] - exp(-pLc + piH[] - piL[]) * (1 + exp(-pLc + piH[] - piL[]))^-2 * E[][UH[1]]) + lambda__LOWREGIME_2[] * (-sigma * (-(1 + exp(-pLc + piH[] - piL[]))^-1 + piH[] * exp(-pLc + piH[] - piL[]) * (1 + exp(-pLc + piH[] - piL[]))^-2 - piL[] * exp(-pLc + piH[] - piL[]) * (1 + exp(-pLc + piH[] - piL[]))^-2) - yH[] * exp(-pLc + piH[] - piL[]) * (1 + exp(-pLc + piH[] - piL[]))^-2 + yL[] * exp(-pLc + piH[] - piL[]) * (1 + exp(-pLc + piH[] - piL[]))^-2) + beta * lambda__LOWREGIME_1[] * ((1 + exp(-pLc + piH[] - piL[]))^-1 + piL[] * exp(-pLc + piH[] - piL[]) * (1 + exp(-pLc + piH[] - piL[]))^-2 - piH[] * exp(-pLc + piH[] - piL[]) * (1 + exp(-pLc + piH[] - piL[]))^-2) - beta * (1 + exp(-pLc + piH[] - piL[]))^-1 * E[][lambda__LOWREGIME_1[1]] = 0",
                 "-beta * sigma * (1 + exp(-pHc - piH[] + piL[]))^-1 * E[][lambda__HIGHREGIME_2[1]] = 0",
                 "-beta * sigma * (1 + exp(-pLc + piH[] - piL[]))^-1 * E[][lambda__LOWREGIME_2[1]] = 0")

# calibrating equations
calibr_equations__ <- character(0)

# variables / equations map
vareqmap__ <- sparseMatrix(i = c(1, 2, 2, 2, 3, 3, 3, 4, 4, 4,
                                 4, 4, 5, 5, 5, 5, 5, 6, 6, 6,
                                 6, 6, 6, 7, 7, 7, 7, 7, 7, 8,
                                 8, 8, 8, 8, 9, 9, 9, 9, 9, 10,
                                 10, 10, 10, 10, 10, 10, 10, 11, 11, 11,
                                 11, 11, 11, 11, 11, 12, 12, 12, 13, 13,
                                 13),
                           j = c(9, 7, 8, 10, 7, 8, 11, 3, 5, 7,
                                 8, 10, 4, 6, 7, 8, 11, 1, 7, 8,
                                 9, 10, 11, 2, 7, 8, 9, 10, 11, 7,
                                 8, 10, 12, 13, 7, 8, 11, 12, 13, 3,
                                 5, 7, 8, 10, 11, 12, 13, 4, 6, 7,
                                 8, 10, 11, 12, 13, 5, 7, 8, 6, 7,
                                 8),
                           x = c(3, 3, 2, 1, 2, 3, 1, 4, 6, 2,
                                 2, 2, 4, 6, 2, 2, 2, 1, 2, 2,
                                 1, 3, 2, 1, 2, 2, 1, 2, 3, 2,
                                 2, 2, 6, 4, 2, 2, 2, 4, 6, 6,
                                 2, 2, 2, 2, 2, 4, 4, 6, 2, 2,
                                 2, 2, 2, 4, 4, 4, 2, 2, 4, 2,
                                 2),
                           dims = c(13, 13))

# variables / calibrating equations map
varcalibreqmap__ <- sparseMatrix(i = NULL, j = NULL, dims = c(0, 13))

# calibrated parameters / equations map
calibrpareqmap__ <- sparseMatrix(i = NULL, j = NULL, dims = c(13, 0))

# calibrated parameters / calibrating equations map
calibrparcalibreqmap__ <- sparseMatrix(i = NULL, j = NULL, dims = c(0, 0))

# free parameters / equations map
freepareqmap__ <- sparseMatrix(i = c(1, 2, 2, 2, 3, 3, 3, 4, 4, 4,
                                     4, 5, 5, 5, 5, 6, 6, 7, 7, 8,
                                     8, 8, 8, 9, 9, 9, 9, 10, 10, 10,
                                     11, 11, 11, 12, 12, 12, 13, 13, 13),
                               j = c(3, 1, 2, 4, 1, 2, 5, 1, 2, 4,
                                     7, 1, 2, 5, 7, 4, 6, 5, 6, 1,
                                     2, 4, 7, 1, 2, 5, 7, 1, 4, 6,
                                     1, 5, 6, 1, 4, 6, 1, 5, 6),
                               x = rep(1, 39), dims = c(13, 7))

# free parameters / calibrating equations map
freeparcalibreqmap__ <- sparseMatrix(i = NULL, j = NULL, dims = c(0, 7))

# shocks / equations map
shockeqmap__ <- sparseMatrix(i = c(1),
                             j = c(1),
                             x = rep(1, 1), dims = c(13, 1))

# steady state equations
ss_eq__ <- function(v, pc, pf)
{
    r <- numeric(13)
    r[1] = -v[9] + exp(pf[3] * log(v[9]))
    r[2] = -v[7] + pf[1] * (v[7] * (1 + exp(-pf[4] - v[7] + v[8]))^-1 + v[8] * (1 - (1 + exp(-pf[4] - v[7] + v[8]))^-1)) + pf[2] * v[10]
    r[3] = -v[8] + pf[1] * (v[7] * (1 - (1 + exp(-pf[5] + v[7] - v[8]))^-1) + v[8] * (1 + exp(-pf[5] + v[7] - v[8]))^-1) + pf[2] * v[11]
    r[4] = v[5] * (1 + exp(-pf[4] - v[7] + v[8]))^-1 + pf[1] * (1 + exp(-pf[4] - v[7] + v[8]))^-1 * (-v[5] + pf[2] * v[3]) - pf[2] * pf[7]^-1 * v[10]
    r[5] = v[6] * (1 + exp(-pf[5] + v[7] - v[8]))^-1 + pf[1] * (1 + exp(-pf[5] + v[7] - v[8]))^-1 * (-v[6] + pf[2] * v[4]) - pf[2] * pf[7]^-1 * v[11]
    r[6] = -v[10] - pf[6] * (v[1] - v[9] - v[7] * (1 + exp(-pf[4] - v[7] + v[8]))^-1 - v[8] * (1 - (1 + exp(-pf[4] - v[7] + v[8]))^-1)) + v[10] * (1 + exp(-pf[4] - v[7] + v[8]))^-1 + v[11] * (1 - (1 + exp(-pf[4] - v[7] + v[8]))^-1)
    r[7] = -v[11] - pf[6] * (v[2] - v[9] - v[7] * (1 - (1 + exp(-pf[5] + v[7] - v[8]))^-1) - v[8] * (1 + exp(-pf[5] + v[7] - v[8]))^-1) + v[10] * (1 - (1 + exp(-pf[5] + v[7] - v[8]))^-1) + v[11] * (1 + exp(-pf[5] + v[7] - v[8]))^-1
    r[8] = v[12] + 0.5 * v[7]^2 - pf[1] * (v[12] * (1 + exp(-pf[4] - v[7] + v[8]))^-1 + v[13] * (1 - (1 + exp(-pf[4] - v[7] + v[8]))^-1)) + 0.5 * pf[2] * pf[7]^-1 * v[10]^2
    r[9] = v[13] + 0.5 * v[8]^2 - pf[1] * (v[12] * (1 - (1 + exp(-pf[5] + v[7] - v[8]))^-1) + v[13] * (1 + exp(-pf[5] + v[7] - v[8]))^-1) + 0.5 * pf[2] * pf[7]^-1 * v[11]^2
    r[10] = -v[7] + pf[1] * (v[12] * exp(-pf[4] - v[7] + v[8]) * (1 + exp(-pf[4] - v[7] + v[8]))^-2 - v[13] * exp(-pf[4] - v[7] + v[8]) * (1 + exp(-pf[4] - v[7] + v[8]))^-2) + v[5] * (-pf[6] * (-(1 + exp(-pf[4] - v[7] + v[8]))^-1 - v[7] * exp(-pf[4] - v[7] + v[8]) * (1 + exp(-pf[4] - v[7] + v[8]))^-2 + v[8] * exp(-pf[4] - v[7] + v[8]) * (1 + exp(-pf[4] - v[7] + v[8]))^-2) + v[10] * exp(-pf[4] - v[7] + v[8]) * (1 + exp(-pf[4] - v[7] + v[8]))^-2 - v[11] * exp(-pf[4] - v[7] + v[8]) * (1 + exp(-pf[4] - v[7] + v[8]))^-2) - pf[1] * v[3] * (1 + exp(-pf[4] - v[7] + v[8]))^-1 + pf[1] * v[3] * ((1 + exp(-pf[4] - v[7] + v[8]))^-1 + v[7] * exp(-pf[4] - v[7] + v[8]) * (1 + exp(-pf[4] - v[7] + v[8]))^-2 - v[8] * exp(-pf[4] - v[7] + v[8]) * (1 + exp(-pf[4] - v[7] + v[8]))^-2)
    r[11] = -v[8] + pf[1] * (v[13] * exp(-pf[5] + v[7] - v[8]) * (1 + exp(-pf[5] + v[7] - v[8]))^-2 - v[12] * exp(-pf[5] + v[7] - v[8]) * (1 + exp(-pf[5] + v[7] - v[8]))^-2) + v[6] * (-pf[6] * (-(1 + exp(-pf[5] + v[7] - v[8]))^-1 + v[7] * exp(-pf[5] + v[7] - v[8]) * (1 + exp(-pf[5] + v[7] - v[8]))^-2 - v[8] * exp(-pf[5] + v[7] - v[8]) * (1 + exp(-pf[5] + v[7] - v[8]))^-2) - v[10] * exp(-pf[5] + v[7] - v[8]) * (1 + exp(-pf[5] + v[7] - v[8]))^-2 + v[11] * exp(-pf[5] + v[7] - v[8]) * (1 + exp(-pf[5] + v[7] - v[8]))^-2) - pf[1] * v[4] * (1 + exp(-pf[5] + v[7] - v[8]))^-1 + pf[1] * v[4] * ((1 + exp(-pf[5] + v[7] - v[8]))^-1 + v[8] * exp(-pf[5] + v[7] - v[8]) * (1 + exp(-pf[5] + v[7] - v[8]))^-2 - v[7] * exp(-pf[5] + v[7] - v[8]) * (1 + exp(-pf[5] + v[7] - v[8]))^-2)
    r[12] = -pf[1] * pf[6] * v[5] * (1 + exp(-pf[4] - v[7] + v[8]))^-1
    r[13] = -pf[1] * pf[6] * v[6] * (1 + exp(-pf[5] + v[7] - v[8]))^-1

    return(r)
}

# calibrating equations
calibr_eq__ <- function(v, pc, pf)
{
    r <- numeric(0)

    return(r)
}

# steady state and calibrating equations Jacobian
ss_calibr_eq_jacob__ <- function(v, pc, pf)
{
    r <- numeric(0)
    jac <- numeric(61)
    jac[1] = -1 + pf[3] * v[9]^-1 * exp(pf[3] * log(v[9]))
    jac[2] = -1 + pf[1] * ((1 + exp(-pf[4] - v[7] + v[8]))^-1 + v[7] * exp(-pf[4] - v[7] + v[8]) * (1 + exp(-pf[4] - v[7] + v[8]))^-2 - v[8] * exp(-pf[4] - v[7] + v[8]) * (1 + exp(-pf[4] - v[7] + v[8]))^-2)
    jac[3] = pf[1] * (1 - (1 + exp(-pf[4] - v[7] + v[8]))^-1 + v[8] * exp(-pf[4] - v[7] + v[8]) * (1 + exp(-pf[4] - v[7] + v[8]))^-2 - v[7] * exp(-pf[4] - v[7] + v[8]) * (1 + exp(-pf[4] - v[7] + v[8]))^-2)
    jac[4] = pf[2]
    jac[5] = pf[1] * (1 - (1 + exp(-pf[5] + v[7] - v[8]))^-1 + v[7] * exp(-pf[5] + v[7] - v[8]) * (1 + exp(-pf[5] + v[7] - v[8]))^-2 - v[8] * exp(-pf[5] + v[7] - v[8]) * (1 + exp(-pf[5] + v[7] - v[8]))^-2)
    jac[6] = -1 + pf[1] * ((1 + exp(-pf[5] + v[7] - v[8]))^-1 + v[8] * exp(-pf[5] + v[7] - v[8]) * (1 + exp(-pf[5] + v[7] - v[8]))^-2 - v[7] * exp(-pf[5] + v[7] - v[8]) * (1 + exp(-pf[5] + v[7] - v[8]))^-2)
    jac[7] = pf[2]
    jac[8] = pf[1] * pf[2] * (1 + exp(-pf[4] - v[7] + v[8]))^-1
    jac[9] = (1 + exp(-pf[4] - v[7] + v[8]))^-1 - pf[1] * (1 + exp(-pf[4] - v[7] + v[8]))^-1
    jac[10] = v[5] * exp(-pf[4] - v[7] + v[8]) * (1 + exp(-pf[4] - v[7] + v[8]))^-2 + pf[1] * exp(-pf[4] - v[7] + v[8]) * (1 + exp(-pf[4] - v[7] + v[8]))^-2 * (-v[5] + pf[2] * v[3])
    jac[11] = -v[5] * exp(-pf[4] - v[7] + v[8]) * (1 + exp(-pf[4] - v[7] + v[8]))^-2 - pf[1] * exp(-pf[4] - v[7] + v[8]) * (1 + exp(-pf[4] - v[7] + v[8]))^-2 * (-v[5] + pf[2] * v[3])
    jac[12] = -pf[2] * pf[7]^-1
    jac[13] = pf[1] * pf[2] * (1 + exp(-pf[5] + v[7] - v[8]))^-1
    jac[14] = (1 + exp(-pf[5] + v[7] - v[8]))^-1 - pf[1] * (1 + exp(-pf[5] + v[7] - v[8]))^-1
    jac[15] = -v[6] * exp(-pf[5] + v[7] - v[8]) * (1 + exp(-pf[5] + v[7] - v[8]))^-2 - pf[1] * exp(-pf[5] + v[7] - v[8]) * (1 + exp(-pf[5] + v[7] - v[8]))^-2 * (-v[6] + pf[2] * v[4])
    jac[16] = v[6] * exp(-pf[5] + v[7] - v[8]) * (1 + exp(-pf[5] + v[7] - v[8]))^-2 + pf[1] * exp(-pf[5] + v[7] - v[8]) * (1 + exp(-pf[5] + v[7] - v[8]))^-2 * (-v[6] + pf[2] * v[4])
    jac[17] = -pf[2] * pf[7]^-1
    jac[18] = -pf[6]
    jac[19] = -pf[6] * (-(1 + exp(-pf[4] - v[7] + v[8]))^-1 - v[7] * exp(-pf[4] - v[7] + v[8]) * (1 + exp(-pf[4] - v[7] + v[8]))^-2 + v[8] * exp(-pf[4] - v[7] + v[8]) * (1 + exp(-pf[4] - v[7] + v[8]))^-2) + v[10] * exp(-pf[4] - v[7] + v[8]) * (1 + exp(-pf[4] - v[7] + v[8]))^-2 - v[11] * exp(-pf[4] - v[7] + v[8]) * (1 + exp(-pf[4] - v[7] + v[8]))^-2
    jac[20] = -pf[6] * (-1 + (1 + exp(-pf[4] - v[7] + v[8]))^-1 + v[7] * exp(-pf[4] - v[7] + v[8]) * (1 + exp(-pf[4] - v[7] + v[8]))^-2 - v[8] * exp(-pf[4] - v[7] + v[8]) * (1 + exp(-pf[4] - v[7] + v[8]))^-2) - v[10] * exp(-pf[4] - v[7] + v[8]) * (1 + exp(-pf[4] - v[7] + v[8]))^-2 + v[11] * exp(-pf[4] - v[7] + v[8]) * (1 + exp(-pf[4] - v[7] + v[8]))^-2
    jac[21] = pf[6]
    jac[22] = -1 + (1 + exp(-pf[4] - v[7] + v[8]))^-1
    jac[23] = 1 - (1 + exp(-pf[4] - v[7] + v[8]))^-1
    jac[24] = -pf[6]
    jac[25] = -pf[6] * (-1 + (1 + exp(-pf[5] + v[7] - v[8]))^-1 - v[7] * exp(-pf[5] + v[7] - v[8]) * (1 + exp(-pf[5] + v[7] - v[8]))^-2 + v[8] * exp(-pf[5] + v[7] - v[8]) * (1 + exp(-pf[5] + v[7] - v[8]))^-2) + v[10] * exp(-pf[5] + v[7] - v[8]) * (1 + exp(-pf[5] + v[7] - v[8]))^-2 - v[11] * exp(-pf[5] + v[7] - v[8]) * (1 + exp(-pf[5] + v[7] - v[8]))^-2
    jac[26] = -pf[6] * (-(1 + exp(-pf[5] + v[7] - v[8]))^-1 + v[7] * exp(-pf[5] + v[7] - v[8]) * (1 + exp(-pf[5] + v[7] - v[8]))^-2 - v[8] * exp(-pf[5] + v[7] - v[8]) * (1 + exp(-pf[5] + v[7] - v[8]))^-2) - v[10] * exp(-pf[5] + v[7] - v[8]) * (1 + exp(-pf[5] + v[7] - v[8]))^-2 + v[11] * exp(-pf[5] + v[7] - v[8]) * (1 + exp(-pf[5] + v[7] - v[8]))^-2
    jac[27] = pf[6]
    jac[28] = 1 - (1 + exp(-pf[5] + v[7] - v[8]))^-1
    jac[29] = -1 + (1 + exp(-pf[5] + v[7] - v[8]))^-1
    jac[30] = v[7] - pf[1] * (v[12] * exp(-pf[4] - v[7] + v[8]) * (1 + exp(-pf[4] - v[7] + v[8]))^-2 - v[13] * exp(-pf[4] - v[7] + v[8]) * (1 + exp(-pf[4] - v[7] + v[8]))^-2)
    jac[31] = -pf[1] * (v[13] * exp(-pf[4] - v[7] + v[8]) * (1 + exp(-pf[4] - v[7] + v[8]))^-2 - v[12] * exp(-pf[4] - v[7] + v[8]) * (1 + exp(-pf[4] - v[7] + v[8]))^-2)
    jac[32] = pf[2] * pf[7]^-1 * v[10]
    jac[33] = 1 - pf[1] * (1 + exp(-pf[4] - v[7] + v[8]))^-1
    jac[34] = -pf[1] * (1 - (1 + exp(-pf[4] - v[7] + v[8]))^-1)
    jac[35] = -pf[1] * (v[12] * exp(-pf[5] + v[7] - v[8]) * (1 + exp(-pf[5] + v[7] - v[8]))^-2 - v[13] * exp(-pf[5] + v[7] - v[8]) * (1 + exp(-pf[5] + v[7] - v[8]))^-2)
    jac[36] = v[8] - pf[1] * (v[13] * exp(-pf[5] + v[7] - v[8]) * (1 + exp(-pf[5] + v[7] - v[8]))^-2 - v[12] * exp(-pf[5] + v[7] - v[8]) * (1 + exp(-pf[5] + v[7] - v[8]))^-2)
    jac[37] = pf[2] * pf[7]^-1 * v[11]
    jac[38] = -pf[1] * (1 - (1 + exp(-pf[5] + v[7] - v[8]))^-1)
    jac[39] = 1 - pf[1] * (1 + exp(-pf[5] + v[7] - v[8]))^-1
    jac[40] = pf[1] * ((1 + exp(-pf[4] - v[7] + v[8]))^-1 + v[7] * exp(-pf[4] - v[7] + v[8]) * (1 + exp(-pf[4] - v[7] + v[8]))^-2 - v[8] * exp(-pf[4] - v[7] + v[8]) * (1 + exp(-pf[4] - v[7] + v[8]))^-2) - pf[1] * (1 + exp(-pf[4] - v[7] + v[8]))^-1
    jac[41] = -pf[6] * (-(1 + exp(-pf[4] - v[7] + v[8]))^-1 - v[7] * exp(-pf[4] - v[7] + v[8]) * (1 + exp(-pf[4] - v[7] + v[8]))^-2 + v[8] * exp(-pf[4] - v[7] + v[8]) * (1 + exp(-pf[4] - v[7] + v[8]))^-2) + v[10] * exp(-pf[4] - v[7] + v[8]) * (1 + exp(-pf[4] - v[7] + v[8]))^-2 - v[11] * exp(-pf[4] - v[7] + v[8]) * (1 + exp(-pf[4] - v[7] + v[8]))^-2
    jac[42] = -1 + pf[1] * (-v[12] * exp(-pf[4] - v[7] + v[8]) * (1 + exp(-pf[4] - v[7] + v[8]))^-2 + 2 * v[12] * exp(-pf[4] - v[7] + v[8])^2 * (1 + exp(-pf[4] - v[7] + v[8]))^-3 + v[13] * exp(-pf[4] - v[7] + v[8]) * (1 + exp(-pf[4] - v[7] + v[8]))^-2 - 2 * v[13] * exp(-pf[4] - v[7] + v[8])^2 * (1 + exp(-pf[4] - v[7] + v[8]))^-3) + v[5] * (-pf[6] * (-2 * exp(-pf[4] - v[7] + v[8]) * (1 + exp(-pf[4] - v[7] + v[8]))^-2 + v[7] * exp(-pf[4] - v[7] + v[8]) * (1 + exp(-pf[4] - v[7] + v[8]))^-2 - 2 * v[7] * exp(-pf[4] - v[7] + v[8])^2 * (1 + exp(-pf[4] - v[7] + v[8]))^-3 - v[8] * exp(-pf[4] - v[7] + v[8]) * (1 + exp(-pf[4] - v[7] + v[8]))^-2 + 2 * v[8] * exp(-pf[4] - v[7] + v[8])^2 * (1 + exp(-pf[4] - v[7] + v[8]))^-3) - v[10] * exp(-pf[4] - v[7] + v[8]) * (1 + exp(-pf[4] - v[7] + v[8]))^-2 + 2 * v[10] * exp(-pf[4] - v[7] + v[8])^2 * (1 + exp(-pf[4] - v[7] + v[8]))^-3 + v[11] * exp(-pf[4] - v[7] + v[8]) * (1 + exp(-pf[4] - v[7] + v[8]))^-2 - 2 * v[11] * exp(-pf[4] - v[7] + v[8])^2 * (1 + exp(-pf[4] - v[7] + v[8]))^-3) + pf[1] * v[3] * (2 * exp(-pf[4] - v[7] + v[8]) * (1 + exp(-pf[4] - v[7] + v[8]))^-2 - v[7] * exp(-pf[4] - v[7] + v[8]) * (1 + exp(-pf[4] - v[7] + v[8]))^-2 + 2 * v[7] * exp(-pf[4] - v[7] + v[8])^2 * (1 + exp(-pf[4] - v[7] + v[8]))^-3 + v[8] * exp(-pf[4] - v[7] + v[8]) * (1 + exp(-pf[4] - v[7] + v[8]))^-2 - 2 * v[8] * exp(-pf[4] - v[7] + v[8])^2 * (1 + exp(-pf[4] - v[7] + v[8]))^-3) - pf[1] * v[3] * exp(-pf[4] - v[7] + v[8]) * (1 + exp(-pf[4] - v[7] + v[8]))^-2
    jac[43] = pf[1] * (v[12] * exp(-pf[4] - v[7] + v[8]) * (1 + exp(-pf[4] - v[7] + v[8]))^-2 - 2 * v[12] * exp(-pf[4] - v[7] + v[8])^2 * (1 + exp(-pf[4] - v[7] + v[8]))^-3 - v[13] * exp(-pf[4] - v[7] + v[8]) * (1 + exp(-pf[4] - v[7] + v[8]))^-2 + 2 * v[13] * exp(-pf[4] - v[7] + v[8])^2 * (1 + exp(-pf[4] - v[7] + v[8]))^-3) + v[5] * (-pf[6] * (2 * exp(-pf[4] - v[7] + v[8]) * (1 + exp(-pf[4] - v[7] + v[8]))^-2 - v[7] * exp(-pf[4] - v[7] + v[8]) * (1 + exp(-pf[4] - v[7] + v[8]))^-2 + 2 * v[7] * exp(-pf[4] - v[7] + v[8])^2 * (1 + exp(-pf[4] - v[7] + v[8]))^-3 + v[8] * exp(-pf[4] - v[7] + v[8]) * (1 + exp(-pf[4] - v[7] + v[8]))^-2 - 2 * v[8] * exp(-pf[4] - v[7] + v[8])^2 * (1 + exp(-pf[4] - v[7] + v[8]))^-3) + v[10] * exp(-pf[4] - v[7] + v[8]) * (1 + exp(-pf[4] - v[7] + v[8]))^-2 - 2 * v[10] * exp(-pf[4] - v[7] + v[8])^2 * (1 + exp(-pf[4] - v[7] + v[8]))^-3 - v[11] * exp(-pf[4] - v[7] + v[8]) * (1 + exp(-pf[4] - v[7] + v[8]))^-2 + 2 * v[11] * exp(-pf[4] - v[7] + v[8])^2 * (1 + exp(-pf[4] - v[7] + v[8]))^-3) + pf[1] * v[3] * (-2 * exp(-pf[4] - v[7] + v[8]) * (1 + exp(-pf[4] - v[7] + v[8]))^-2 + v[7] * exp(-pf[4] - v[7] + v[8]) * (1 + exp(-pf[4] - v[7] + v[8]))^-2 - 2 * v[7] * exp(-pf[4] - v[7] + v[8])^2 * (1 + exp(-pf[4] - v[7] + v[8]))^-3 - v[8] * exp(-pf[4] - v[7] + v[8]) * (1 + exp(-pf[4] - v[7] + v[8]))^-2 + 2 * v[8] * exp(-pf[4] - v[7] + v[8])^2 * (1 + exp(-pf[4] - v[7] + v[8]))^-3) + pf[1] * v[3] * exp(-pf[4] - v[7] + v[8]) * (1 + exp(-pf[4] - v[7] + v[8]))^-2
    jac[44] = v[5] * exp(-pf[4] - v[7] + v[8]) * (1 + exp(-pf[4] - v[7] + v[8]))^-2
    jac[45] = -v[5] * exp(-pf[4] - v[7] + v[8]) * (1 + exp(-pf[4] - v[7] + v[8]))^-2
    jac[46] = pf[1] * exp(-pf[4] - v[7] + v[8]) * (1 + exp(-pf[4] - v[7] + v[8]))^-2
    jac[47] = -pf[1] * exp(-pf[4] - v[7] + v[8]) * (1 + exp(-pf[4] - v[7] + v[8]))^-2
    jac[48] = pf[1] * ((1 + exp(-pf[5] + v[7] - v[8]))^-1 + v[8] * exp(-pf[5] + v[7] - v[8]) * (1 + exp(-pf[5] + v[7] - v[8]))^-2 - v[7] * exp(-pf[5] + v[7] - v[8]) * (1 + exp(-pf[5] + v[7] - v[8]))^-2) - pf[1] * (1 + exp(-pf[5] + v[7] - v[8]))^-1
    jac[49] = -pf[6] * (-(1 + exp(-pf[5] + v[7] - v[8]))^-1 + v[7] * exp(-pf[5] + v[7] - v[8]) * (1 + exp(-pf[5] + v[7] - v[8]))^-2 - v[8] * exp(-pf[5] + v[7] - v[8]) * (1 + exp(-pf[5] + v[7] - v[8]))^-2) - v[10] * exp(-pf[5] + v[7] - v[8]) * (1 + exp(-pf[5] + v[7] - v[8]))^-2 + v[11] * exp(-pf[5] + v[7] - v[8]) * (1 + exp(-pf[5] + v[7] - v[8]))^-2
    jac[50] = pf[1] * (-v[12] * exp(-pf[5] + v[7] - v[8]) * (1 + exp(-pf[5] + v[7] - v[8]))^-2 + 2 * v[12] * exp(-pf[5] + v[7] - v[8])^2 * (1 + exp(-pf[5] + v[7] - v[8]))^-3 + v[13] * exp(-pf[5] + v[7] - v[8]) * (1 + exp(-pf[5] + v[7] - v[8]))^-2 - 2 * v[13] * exp(-pf[5] + v[7] - v[8])^2 * (1 + exp(-pf[5] + v[7] - v[8]))^-3) + v[6] * (-pf[6] * (2 * exp(-pf[5] + v[7] - v[8]) * (1 + exp(-pf[5] + v[7] - v[8]))^-2 + v[7] * exp(-pf[5] + v[7] - v[8]) * (1 + exp(-pf[5] + v[7] - v[8]))^-2 - 2 * v[7] * exp(-pf[5] + v[7] - v[8])^2 * (1 + exp(-pf[5] + v[7] - v[8]))^-3 - v[8] * exp(-pf[5] + v[7] - v[8]) * (1 + exp(-pf[5] + v[7] - v[8]))^-2 + 2 * v[8] * exp(-pf[5] + v[7] - v[8])^2 * (1 + exp(-pf[5] + v[7] - v[8]))^-3) - v[10] * exp(-pf[5] + v[7] - v[8]) * (1 + exp(-pf[5] + v[7] - v[8]))^-2 + 2 * v[10] * exp(-pf[5] + v[7] - v[8])^2 * (1 + exp(-pf[5] + v[7] - v[8]))^-3 + v[11] * exp(-pf[5] + v[7] - v[8]) * (1 + exp(-pf[5] + v[7] - v[8]))^-2 - 2 * v[11] * exp(-pf[5] + v[7] - v[8])^2 * (1 + exp(-pf[5] + v[7] - v[8]))^-3) + pf[1] * v[4] * (-2 * exp(-pf[5] + v[7] - v[8]) * (1 + exp(-pf[5] + v[7] - v[8]))^-2 - v[7] * exp(-pf[5] + v[7] - v[8]) * (1 + exp(-pf[5] + v[7] - v[8]))^-2 + 2 * v[7] * exp(-pf[5] + v[7] - v[8])^2 * (1 + exp(-pf[5] + v[7] - v[8]))^-3 + v[8] * exp(-pf[5] + v[7] - v[8]) * (1 + exp(-pf[5] + v[7] - v[8]))^-2 - 2 * v[8] * exp(-pf[5] + v[7] - v[8])^2 * (1 + exp(-pf[5] + v[7] - v[8]))^-3) + pf[1] * v[4] * exp(-pf[5] + v[7] - v[8]) * (1 + exp(-pf[5] + v[7] - v[8]))^-2
    jac[51] = -1 + pf[1] * (v[12] * exp(-pf[5] + v[7] - v[8]) * (1 + exp(-pf[5] + v[7] - v[8]))^-2 - 2 * v[12] * exp(-pf[5] + v[7] - v[8])^2 * (1 + exp(-pf[5] + v[7] - v[8]))^-3 - v[13] * exp(-pf[5] + v[7] - v[8]) * (1 + exp(-pf[5] + v[7] - v[8]))^-2 + 2 * v[13] * exp(-pf[5] + v[7] - v[8])^2 * (1 + exp(-pf[5] + v[7] - v[8]))^-3) + v[6] * (-pf[6] * (-2 * exp(-pf[5] + v[7] - v[8]) * (1 + exp(-pf[5] + v[7] - v[8]))^-2 - v[7] * exp(-pf[5] + v[7] - v[8]) * (1 + exp(-pf[5] + v[7] - v[8]))^-2 + 2 * v[7] * exp(-pf[5] + v[7] - v[8])^2 * (1 + exp(-pf[5] + v[7] - v[8]))^-3 + v[8] * exp(-pf[5] + v[7] - v[8]) * (1 + exp(-pf[5] + v[7] - v[8]))^-2 - 2 * v[8] * exp(-pf[5] + v[7] - v[8])^2 * (1 + exp(-pf[5] + v[7] - v[8]))^-3) + v[10] * exp(-pf[5] + v[7] - v[8]) * (1 + exp(-pf[5] + v[7] - v[8]))^-2 - 2 * v[10] * exp(-pf[5] + v[7] - v[8])^2 * (1 + exp(-pf[5] + v[7] - v[8]))^-3 - v[11] * exp(-pf[5] + v[7] - v[8]) * (1 + exp(-pf[5] + v[7] - v[8]))^-2 + 2 * v[11] * exp(-pf[5] + v[7] - v[8])^2 * (1 + exp(-pf[5] + v[7] - v[8]))^-3) + pf[1] * v[4] * (2 * exp(-pf[5] + v[7] - v[8]) * (1 + exp(-pf[5] + v[7] - v[8]))^-2 + v[7] * exp(-pf[5] + v[7] - v[8]) * (1 + exp(-pf[5] + v[7] - v[8]))^-2 - 2 * v[7] * exp(-pf[5] + v[7] - v[8])^2 * (1 + exp(-pf[5] + v[7] - v[8]))^-3 - v[8] * exp(-pf[5] + v[7] - v[8]) * (1 + exp(-pf[5] + v[7] - v[8]))^-2 + 2 * v[8] * exp(-pf[5] + v[7] - v[8])^2 * (1 + exp(-pf[5] + v[7] - v[8]))^-3) - pf[1] * v[4] * exp(-pf[5] + v[7] - v[8]) * (1 + exp(-pf[5] + v[7] - v[8]))^-2
    jac[52] = -v[6] * exp(-pf[5] + v[7] - v[8]) * (1 + exp(-pf[5] + v[7] - v[8]))^-2
    jac[53] = v[6] * exp(-pf[5] + v[7] - v[8]) * (1 + exp(-pf[5] + v[7] - v[8]))^-2
    jac[54] = -pf[1] * exp(-pf[5] + v[7] - v[8]) * (1 + exp(-pf[5] + v[7] - v[8]))^-2
    jac[55] = pf[1] * exp(-pf[5] + v[7] - v[8]) * (1 + exp(-pf[5] + v[7] - v[8]))^-2
    jac[56] = -pf[1] * pf[6] * (1 + exp(-pf[4] - v[7] + v[8]))^-1
    jac[57] = -pf[1] * pf[6] * v[5] * exp(-pf[4] - v[7] + v[8]) * (1 + exp(-pf[4] - v[7] + v[8]))^-2
    jac[58] = pf[1] * pf[6] * v[5] * exp(-pf[4] - v[7] + v[8]) * (1 + exp(-pf[4] - v[7] + v[8]))^-2
    jac[59] = -pf[1] * pf[6] * (1 + exp(-pf[5] + v[7] - v[8]))^-1
    jac[60] = pf[1] * pf[6] * v[6] * exp(-pf[5] + v[7] - v[8]) * (1 + exp(-pf[5] + v[7] - v[8]))^-2
    jac[61] = -pf[1] * pf[6] * v[6] * exp(-pf[5] + v[7] - v[8]) * (1 + exp(-pf[5] + v[7] - v[8]))^-2
    jacob <- sparseMatrix(i = c(1, 2, 2, 2, 3, 3, 3, 4, 4, 4,
                                4, 4, 5, 5, 5, 5, 5, 6, 6, 6,
                                6, 6, 6, 7, 7, 7, 7, 7, 7, 8,
                                8, 8, 8, 8, 9, 9, 9, 9, 9, 10,
                                10, 10, 10, 10, 10, 10, 10, 11, 11, 11,
                                11, 11, 11, 11, 11, 12, 12, 12, 13, 13,
                                13),
                          j = c(9, 7, 8, 10, 7, 8, 11, 3, 5, 7,
                                8, 10, 4, 6, 7, 8, 11, 1, 7, 8,
                                9, 10, 11, 2, 7, 8, 9, 10, 11, 7,
                                8, 10, 12, 13, 7, 8, 11, 12, 13, 3,
                                5, 7, 8, 10, 11, 12, 13, 4, 6, 7,
                                8, 10, 11, 12, 13, 5, 7, 8, 6, 7,
                                8),
                          x = jac, dims = c(13, 13))

    return(jacob)
}

# 1st order perturbation
pert1__ <- function(v, pc, pf)
{
    Atm1x <- numeric(11)
    Atm1x[1] = pf[3] * v[9]^-1 * exp(pf[3] * log(v[9]))
    Atm1x[2] = -1
    Atm1x[3] = pf[2]
    Atm1x[4] = -1
    Atm1x[5] = pf[2]
    Atm1x[6] = -pf[6]
    Atm1x[7] = pf[6]
    Atm1x[8] = -1
    Atm1x[9] = -pf[6]
    Atm1x[10] = pf[6]
    Atm1x[11] = -1
    Atm1 <- sparseMatrix(i = c(1, 2, 2, 3, 3, 6, 6, 6, 7, 7,
                               7),
                         j = c(9, 7, 10, 8, 11, 1, 9, 10, 2, 9,
                               11),
                         x = Atm1x, dims = c(13, 13))

    Atx <- numeric(45)
    Atx[1] = -1
    Atx[2] = pf[1] * ((1 + exp(-pf[4] - v[7] + v[8]))^-1 + v[7] * exp(-pf[4] - v[7] + v[8]) * (1 + exp(-pf[4] - v[7] + v[8]))^-2 - v[8] * exp(-pf[4] - v[7] + v[8]) * (1 + exp(-pf[4] - v[7] + v[8]))^-2)
    Atx[3] = pf[1] * (1 - (1 + exp(-pf[4] - v[7] + v[8]))^-1 - v[7] * exp(-pf[4] - v[7] + v[8]) * (1 + exp(-pf[4] - v[7] + v[8]))^-2 + v[8] * exp(-pf[4] - v[7] + v[8]) * (1 + exp(-pf[4] - v[7] + v[8]))^-2)
    Atx[4] = pf[1] * (1 - (1 + exp(-pf[5] + v[7] - v[8]))^-1 + v[7] * exp(-pf[5] + v[7] - v[8]) * (1 + exp(-pf[5] + v[7] - v[8]))^-2 - v[8] * exp(-pf[5] + v[7] - v[8]) * (1 + exp(-pf[5] + v[7] - v[8]))^-2)
    Atx[5] = pf[1] * ((1 + exp(-pf[5] + v[7] - v[8]))^-1 + v[8] * exp(-pf[5] + v[7] - v[8]) * (1 + exp(-pf[5] + v[7] - v[8]))^-2 - v[7] * exp(-pf[5] + v[7] - v[8]) * (1 + exp(-pf[5] + v[7] - v[8]))^-2)
    Atx[6] = (1 + exp(-pf[4] - v[7] + v[8]))^-1
    Atx[7] = v[5] * exp(-pf[4] - v[7] + v[8]) * (1 + exp(-pf[4] - v[7] + v[8]))^-2 + pf[1] * exp(-pf[4] - v[7] + v[8]) * (1 + exp(-pf[4] - v[7] + v[8]))^-2 * (-v[5] + pf[2] * v[3])
    Atx[8] = -v[5] * exp(-pf[4] - v[7] + v[8]) * (1 + exp(-pf[4] - v[7] + v[8]))^-2 - pf[1] * exp(-pf[4] - v[7] + v[8]) * (1 + exp(-pf[4] - v[7] + v[8]))^-2 * (-v[5] + pf[2] * v[3])
    Atx[9] = -pf[2] * pf[7]^-1
    Atx[10] = (1 + exp(-pf[5] + v[7] - v[8]))^-1
    Atx[11] = -v[6] * exp(-pf[5] + v[7] - v[8]) * (1 + exp(-pf[5] + v[7] - v[8]))^-2 - pf[1] * exp(-pf[5] + v[7] - v[8]) * (1 + exp(-pf[5] + v[7] - v[8]))^-2 * (-v[6] + pf[2] * v[4])
    Atx[12] = v[6] * exp(-pf[5] + v[7] - v[8]) * (1 + exp(-pf[5] + v[7] - v[8]))^-2 + pf[1] * exp(-pf[5] + v[7] - v[8]) * (1 + exp(-pf[5] + v[7] - v[8]))^-2 * (-v[6] + pf[2] * v[4])
    Atx[13] = -pf[2] * pf[7]^-1
    Atx[14] = -pf[6] * (-(1 + exp(-pf[4] - v[7] + v[8]))^-1 - v[7] * exp(-pf[4] - v[7] + v[8]) * (1 + exp(-pf[4] - v[7] + v[8]))^-2 + v[8] * exp(-pf[4] - v[7] + v[8]) * (1 + exp(-pf[4] - v[7] + v[8]))^-2) + v[10] * exp(-pf[4] - v[7] + v[8]) * (1 + exp(-pf[4] - v[7] + v[8]))^-2 - v[11] * exp(-pf[4] - v[7] + v[8]) * (1 + exp(-pf[4] - v[7] + v[8]))^-2
    Atx[15] = -pf[6] * (-1 + (1 + exp(-pf[4] - v[7] + v[8]))^-1 + v[7] * exp(-pf[4] - v[7] + v[8]) * (1 + exp(-pf[4] - v[7] + v[8]))^-2 - v[8] * exp(-pf[4] - v[7] + v[8]) * (1 + exp(-pf[4] - v[7] + v[8]))^-2) - v[10] * exp(-pf[4] - v[7] + v[8]) * (1 + exp(-pf[4] - v[7] + v[8]))^-2 + v[11] * exp(-pf[4] - v[7] + v[8]) * (1 + exp(-pf[4] - v[7] + v[8]))^-2
    Atx[16] = (1 + exp(-pf[4] - v[7] + v[8]))^-1
    Atx[17] = 1 - (1 + exp(-pf[4] - v[7] + v[8]))^-1
    Atx[18] = -pf[6] * (-1 + (1 + exp(-pf[5] + v[7] - v[8]))^-1 - v[7] * exp(-pf[5] + v[7] - v[8]) * (1 + exp(-pf[5] + v[7] - v[8]))^-2 + v[8] * exp(-pf[5] + v[7] - v[8]) * (1 + exp(-pf[5] + v[7] - v[8]))^-2) + v[10] * exp(-pf[5] + v[7] - v[8]) * (1 + exp(-pf[5] + v[7] - v[8]))^-2 - v[11] * exp(-pf[5] + v[7] - v[8]) * (1 + exp(-pf[5] + v[7] - v[8]))^-2
    Atx[19] = -pf[6] * (-(1 + exp(-pf[5] + v[7] - v[8]))^-1 + v[7] * exp(-pf[5] + v[7] - v[8]) * (1 + exp(-pf[5] + v[7] - v[8]))^-2 - v[8] * exp(-pf[5] + v[7] - v[8]) * (1 + exp(-pf[5] + v[7] - v[8]))^-2) - v[10] * exp(-pf[5] + v[7] - v[8]) * (1 + exp(-pf[5] + v[7] - v[8]))^-2 + v[11] * exp(-pf[5] + v[7] - v[8]) * (1 + exp(-pf[5] + v[7] - v[8]))^-2
    Atx[20] = 1 - (1 + exp(-pf[5] + v[7] - v[8]))^-1
    Atx[21] = (1 + exp(-pf[5] + v[7] - v[8]))^-1
    Atx[22] = v[7] - pf[1] * (v[12] * exp(-pf[4] - v[7] + v[8]) * (1 + exp(-pf[4] - v[7] + v[8]))^-2 - v[13] * exp(-pf[4] - v[7] + v[8]) * (1 + exp(-pf[4] - v[7] + v[8]))^-2)
    Atx[23] = -pf[1] * (v[13] * exp(-pf[4] - v[7] + v[8]) * (1 + exp(-pf[4] - v[7] + v[8]))^-2 - v[12] * exp(-pf[4] - v[7] + v[8]) * (1 + exp(-pf[4] - v[7] + v[8]))^-2)
    Atx[24] = pf[2] * pf[7]^-1 * v[10]
    Atx[25] = 1
    Atx[26] = -pf[1] * (v[12] * exp(-pf[5] + v[7] - v[8]) * (1 + exp(-pf[5] + v[7] - v[8]))^-2 - v[13] * exp(-pf[5] + v[7] - v[8]) * (1 + exp(-pf[5] + v[7] - v[8]))^-2)
    Atx[27] = v[8] - pf[1] * (v[13] * exp(-pf[5] + v[7] - v[8]) * (1 + exp(-pf[5] + v[7] - v[8]))^-2 - v[12] * exp(-pf[5] + v[7] - v[8]) * (1 + exp(-pf[5] + v[7] - v[8]))^-2)
    Atx[28] = pf[2] * pf[7]^-1 * v[11]
    Atx[29] = 1
    Atx[30] = pf[1] * ((1 + exp(-pf[4] - v[7] + v[8]))^-1 + v[7] * exp(-pf[4] - v[7] + v[8]) * (1 + exp(-pf[4] - v[7] + v[8]))^-2 - v[8] * exp(-pf[4] - v[7] + v[8]) * (1 + exp(-pf[4] - v[7] + v[8]))^-2)
    Atx[31] = -pf[6] * (-(1 + exp(-pf[4] - v[7] + v[8]))^-1 - v[7] * exp(-pf[4] - v[7] + v[8]) * (1 + exp(-pf[4] - v[7] + v[8]))^-2 + v[8] * exp(-pf[4] - v[7] + v[8]) * (1 + exp(-pf[4] - v[7] + v[8]))^-2) + v[10] * exp(-pf[4] - v[7] + v[8]) * (1 + exp(-pf[4] - v[7] + v[8]))^-2 - v[11] * exp(-pf[4] - v[7] + v[8]) * (1 + exp(-pf[4] - v[7] + v[8]))^-2
    Atx[32] = -1 + pf[1] * (-v[12] * exp(-pf[4] - v[7] + v[8]) * (1 + exp(-pf[4] - v[7] + v[8]))^-2 + 2 * v[12] * exp(-pf[4] - v[7] + v[8])^2 * (1 + exp(-pf[4] - v[7] + v[8]))^-3 + v[13] * exp(-pf[4] - v[7] + v[8]) * (1 + exp(-pf[4] - v[7] + v[8]))^-2 - 2 * v[13] * exp(-pf[4] - v[7] + v[8])^2 * (1 + exp(-pf[4] - v[7] + v[8]))^-3) + v[5] * (-pf[6] * (-2 * exp(-pf[4] - v[7] + v[8]) * (1 + exp(-pf[4] - v[7] + v[8]))^-2 + v[7] * exp(-pf[4] - v[7] + v[8]) * (1 + exp(-pf[4] - v[7] + v[8]))^-2 - 2 * v[7] * exp(-pf[4] - v[7] + v[8])^2 * (1 + exp(-pf[4] - v[7] + v[8]))^-3 - v[8] * exp(-pf[4] - v[7] + v[8]) * (1 + exp(-pf[4] - v[7] + v[8]))^-2 + 2 * v[8] * exp(-pf[4] - v[7] + v[8])^2 * (1 + exp(-pf[4] - v[7] + v[8]))^-3) - v[10] * exp(-pf[4] - v[7] + v[8]) * (1 + exp(-pf[4] - v[7] + v[8]))^-2 + 2 * v[10] * exp(-pf[4] - v[7] + v[8])^2 * (1 + exp(-pf[4] - v[7] + v[8]))^-3 + v[11] * exp(-pf[4] - v[7] + v[8]) * (1 + exp(-pf[4] - v[7] + v[8]))^-2 - 2 * v[11] * exp(-pf[4] - v[7] + v[8])^2 * (1 + exp(-pf[4] - v[7] + v[8]))^-3) + pf[1] * v[3] * (2 * exp(-pf[4] - v[7] + v[8]) * (1 + exp(-pf[4] - v[7] + v[8]))^-2 - v[7] * exp(-pf[4] - v[7] + v[8]) * (1 + exp(-pf[4] - v[7] + v[8]))^-2 + 2 * v[7] * exp(-pf[4] - v[7] + v[8])^2 * (1 + exp(-pf[4] - v[7] + v[8]))^-3 + v[8] * exp(-pf[4] - v[7] + v[8]) * (1 + exp(-pf[4] - v[7] + v[8]))^-2 - 2 * v[8] * exp(-pf[4] - v[7] + v[8])^2 * (1 + exp(-pf[4] - v[7] + v[8]))^-3) - pf[1] * v[3] * exp(-pf[4] - v[7] + v[8]) * (1 + exp(-pf[4] - v[7] + v[8]))^-2
    Atx[33] = pf[1] * (v[12] * exp(-pf[4] - v[7] + v[8]) * (1 + exp(-pf[4] - v[7] + v[8]))^-2 - 2 * v[12] * exp(-pf[4] - v[7] + v[8])^2 * (1 + exp(-pf[4] - v[7] + v[8]))^-3 - v[13] * exp(-pf[4] - v[7] + v[8]) * (1 + exp(-pf[4] - v[7] + v[8]))^-2 + 2 * v[13] * exp(-pf[4] - v[7] + v[8])^2 * (1 + exp(-pf[4] - v[7] + v[8]))^-3) + v[5] * (-pf[6] * (2 * exp(-pf[4] - v[7] + v[8]) * (1 + exp(-pf[4] - v[7] + v[8]))^-2 - v[7] * exp(-pf[4] - v[7] + v[8]) * (1 + exp(-pf[4] - v[7] + v[8]))^-2 + 2 * v[7] * exp(-pf[4] - v[7] + v[8])^2 * (1 + exp(-pf[4] - v[7] + v[8]))^-3 + v[8] * exp(-pf[4] - v[7] + v[8]) * (1 + exp(-pf[4] - v[7] + v[8]))^-2 - 2 * v[8] * exp(-pf[4] - v[7] + v[8])^2 * (1 + exp(-pf[4] - v[7] + v[8]))^-3) + v[10] * exp(-pf[4] - v[7] + v[8]) * (1 + exp(-pf[4] - v[7] + v[8]))^-2 - 2 * v[10] * exp(-pf[4] - v[7] + v[8])^2 * (1 + exp(-pf[4] - v[7] + v[8]))^-3 - v[11] * exp(-pf[4] - v[7] + v[8]) * (1 + exp(-pf[4] - v[7] + v[8]))^-2 + 2 * v[11] * exp(-pf[4] - v[7] + v[8])^2 * (1 + exp(-pf[4] - v[7] + v[8]))^-3) + pf[1] * v[3] * (-2 * exp(-pf[4] - v[7] + v[8]) * (1 + exp(-pf[4] - v[7] + v[8]))^-2 + v[7] * exp(-pf[4] - v[7] + v[8]) * (1 + exp(-pf[4] - v[7] + v[8]))^-2 - 2 * v[7] * exp(-pf[4] - v[7] + v[8])^2 * (1 + exp(-pf[4] - v[7] + v[8]))^-3 - v[8] * exp(-pf[4] - v[7] + v[8]) * (1 + exp(-pf[4] - v[7] + v[8]))^-2 + 2 * v[8] * exp(-pf[4] - v[7] + v[8])^2 * (1 + exp(-pf[4] - v[7] + v[8]))^-3) + pf[1] * v[3] * exp(-pf[4] - v[7] + v[8]) * (1 + exp(-pf[4] - v[7] + v[8]))^-2
    Atx[34] = v[5] * exp(-pf[4] - v[7] + v[8]) * (1 + exp(-pf[4] - v[7] + v[8]))^-2
    Atx[35] = -v[5] * exp(-pf[4] - v[7] + v[8]) * (1 + exp(-pf[4] - v[7] + v[8]))^-2
    Atx[36] = pf[1] * ((1 + exp(-pf[5] + v[7] - v[8]))^-1 + v[8] * exp(-pf[5] + v[7] - v[8]) * (1 + exp(-pf[5] + v[7] - v[8]))^-2 - v[7] * exp(-pf[5] + v[7] - v[8]) * (1 + exp(-pf[5] + v[7] - v[8]))^-2)
    Atx[37] = -pf[6] * (-(1 + exp(-pf[5] + v[7] - v[8]))^-1 + v[7] * exp(-pf[5] + v[7] - v[8]) * (1 + exp(-pf[5] + v[7] - v[8]))^-2 - v[8] * exp(-pf[5] + v[7] - v[8]) * (1 + exp(-pf[5] + v[7] - v[8]))^-2) - v[10] * exp(-pf[5] + v[7] - v[8]) * (1 + exp(-pf[5] + v[7] - v[8]))^-2 + v[11] * exp(-pf[5] + v[7] - v[8]) * (1 + exp(-pf[5] + v[7] - v[8]))^-2
    Atx[38] = pf[1] * (-v[12] * exp(-pf[5] + v[7] - v[8]) * (1 + exp(-pf[5] + v[7] - v[8]))^-2 + 2 * v[12] * exp(-pf[5] + v[7] - v[8])^2 * (1 + exp(-pf[5] + v[7] - v[8]))^-3 + v[13] * exp(-pf[5] + v[7] - v[8]) * (1 + exp(-pf[5] + v[7] - v[8]))^-2 - 2 * v[13] * exp(-pf[5] + v[7] - v[8])^2 * (1 + exp(-pf[5] + v[7] - v[8]))^-3) + v[6] * (-pf[6] * (2 * exp(-pf[5] + v[7] - v[8]) * (1 + exp(-pf[5] + v[7] - v[8]))^-2 + v[7] * exp(-pf[5] + v[7] - v[8]) * (1 + exp(-pf[5] + v[7] - v[8]))^-2 - 2 * v[7] * exp(-pf[5] + v[7] - v[8])^2 * (1 + exp(-pf[5] + v[7] - v[8]))^-3 - v[8] * exp(-pf[5] + v[7] - v[8]) * (1 + exp(-pf[5] + v[7] - v[8]))^-2 + 2 * v[8] * exp(-pf[5] + v[7] - v[8])^2 * (1 + exp(-pf[5] + v[7] - v[8]))^-3) - v[10] * exp(-pf[5] + v[7] - v[8]) * (1 + exp(-pf[5] + v[7] - v[8]))^-2 + 2 * v[10] * exp(-pf[5] + v[7] - v[8])^2 * (1 + exp(-pf[5] + v[7] - v[8]))^-3 + v[11] * exp(-pf[5] + v[7] - v[8]) * (1 + exp(-pf[5] + v[7] - v[8]))^-2 - 2 * v[11] * exp(-pf[5] + v[7] - v[8])^2 * (1 + exp(-pf[5] + v[7] - v[8]))^-3) + pf[1] * v[4] * (-2 * exp(-pf[5] + v[7] - v[8]) * (1 + exp(-pf[5] + v[7] - v[8]))^-2 - v[7] * exp(-pf[5] + v[7] - v[8]) * (1 + exp(-pf[5] + v[7] - v[8]))^-2 + 2 * v[7] * exp(-pf[5] + v[7] - v[8])^2 * (1 + exp(-pf[5] + v[7] - v[8]))^-3 + v[8] * exp(-pf[5] + v[7] - v[8]) * (1 + exp(-pf[5] + v[7] - v[8]))^-2 - 2 * v[8] * exp(-pf[5] + v[7] - v[8])^2 * (1 + exp(-pf[5] + v[7] - v[8]))^-3) + pf[1] * v[4] * exp(-pf[5] + v[7] - v[8]) * (1 + exp(-pf[5] + v[7] - v[8]))^-2
    Atx[39] = -1 + pf[1] * (v[12] * exp(-pf[5] + v[7] - v[8]) * (1 + exp(-pf[5] + v[7] - v[8]))^-2 - 2 * v[12] * exp(-pf[5] + v[7] - v[8])^2 * (1 + exp(-pf[5] + v[7] - v[8]))^-3 - v[13] * exp(-pf[5] + v[7] - v[8]) * (1 + exp(-pf[5] + v[7] - v[8]))^-2 + 2 * v[13] * exp(-pf[5] + v[7] - v[8])^2 * (1 + exp(-pf[5] + v[7] - v[8]))^-3) + v[6] * (-pf[6] * (-2 * exp(-pf[5] + v[7] - v[8]) * (1 + exp(-pf[5] + v[7] - v[8]))^-2 - v[7] * exp(-pf[5] + v[7] - v[8]) * (1 + exp(-pf[5] + v[7] - v[8]))^-2 + 2 * v[7] * exp(-pf[5] + v[7] - v[8])^2 * (1 + exp(-pf[5] + v[7] - v[8]))^-3 + v[8] * exp(-pf[5] + v[7] - v[8]) * (1 + exp(-pf[5] + v[7] - v[8]))^-2 - 2 * v[8] * exp(-pf[5] + v[7] - v[8])^2 * (1 + exp(-pf[5] + v[7] - v[8]))^-3) + v[10] * exp(-pf[5] + v[7] - v[8]) * (1 + exp(-pf[5] + v[7] - v[8]))^-2 - 2 * v[10] * exp(-pf[5] + v[7] - v[8])^2 * (1 + exp(-pf[5] + v[7] - v[8]))^-3 - v[11] * exp(-pf[5] + v[7] - v[8]) * (1 + exp(-pf[5] + v[7] - v[8]))^-2 + 2 * v[11] * exp(-pf[5] + v[7] - v[8])^2 * (1 + exp(-pf[5] + v[7] - v[8]))^-3) + pf[1] * v[4] * (2 * exp(-pf[5] + v[7] - v[8]) * (1 + exp(-pf[5] + v[7] - v[8]))^-2 + v[7] * exp(-pf[5] + v[7] - v[8]) * (1 + exp(-pf[5] + v[7] - v[8]))^-2 - 2 * v[7] * exp(-pf[5] + v[7] - v[8])^2 * (1 + exp(-pf[5] + v[7] - v[8]))^-3 - v[8] * exp(-pf[5] + v[7] - v[8]) * (1 + exp(-pf[5] + v[7] - v[8]))^-2 + 2 * v[8] * exp(-pf[5] + v[7] - v[8])^2 * (1 + exp(-pf[5] + v[7] - v[8]))^-3) - pf[1] * v[4] * exp(-pf[5] + v[7] - v[8]) * (1 + exp(-pf[5] + v[7] - v[8]))^-2
    Atx[40] = -v[6] * exp(-pf[5] + v[7] - v[8]) * (1 + exp(-pf[5] + v[7] - v[8]))^-2
    Atx[41] = v[6] * exp(-pf[5] + v[7] - v[8]) * (1 + exp(-pf[5] + v[7] - v[8]))^-2
    Atx[42] = -pf[1] * pf[6] * v[5] * exp(-pf[4] - v[7] + v[8]) * (1 + exp(-pf[4] - v[7] + v[8]))^-2
    Atx[43] = pf[1] * pf[6] * v[5] * exp(-pf[4] - v[7] + v[8]) * (1 + exp(-pf[4] - v[7] + v[8]))^-2
    Atx[44] = pf[1] * pf[6] * v[6] * exp(-pf[5] + v[7] - v[8]) * (1 + exp(-pf[5] + v[7] - v[8]))^-2
    Atx[45] = -pf[1] * pf[6] * v[6] * exp(-pf[5] + v[7] - v[8]) * (1 + exp(-pf[5] + v[7] - v[8]))^-2
    At <- sparseMatrix(i = c(1, 2, 2, 3, 3, 4, 4, 4, 4, 5,
                             5, 5, 5, 6, 6, 6, 6, 7, 7, 7,
                             7, 8, 8, 8, 8, 9, 9, 9, 9, 10,
                             10, 10, 10, 10, 10, 11, 11, 11, 11, 11,
                             11, 12, 12, 13, 13),
                       j = c(9, 7, 8, 7, 8, 5, 7, 8, 10, 6,
                             7, 8, 11, 7, 8, 10, 11, 7, 8, 10,
                             11, 7, 8, 10, 12, 7, 8, 11, 13, 3,
                             5, 7, 8, 10, 11, 4, 6, 7, 8, 10,
                             11, 7, 8, 7, 8),
                         x = Atx, dims = c(13, 13))

    Atp1x <- numeric(16)
    Atp1x[1] = pf[1] * pf[2] * (1 + exp(-pf[4] - v[7] + v[8]))^-1
    Atp1x[2] = -pf[1] * (1 + exp(-pf[4] - v[7] + v[8]))^-1
    Atp1x[3] = pf[1] * pf[2] * (1 + exp(-pf[5] + v[7] - v[8]))^-1
    Atp1x[4] = -pf[1] * (1 + exp(-pf[5] + v[7] - v[8]))^-1
    Atp1x[5] = -pf[1] * (1 + exp(-pf[4] - v[7] + v[8]))^-1
    Atp1x[6] = -pf[1] * (1 - (1 + exp(-pf[4] - v[7] + v[8]))^-1)
    Atp1x[7] = -pf[1] * (1 - (1 + exp(-pf[5] + v[7] - v[8]))^-1)
    Atp1x[8] = -pf[1] * (1 + exp(-pf[5] + v[7] - v[8]))^-1
    Atp1x[9] = -pf[1] * (1 + exp(-pf[4] - v[7] + v[8]))^-1
    Atp1x[10] = pf[1] * exp(-pf[4] - v[7] + v[8]) * (1 + exp(-pf[4] - v[7] + v[8]))^-2
    Atp1x[11] = -pf[1] * exp(-pf[4] - v[7] + v[8]) * (1 + exp(-pf[4] - v[7] + v[8]))^-2
    Atp1x[12] = -pf[1] * (1 + exp(-pf[5] + v[7] - v[8]))^-1
    Atp1x[13] = -pf[1] * exp(-pf[5] + v[7] - v[8]) * (1 + exp(-pf[5] + v[7] - v[8]))^-2
    Atp1x[14] = pf[1] * exp(-pf[5] + v[7] - v[8]) * (1 + exp(-pf[5] + v[7] - v[8]))^-2
    Atp1x[15] = -pf[1] * pf[6] * (1 + exp(-pf[4] - v[7] + v[8]))^-1
    Atp1x[16] = -pf[1] * pf[6] * (1 + exp(-pf[5] + v[7] - v[8]))^-1
    Atp1 <- sparseMatrix(i = c(4, 4, 5, 5, 8, 8, 9, 9, 10, 10,
                               10, 11, 11, 11, 12, 13),
                         j = c(3, 5, 4, 6, 12, 13, 12, 13, 3, 12,
                               13, 4, 12, 13, 5, 6),
                         x = Atp1x, dims = c(13, 13))

    Aepsx <- numeric(1)
    Aepsx[1] = exp(pf[3] * log(v[9]))
    Aeps <- sparseMatrix(i = c(1),
                         j = c(1),
                         x = Aepsx, dims = c(13, 1))

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
