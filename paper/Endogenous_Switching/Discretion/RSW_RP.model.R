# Generated on 2024-10-28 10:32:44 by gEcon ver. 1.2.1 (2023-01-18)
# http://gecon.r-forge.r-project.org/

# Model name: RSW_RP

# info
info__ <- c("RSW_RP", "C:/Users/fm007/Documents/GitHub/OptimalRSW/paper/Endogenous_Switching/Discretion/RSW_RP.gcn", "2024-10-28 10:32:44", "false")

# index sets
index_sets__ <- list()

# variables
variables__ <- c("etapi",
                 "iH",
                 "iL",
                 "lambda__HIGHREGIME_1",
                 "lambda__HIGHREGIME_2",
                 "lambda__LOWREGIME_1",
                 "lambda__LOWREGIME_2",
                 "piH",
                 "piL",
                 "yH",
                 "yL",
                 "UH",
                 "UL")

variables_tex__ <- c("{e\\!t\\!a\\!p\\!i}",
                     "{i\\!H}",
                     "{i\\!L}",
                     "\\lambda^{\\mathrm{HIGHREGIME}^{\\mathrm{1}}}",
                     "\\lambda^{\\mathrm{HIGHREGIME}^{\\mathrm{2}}}",
                     "\\lambda^{\\mathrm{LOWREGIME}^{\\mathrm{1}}}",
                     "\\lambda^{\\mathrm{LOWREGIME}^{\\mathrm{2}}}",
                     "{p\\!i\\!H}",
                     "{p\\!i\\!L}",
                     "{y\\!H}",
                     "{y\\!L}",
                     "{U\\!H}",
                     "{U\\!L}")

# shocks
shocks__ <- c("epsilon_pi")

shocks_tex__ <- c("\\epsilon^{\\pi}")

# parameters
parameters__ <- c("beta",
                  "kappa",
                  "phi",
                  "pitH",
                  "pitCB",
                  "pitL",
                  "pHss",
                  "pLss",
                  "sigma",
                  "tau",
                  "theta")

parameters_tex__ <- c("\\beta",
                     "\\kappa",
                     "\\phi",
                     "{p\\!i\\!t\\!H}",
                     "{p\\!i\\!t\\!C\\!B}",
                     "{p\\!i\\!t\\!L}",
                     "{p\\!H\\!s\\!s}",
                     "{p\\!L\\!s\\!s}",
                     "\\sigma",
                     "\\tau",
                     "\\theta")

# free parameters
parameters_free__ <- c("beta",
                       "kappa",
                       "phi",
                       "pitH",
                       "pitCB",
                       "pitL",
                       "pHss",
                       "pLss",
                       "sigma",
                       "tau",
                       "theta")

# free parameters' values
parameters_free_val__ <- c(0.99,
                           0.2465,
                           0.95,
                           0,
                           0,
                           2,
                           0.99,
                           0.99,
                           1,
                           0.001,
                           6)

# equations
equations__ <- c("-etapi[] + exp(epsilon_pi[] + phi * log(etapi[-1])) = 0",
                 "lambda__HIGHREGIME_2[] * (pHss + tau * (-pitCB + piH[])^2) + (beta - beta * (1 - pHss - tau * (-pitCB + piH[])^2)) * E[][-lambda__HIGHREGIME_2[1] + kappa * lambda__HIGHREGIME_1[1]] - kappa * theta^-1 * yH[] = 0",
                 "lambda__LOWREGIME_2[] * (pLss - tau * (-pitCB + piL[])^2) + (beta - beta * (1 - pLss + tau * (-pitCB + piL[])^2)) * E[][-lambda__LOWREGIME_2[1] + kappa * lambda__LOWREGIME_1[1]] - kappa * theta^-1 * yL[] = 0",
                 "-piH[-1] + log(etapi[-1]) + kappa * yH[-1] + beta * piH[] * (pHss + tau * (-pitCB + piH[])^2) + beta * (-piH[] + piL[]) * (1 - pHss - tau * (-pitCB + piH[])^2) = 0",
                 "-piL[-1] + log(etapi[-1]) + kappa * yL[-1] + beta * piL[] * (pLss - tau * (-pitCB + piL[])^2) + beta * (piH[] - piL[]) * (1 - pLss + tau * (-pitCB + piL[])^2) = 0",
                 "-yH[-1] + yH[] - sigma * (iH[-1] - piH[]) + (-yH[] + yL[]) * (1 - pHss - tau * (-pitCB + piH[])^2) + sigma * (-piH[] + piL[]) * (1 - pHss - tau * (-pitCB + piH[])^2) = 0",
                 "-yL[-1] + yL[] - sigma * (iL[-1] - piL[]) + (yH[] - yL[]) * (1 - pLss + tau * (-pitCB + piL[])^2) + sigma * (piH[] - piL[]) * (1 - pLss + tau * (-pitCB + piL[])^2) = 0",
                 "UH[] + 0.5 * (pitH - pitCB + piH[])^2 - beta * E[][UH[1]] - beta * (-E[][UH[1]] + E[][UL[1]]) * (1 - pHss - tau * (-pitCB + piH[])^2) + 0.5 * kappa * theta^-1 * yH[]^2 = 0",
                 "UL[] + 0.5 * (-pitCB + pitL + piL[])^2 - beta * E[][UL[1]] - beta * (E[][UH[1]] - E[][UL[1]]) * (1 - pLss + tau * (-pitCB + piL[])^2) + 0.5 * kappa * theta^-1 * yL[]^2 = 0",
                 "-pitH + pitCB - piH[] + lambda__HIGHREGIME_1[] * (beta * (pHss + tau * (-pitCB + piH[])^2) - beta * (1 - pHss - tau * (-pitCB + piH[])^2) + 2 * beta * tau * piH[] * (-pitCB + piH[]) - 2 * beta * tau * (-pitCB + piH[]) * (-piH[] + piL[])) + lambda__HIGHREGIME_2[] * (sigma - sigma * (1 - pHss - tau * (-pitCB + piH[])^2) - 2 * tau * (-pitCB + piH[]) * (-yH[] + yL[]) - 2 * sigma * tau * (-pitCB + piH[]) * (-piH[] + piL[])) - (beta - beta * (1 - pHss - tau * (-pitCB + piH[])^2)) * E[][lambda__HIGHREGIME_1[1]] - 2 * beta * tau * (-pitCB + piH[]) * (-E[][UH[1]] + E[][UL[1]]) = 0",
                 "pitCB - pitL - piL[] + lambda__LOWREGIME_1[] * (beta * (pLss - tau * (-pitCB + piL[])^2) - beta * (1 - pLss + tau * (-pitCB + piL[])^2) - 2 * beta * tau * piL[] * (-pitCB + piL[]) + 2 * beta * tau * (-pitCB + piL[]) * (piH[] - piL[])) + lambda__LOWREGIME_2[] * (sigma - sigma * (1 - pLss + tau * (-pitCB + piL[])^2) + 2 * tau * (-pitCB + piL[]) * (yH[] - yL[]) + 2 * sigma * tau * (-pitCB + piL[]) * (piH[] - piL[])) - (beta - beta * (1 - pLss + tau * (-pitCB + piL[])^2)) * E[][lambda__LOWREGIME_1[1]] + 2 * beta * tau * (-pitCB + piL[]) * (E[][UH[1]] - E[][UL[1]]) = 0",
                 "-sigma * (beta - beta * (1 - pHss - tau * (-pitCB + piH[])^2)) * E[][lambda__HIGHREGIME_2[1]] = 0",
                 "-sigma * (beta - beta * (1 - pLss + tau * (-pitCB + piL[])^2)) * E[][lambda__LOWREGIME_2[1]] = 0")

# calibrating equations
calibr_equations__ <- character(0)

# variables / equations map
vareqmap__ <- sparseMatrix(i = c(1, 2, 2, 2, 2, 3, 3, 3, 3, 4,
                                 4, 4, 4, 5, 5, 5, 5, 6, 6, 6,
                                 6, 6, 7, 7, 7, 7, 7, 8, 8, 8,
                                 8, 9, 9, 9, 9, 10, 10, 10, 10, 10,
                                 10, 10, 10, 11, 11, 11, 11, 11, 11, 11,
                                 11, 12, 12, 13, 13),
                           j = c(1, 4, 5, 8, 10, 6, 7, 9, 11, 1,
                                 8, 9, 10, 1, 8, 9, 11, 2, 8, 9,
                                 10, 11, 3, 8, 9, 10, 11, 8, 10, 12,
                                 13, 9, 11, 12, 13, 4, 5, 8, 9, 10,
                                 11, 12, 13, 6, 7, 8, 9, 10, 11, 12,
                                 13, 5, 8, 7, 9),
                           x = c(3, 4, 6, 2, 2, 4, 6, 2, 2, 1,
                                 3, 2, 1, 1, 2, 3, 1, 1, 2, 2,
                                 3, 2, 1, 2, 2, 2, 3, 2, 2, 6,
                                 4, 2, 2, 4, 6, 6, 2, 2, 2, 2,
                                 2, 4, 4, 6, 2, 2, 2, 2, 2, 4,
                                 4, 4, 2, 4, 2),
                           dims = c(13, 13))

# variables / calibrating equations map
varcalibreqmap__ <- sparseMatrix(i = NULL, j = NULL, dims = c(0, 13))

# calibrated parameters / equations map
calibrpareqmap__ <- sparseMatrix(i = NULL, j = NULL, dims = c(13, 0))

# calibrated parameters / calibrating equations map
calibrparcalibreqmap__ <- sparseMatrix(i = NULL, j = NULL, dims = c(0, 0))

# free parameters / equations map
freepareqmap__ <- sparseMatrix(i = c(1, 2, 2, 2, 2, 2, 2, 3, 3, 3,
                                     3, 3, 3, 4, 4, 4, 4, 4, 5, 5,
                                     5, 5, 5, 6, 6, 6, 6, 7, 7, 7,
                                     7, 8, 8, 8, 8, 8, 8, 8, 9, 9,
                                     9, 9, 9, 9, 9, 10, 10, 10, 10, 10,
                                     10, 11, 11, 11, 11, 11, 11, 12, 12, 12,
                                     12, 12, 13, 13, 13, 13, 13),
                               j = c(3, 1, 2, 5, 7, 10, 11, 1, 2, 5,
                                     8, 10, 11, 1, 2, 5, 7, 10, 1, 2,
                                     5, 8, 10, 5, 7, 9, 10, 5, 8, 9,
                                     10, 1, 2, 4, 5, 7, 10, 11, 1, 2,
                                     5, 6, 8, 10, 11, 1, 4, 5, 7, 9,
                                     10, 1, 5, 6, 8, 9, 10, 1, 5, 7,
                                     9, 10, 1, 5, 8, 9, 10),
                               x = rep(1, 67), dims = c(13, 11))

# free parameters / calibrating equations map
freeparcalibreqmap__ <- sparseMatrix(i = NULL, j = NULL, dims = c(0, 11))

# shocks / equations map
shockeqmap__ <- sparseMatrix(i = c(1),
                             j = c(1),
                             x = rep(1, 1), dims = c(13, 1))

# steady state equations
ss_eq__ <- function(v, pc, pf)
{
    r <- numeric(13)
    r[1] = -v[1] + exp(pf[3] * log(v[1]))
    r[2] = v[5] * (pf[7] + pf[10] * (-pf[5] + v[8])^2) + (pf[1] - pf[1] * (1 - pf[7] - pf[10] * (-pf[5] + v[8])^2)) * (-v[5] + pf[2] * v[4]) - pf[2] * pf[11]^-1 * v[10]
    r[3] = v[7] * (pf[8] - pf[10] * (-pf[5] + v[9])^2) + (pf[1] - pf[1] * (1 - pf[8] + pf[10] * (-pf[5] + v[9])^2)) * (-v[7] + pf[2] * v[6]) - pf[2] * pf[11]^-1 * v[11]
    r[4] = -v[8] + log(v[1]) + pf[2] * v[10] + pf[1] * v[8] * (pf[7] + pf[10] * (-pf[5] + v[8])^2) + pf[1] * (-v[8] + v[9]) * (1 - pf[7] - pf[10] * (-pf[5] + v[8])^2)
    r[5] = -v[9] + log(v[1]) + pf[2] * v[11] + pf[1] * v[9] * (pf[8] - pf[10] * (-pf[5] + v[9])^2) + pf[1] * (v[8] - v[9]) * (1 - pf[8] + pf[10] * (-pf[5] + v[9])^2)
    r[6] = (-v[10] + v[11]) * (1 - pf[7] - pf[10] * (-pf[5] + v[8])^2) - pf[9] * (v[2] - v[8]) + pf[9] * (-v[8] + v[9]) * (1 - pf[7] - pf[10] * (-pf[5] + v[8])^2)
    r[7] = (v[10] - v[11]) * (1 - pf[8] + pf[10] * (-pf[5] + v[9])^2) - pf[9] * (v[3] - v[9]) + pf[9] * (v[8] - v[9]) * (1 - pf[8] + pf[10] * (-pf[5] + v[9])^2)
    r[8] = v[12] + 0.5 * (pf[4] - pf[5] + v[8])^2 - pf[1] * v[12] - pf[1] * (-v[12] + v[13]) * (1 - pf[7] - pf[10] * (-pf[5] + v[8])^2) + 0.5 * pf[2] * pf[11]^-1 * v[10]^2
    r[9] = v[13] + 0.5 * (-pf[5] + pf[6] + v[9])^2 - pf[1] * v[13] - pf[1] * (v[12] - v[13]) * (1 - pf[8] + pf[10] * (-pf[5] + v[9])^2) + 0.5 * pf[2] * pf[11]^-1 * v[11]^2
    r[10] = -pf[4] + pf[5] - v[8] - v[4] * (pf[1] - pf[1] * (1 - pf[7] - pf[10] * (-pf[5] + v[8])^2)) + v[4] * (pf[1] * (pf[7] + pf[10] * (-pf[5] + v[8])^2) - pf[1] * (1 - pf[7] - pf[10] * (-pf[5] + v[8])^2) + 2 * pf[1] * pf[10] * v[8] * (-pf[5] + v[8]) - 2 * pf[1] * pf[10] * (-pf[5] + v[8]) * (-v[8] + v[9])) + v[5] * (pf[9] - pf[9] * (1 - pf[7] - pf[10] * (-pf[5] + v[8])^2) - 2 * pf[10] * (-pf[5] + v[8]) * (-v[10] + v[11]) - 2 * pf[9] * pf[10] * (-pf[5] + v[8]) * (-v[8] + v[9])) - 2 * pf[1] * pf[10] * (-pf[5] + v[8]) * (-v[12] + v[13])
    r[11] = pf[5] - pf[6] - v[9] - v[6] * (pf[1] - pf[1] * (1 - pf[8] + pf[10] * (-pf[5] + v[9])^2)) + v[6] * (pf[1] * (pf[8] - pf[10] * (-pf[5] + v[9])^2) - pf[1] * (1 - pf[8] + pf[10] * (-pf[5] + v[9])^2) - 2 * pf[1] * pf[10] * v[9] * (-pf[5] + v[9]) + 2 * pf[1] * pf[10] * (-pf[5] + v[9]) * (v[8] - v[9])) + v[7] * (pf[9] - pf[9] * (1 - pf[8] + pf[10] * (-pf[5] + v[9])^2) + 2 * pf[10] * (-pf[5] + v[9]) * (v[10] - v[11]) + 2 * pf[9] * pf[10] * (-pf[5] + v[9]) * (v[8] - v[9])) + 2 * pf[1] * pf[10] * (-pf[5] + v[9]) * (v[12] - v[13])
    r[12] = -pf[9] * v[5] * (pf[1] - pf[1] * (1 - pf[7] - pf[10] * (-pf[5] + v[8])^2))
    r[13] = -pf[9] * v[7] * (pf[1] - pf[1] * (1 - pf[8] + pf[10] * (-pf[5] + v[9])^2))

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
    jac <- numeric(55)
    jac[1] = -1 + pf[3] * v[1]^-1 * exp(pf[3] * log(v[1]))
    jac[2] = pf[2] * (pf[1] - pf[1] * (1 - pf[7] - pf[10] * (-pf[5] + v[8])^2))
    jac[3] = -pf[1] + pf[7] + pf[1] * (1 - pf[7] - pf[10] * (-pf[5] + v[8])^2) + pf[10] * (-pf[5] + v[8])^2
    jac[4] = 2 * pf[10] * v[5] * (-pf[5] + v[8]) + 2 * pf[1] * pf[10] * (-pf[5] + v[8]) * (-v[5] + pf[2] * v[4])
    jac[5] = -pf[2] * pf[11]^-1
    jac[6] = pf[2] * (pf[1] - pf[1] * (1 - pf[8] + pf[10] * (-pf[5] + v[9])^2))
    jac[7] = -pf[1] + pf[8] + pf[1] * (1 - pf[8] + pf[10] * (-pf[5] + v[9])^2) - pf[10] * (-pf[5] + v[9])^2
    jac[8] = -2 * pf[10] * v[7] * (-pf[5] + v[9]) - 2 * pf[1] * pf[10] * (-pf[5] + v[9]) * (-v[7] + pf[2] * v[6])
    jac[9] = -pf[2] * pf[11]^-1
    jac[10] = v[1]^-1
    jac[11] = -1 + pf[1] * (pf[7] + pf[10] * (-pf[5] + v[8])^2) - pf[1] * (1 - pf[7] - pf[10] * (-pf[5] + v[8])^2) + 2 * pf[1] * pf[10] * v[8] * (-pf[5] + v[8]) - 2 * pf[1] * pf[10] * (-pf[5] + v[8]) * (-v[8] + v[9])
    jac[12] = pf[1] * (1 - pf[7] - pf[10] * (-pf[5] + v[8])^2)
    jac[13] = pf[2]
    jac[14] = v[1]^-1
    jac[15] = pf[1] * (1 - pf[8] + pf[10] * (-pf[5] + v[9])^2)
    jac[16] = -1 + pf[1] * (pf[8] - pf[10] * (-pf[5] + v[9])^2) - pf[1] * (1 - pf[8] + pf[10] * (-pf[5] + v[9])^2) - 2 * pf[1] * pf[10] * v[9] * (-pf[5] + v[9]) + 2 * pf[1] * pf[10] * (-pf[5] + v[9]) * (v[8] - v[9])
    jac[17] = pf[2]
    jac[18] = -pf[9]
    jac[19] = pf[9] - pf[9] * (1 - pf[7] - pf[10] * (-pf[5] + v[8])^2) - 2 * pf[10] * (-pf[5] + v[8]) * (-v[10] + v[11]) - 2 * pf[9] * pf[10] * (-pf[5] + v[8]) * (-v[8] + v[9])
    jac[20] = pf[9] * (1 - pf[7] - pf[10] * (-pf[5] + v[8])^2)
    jac[21] = -1 + pf[7] + pf[10] * (-pf[5] + v[8])^2
    jac[22] = 1 - pf[7] - pf[10] * (-pf[5] + v[8])^2
    jac[23] = -pf[9]
    jac[24] = pf[9] * (1 - pf[8] + pf[10] * (-pf[5] + v[9])^2)
    jac[25] = pf[9] - pf[9] * (1 - pf[8] + pf[10] * (-pf[5] + v[9])^2) + 2 * pf[10] * (-pf[5] + v[9]) * (v[10] - v[11]) + 2 * pf[9] * pf[10] * (-pf[5] + v[9]) * (v[8] - v[9])
    jac[26] = 1 - pf[8] + pf[10] * (-pf[5] + v[9])^2
    jac[27] = -1 + pf[8] - pf[10] * (-pf[5] + v[9])^2
    jac[28] = pf[4] - pf[5] + v[8] + 2 * pf[1] * pf[10] * (-pf[5] + v[8]) * (-v[12] + v[13])
    jac[29] = pf[2] * pf[11]^-1 * v[10]
    jac[30] = 1 - pf[1] + pf[1] * (1 - pf[7] - pf[10] * (-pf[5] + v[8])^2)
    jac[31] = -pf[1] * (1 - pf[7] - pf[10] * (-pf[5] + v[8])^2)
    jac[32] = -pf[5] + pf[6] + v[9] - 2 * pf[1] * pf[10] * (-pf[5] + v[9]) * (v[12] - v[13])
    jac[33] = pf[2] * pf[11]^-1 * v[11]
    jac[34] = -pf[1] * (1 - pf[8] + pf[10] * (-pf[5] + v[9])^2)
    jac[35] = 1 - pf[1] + pf[1] * (1 - pf[8] + pf[10] * (-pf[5] + v[9])^2)
    jac[36] = -pf[1] + pf[1] * (pf[7] + pf[10] * (-pf[5] + v[8])^2) + 2 * pf[1] * pf[10] * v[8] * (-pf[5] + v[8]) - 2 * pf[1] * pf[10] * (-pf[5] + v[8]) * (-v[8] + v[9])
    jac[37] = pf[9] - pf[9] * (1 - pf[7] - pf[10] * (-pf[5] + v[8])^2) - 2 * pf[10] * (-pf[5] + v[8]) * (-v[10] + v[11]) - 2 * pf[9] * pf[10] * (-pf[5] + v[8]) * (-v[8] + v[9])
    jac[38] = -1 + v[4] * (2 * pf[1] * pf[10] * v[8] + 8 * pf[1] * pf[10] * (-pf[5] + v[8]) - 2 * pf[1] * pf[10] * (-v[8] + v[9])) + v[5] * (-2 * pf[10] * (-v[10] + v[11]) + 2 * pf[9] * pf[10] * (-pf[5] + v[8]) - 2 * pf[9] * pf[10] * (-v[8] + v[9]) + 2 * pf[9] * pf[10] * (-pf[5] + v[8])) - 2 * pf[1] * pf[10] * (-v[12] + v[13]) - 2 * pf[1] * pf[10] * v[4] * (-pf[5] + v[8])
    jac[39] = -2 * pf[1] * pf[10] * v[4] * (-pf[5] + v[8]) - 2 * pf[9] * pf[10] * v[5] * (-pf[5] + v[8])
    jac[40] = 2 * pf[10] * v[5] * (-pf[5] + v[8])
    jac[41] = -2 * pf[10] * v[5] * (-pf[5] + v[8])
    jac[42] = 2 * pf[1] * pf[10] * (-pf[5] + v[8])
    jac[43] = -2 * pf[1] * pf[10] * (-pf[5] + v[8])
    jac[44] = -pf[1] + pf[1] * (pf[8] - pf[10] * (-pf[5] + v[9])^2) - 2 * pf[1] * pf[10] * v[9] * (-pf[5] + v[9]) + 2 * pf[1] * pf[10] * (-pf[5] + v[9]) * (v[8] - v[9])
    jac[45] = pf[9] - pf[9] * (1 - pf[8] + pf[10] * (-pf[5] + v[9])^2) + 2 * pf[10] * (-pf[5] + v[9]) * (v[10] - v[11]) + 2 * pf[9] * pf[10] * (-pf[5] + v[9]) * (v[8] - v[9])
    jac[46] = 2 * pf[1] * pf[10] * v[6] * (-pf[5] + v[9]) + 2 * pf[9] * pf[10] * v[7] * (-pf[5] + v[9])
    jac[47] = -1 + v[6] * (-2 * pf[1] * pf[10] * v[9] - 8 * pf[1] * pf[10] * (-pf[5] + v[9]) + 2 * pf[1] * pf[10] * (v[8] - v[9])) + v[7] * (2 * pf[10] * (v[10] - v[11]) - 2 * pf[9] * pf[10] * (-pf[5] + v[9]) + 2 * pf[9] * pf[10] * (v[8] - v[9]) - 2 * pf[9] * pf[10] * (-pf[5] + v[9])) + 2 * pf[1] * pf[10] * (v[12] - v[13]) + 2 * pf[1] * pf[10] * v[6] * (-pf[5] + v[9])
    jac[48] = 2 * pf[10] * v[7] * (-pf[5] + v[9])
    jac[49] = -2 * pf[10] * v[7] * (-pf[5] + v[9])
    jac[50] = 2 * pf[1] * pf[10] * (-pf[5] + v[9])
    jac[51] = -2 * pf[1] * pf[10] * (-pf[5] + v[9])
    jac[52] = -pf[9] * (pf[1] - pf[1] * (1 - pf[7] - pf[10] * (-pf[5] + v[8])^2))
    jac[53] = -2 * pf[1] * pf[9] * pf[10] * v[5] * (-pf[5] + v[8])
    jac[54] = -pf[9] * (pf[1] - pf[1] * (1 - pf[8] + pf[10] * (-pf[5] + v[9])^2))
    jac[55] = 2 * pf[1] * pf[9] * pf[10] * v[7] * (-pf[5] + v[9])
    jacob <- sparseMatrix(i = c(1, 2, 2, 2, 2, 3, 3, 3, 3, 4,
                                4, 4, 4, 5, 5, 5, 5, 6, 6, 6,
                                6, 6, 7, 7, 7, 7, 7, 8, 8, 8,
                                8, 9, 9, 9, 9, 10, 10, 10, 10, 10,
                                10, 10, 10, 11, 11, 11, 11, 11, 11, 11,
                                11, 12, 12, 13, 13),
                          j = c(1, 4, 5, 8, 10, 6, 7, 9, 11, 1,
                                8, 9, 10, 1, 8, 9, 11, 2, 8, 9,
                                10, 11, 3, 8, 9, 10, 11, 8, 10, 12,
                                13, 9, 11, 12, 13, 4, 5, 8, 9, 10,
                                11, 12, 13, 6, 7, 8, 9, 10, 11, 12,
                                13, 5, 8, 7, 9),
                          x = jac, dims = c(13, 13))

    return(jacob)
}

# 1st order perturbation
pert1__ <- function(v, pc, pf)
{
    Atm1x <- numeric(11)
    Atm1x[1] = pf[3] * v[1]^-1 * exp(pf[3] * log(v[1]))
    Atm1x[2] = v[1]^-1
    Atm1x[3] = -1
    Atm1x[4] = pf[2]
    Atm1x[5] = v[1]^-1
    Atm1x[6] = -1
    Atm1x[7] = pf[2]
    Atm1x[8] = -pf[9]
    Atm1x[9] = -1
    Atm1x[10] = -pf[9]
    Atm1x[11] = -1
    Atm1 <- sparseMatrix(i = c(1, 4, 4, 4, 5, 5, 5, 6, 6, 7,
                               7),
                         j = c(1, 1, 8, 10, 1, 9, 11, 2, 10, 3,
                               11),
                         x = Atm1x, dims = c(13, 13))

    Atx <- numeric(39)
    Atx[1] = -1
    Atx[2] = pf[7] + pf[10] * (-pf[5] + v[8])^2
    Atx[3] = 2 * pf[10] * v[5] * (-pf[5] + v[8]) + 2 * pf[1] * pf[10] * (-pf[5] + v[8]) * (-v[5] + pf[2] * v[4])
    Atx[4] = -pf[2] * pf[11]^-1
    Atx[5] = pf[8] - pf[10] * (-pf[5] + v[9])^2
    Atx[6] = -2 * pf[10] * v[7] * (-pf[5] + v[9]) - 2 * pf[1] * pf[10] * (-pf[5] + v[9]) * (-v[7] + pf[2] * v[6])
    Atx[7] = -pf[2] * pf[11]^-1
    Atx[8] = pf[1] * (pf[7] + pf[10] * (-pf[5] + v[8])^2) - pf[1] * (1 - pf[7] - pf[10] * (-pf[5] + v[8])^2) + 2 * pf[1] * pf[10] * v[8] * (-pf[5] + v[8]) - 2 * pf[1] * pf[10] * (-pf[5] + v[8]) * (-v[8] + v[9])
    Atx[9] = pf[1] * (1 - pf[7] - pf[10] * (-pf[5] + v[8])^2)
    Atx[10] = pf[1] * (1 - pf[8] + pf[10] * (-pf[5] + v[9])^2)
    Atx[11] = pf[1] * (pf[8] - pf[10] * (-pf[5] + v[9])^2) - pf[1] * (1 - pf[8] + pf[10] * (-pf[5] + v[9])^2) - 2 * pf[1] * pf[10] * v[9] * (-pf[5] + v[9]) + 2 * pf[1] * pf[10] * (-pf[5] + v[9]) * (v[8] - v[9])
    Atx[12] = pf[9] - pf[9] * (1 - pf[7] - pf[10] * (-pf[5] + v[8])^2) - 2 * pf[10] * (-pf[5] + v[8]) * (-v[10] + v[11]) - 2 * pf[9] * pf[10] * (-pf[5] + v[8]) * (-v[8] + v[9])
    Atx[13] = pf[9] * (1 - pf[7] - pf[10] * (-pf[5] + v[8])^2)
    Atx[14] = pf[7] + pf[10] * (-pf[5] + v[8])^2
    Atx[15] = 1 - pf[7] - pf[10] * (-pf[5] + v[8])^2
    Atx[16] = pf[9] * (1 - pf[8] + pf[10] * (-pf[5] + v[9])^2)
    Atx[17] = pf[9] - pf[9] * (1 - pf[8] + pf[10] * (-pf[5] + v[9])^2) + 2 * pf[10] * (-pf[5] + v[9]) * (v[10] - v[11]) + 2 * pf[9] * pf[10] * (-pf[5] + v[9]) * (v[8] - v[9])
    Atx[18] = 1 - pf[8] + pf[10] * (-pf[5] + v[9])^2
    Atx[19] = pf[8] - pf[10] * (-pf[5] + v[9])^2
    Atx[20] = pf[4] - pf[5] + v[8] + 2 * pf[1] * pf[10] * (-pf[5] + v[8]) * (-v[12] + v[13])
    Atx[21] = pf[2] * pf[11]^-1 * v[10]
    Atx[22] = 1
    Atx[23] = -pf[5] + pf[6] + v[9] - 2 * pf[1] * pf[10] * (-pf[5] + v[9]) * (v[12] - v[13])
    Atx[24] = pf[2] * pf[11]^-1 * v[11]
    Atx[25] = 1
    Atx[26] = pf[1] * (pf[7] + pf[10] * (-pf[5] + v[8])^2) - pf[1] * (1 - pf[7] - pf[10] * (-pf[5] + v[8])^2) + 2 * pf[1] * pf[10] * v[8] * (-pf[5] + v[8]) - 2 * pf[1] * pf[10] * (-pf[5] + v[8]) * (-v[8] + v[9])
    Atx[27] = pf[9] - pf[9] * (1 - pf[7] - pf[10] * (-pf[5] + v[8])^2) - 2 * pf[10] * (-pf[5] + v[8]) * (-v[10] + v[11]) - 2 * pf[9] * pf[10] * (-pf[5] + v[8]) * (-v[8] + v[9])
    Atx[28] = -1 + v[4] * (2 * pf[1] * pf[10] * v[8] + 8 * pf[1] * pf[10] * (-pf[5] + v[8]) - 2 * pf[1] * pf[10] * (-v[8] + v[9])) + v[5] * (-2 * pf[10] * (-v[10] + v[11]) + 4 * pf[9] * pf[10] * (-pf[5] + v[8]) - 2 * pf[9] * pf[10] * (-v[8] + v[9])) - 2 * pf[1] * pf[10] * (-v[12] + v[13]) - 2 * pf[1] * pf[10] * v[4] * (-pf[5] + v[8])
    Atx[29] = -2 * pf[1] * pf[10] * v[4] * (-pf[5] + v[8]) - 2 * pf[9] * pf[10] * v[5] * (-pf[5] + v[8])
    Atx[30] = 2 * pf[10] * v[5] * (-pf[5] + v[8])
    Atx[31] = -2 * pf[10] * v[5] * (-pf[5] + v[8])
    Atx[32] = pf[1] * (pf[8] - pf[10] * (-pf[5] + v[9])^2) - pf[1] * (1 - pf[8] + pf[10] * (-pf[5] + v[9])^2) - 2 * pf[1] * pf[10] * v[9] * (-pf[5] + v[9]) + 2 * pf[1] * pf[10] * (-pf[5] + v[9]) * (v[8] - v[9])
    Atx[33] = pf[9] - pf[9] * (1 - pf[8] + pf[10] * (-pf[5] + v[9])^2) + 2 * pf[10] * (-pf[5] + v[9]) * (v[10] - v[11]) + 2 * pf[9] * pf[10] * (-pf[5] + v[9]) * (v[8] - v[9])
    Atx[34] = 2 * pf[1] * pf[10] * v[6] * (-pf[5] + v[9]) + 2 * pf[9] * pf[10] * v[7] * (-pf[5] + v[9])
    Atx[35] = -1 + v[6] * (-2 * pf[1] * pf[10] * v[9] - 8 * pf[1] * pf[10] * (-pf[5] + v[9]) + 2 * pf[1] * pf[10] * (v[8] - v[9])) + v[7] * (2 * pf[10] * (v[10] - v[11]) - 4 * pf[9] * pf[10] * (-pf[5] + v[9]) + 2 * pf[9] * pf[10] * (v[8] - v[9])) + 2 * pf[1] * pf[10] * (v[12] - v[13]) + 2 * pf[1] * pf[10] * v[6] * (-pf[5] + v[9])
    Atx[36] = 2 * pf[10] * v[7] * (-pf[5] + v[9])
    Atx[37] = -2 * pf[10] * v[7] * (-pf[5] + v[9])
    Atx[38] = -2 * pf[1] * pf[9] * pf[10] * v[5] * (-pf[5] + v[8])
    Atx[39] = 2 * pf[1] * pf[9] * pf[10] * v[7] * (-pf[5] + v[9])
    At <- sparseMatrix(i = c(1, 2, 2, 2, 3, 3, 3, 4, 4, 5,
                             5, 6, 6, 6, 6, 7, 7, 7, 7, 8,
                             8, 8, 9, 9, 9, 10, 10, 10, 10, 10,
                             10, 11, 11, 11, 11, 11, 11, 12, 13),
                       j = c(1, 5, 8, 10, 7, 9, 11, 8, 9, 8,
                             9, 8, 9, 10, 11, 8, 9, 10, 11, 8,
                             10, 12, 9, 11, 13, 4, 5, 8, 9, 10,
                             11, 6, 7, 8, 9, 10, 11, 8, 9),
                         x = Atx, dims = c(13, 13))

    Atp1x <- numeric(16)
    Atp1x[1] = pf[2] * (pf[1] - pf[1] * (1 - pf[7] - pf[10] * (-pf[5] + v[8])^2))
    Atp1x[2] = -pf[1] + pf[1] * (1 - pf[7] - pf[10] * (-pf[5] + v[8])^2)
    Atp1x[3] = pf[2] * (pf[1] - pf[1] * (1 - pf[8] + pf[10] * (-pf[5] + v[9])^2))
    Atp1x[4] = -pf[1] + pf[1] * (1 - pf[8] + pf[10] * (-pf[5] + v[9])^2)
    Atp1x[5] = -pf[1] + pf[1] * (1 - pf[7] - pf[10] * (-pf[5] + v[8])^2)
    Atp1x[6] = -pf[1] * (1 - pf[7] - pf[10] * (-pf[5] + v[8])^2)
    Atp1x[7] = -pf[1] * (1 - pf[8] + pf[10] * (-pf[5] + v[9])^2)
    Atp1x[8] = -pf[1] + pf[1] * (1 - pf[8] + pf[10] * (-pf[5] + v[9])^2)
    Atp1x[9] = -pf[1] + pf[1] * (1 - pf[7] - pf[10] * (-pf[5] + v[8])^2)
    Atp1x[10] = 2 * pf[1] * pf[10] * (-pf[5] + v[8])
    Atp1x[11] = -2 * pf[1] * pf[10] * (-pf[5] + v[8])
    Atp1x[12] = -pf[1] + pf[1] * (1 - pf[8] + pf[10] * (-pf[5] + v[9])^2)
    Atp1x[13] = 2 * pf[1] * pf[10] * (-pf[5] + v[9])
    Atp1x[14] = -2 * pf[1] * pf[10] * (-pf[5] + v[9])
    Atp1x[15] = -pf[9] * (pf[1] - pf[1] * (1 - pf[7] - pf[10] * (-pf[5] + v[8])^2))
    Atp1x[16] = -pf[9] * (pf[1] - pf[1] * (1 - pf[8] + pf[10] * (-pf[5] + v[9])^2))
    Atp1 <- sparseMatrix(i = c(2, 2, 3, 3, 8, 8, 9, 9, 10, 10,
                               10, 11, 11, 11, 12, 13),
                         j = c(4, 5, 6, 7, 12, 13, 12, 13, 4, 12,
                               13, 6, 12, 13, 5, 7),
                         x = Atp1x, dims = c(13, 13))

    Aepsx <- numeric(1)
    Aepsx[1] = exp(pf[3] * log(v[1]))
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
