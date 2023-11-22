# Generated on 2023-11-22 18:10:38 by gEcon ver. 1.2.1 (2023-01-18)
# http://gecon.r-forge.r-project.org/

# Model name: RSW_RP_PERS

# info
info__ <- c("RSW_RP_PERS", "C:/Users/fm007/Documents/GitHub/OptimalRSW/RSW_RP_PERS.gcn", "2023-11-22 18:10:38", "false")

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
                 "lambda__HIGHREGIME_piH__lag_1",
                 "lambda__LOWREGIME_piL__lag_1",
                 "piH",
                 "piL",
                 "piH__lag_1",
                 "piL__lag_1",
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
                     "\\lambda^{\\mathrm{HIGHREGIME}^{\\mathrm{piH}^{\\mathrm{lag}^{\\mathrm{1}}}}}",
                     "\\lambda^{\\mathrm{LOWREGIME}^{\\mathrm{piL}^{\\mathrm{lag}^{\\mathrm{1}}}}}",
                     "{p\\!i\\!H}",
                     "{p\\!i\\!L}",
                     "{p\\!i\\!H}^{\\mathrm{lag}^{\\mathrm{1}}}",
                     "{p\\!i\\!L}^{\\mathrm{lag}^{\\mathrm{1}}}",
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
                  "phipi",
                  "pitH",
                  "pitCB",
                  "pitL",
                  "pH",
                  "pL",
                  "sigma",
                  "theta")

parameters_tex__ <- c("\\beta",
                     "\\kappa",
                     "\\phi",
                     "{p\\!h\\!i\\!p\\!i}",
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
                       "phi",
                       "phipi",
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
                           0.95,
                           0.8,
                           0,
                           0,
                           2,
                           0.99,
                           0.99,
                           1,
                           6)

# equations
equations__ <- c("piH[-1] - piH__lag_1[] = 0",
                 "piL[-1] - piL__lag_1[] = 0",
                 "-etapi[] + exp(epsilon_pi[] + phi * log(etapi[-1])) = 0",
                 "-lambda__HIGHREGIME_piH__lag_1[] + phipi * (beta - beta * (1 - pH)) * E[][lambda__HIGHREGIME_1[1]] = 0",
                 "-lambda__LOWREGIME_piL__lag_1[] + phipi * (beta - beta * (1 - pL)) * E[][lambda__LOWREGIME_1[1]] = 0",
                 "pH * lambda__HIGHREGIME_2[] + (beta - beta * (1 - pH)) * (kappa * E[][lambda__HIGHREGIME_1[1]] - E[][lambda__HIGHREGIME_2[1]]) - kappa * theta^-1 * yH[] = 0",
                 "pL * lambda__LOWREGIME_2[] + (beta - beta * (1 - pL)) * (kappa * E[][lambda__LOWREGIME_1[1]] - E[][lambda__LOWREGIME_2[1]]) - kappa * theta^-1 * yL[] = 0",
                 "-yH[-1] + yH[] - sigma * (iH[-1] - piH[]) + (1 - pH) * (-yH[] + yL[]) + sigma * (1 - pH) * (-piH[] + piL[]) = 0",
                 "-yL[-1] + yL[] - sigma * (iL[-1] - piL[]) + (1 - pL) * (yH[] - yL[]) + sigma * (1 - pL) * (piH[] - piL[]) = 0",
                 "UH[] + 0.5 * (pitH - pitCB + piH[])^2 - beta * E[][UH[1]] - beta * (1 - pH) * (-E[][UH[1]] + E[][UL[1]]) + 0.5 * kappa * theta^-1 * yH[]^2 = 0",
                 "UL[] + 0.5 * (-pitCB + pitL + piL[])^2 - beta * E[][UL[1]] - beta * (1 - pL) * (E[][UH[1]] - E[][UL[1]]) + 0.5 * kappa * theta^-1 * yL[]^2 = 0",
                 "-pitH + pitCB - piH[] + lambda__HIGHREGIME_1[] * (beta * (1 - phipi) - beta * (1 - phipi) * (1 - pH)) + lambda__HIGHREGIME_2[] * (sigma - sigma * (1 - pH)) + (beta - beta * (1 - pH)) * (-E[][lambda__HIGHREGIME_1[1]] + E[][lambda__HIGHREGIME_piH__lag_1[1]]) = 0",
                 "pitCB - pitL - piL[] + lambda__LOWREGIME_1[] * (beta * pL * (1 - phipi) - beta * (1 - phipi) * (1 - pL)) + lambda__LOWREGIME_2[] * (sigma - sigma * (1 - pL)) + (beta - beta * (1 - pL)) * (-E[][lambda__LOWREGIME_1[1]] + E[][lambda__LOWREGIME_piL__lag_1[1]]) = 0",
                 "-piH[-1] + log(etapi[-1]) + kappa * yH[-1] + phipi * piH__lag_1[-1] + beta * piH[] * (1 - phipi) + beta * (1 - phipi) * (1 - pH) * (-piH[] + piL[]) = 0",
                 "-piL[-1] + log(etapi[-1]) + kappa * yL[-1] + phipi * piL__lag_1[-1] + beta * pL * piL[] * (1 - phipi) + beta * (1 - phipi) * (1 - pL) * (piH[] - piL[]) = 0",
                 "-sigma * (beta - beta * (1 - pH)) * E[][lambda__HIGHREGIME_2[1]] = 0",
                 "-sigma * (beta - beta * (1 - pL)) * E[][lambda__LOWREGIME_2[1]] = 0")

# calibrating equations
calibr_equations__ <- character(0)

# variables / equations map
vareqmap__ <- sparseMatrix(i = c(1, 1, 2, 2, 3, 4, 4, 5, 5, 6,
                                 6, 6, 7, 7, 7, 8, 8, 8, 8, 8,
                                 9, 9, 9, 9, 9, 10, 10, 10, 10, 11,
                                 11, 11, 11, 12, 12, 12, 12, 13, 13, 13,
                                 13, 14, 14, 14, 14, 14, 15, 15, 15, 15,
                                 15, 16, 17),
                           j = c(10, 12, 11, 13, 1, 4, 8, 6, 9, 4,
                                 5, 14, 6, 7, 15, 2, 10, 11, 14, 15,
                                 3, 10, 11, 14, 15, 10, 14, 16, 17, 11,
                                 15, 16, 17, 4, 5, 8, 10, 6, 7, 9,
                                 11, 1, 10, 11, 12, 14, 1, 10, 11, 13,
                                 15, 5, 7),
                           x = c(1, 2, 1, 2, 3, 4, 2, 4, 2, 4,
                                 6, 2, 4, 6, 2, 1, 2, 2, 3, 2,
                                 1, 2, 2, 2, 3, 2, 2, 6, 4, 2,
                                 2, 4, 6, 6, 2, 4, 2, 6, 2, 4,
                                 2, 1, 3, 2, 1, 1, 1, 2, 3, 1,
                                 1, 4, 4),
                           dims = c(17, 17))

# variables / calibrating equations map
varcalibreqmap__ <- sparseMatrix(i = NULL, j = NULL, dims = c(0, 17))

# calibrated parameters / equations map
calibrpareqmap__ <- sparseMatrix(i = NULL, j = NULL, dims = c(17, 0))

# calibrated parameters / calibrating equations map
calibrparcalibreqmap__ <- sparseMatrix(i = NULL, j = NULL, dims = c(0, 0))

# free parameters / equations map
freepareqmap__ <- sparseMatrix(i = c(3, 4, 4, 4, 5, 5, 5, 6, 6, 6,
                                     6, 7, 7, 7, 7, 8, 8, 9, 9, 10,
                                     10, 10, 10, 10, 10, 11, 11, 11, 11, 11,
                                     11, 12, 12, 12, 12, 12, 12, 13, 13, 13,
                                     13, 13, 13, 14, 14, 14, 14, 15, 15, 15,
                                     15, 16, 16, 16, 17, 17, 17),
                               j = c(3, 1, 4, 8, 1, 4, 9, 1, 2, 8,
                                     11, 1, 2, 9, 11, 8, 10, 9, 10, 1,
                                     2, 5, 6, 8, 11, 1, 2, 6, 7, 9,
                                     11, 1, 4, 5, 6, 8, 10, 1, 4, 6,
                                     7, 9, 10, 1, 2, 4, 8, 1, 2, 4,
                                     9, 1, 8, 10, 1, 9, 10),
                               x = rep(1, 57), dims = c(17, 11))

# free parameters / calibrating equations map
freeparcalibreqmap__ <- sparseMatrix(i = NULL, j = NULL, dims = c(0, 11))

# shocks / equations map
shockeqmap__ <- sparseMatrix(i = c(3),
                             j = c(1),
                             x = rep(1, 1), dims = c(17, 1))

# steady state equations
ss_eq__ <- function(v, pc, pf)
{
    r <- numeric(17)
    r[1] = v[10] - v[12]
    r[2] = v[11] - v[13]
    r[3] = -v[1] + exp(pf[3] * log(v[1]))
    r[4] = -v[8] + pf[4] * v[4] * (pf[1] - pf[1] * (1 - pf[8]))
    r[5] = -v[9] + pf[4] * v[6] * (pf[1] - pf[1] * (1 - pf[9]))
    r[6] = pf[8] * v[5] + (pf[1] - pf[1] * (1 - pf[8])) * (-v[5] + pf[2] * v[4]) - pf[2] * pf[11]^-1 * v[14]
    r[7] = pf[9] * v[7] + (pf[1] - pf[1] * (1 - pf[9])) * (-v[7] + pf[2] * v[6]) - pf[2] * pf[11]^-1 * v[15]
    r[8] = (1 - pf[8]) * (-v[14] + v[15]) - pf[10] * (v[2] - v[10]) + pf[10] * (1 - pf[8]) * (-v[10] + v[11])
    r[9] = (1 - pf[9]) * (v[14] - v[15]) - pf[10] * (v[3] - v[11]) + pf[10] * (1 - pf[9]) * (v[10] - v[11])
    r[10] = v[16] + 0.5 * (pf[5] - pf[6] + v[10])^2 - pf[1] * v[16] - pf[1] * (1 - pf[8]) * (-v[16] + v[17]) + 0.5 * pf[2] * pf[11]^-1 * v[14]^2
    r[11] = v[17] + 0.5 * (-pf[6] + pf[7] + v[11])^2 - pf[1] * v[17] - pf[1] * (1 - pf[9]) * (v[16] - v[17]) + 0.5 * pf[2] * pf[11]^-1 * v[15]^2
    r[12] = -pf[5] + pf[6] - v[10] + v[4] * (pf[1] * (1 - pf[4]) - pf[1] * (1 - pf[4]) * (1 - pf[8])) + v[5] * (pf[10] - pf[10] * (1 - pf[8])) + (pf[1] - pf[1] * (1 - pf[8])) * (-v[4] + v[8])
    r[13] = pf[6] - pf[7] - v[11] + v[6] * (pf[1] * pf[9] * (1 - pf[4]) - pf[1] * (1 - pf[4]) * (1 - pf[9])) + v[7] * (pf[10] - pf[10] * (1 - pf[9])) + (pf[1] - pf[1] * (1 - pf[9])) * (-v[6] + v[9])
    r[14] = -v[10] + log(v[1]) + pf[2] * v[14] + pf[4] * v[12] + pf[1] * v[10] * (1 - pf[4]) + pf[1] * (1 - pf[4]) * (1 - pf[8]) * (-v[10] + v[11])
    r[15] = -v[11] + log(v[1]) + pf[2] * v[15] + pf[4] * v[13] + pf[1] * pf[9] * v[11] * (1 - pf[4]) + pf[1] * (1 - pf[4]) * (1 - pf[9]) * (v[10] - v[11])
    r[16] = -pf[10] * v[5] * (pf[1] - pf[1] * (1 - pf[8]))
    r[17] = -pf[10] * v[7] * (pf[1] - pf[1] * (1 - pf[9]))

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
    jac <- numeric(53)
    jac[1] = 1
    jac[2] = -1
    jac[3] = 1
    jac[4] = -1
    jac[5] = -1 + pf[3] * v[1]^-1 * exp(pf[3] * log(v[1]))
    jac[6] = pf[4] * (pf[1] - pf[1] * (1 - pf[8]))
    jac[7] = -1
    jac[8] = pf[4] * (pf[1] - pf[1] * (1 - pf[9]))
    jac[9] = -1
    jac[10] = pf[2] * (pf[1] - pf[1] * (1 - pf[8]))
    jac[11] = -pf[1] + pf[8] + pf[1] * (1 - pf[8])
    jac[12] = -pf[2] * pf[11]^-1
    jac[13] = pf[2] * (pf[1] - pf[1] * (1 - pf[9]))
    jac[14] = -pf[1] + pf[9] + pf[1] * (1 - pf[9])
    jac[15] = -pf[2] * pf[11]^-1
    jac[16] = -pf[10]
    jac[17] = pf[10] - pf[10] * (1 - pf[8])
    jac[18] = pf[10] * (1 - pf[8])
    jac[19] = -1 + pf[8]
    jac[20] = 1 - pf[8]
    jac[21] = -pf[10]
    jac[22] = pf[10] * (1 - pf[9])
    jac[23] = pf[10] - pf[10] * (1 - pf[9])
    jac[24] = 1 - pf[9]
    jac[25] = -1 + pf[9]
    jac[26] = pf[5] - pf[6] + v[10]
    jac[27] = pf[2] * pf[11]^-1 * v[14]
    jac[28] = 1 - pf[1] + pf[1] * (1 - pf[8])
    jac[29] = -pf[1] * (1 - pf[8])
    jac[30] = -pf[6] + pf[7] + v[11]
    jac[31] = pf[2] * pf[11]^-1 * v[15]
    jac[32] = -pf[1] * (1 - pf[9])
    jac[33] = 1 - pf[1] + pf[1] * (1 - pf[9])
    jac[34] = -pf[1] + pf[1] * (1 - pf[4]) + pf[1] * (1 - pf[8]) - pf[1] * (1 - pf[4]) * (1 - pf[8])
    jac[35] = pf[10] - pf[10] * (1 - pf[8])
    jac[36] = pf[1] - pf[1] * (1 - pf[8])
    jac[37] = -1
    jac[38] = -pf[1] + pf[1] * (1 - pf[9]) + pf[1] * pf[9] * (1 - pf[4]) - pf[1] * (1 - pf[4]) * (1 - pf[9])
    jac[39] = pf[10] - pf[10] * (1 - pf[9])
    jac[40] = pf[1] - pf[1] * (1 - pf[9])
    jac[41] = -1
    jac[42] = v[1]^-1
    jac[43] = -1 + pf[1] * (1 - pf[4]) - pf[1] * (1 - pf[4]) * (1 - pf[8])
    jac[44] = pf[1] * (1 - pf[4]) * (1 - pf[8])
    jac[45] = pf[4]
    jac[46] = pf[2]
    jac[47] = v[1]^-1
    jac[48] = pf[1] * (1 - pf[4]) * (1 - pf[9])
    jac[49] = -1 + pf[1] * pf[9] * (1 - pf[4]) - pf[1] * (1 - pf[4]) * (1 - pf[9])
    jac[50] = pf[4]
    jac[51] = pf[2]
    jac[52] = -pf[10] * (pf[1] - pf[1] * (1 - pf[8]))
    jac[53] = -pf[10] * (pf[1] - pf[1] * (1 - pf[9]))
    jacob <- sparseMatrix(i = c(1, 1, 2, 2, 3, 4, 4, 5, 5, 6,
                                6, 6, 7, 7, 7, 8, 8, 8, 8, 8,
                                9, 9, 9, 9, 9, 10, 10, 10, 10, 11,
                                11, 11, 11, 12, 12, 12, 12, 13, 13, 13,
                                13, 14, 14, 14, 14, 14, 15, 15, 15, 15,
                                15, 16, 17),
                          j = c(10, 12, 11, 13, 1, 4, 8, 6, 9, 4,
                                5, 14, 6, 7, 15, 2, 10, 11, 14, 15,
                                3, 10, 11, 14, 15, 10, 14, 16, 17, 11,
                                15, 16, 17, 4, 5, 8, 10, 6, 7, 9,
                                11, 1, 10, 11, 12, 14, 1, 10, 11, 13,
                                15, 5, 7),
                          x = jac, dims = c(17, 17))

    return(jacob)
}

# 1st order perturbation
pert1__ <- function(v, pc, pf)
{
    Atm1x <- numeric(15)
    Atm1x[1] = 1
    Atm1x[2] = 1
    Atm1x[3] = pf[3] * v[1]^-1 * exp(pf[3] * log(v[1]))
    Atm1x[4] = -pf[10]
    Atm1x[5] = -1
    Atm1x[6] = -pf[10]
    Atm1x[7] = -1
    Atm1x[8] = v[1]^-1
    Atm1x[9] = -1
    Atm1x[10] = pf[4]
    Atm1x[11] = pf[2]
    Atm1x[12] = v[1]^-1
    Atm1x[13] = -1
    Atm1x[14] = pf[4]
    Atm1x[15] = pf[2]
    Atm1 <- sparseMatrix(i = c(1, 2, 3, 8, 8, 9, 9, 14, 14, 14,
                               14, 15, 15, 15, 15),
                         j = c(10, 11, 1, 2, 14, 3, 15, 1, 10, 12,
                               14, 1, 11, 13, 15),
                         x = Atm1x, dims = c(17, 17))

    Atx <- numeric(33)
    Atx[1] = -1
    Atx[2] = -1
    Atx[3] = -1
    Atx[4] = -1
    Atx[5] = -1
    Atx[6] = pf[8]
    Atx[7] = -pf[2] * pf[11]^-1
    Atx[8] = pf[9]
    Atx[9] = -pf[2] * pf[11]^-1
    Atx[10] = pf[10] - pf[10] * (1 - pf[8])
    Atx[11] = pf[10] * (1 - pf[8])
    Atx[12] = pf[8]
    Atx[13] = 1 - pf[8]
    Atx[14] = pf[10] * (1 - pf[9])
    Atx[15] = pf[10] - pf[10] * (1 - pf[9])
    Atx[16] = 1 - pf[9]
    Atx[17] = pf[9]
    Atx[18] = pf[5] - pf[6] + v[10]
    Atx[19] = pf[2] * pf[11]^-1 * v[14]
    Atx[20] = 1
    Atx[21] = -pf[6] + pf[7] + v[11]
    Atx[22] = pf[2] * pf[11]^-1 * v[15]
    Atx[23] = 1
    Atx[24] = pf[1] * (1 - pf[4]) - pf[1] * (1 - pf[4]) * (1 - pf[8])
    Atx[25] = pf[10] - pf[10] * (1 - pf[8])
    Atx[26] = -1
    Atx[27] = pf[1] * pf[9] * (1 - pf[4]) - pf[1] * (1 - pf[4]) * (1 - pf[9])
    Atx[28] = pf[10] - pf[10] * (1 - pf[9])
    Atx[29] = -1
    Atx[30] = pf[1] * (1 - pf[4]) - pf[1] * (1 - pf[4]) * (1 - pf[8])
    Atx[31] = pf[1] * (1 - pf[4]) * (1 - pf[8])
    Atx[32] = pf[1] * (1 - pf[4]) * (1 - pf[9])
    Atx[33] = pf[1] * pf[9] * (1 - pf[4]) - pf[1] * (1 - pf[4]) * (1 - pf[9])
    At <- sparseMatrix(i = c(1, 2, 3, 4, 5, 6, 6, 7, 7, 8,
                             8, 8, 8, 9, 9, 9, 9, 10, 10, 10,
                             11, 11, 11, 12, 12, 12, 13, 13, 13, 14,
                             14, 15, 15),
                       j = c(12, 13, 1, 8, 9, 5, 14, 7, 15, 10,
                             11, 14, 15, 10, 11, 14, 15, 10, 14, 16,
                             11, 15, 17, 4, 5, 10, 6, 7, 11, 10,
                             11, 10, 11),
                         x = Atx, dims = c(17, 17))

    Atp1x <- numeric(16)
    Atp1x[1] = pf[4] * (pf[1] - pf[1] * (1 - pf[8]))
    Atp1x[2] = pf[4] * (pf[1] - pf[1] * (1 - pf[9]))
    Atp1x[3] = pf[2] * (pf[1] - pf[1] * (1 - pf[8]))
    Atp1x[4] = -pf[1] + pf[1] * (1 - pf[8])
    Atp1x[5] = pf[2] * (pf[1] - pf[1] * (1 - pf[9]))
    Atp1x[6] = -pf[1] + pf[1] * (1 - pf[9])
    Atp1x[7] = -pf[1] + pf[1] * (1 - pf[8])
    Atp1x[8] = -pf[1] * (1 - pf[8])
    Atp1x[9] = -pf[1] * (1 - pf[9])
    Atp1x[10] = -pf[1] + pf[1] * (1 - pf[9])
    Atp1x[11] = -pf[1] + pf[1] * (1 - pf[8])
    Atp1x[12] = pf[1] - pf[1] * (1 - pf[8])
    Atp1x[13] = -pf[1] + pf[1] * (1 - pf[9])
    Atp1x[14] = pf[1] - pf[1] * (1 - pf[9])
    Atp1x[15] = -pf[10] * (pf[1] - pf[1] * (1 - pf[8]))
    Atp1x[16] = -pf[10] * (pf[1] - pf[1] * (1 - pf[9]))
    Atp1 <- sparseMatrix(i = c(4, 5, 6, 6, 7, 7, 10, 10, 11, 11,
                               12, 12, 13, 13, 16, 17),
                         j = c(4, 6, 4, 5, 6, 7, 16, 17, 16, 17,
                               4, 8, 6, 9, 5, 7),
                         x = Atp1x, dims = c(17, 17))

    Aepsx <- numeric(1)
    Aepsx[1] = exp(pf[3] * log(v[1]))
    Aeps <- sparseMatrix(i = c(3),
                         j = c(1),
                         x = Aepsx, dims = c(17, 1))

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
