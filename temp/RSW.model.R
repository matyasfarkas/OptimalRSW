# Generated on 2025-02-07 10:12:54 by gEcon ver. 1.2.1 (2023-01-18)
# http://gecon.r-forge.r-project.org/

# Model name: RSW

# info
info__ <- c("RSW", "C:/Users/fm007/Documents/GitHub/OptimalRSW/temp/RSW.gcn", "2025-02-07 10:12:54", "false")

# index sets
index_sets__ <- list()

# variables
variables__ <- c("etapi",
                 "etag",
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
                     "{e\\!t\\!a\\!g}",
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
shocks__ <- c("epsilon_pi",
              "epsilon_g")

shocks_tex__ <- c("\\epsilon^{\\pi}",
                  "\\epsilon^{\\mathrm{g}}")

# parameters
parameters__ <- c("beta",
                  "kappa",
                  "phi",
                  "phig",
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
                     "{p\\!h\\!i\\!g}",
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
                       "phig",
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
                           0.99,
                           0,
                           0,
                           2,
                           0.99,
                           0.99,
                           1,
                           6)

# equations
equations__ <- c("-etapi[] + exp(epsilon_pi[] + phi * log(etapi[-1])) = 0",
                 "-etag[] + exp(epsilon_g[] + phig * log(etag[-1])) = 0",
                 "pH * lambda__HIGHREGIME_2[] + beta * pH * (kappa * E[][lambda__HIGHREGIME_1[1]] - E[][lambda__HIGHREGIME_2[1]]) - kappa * theta^-1 * yH[] = 0",
                 "pL * lambda__LOWREGIME_2[] + beta * pL * (kappa * E[][lambda__LOWREGIME_1[1]] - E[][lambda__LOWREGIME_2[1]]) - kappa * theta^-1 * yL[] = 0",
                 "-piH[-1] + log(etapi[-1]) + beta * (pH * piH[] + piL[] * (1 - pH)) + kappa * yH[-1] = 0",
                 "-piL[-1] + log(etapi[-1]) + beta * (pL * piL[] + piH[] * (1 - pL)) + kappa * yL[-1] = 0",
                 "UH[] + 0.5 * (pitH - pitCB + piH[])^2 - beta * (pH * E[][UH[1]] + (1 - pH) * E[][UL[1]]) + 0.5 * kappa * theta^-1 * yH[]^2 = 0",
                 "UL[] + 0.5 * (-pitCB + pitL + piL[])^2 - beta * (pL * E[][UL[1]] + (1 - pL) * E[][UH[1]]) + 0.5 * kappa * theta^-1 * yL[]^2 = 0",
                 "-yH[-1] + log(etag[-1]) + pH * yH[] - sigma * (iH[-1] - pH * piH[] - piL[] * (1 - pH)) + yL[] * (1 - pH) = 0",
                 "-yL[-1] + log(etag[-1]) + pL * yL[] - sigma * (iL[-1] - pL * piL[] - piH[] * (1 - pL)) + yH[] * (1 - pL) = 0",
                 "-pitH + pitCB - piH[] + beta * pH * lambda__HIGHREGIME_1[] - beta * pH * E[][lambda__HIGHREGIME_1[1]] + pH * sigma * lambda__HIGHREGIME_2[] = 0",
                 "pitCB - pitL - piL[] + beta * pL * lambda__LOWREGIME_1[] - beta * pL * E[][lambda__LOWREGIME_1[1]] + pL * sigma * lambda__LOWREGIME_2[] = 0",
                 "-beta * pH * sigma * E[][lambda__HIGHREGIME_2[1]] = 0",
                 "-beta * pL * sigma * E[][lambda__LOWREGIME_2[1]] = 0")

# calibrating equations
calibr_equations__ <- character(0)

# variables / equations map
vareqmap__ <- sparseMatrix(i = c(1, 2, 3, 3, 3, 4, 4, 4, 5, 5,
                                 5, 5, 6, 6, 6, 6, 7, 7, 7, 7,
                                 8, 8, 8, 8, 9, 9, 9, 9, 9, 9,
                                 10, 10, 10, 10, 10, 10, 11, 11, 11, 12,
                                 12, 12, 13, 14),
                           j = c(1, 2, 5, 6, 11, 7, 8, 12, 1, 9,
                                 10, 11, 1, 9, 10, 12, 9, 11, 13, 14,
                                 10, 12, 13, 14, 2, 3, 9, 10, 11, 12,
                                 2, 4, 9, 10, 11, 12, 5, 6, 9, 7,
                                 8, 10, 6, 8),
                           x = c(3, 3, 4, 6, 2, 4, 6, 2, 1, 3,
                                 2, 1, 1, 2, 3, 1, 2, 2, 6, 4,
                                 2, 2, 4, 6, 1, 1, 2, 2, 3, 2,
                                 1, 1, 2, 2, 2, 3, 6, 2, 2, 6,
                                 2, 2, 4, 4),
                           dims = c(14, 14))

# variables / calibrating equations map
varcalibreqmap__ <- sparseMatrix(i = NULL, j = NULL, dims = c(0, 14))

# calibrated parameters / equations map
calibrpareqmap__ <- sparseMatrix(i = NULL, j = NULL, dims = c(14, 0))

# calibrated parameters / calibrating equations map
calibrparcalibreqmap__ <- sparseMatrix(i = NULL, j = NULL, dims = c(0, 0))

# free parameters / equations map
freepareqmap__ <- sparseMatrix(i = c(1, 2, 3, 3, 3, 3, 4, 4, 4, 4,
                                     5, 5, 5, 6, 6, 6, 7, 7, 7, 7,
                                     7, 7, 8, 8, 8, 8, 8, 8, 9, 9,
                                     10, 10, 11, 11, 11, 11, 11, 12, 12, 12,
                                     12, 12, 13, 13, 13, 14, 14, 14),
                               j = c(3, 4, 1, 2, 8, 11, 1, 2, 9, 11,
                                     1, 2, 8, 1, 2, 9, 1, 2, 5, 6,
                                     8, 11, 1, 2, 6, 7, 9, 11, 8, 10,
                                     9, 10, 1, 5, 6, 8, 10, 1, 6, 7,
                                     9, 10, 1, 8, 10, 1, 9, 10),
                               x = rep(1, 48), dims = c(14, 11))

# free parameters / calibrating equations map
freeparcalibreqmap__ <- sparseMatrix(i = NULL, j = NULL, dims = c(0, 11))

# shocks / equations map
shockeqmap__ <- sparseMatrix(i = c(1, 2),
                             j = c(1, 2),
                             x = rep(1, 2), dims = c(14, 2))

# steady state equations
ss_eq__ <- function(v, pc, pf)
{
    r <- numeric(14)
    r[1] = -v[1] + exp(pf[3] * log(v[1]))
    r[2] = -v[2] + exp(pf[4] * log(v[2]))
    r[3] = pf[8] * v[6] + pf[1] * pf[8] * (-v[6] + pf[2] * v[5]) - pf[2] * pf[11]^-1 * v[11]
    r[4] = pf[9] * v[8] + pf[1] * pf[9] * (-v[8] + pf[2] * v[7]) - pf[2] * pf[11]^-1 * v[12]
    r[5] = -v[9] + log(v[1]) + pf[1] * (pf[8] * v[9] + v[10] * (1 - pf[8])) + pf[2] * v[11]
    r[6] = -v[10] + log(v[1]) + pf[1] * (pf[9] * v[10] + v[9] * (1 - pf[9])) + pf[2] * v[12]
    r[7] = v[13] + 0.5 * (pf[5] - pf[6] + v[9])^2 - pf[1] * (pf[8] * v[13] + v[14] * (1 - pf[8])) + 0.5 * pf[2] * pf[11]^-1 * v[11]^2
    r[8] = v[14] + 0.5 * (-pf[6] + pf[7] + v[10])^2 - pf[1] * (pf[9] * v[14] + v[13] * (1 - pf[9])) + 0.5 * pf[2] * pf[11]^-1 * v[12]^2
    r[9] = -v[11] + log(v[2]) + pf[8] * v[11] - pf[10] * (v[3] - pf[8] * v[9] - v[10] * (1 - pf[8])) + v[12] * (1 - pf[8])
    r[10] = -v[12] + log(v[2]) + pf[9] * v[12] - pf[10] * (v[4] - pf[9] * v[10] - v[9] * (1 - pf[9])) + v[11] * (1 - pf[9])
    r[11] = -pf[5] + pf[6] - v[9] + pf[8] * pf[10] * v[6]
    r[12] = pf[6] - pf[7] - v[10] + pf[9] * pf[10] * v[8]
    r[13] = -pf[1] * pf[8] * pf[10] * v[6]
    r[14] = -pf[1] * pf[9] * pf[10] * v[8]

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
    jac <- numeric(42)
    jac[1] = -1 + pf[3] * v[1]^-1 * exp(pf[3] * log(v[1]))
    jac[2] = -1 + pf[4] * v[2]^-1 * exp(pf[4] * log(v[2]))
    jac[3] = pf[1] * pf[2] * pf[8]
    jac[4] = pf[8] - pf[1] * pf[8]
    jac[5] = -pf[2] * pf[11]^-1
    jac[6] = pf[1] * pf[2] * pf[9]
    jac[7] = pf[9] - pf[1] * pf[9]
    jac[8] = -pf[2] * pf[11]^-1
    jac[9] = v[1]^-1
    jac[10] = -1 + pf[1] * pf[8]
    jac[11] = pf[1] * (1 - pf[8])
    jac[12] = pf[2]
    jac[13] = v[1]^-1
    jac[14] = pf[1] * (1 - pf[9])
    jac[15] = -1 + pf[1] * pf[9]
    jac[16] = pf[2]
    jac[17] = pf[5] - pf[6] + v[9]
    jac[18] = pf[2] * pf[11]^-1 * v[11]
    jac[19] = 1 - pf[1] * pf[8]
    jac[20] = -pf[1] * (1 - pf[8])
    jac[21] = -pf[6] + pf[7] + v[10]
    jac[22] = pf[2] * pf[11]^-1 * v[12]
    jac[23] = -pf[1] * (1 - pf[9])
    jac[24] = 1 - pf[1] * pf[9]
    jac[25] = v[2]^-1
    jac[26] = -pf[10]
    jac[27] = pf[8] * pf[10]
    jac[28] = -pf[10] * (-1 + pf[8])
    jac[29] = -1 + pf[8]
    jac[30] = 1 - pf[8]
    jac[31] = v[2]^-1
    jac[32] = -pf[10]
    jac[33] = -pf[10] * (-1 + pf[9])
    jac[34] = pf[9] * pf[10]
    jac[35] = 1 - pf[9]
    jac[36] = -1 + pf[9]
    jac[37] = pf[8] * pf[10]
    jac[38] = -1
    jac[39] = pf[9] * pf[10]
    jac[40] = -1
    jac[41] = -pf[1] * pf[8] * pf[10]
    jac[42] = -pf[1] * pf[9] * pf[10]
    jacob <- sparseMatrix(i = c(1, 2, 3, 3, 3, 4, 4, 4, 5, 5,
                                5, 5, 6, 6, 6, 6, 7, 7, 7, 7,
                                8, 8, 8, 8, 9, 9, 9, 9, 9, 9,
                                10, 10, 10, 10, 10, 10, 11, 11, 12, 12,
                                13, 14),
                          j = c(1, 2, 5, 6, 11, 7, 8, 12, 1, 9,
                                10, 11, 1, 9, 10, 12, 9, 11, 13, 14,
                                10, 12, 13, 14, 2, 3, 9, 10, 11, 12,
                                2, 4, 9, 10, 11, 12, 6, 9, 8, 10,
                                6, 8),
                          x = jac, dims = c(14, 14))

    return(jacob)
}

# 1st order perturbation
pert1__ <- function(v, pc, pf)
{
    Atm1x <- numeric(14)
    Atm1x[1] = pf[3] * v[1]^-1 * exp(pf[3] * log(v[1]))
    Atm1x[2] = pf[4] * v[2]^-1 * exp(pf[4] * log(v[2]))
    Atm1x[3] = v[1]^-1
    Atm1x[4] = -1
    Atm1x[5] = pf[2]
    Atm1x[6] = v[1]^-1
    Atm1x[7] = -1
    Atm1x[8] = pf[2]
    Atm1x[9] = v[2]^-1
    Atm1x[10] = -pf[10]
    Atm1x[11] = -1
    Atm1x[12] = v[2]^-1
    Atm1x[13] = -pf[10]
    Atm1x[14] = -1
    Atm1 <- sparseMatrix(i = c(1, 2, 5, 5, 5, 6, 6, 6, 9, 9,
                               9, 10, 10, 10),
                         j = c(1, 2, 1, 9, 11, 1, 10, 12, 2, 3,
                               11, 2, 4, 12),
                         x = Atm1x, dims = c(14, 14))

    Atx <- numeric(30)
    Atx[1] = -1
    Atx[2] = -1
    Atx[3] = pf[8]
    Atx[4] = -pf[2] * pf[11]^-1
    Atx[5] = pf[9]
    Atx[6] = -pf[2] * pf[11]^-1
    Atx[7] = pf[1] * pf[8]
    Atx[8] = pf[1] * (1 - pf[8])
    Atx[9] = pf[1] * (1 - pf[9])
    Atx[10] = pf[1] * pf[9]
    Atx[11] = pf[5] - pf[6] + v[9]
    Atx[12] = pf[2] * pf[11]^-1 * v[11]
    Atx[13] = 1
    Atx[14] = -pf[6] + pf[7] + v[10]
    Atx[15] = pf[2] * pf[11]^-1 * v[12]
    Atx[16] = 1
    Atx[17] = pf[8] * pf[10]
    Atx[18] = -pf[10] * (-1 + pf[8])
    Atx[19] = pf[8]
    Atx[20] = 1 - pf[8]
    Atx[21] = -pf[10] * (-1 + pf[9])
    Atx[22] = pf[9] * pf[10]
    Atx[23] = 1 - pf[9]
    Atx[24] = pf[9]
    Atx[25] = pf[1] * pf[8]
    Atx[26] = pf[8] * pf[10]
    Atx[27] = -1
    Atx[28] = pf[1] * pf[9]
    Atx[29] = pf[9] * pf[10]
    Atx[30] = -1
    At <- sparseMatrix(i = c(1, 2, 3, 3, 4, 4, 5, 5, 6, 6,
                             7, 7, 7, 8, 8, 8, 9, 9, 9, 9,
                             10, 10, 10, 10, 11, 11, 11, 12, 12, 12),
                       j = c(1, 2, 6, 11, 8, 12, 9, 10, 9, 10,
                             9, 11, 13, 10, 12, 14, 9, 10, 11, 12,
                             9, 10, 11, 12, 5, 6, 9, 7, 8, 10),
                         x = Atx, dims = c(14, 14))

    Atp1x <- numeric(12)
    Atp1x[1] = pf[1] * pf[2] * pf[8]
    Atp1x[2] = -pf[1] * pf[8]
    Atp1x[3] = pf[1] * pf[2] * pf[9]
    Atp1x[4] = -pf[1] * pf[9]
    Atp1x[5] = -pf[1] * pf[8]
    Atp1x[6] = -pf[1] * (1 - pf[8])
    Atp1x[7] = -pf[1] * (1 - pf[9])
    Atp1x[8] = -pf[1] * pf[9]
    Atp1x[9] = -pf[1] * pf[8]
    Atp1x[10] = -pf[1] * pf[9]
    Atp1x[11] = -pf[1] * pf[8] * pf[10]
    Atp1x[12] = -pf[1] * pf[9] * pf[10]
    Atp1 <- sparseMatrix(i = c(3, 3, 4, 4, 7, 7, 8, 8, 11, 12,
                               13, 14),
                         j = c(5, 6, 7, 8, 13, 14, 13, 14, 5, 7,
                               6, 8),
                         x = Atp1x, dims = c(14, 14))

    Aepsx <- numeric(2)
    Aepsx[1] = exp(pf[3] * log(v[1]))
    Aepsx[2] = exp(pf[4] * log(v[2]))
    Aeps <- sparseMatrix(i = c(1, 2),
                         j = c(1, 2),
                         x = Aepsx, dims = c(14, 2))

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
