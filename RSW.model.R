# Generated on 2023-07-17 16:52:04 by gEcon ver. 1.2.1 (2023-01-18)
# http://gecon.r-forge.r-project.org/

# Model name: RSW

# info
info__ <- c("RSW", "C:/Users/fm007/Documents/GitHub/OptimalRSW/RSW.gcn", "2023-07-17 16:52:04", "false")

# index sets
index_sets__ <- list()

# variables
variables__ <- c("iH",
                 "iL",
                 "lambda__HIGHREGIME_1",
                 "lambda__HIGHREGIME_2",
                 "lambda__LOWREGIME_1",
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
                     "\\lambda^{\\mathrm{HIGHREGIME}^{\\mathrm{2}}}",
                     "\\lambda^{\\mathrm{LOWREGIME}^{\\mathrm{1}}}",
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
                  "pH",
                  "pL",
                  "sigma",
                  "theta")

parameters_tex__ <- c("\\beta",
                     "\\kappa",
                     "\\phi",
                     "{p\\!H}",
                     "{p\\!L}",
                     "\\sigma",
                     "\\theta")

# free parameters
parameters_free__ <- c("beta",
                       "kappa",
                       "phi",
                       "pH",
                       "pL",
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
                 "-piH[-1] + beta * (pH * piH[] + piL[] * (1 - pH)) + kappa * yH[-1] = 0",
                 "-piL[-1] + beta * (pL * piL[] + piH[] * (1 - pL)) + kappa * yL[-1] = 0",
                 "pH * lambda__HIGHREGIME_2[] + beta * pH * (kappa * E[][lambda__HIGHREGIME_1[1]] - E[][lambda__HIGHREGIME_2[1]]) - kappa * theta^-1 * yH[] = 0",
                 "pL * lambda__LOWREGIME_2[] + beta * pL * (kappa * E[][lambda__LOWREGIME_1[1]] - E[][lambda__LOWREGIME_2[1]]) - kappa * theta^-1 * yL[] = 0",
                 "-yH[-1] + pH * yH[] - sigma * (iH[-1] - rn[-1] - pH * piH[] - piL[] * (1 - pH)) + yL[] * (1 - pH) = 0",
                 "-yL[-1] + pL * yL[] - sigma * (iL[-1] - rn[-1] - pL * piL[] - piH[] * (1 - pL)) + yH[] * (1 - pL) = 0",
                 "-piH[] + beta * pH * lambda__HIGHREGIME_1[] - beta * pH * E[][lambda__HIGHREGIME_1[1]] + pH * sigma * lambda__HIGHREGIME_2[] = 0",
                 "-piL[] + beta * pL * lambda__LOWREGIME_1[] - beta * pL * E[][lambda__LOWREGIME_1[1]] + pL * sigma * lambda__LOWREGIME_2[] = 0",
                 "UH[] + 0.5 * piH[]^2 - beta * (pH * E[][UH[1]] + (1 - pH) * E[][UL[1]]) + 0.5 * kappa * theta^-1 * yH[]^2 = 0",
                 "UL[] + 0.5 * piL[]^2 - beta * (pL * E[][UL[1]] + (1 - pL) * E[][UH[1]]) + 0.5 * kappa * theta^-1 * yL[]^2 = 0",
                 "-beta * pH * sigma * E[][lambda__HIGHREGIME_2[1]] = 0",
                 "-beta * pL * sigma * E[][lambda__LOWREGIME_2[1]] = 0")

# calibrating equations
calibr_equations__ <- character(0)

# variables / equations map
vareqmap__ <- sparseMatrix(i = c(1, 2, 2, 2, 3, 3, 3, 4, 4, 4,
                                 5, 5, 5, 6, 6, 6, 6, 6, 6, 7,
                                 7, 7, 7, 7, 7, 8, 8, 8, 9, 9,
                                 9, 10, 10, 10, 10, 11, 11, 11, 11, 12,
                                 13),
                           j = c(9, 7, 8, 10, 7, 8, 11, 3, 4, 10,
                                 5, 6, 11, 1, 7, 8, 9, 10, 11, 2,
                                 7, 8, 9, 10, 11, 3, 4, 7, 5, 6,
                                 8, 7, 10, 12, 13, 8, 11, 12, 13, 4,
                                 6),
                           x = c(3, 3, 2, 1, 2, 3, 1, 4, 6, 2,
                                 4, 6, 2, 1, 2, 2, 1, 3, 2, 1,
                                 2, 2, 1, 2, 3, 6, 2, 2, 6, 2,
                                 2, 2, 2, 6, 4, 2, 2, 4, 6, 4,
                                 4),
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
                                     8, 8, 9, 9, 9, 10, 10, 10, 10, 11,
                                     11, 11, 11, 12, 12, 12, 13, 13, 13),
                               j = c(3, 1, 2, 4, 1, 2, 5, 1, 2, 4,
                                     7, 1, 2, 5, 7, 4, 6, 5, 6, 1,
                                     4, 6, 1, 5, 6, 1, 2, 4, 7, 1,
                                     2, 5, 7, 1, 4, 6, 1, 5, 6),
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
    r[2] = -v[7] + pf[1] * (pf[4] * v[7] + v[8] * (1 - pf[4])) + pf[2] * v[10]
    r[3] = -v[8] + pf[1] * (pf[5] * v[8] + v[7] * (1 - pf[5])) + pf[2] * v[11]
    r[4] = pf[4] * v[4] + pf[1] * pf[4] * (-v[4] + pf[2] * v[3]) - pf[2] * pf[7]^-1 * v[10]
    r[5] = pf[5] * v[6] + pf[1] * pf[5] * (-v[6] + pf[2] * v[5]) - pf[2] * pf[7]^-1 * v[11]
    r[6] = -v[10] + pf[4] * v[10] - pf[6] * (v[1] - v[9] - pf[4] * v[7] - v[8] * (1 - pf[4])) + v[11] * (1 - pf[4])
    r[7] = -v[11] + pf[5] * v[11] - pf[6] * (v[2] - v[9] - pf[5] * v[8] - v[7] * (1 - pf[5])) + v[10] * (1 - pf[5])
    r[8] = -v[7] + pf[4] * pf[6] * v[4]
    r[9] = -v[8] + pf[5] * pf[6] * v[6]
    r[10] = v[12] + 0.5 * v[7]^2 - pf[1] * (pf[4] * v[12] + v[13] * (1 - pf[4])) + 0.5 * pf[2] * pf[7]^-1 * v[10]^2
    r[11] = v[13] + 0.5 * v[8]^2 - pf[1] * (pf[5] * v[13] + v[12] * (1 - pf[5])) + 0.5 * pf[2] * pf[7]^-1 * v[11]^2
    r[12] = -pf[1] * pf[4] * pf[6] * v[4]
    r[13] = -pf[1] * pf[5] * pf[6] * v[6]

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
    jac <- numeric(39)
    jac[1] = -1 + pf[3] * v[9]^-1 * exp(pf[3] * log(v[9]))
    jac[2] = -1 + pf[1] * pf[4]
    jac[3] = pf[1] * (1 - pf[4])
    jac[4] = pf[2]
    jac[5] = pf[1] * (1 - pf[5])
    jac[6] = -1 + pf[1] * pf[5]
    jac[7] = pf[2]
    jac[8] = pf[1] * pf[2] * pf[4]
    jac[9] = pf[4] - pf[1] * pf[4]
    jac[10] = -pf[2] * pf[7]^-1
    jac[11] = pf[1] * pf[2] * pf[5]
    jac[12] = pf[5] - pf[1] * pf[5]
    jac[13] = -pf[2] * pf[7]^-1
    jac[14] = -pf[6]
    jac[15] = pf[4] * pf[6]
    jac[16] = -pf[6] * (-1 + pf[4])
    jac[17] = pf[6]
    jac[18] = -1 + pf[4]
    jac[19] = 1 - pf[4]
    jac[20] = -pf[6]
    jac[21] = -pf[6] * (-1 + pf[5])
    jac[22] = pf[5] * pf[6]
    jac[23] = pf[6]
    jac[24] = 1 - pf[5]
    jac[25] = -1 + pf[5]
    jac[26] = pf[4] * pf[6]
    jac[27] = -1
    jac[28] = pf[5] * pf[6]
    jac[29] = -1
    jac[30] = v[7]
    jac[31] = pf[2] * pf[7]^-1 * v[10]
    jac[32] = 1 - pf[1] * pf[4]
    jac[33] = -pf[1] * (1 - pf[4])
    jac[34] = v[8]
    jac[35] = pf[2] * pf[7]^-1 * v[11]
    jac[36] = -pf[1] * (1 - pf[5])
    jac[37] = 1 - pf[1] * pf[5]
    jac[38] = -pf[1] * pf[4] * pf[6]
    jac[39] = -pf[1] * pf[5] * pf[6]
    jacob <- sparseMatrix(i = c(1, 2, 2, 2, 3, 3, 3, 4, 4, 4,
                                5, 5, 5, 6, 6, 6, 6, 6, 6, 7,
                                7, 7, 7, 7, 7, 8, 8, 9, 9, 10,
                                10, 10, 10, 11, 11, 11, 11, 12, 13),
                          j = c(9, 7, 8, 10, 7, 8, 11, 3, 4, 10,
                                5, 6, 11, 1, 7, 8, 9, 10, 11, 2,
                                7, 8, 9, 10, 11, 4, 7, 6, 8, 7,
                                10, 12, 13, 8, 11, 12, 13, 4, 6),
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

    Atx <- numeric(29)
    Atx[1] = -1
    Atx[2] = pf[1] * pf[4]
    Atx[3] = pf[1] * (1 - pf[4])
    Atx[4] = pf[1] * (1 - pf[5])
    Atx[5] = pf[1] * pf[5]
    Atx[6] = pf[4]
    Atx[7] = -pf[2] * pf[7]^-1
    Atx[8] = pf[5]
    Atx[9] = -pf[2] * pf[7]^-1
    Atx[10] = pf[4] * pf[6]
    Atx[11] = -pf[6] * (-1 + pf[4])
    Atx[12] = pf[4]
    Atx[13] = 1 - pf[4]
    Atx[14] = -pf[6] * (-1 + pf[5])
    Atx[15] = pf[5] * pf[6]
    Atx[16] = 1 - pf[5]
    Atx[17] = pf[5]
    Atx[18] = pf[1] * pf[4]
    Atx[19] = pf[4] * pf[6]
    Atx[20] = -1
    Atx[21] = pf[1] * pf[5]
    Atx[22] = pf[5] * pf[6]
    Atx[23] = -1
    Atx[24] = v[7]
    Atx[25] = pf[2] * pf[7]^-1 * v[10]
    Atx[26] = 1
    Atx[27] = v[8]
    Atx[28] = pf[2] * pf[7]^-1 * v[11]
    Atx[29] = 1
    At <- sparseMatrix(i = c(1, 2, 2, 3, 3, 4, 4, 5, 5, 6,
                             6, 6, 6, 7, 7, 7, 7, 8, 8, 8,
                             9, 9, 9, 10, 10, 10, 11, 11, 11),
                       j = c(9, 7, 8, 7, 8, 4, 10, 6, 11, 7,
                             8, 10, 11, 7, 8, 10, 11, 3, 4, 7,
                             5, 6, 8, 7, 10, 12, 8, 11, 13),
                         x = Atx, dims = c(13, 13))

    Atp1x <- numeric(12)
    Atp1x[1] = pf[1] * pf[2] * pf[4]
    Atp1x[2] = -pf[1] * pf[4]
    Atp1x[3] = pf[1] * pf[2] * pf[5]
    Atp1x[4] = -pf[1] * pf[5]
    Atp1x[5] = -pf[1] * pf[4]
    Atp1x[6] = -pf[1] * pf[5]
    Atp1x[7] = -pf[1] * pf[4]
    Atp1x[8] = -pf[1] * (1 - pf[4])
    Atp1x[9] = -pf[1] * (1 - pf[5])
    Atp1x[10] = -pf[1] * pf[5]
    Atp1x[11] = -pf[1] * pf[4] * pf[6]
    Atp1x[12] = -pf[1] * pf[5] * pf[6]
    Atp1 <- sparseMatrix(i = c(4, 4, 5, 5, 8, 9, 10, 10, 11, 11,
                               12, 13),
                         j = c(3, 4, 5, 6, 3, 5, 12, 13, 12, 13,
                               4, 6),
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
