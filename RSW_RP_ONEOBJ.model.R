# Generated on 2023-12-01 17:48:04 by gEcon ver. 1.2.1 (2023-01-18)
# http://gecon.r-forge.r-project.org/

# Model name: RSW_RP_ONEOBJ

# info
info__ <- c("RSW_RP_ONEOBJ", "C:/Users/fm007/Documents/GitHub/OptimalRSW/RSW_RP_ONEOBJ.gcn", "2023-12-01 17:48:04", "false")

# index sets
index_sets__ <- list()

# variables
variables__ <- c("etapi",
                 "i",
                 "lambda__OPTIMALMP_1",
                 "lambda__OPTIMALMP_2",
                 "lambda__OPTIMALMP_3",
                 "lambda__OPTIMALMP_4",
                 "piH",
                 "piL",
                 "yH",
                 "yL",
                 "U")

variables_tex__ <- c("{e\\!t\\!a\\!p\\!i}",
                     "i",
                     "\\lambda^{\\mathrm{OPTIMALMP}^{\\mathrm{1}}}",
                     "\\lambda^{\\mathrm{OPTIMALMP}^{\\mathrm{2}}}",
                     "\\lambda^{\\mathrm{OPTIMALMP}^{\\mathrm{3}}}",
                     "\\lambda^{\\mathrm{OPTIMALMP}^{\\mathrm{4}}}",
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
                 "beta * (kappa * E[][lambda__OPTIMALMP_1[1]] - E[][lambda__OPTIMALMP_2[1]]) + pH * lambda__OPTIMALMP_2[] + lambda__OPTIMALMP_4[] * (1 - pL) - 0.5 * kappa * theta^-1 * yH[] = 0",
                 "beta * (kappa * E[][lambda__OPTIMALMP_3[1]] - E[][lambda__OPTIMALMP_4[1]]) + pL * lambda__OPTIMALMP_4[] + lambda__OPTIMALMP_2[] * (1 - pH) - 0.5 * kappa * theta^-1 * yL[] = 0",
                 "-piH[-1] + log(etapi[-1]) + beta * piH[] + kappa * yH[-1] + beta * (1 - pH) * (-piH[] + piL[]) = 0",
                 "-piL[-1] + log(etapi[-1]) + kappa * yL[-1] + beta * pL * piL[] + beta * (1 - pL) * (piH[] - piL[]) = 0",
                 "-yH[-1] + yH[] - sigma * (i[-1] - piH[]) + (1 - pH) * (-yH[] + yL[]) + sigma * (1 - pH) * (-piH[] + piL[]) = 0",
                 "-yL[-1] + yL[] - sigma * (i[-1] - piL[]) + (1 - pL) * (yH[] - yL[]) + sigma * (1 - pL) * (piH[] - piL[]) = 0",
                 "U[] + 0.25 * (pitH - pitCB + piH[])^2 + 0.25 * (-pitCB + pitL + piL[])^2 - beta * E[][U[1]] + 0.25 * kappa * theta^-1 * yH[]^2 + 0.25 * kappa * theta^-1 * yL[]^2 = 0",
                 "-0.5 * pitH + 0.5 * pitCB - 0.5 * piH[] - beta * E[][lambda__OPTIMALMP_1[1]] + lambda__OPTIMALMP_1[] * (beta - beta * (1 - pH)) + lambda__OPTIMALMP_2[] * (sigma - sigma * (1 - pH)) + beta * lambda__OPTIMALMP_3[] * (1 - pL) + sigma * lambda__OPTIMALMP_4[] * (1 - pL) = 0",
                 "0.5 * pitCB - 0.5 * pitL - 0.5 * piL[] - beta * E[][lambda__OPTIMALMP_3[1]] + lambda__OPTIMALMP_3[] * (beta * pL - beta * (1 - pL)) + lambda__OPTIMALMP_4[] * (sigma - sigma * (1 - pL)) + beta * lambda__OPTIMALMP_1[] * (1 - pH) + sigma * lambda__OPTIMALMP_2[] * (1 - pH) = 0",
                 "beta * (-sigma * E[][lambda__OPTIMALMP_2[1]] - sigma * E[][lambda__OPTIMALMP_4[1]]) = 0")

# calibrating equations
calibr_equations__ <- character(0)

# variables / equations map
vareqmap__ <- sparseMatrix(i = c(1, 2, 2, 2, 2, 3, 3, 3, 3, 4,
                                 4, 4, 4, 5, 5, 5, 5, 6, 6, 6,
                                 6, 6, 7, 7, 7, 7, 7, 8, 8, 8,
                                 8, 8, 9, 9, 9, 9, 9, 10, 10, 10,
                                 10, 10, 11, 11),
                           j = c(1, 3, 4, 6, 9, 4, 5, 6, 10, 1,
                                 7, 8, 9, 1, 7, 8, 10, 2, 7, 8,
                                 9, 10, 2, 7, 8, 9, 10, 7, 8, 9,
                                 10, 11, 3, 4, 5, 6, 7, 3, 4, 5,
                                 6, 8, 4, 6),
                           x = c(3, 4, 6, 2, 2, 2, 4, 6, 2, 1,
                                 3, 2, 1, 1, 2, 3, 1, 1, 2, 2,
                                 3, 2, 1, 2, 2, 2, 3, 2, 2, 2,
                                 2, 6, 6, 2, 2, 2, 2, 2, 2, 6,
                                 2, 2, 4, 4),
                           dims = c(11, 11))

# variables / calibrating equations map
varcalibreqmap__ <- sparseMatrix(i = NULL, j = NULL, dims = c(0, 11))

# calibrated parameters / equations map
calibrpareqmap__ <- sparseMatrix(i = NULL, j = NULL, dims = c(11, 0))

# calibrated parameters / calibrating equations map
calibrparcalibreqmap__ <- sparseMatrix(i = NULL, j = NULL, dims = c(0, 0))

# free parameters / equations map
freepareqmap__ <- sparseMatrix(i = c(1, 2, 2, 2, 2, 2, 3, 3, 3, 3,
                                     3, 4, 4, 4, 5, 5, 5, 6, 6, 7,
                                     7, 8, 8, 8, 8, 8, 8, 9, 9, 9,
                                     9, 9, 9, 10, 10, 10, 10, 10, 10, 11,
                                     11),
                               j = c(3, 1, 2, 7, 8, 10, 1, 2, 7, 8,
                                     10, 1, 2, 7, 1, 2, 8, 7, 9, 8,
                                     9, 1, 2, 4, 5, 6, 10, 1, 4, 5,
                                     7, 8, 9, 1, 5, 6, 7, 8, 9, 1,
                                     9),
                               x = rep(1, 41), dims = c(11, 10))

# free parameters / calibrating equations map
freeparcalibreqmap__ <- sparseMatrix(i = NULL, j = NULL, dims = c(0, 10))

# shocks / equations map
shockeqmap__ <- sparseMatrix(i = c(1),
                             j = c(1),
                             x = rep(1, 1), dims = c(11, 1))

# steady state equations
ss_eq__ <- function(v, pc, pf)
{
    r <- numeric(11)
    r[1] = -v[1] + exp(pf[3] * log(v[1]))
    r[2] = pf[1] * (-v[4] + pf[2] * v[3]) + pf[7] * v[4] + v[6] * (1 - pf[8]) - 0.5 * pf[2] * pf[10]^-1 * v[9]
    r[3] = pf[1] * (-v[6] + pf[2] * v[5]) + pf[8] * v[6] + v[4] * (1 - pf[7]) - 0.5 * pf[2] * pf[10]^-1 * v[10]
    r[4] = -v[7] + log(v[1]) + pf[1] * v[7] + pf[2] * v[9] + pf[1] * (1 - pf[7]) * (-v[7] + v[8])
    r[5] = -v[8] + log(v[1]) + pf[2] * v[10] + pf[1] * pf[8] * v[8] + pf[1] * (1 - pf[8]) * (v[7] - v[8])
    r[6] = (1 - pf[7]) * (-v[9] + v[10]) - pf[9] * (v[2] - v[7]) + pf[9] * (1 - pf[7]) * (-v[7] + v[8])
    r[7] = (1 - pf[8]) * (v[9] - v[10]) - pf[9] * (v[2] - v[8]) + pf[9] * (1 - pf[8]) * (v[7] - v[8])
    r[8] = v[11] + 0.25 * (pf[4] - pf[5] + v[7])^2 + 0.25 * (-pf[5] + pf[6] + v[8])^2 - pf[1] * v[11] + 0.25 * pf[2] * pf[10]^-1 * v[9]^2 + 0.25 * pf[2] * pf[10]^-1 * v[10]^2
    r[9] = -0.5 * pf[4] + 0.5 * pf[5] - 0.5 * v[7] - pf[1] * v[3] + v[3] * (pf[1] - pf[1] * (1 - pf[7])) + v[4] * (pf[9] - pf[9] * (1 - pf[7])) + pf[1] * v[5] * (1 - pf[8]) + pf[9] * v[6] * (1 - pf[8])
    r[10] = 0.5 * pf[5] - 0.5 * pf[6] - 0.5 * v[8] - pf[1] * v[5] + v[5] * (pf[1] * pf[8] - pf[1] * (1 - pf[8])) + v[6] * (pf[9] - pf[9] * (1 - pf[8])) + pf[1] * v[3] * (1 - pf[7]) + pf[9] * v[4] * (1 - pf[7])
    r[11] = pf[1] * (-pf[9] * v[4] - pf[9] * v[6])

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
    jac <- numeric(44)
    jac[1] = -1 + pf[3] * v[1]^-1 * exp(pf[3] * log(v[1]))
    jac[2] = pf[1] * pf[2]
    jac[3] = -pf[1] + pf[7]
    jac[4] = 1 - pf[8]
    jac[5] = -0.5 * pf[2] * pf[10]^-1
    jac[6] = 1 - pf[7]
    jac[7] = pf[1] * pf[2]
    jac[8] = -pf[1] + pf[8]
    jac[9] = -0.5 * pf[2] * pf[10]^-1
    jac[10] = v[1]^-1
    jac[11] = -1 + pf[1] - pf[1] * (1 - pf[7])
    jac[12] = pf[1] * (1 - pf[7])
    jac[13] = pf[2]
    jac[14] = v[1]^-1
    jac[15] = pf[1] * (1 - pf[8])
    jac[16] = -1 + pf[1] * pf[8] - pf[1] * (1 - pf[8])
    jac[17] = pf[2]
    jac[18] = -pf[9]
    jac[19] = pf[9] - pf[9] * (1 - pf[7])
    jac[20] = pf[9] * (1 - pf[7])
    jac[21] = -1 + pf[7]
    jac[22] = 1 - pf[7]
    jac[23] = -pf[9]
    jac[24] = pf[9] * (1 - pf[8])
    jac[25] = pf[9] - pf[9] * (1 - pf[8])
    jac[26] = 1 - pf[8]
    jac[27] = -1 + pf[8]
    jac[28] = 0.5 * pf[4] - 0.5 * pf[5] + 0.5 * v[7]
    jac[29] = -0.5 * pf[5] + 0.5 * pf[6] + 0.5 * v[8]
    jac[30] = 0.5 * pf[2] * pf[10]^-1 * v[9]
    jac[31] = 0.5 * pf[2] * pf[10]^-1 * v[10]
    jac[32] = 1 - pf[1]
    jac[33] = -pf[1] * (1 - pf[7])
    jac[34] = pf[9] - pf[9] * (1 - pf[7])
    jac[35] = pf[1] * (1 - pf[8])
    jac[36] = pf[9] * (1 - pf[8])
    jac[37] = -0.5
    jac[38] = pf[1] * (1 - pf[7])
    jac[39] = pf[9] * (1 - pf[7])
    jac[40] = -pf[1] + pf[1] * pf[8] - pf[1] * (1 - pf[8])
    jac[41] = pf[9] - pf[9] * (1 - pf[8])
    jac[42] = -0.5
    jac[43] = -pf[1] * pf[9]
    jac[44] = -pf[1] * pf[9]
    jacob <- sparseMatrix(i = c(1, 2, 2, 2, 2, 3, 3, 3, 3, 4,
                                4, 4, 4, 5, 5, 5, 5, 6, 6, 6,
                                6, 6, 7, 7, 7, 7, 7, 8, 8, 8,
                                8, 8, 9, 9, 9, 9, 9, 10, 10, 10,
                                10, 10, 11, 11),
                          j = c(1, 3, 4, 6, 9, 4, 5, 6, 10, 1,
                                7, 8, 9, 1, 7, 8, 10, 2, 7, 8,
                                9, 10, 2, 7, 8, 9, 10, 7, 8, 9,
                                10, 11, 3, 4, 5, 6, 7, 3, 4, 5,
                                6, 8, 4, 6),
                          x = jac, dims = c(11, 11))

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
                         j = c(1, 1, 7, 9, 1, 8, 10, 2, 9, 2,
                               10),
                         x = Atm1x, dims = c(11, 11))

    Atx <- numeric(34)
    Atx[1] = -1
    Atx[2] = pf[7]
    Atx[3] = 1 - pf[8]
    Atx[4] = -0.5 * pf[2] * pf[10]^-1
    Atx[5] = 1 - pf[7]
    Atx[6] = pf[8]
    Atx[7] = -0.5 * pf[2] * pf[10]^-1
    Atx[8] = pf[1] - pf[1] * (1 - pf[7])
    Atx[9] = pf[1] * (1 - pf[7])
    Atx[10] = pf[1] * (1 - pf[8])
    Atx[11] = pf[1] * pf[8] - pf[1] * (1 - pf[8])
    Atx[12] = pf[9] - pf[9] * (1 - pf[7])
    Atx[13] = pf[9] * (1 - pf[7])
    Atx[14] = pf[7]
    Atx[15] = 1 - pf[7]
    Atx[16] = pf[9] * (1 - pf[8])
    Atx[17] = pf[9] - pf[9] * (1 - pf[8])
    Atx[18] = 1 - pf[8]
    Atx[19] = pf[8]
    Atx[20] = 0.5 * pf[4] - 0.5 * pf[5] + 0.5 * v[7]
    Atx[21] = -0.5 * pf[5] + 0.5 * pf[6] + 0.5 * v[8]
    Atx[22] = 0.5 * pf[2] * pf[10]^-1 * v[9]
    Atx[23] = 0.5 * pf[2] * pf[10]^-1 * v[10]
    Atx[24] = 1
    Atx[25] = pf[1] - pf[1] * (1 - pf[7])
    Atx[26] = pf[9] - pf[9] * (1 - pf[7])
    Atx[27] = pf[1] * (1 - pf[8])
    Atx[28] = pf[9] * (1 - pf[8])
    Atx[29] = -0.5
    Atx[30] = pf[1] * (1 - pf[7])
    Atx[31] = pf[9] * (1 - pf[7])
    Atx[32] = pf[1] * pf[8] - pf[1] * (1 - pf[8])
    Atx[33] = pf[9] - pf[9] * (1 - pf[8])
    Atx[34] = -0.5
    At <- sparseMatrix(i = c(1, 2, 2, 2, 3, 3, 3, 4, 4, 5,
                             5, 6, 6, 6, 6, 7, 7, 7, 7, 8,
                             8, 8, 8, 8, 9, 9, 9, 9, 9, 10,
                             10, 10, 10, 10),
                       j = c(1, 4, 6, 9, 4, 6, 10, 7, 8, 7,
                             8, 7, 8, 9, 10, 7, 8, 9, 10, 7,
                             8, 9, 10, 11, 3, 4, 5, 6, 7, 3,
                             4, 5, 6, 8),
                         x = Atx, dims = c(11, 11))

    Atp1x <- numeric(9)
    Atp1x[1] = pf[1] * pf[2]
    Atp1x[2] = -pf[1]
    Atp1x[3] = pf[1] * pf[2]
    Atp1x[4] = -pf[1]
    Atp1x[5] = -pf[1]
    Atp1x[6] = -pf[1]
    Atp1x[7] = -pf[1]
    Atp1x[8] = -pf[1] * pf[9]
    Atp1x[9] = -pf[1] * pf[9]
    Atp1 <- sparseMatrix(i = c(2, 2, 3, 3, 8, 9, 10, 11, 11),
                         j = c(3, 4, 5, 6, 11, 3, 5, 4, 6),
                         x = Atp1x, dims = c(11, 11))

    Aepsx <- numeric(1)
    Aepsx[1] = exp(pf[3] * log(v[1]))
    Aeps <- sparseMatrix(i = c(1),
                         j = c(1),
                         x = Aepsx, dims = c(11, 1))

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
