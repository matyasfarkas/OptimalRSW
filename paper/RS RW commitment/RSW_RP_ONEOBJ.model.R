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
                               j = c(4, 1, 2, 3, 1, 2, 3, 1, 2, 8,
                                     1, 2, 9, 1, 5, 6, 8, 9, 1, 6,
                                     7, 8, 9, 1, 3, 5, 6, 7),
                               x = rep(1, 28), dims = c(8, 11))

# free parameters / calibrating equations map
freeparcalibreqmap__ <- sparseMatrix(i = NULL, j = NULL, dims = c(0, 11))

# shocks / equations map
shockeqmap__ <- sparseMatrix(i = c(1),
                             j = c(1),
                             x = rep(1, 1), dims = c(8, 1))

# steady state equations
ss_eq__ <- function(v, pc, pf)
{
    r <- numeric(8)
    r[1] = -v[1] + exp(pf[4] * log(v[1]))
    r[2] = -0.5 * pf[3] * v[6] + pf[1] * pf[2] * v[2]
    r[3] = -0.5 * pf[3] * v[7] + pf[1] * pf[2] * v[3]
    r[4] = -v[4] + log(v[1]) + pf[1] * v[4] + pf[2] * v[6] + pf[1] * (1 - pf[8]) * (-v[4] + v[5])
    r[5] = -v[5] + log(v[1]) + pf[2] * v[7] + pf[1] * pf[9] * v[5] + pf[1] * (1 - pf[9]) * (v[4] - v[5])
    r[6] = -0.5 * pf[5] + 0.5 * pf[6] - 0.5 * v[4] - pf[1] * v[2] + v[2] * (pf[1] - pf[1] * (1 - pf[8])) + pf[1] * v[3] * (1 - pf[9])
    r[7] = 0.5 * pf[6] - 0.5 * pf[7] - 0.5 * v[5] - pf[1] * v[3] + v[3] * (pf[1] * pf[9] - pf[1] * (1 - pf[9])) + pf[1] * v[2] * (1 - pf[8])
    r[8] = v[8] + 0.25 * (pf[5] - pf[6] + v[4])^2 + 0.25 * (-pf[6] + pf[7] + v[5])^2 - pf[1] * v[8] + 0.25 * pf[3] * v[6]^2 + 0.25 * pf[3] * v[7]^2

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
    jac <- numeric(24)
    jac[1] = -1 + pf[4] * v[1]^-1 * exp(pf[4] * log(v[1]))
    jac[2] = pf[1] * pf[2]
    jac[3] = -0.5 * pf[3]
    jac[4] = pf[1] * pf[2]
    jac[5] = -0.5 * pf[3]
    jac[6] = v[1]^-1
    jac[7] = -1 + pf[1] - pf[1] * (1 - pf[8])
    jac[8] = pf[1] * (1 - pf[8])
    jac[9] = pf[2]
    jac[10] = v[1]^-1
    jac[11] = pf[1] * (1 - pf[9])
    jac[12] = -1 + pf[1] * pf[9] - pf[1] * (1 - pf[9])
    jac[13] = pf[2]
    jac[14] = -pf[1] * (1 - pf[8])
    jac[15] = pf[1] * (1 - pf[9])
    jac[16] = -0.5
    jac[17] = pf[1] * (1 - pf[8])
    jac[18] = -pf[1] + pf[1] * pf[9] - pf[1] * (1 - pf[9])
    jac[19] = -0.5
    jac[20] = 0.5 * pf[5] - 0.5 * pf[6] + 0.5 * v[4]
    jac[21] = -0.5 * pf[6] + 0.5 * pf[7] + 0.5 * v[5]
    jac[22] = 0.5 * pf[3] * v[6]
    jac[23] = 0.5 * pf[3] * v[7]
    jac[24] = 1 - pf[1]
    jacob <- sparseMatrix(i = c(1, 2, 2, 3, 3, 4, 4, 4, 4, 5,
                                5, 5, 5, 6, 6, 6, 7, 7, 7, 8,
                                8, 8, 8, 8),
                          j = c(1, 2, 6, 3, 7, 1, 4, 5, 6, 1,
                                4, 5, 7, 2, 3, 4, 2, 3, 5, 4,
                                5, 6, 7, 8),
                          x = jac, dims = c(8, 8))

    return(jacob)
}

# 1st order perturbation
pert1__ <- function(v, pc, pf)
{
    Atm1x <- numeric(7)
    Atm1x[1] = pf[4] * v[1]^-1 * exp(pf[4] * log(v[1]))
    Atm1x[2] = v[1]^-1
    Atm1x[3] = -1
    Atm1x[4] = pf[2]
    Atm1x[5] = v[1]^-1
    Atm1x[6] = -1
    Atm1x[7] = pf[2]
    Atm1 <- sparseMatrix(i = c(1, 4, 4, 4, 5, 5, 5),
                         j = c(1, 1, 4, 6, 1, 5, 7),
                         x = Atm1x, dims = c(8, 8))

    Atx <- numeric(18)
    Atx[1] = -1
    Atx[2] = -0.5 * pf[3]
    Atx[3] = -0.5 * pf[3]
    Atx[4] = pf[1] - pf[1] * (1 - pf[8])
    Atx[5] = pf[1] * (1 - pf[8])
    Atx[6] = pf[1] * (1 - pf[9])
    Atx[7] = pf[1] * pf[9] - pf[1] * (1 - pf[9])
    Atx[8] = pf[1] - pf[1] * (1 - pf[8])
    Atx[9] = pf[1] * (1 - pf[9])
    Atx[10] = -0.5
    Atx[11] = pf[1] * (1 - pf[8])
    Atx[12] = pf[1] * pf[9] - pf[1] * (1 - pf[9])
    Atx[13] = -0.5
    Atx[14] = 0.5 * pf[5] - 0.5 * pf[6] + 0.5 * v[4]
    Atx[15] = -0.5 * pf[6] + 0.5 * pf[7] + 0.5 * v[5]
    Atx[16] = 0.5 * pf[3] * v[6]
    Atx[17] = 0.5 * pf[3] * v[7]
    Atx[18] = 1
    At <- sparseMatrix(i = c(1, 2, 3, 4, 4, 5, 5, 6, 6, 6,
                             7, 7, 7, 8, 8, 8, 8, 8),
                       j = c(1, 6, 7, 4, 5, 4, 5, 2, 3, 4,
                             2, 3, 5, 4, 5, 6, 7, 8),
                         x = Atx, dims = c(8, 8))

    Atp1x <- numeric(5)
    Atp1x[1] = pf[1] * pf[2]
    Atp1x[2] = pf[1] * pf[2]
    Atp1x[3] = -pf[1]
    Atp1x[4] = -pf[1]
    Atp1x[5] = -pf[1]
    Atp1 <- sparseMatrix(i = c(2, 3, 6, 7, 8),
                         j = c(2, 3, 2, 3, 8),
                         x = Atp1x, dims = c(8, 8))

    Aepsx <- numeric(1)
    Aepsx[1] = exp(pf[4] * log(v[1]))
    Aeps <- sparseMatrix(i = c(1),
                         j = c(1),
                         x = Aepsx, dims = c(8, 1))

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
