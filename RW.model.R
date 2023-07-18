# Generated on 2023-07-18 11:21:01 by gEcon ver. 1.2.1 (2023-01-18)
# http://gecon.r-forge.r-project.org/

# Model name: RW

# info
info__ <- c("RW", "C:/Users/fm007/Documents/GitHub/OptimalRSW/RW.gcn", "2023-07-18 11:21:01", "false")

# index sets
index_sets__ <- list()

# variables
variables__ <- c("i",
                 "lambda__RW_1",
                 "lambda__RW_2",
                 "pi",
                 "rn",
                 "y",
                 "U")

variables_tex__ <- c("i",
                     "\\lambda^{\\mathrm{RW}^{\\mathrm{1}}}",
                     "\\lambda^{\\mathrm{RW}^{\\mathrm{2}}}",
                     "\\pi",
                     "{r\\!n}",
                     "y",
                     "U")

# shocks
shocks__ <- c("epsilon_Z")

shocks_tex__ <- c("\\epsilon^{\\mathrm{Z}}")

# parameters
parameters__ <- c("beta",
                  "kappa",
                  "phi",
                  "pistar",
                  "sigma",
                  "theta")

parameters_tex__ <- c("\\beta",
                     "\\kappa",
                     "\\phi",
                     "{p\\!i\\!s\\!t\\!a\\!r}",
                     "\\sigma",
                     "\\theta")

# free parameters
parameters_free__ <- c("beta",
                       "kappa",
                       "phi",
                       "pistar",
                       "sigma",
                       "theta")

# free parameters' values
parameters_free_val__ <- c(0.99,
                           0.2465,
                           0.95,
                           0,
                           1,
                           6)

# equations
equations__ <- c("-rn[] + exp(epsilon_Z[] + phi * log(rn[-1])) = 0",
                 "-pi[-1] + beta * pi[] + kappa * y[-1] = 0",
                 "-y[-1] + y[] - sigma * (i[-1] - rn[-1] - pi[]) = 0",
                 "lambda__RW_2[] + beta * (kappa * E[][lambda__RW_1[1]] - E[][lambda__RW_2[1]]) - kappa * theta^-1 * y[] = 0",
                 "U[] + 0.5 * (-pistar + pi[])^2 - beta * E[][U[1]] + 0.5 * kappa * theta^-1 * y[]^2 = 0",
                 "pistar - pi[] + beta * lambda__RW_1[] - beta * E[][lambda__RW_1[1]] + sigma * lambda__RW_2[] = 0",
                 "-beta * sigma * E[][lambda__RW_2[1]] = 0")

# calibrating equations
calibr_equations__ <- character(0)

# variables / equations map
vareqmap__ <- sparseMatrix(i = c(1, 2, 2, 3, 3, 3, 3, 4, 4, 4,
                                 5, 5, 5, 6, 6, 6, 7),
                           j = c(5, 4, 6, 1, 4, 5, 6, 2, 3, 6,
                                 4, 6, 7, 2, 3, 4, 3),
                           x = c(3, 3, 1, 1, 2, 1, 3, 4, 6, 2,
                                 2, 2, 6, 6, 2, 2, 4),
                           dims = c(7, 7))

# variables / calibrating equations map
varcalibreqmap__ <- sparseMatrix(i = NULL, j = NULL, dims = c(0, 7))

# calibrated parameters / equations map
calibrpareqmap__ <- sparseMatrix(i = NULL, j = NULL, dims = c(7, 0))

# calibrated parameters / calibrating equations map
calibrparcalibreqmap__ <- sparseMatrix(i = NULL, j = NULL, dims = c(0, 0))

# free parameters / equations map
freepareqmap__ <- sparseMatrix(i = c(1, 2, 2, 3, 4, 4, 4, 5, 5, 5,
                                     5, 6, 6, 6, 7, 7),
                               j = c(3, 1, 2, 5, 1, 2, 6, 1, 2, 4,
                                     6, 1, 4, 5, 1, 5),
                               x = rep(1, 16), dims = c(7, 6))

# free parameters / calibrating equations map
freeparcalibreqmap__ <- sparseMatrix(i = NULL, j = NULL, dims = c(0, 6))

# shocks / equations map
shockeqmap__ <- sparseMatrix(i = c(1),
                             j = c(1),
                             x = rep(1, 1), dims = c(7, 1))

# steady state equations
ss_eq__ <- function(v, pc, pf)
{
    r <- numeric(7)
    r[1] = -v[5] + exp(pf[3] * log(v[5]))
    r[2] = -v[4] + pf[1] * v[4] + pf[2] * v[6]
    r[3] = -pf[5] * (v[1] - v[4] - v[5])
    r[4] = v[3] + pf[1] * (-v[3] + pf[2] * v[2]) - pf[2] * pf[6]^-1 * v[6]
    r[5] = v[7] + 0.5 * (-pf[4] + v[4])^2 - pf[1] * v[7] + 0.5 * pf[2] * pf[6]^-1 * v[6]^2
    r[6] = pf[4] - v[4] + pf[5] * v[3]
    r[7] = -pf[1] * pf[5] * v[3]

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
    jac <- numeric(15)
    jac[1] = -1 + pf[3] * v[5]^-1 * exp(pf[3] * log(v[5]))
    jac[2] = -1 + pf[1]
    jac[3] = pf[2]
    jac[4] = -pf[5]
    jac[5] = pf[5]
    jac[6] = pf[5]
    jac[7] = pf[1] * pf[2]
    jac[8] = 1 - pf[1]
    jac[9] = -pf[2] * pf[6]^-1
    jac[10] = -pf[4] + v[4]
    jac[11] = pf[2] * pf[6]^-1 * v[6]
    jac[12] = 1 - pf[1]
    jac[13] = pf[5]
    jac[14] = -1
    jac[15] = -pf[1] * pf[5]
    jacob <- sparseMatrix(i = c(1, 2, 2, 3, 3, 3, 4, 4, 4, 5,
                                5, 5, 6, 6, 7),
                          j = c(5, 4, 6, 1, 4, 5, 2, 3, 6, 4,
                                6, 7, 3, 4, 3),
                          x = jac, dims = c(7, 7))

    return(jacob)
}

# 1st order perturbation
pert1__ <- function(v, pc, pf)
{
    Atm1x <- numeric(6)
    Atm1x[1] = pf[3] * v[5]^-1 * exp(pf[3] * log(v[5]))
    Atm1x[2] = -1
    Atm1x[3] = pf[2]
    Atm1x[4] = -pf[5]
    Atm1x[5] = pf[5]
    Atm1x[6] = -1
    Atm1 <- sparseMatrix(i = c(1, 2, 2, 3, 3, 3),
                         j = c(5, 4, 6, 1, 5, 6),
                         x = Atm1x, dims = c(7, 7))

    Atx <- numeric(12)
    Atx[1] = -1
    Atx[2] = pf[1]
    Atx[3] = pf[5]
    Atx[4] = 1
    Atx[5] = 1
    Atx[6] = -pf[2] * pf[6]^-1
    Atx[7] = -pf[4] + v[4]
    Atx[8] = pf[2] * pf[6]^-1 * v[6]
    Atx[9] = 1
    Atx[10] = pf[1]
    Atx[11] = pf[5]
    Atx[12] = -1
    At <- sparseMatrix(i = c(1, 2, 3, 3, 4, 4, 5, 5, 5, 6,
                             6, 6),
                       j = c(5, 4, 4, 6, 3, 6, 4, 6, 7, 2,
                             3, 4),
                         x = Atx, dims = c(7, 7))

    Atp1x <- numeric(5)
    Atp1x[1] = pf[1] * pf[2]
    Atp1x[2] = -pf[1]
    Atp1x[3] = -pf[1]
    Atp1x[4] = -pf[1]
    Atp1x[5] = -pf[1] * pf[5]
    Atp1 <- sparseMatrix(i = c(4, 4, 5, 6, 7),
                         j = c(2, 3, 7, 2, 3),
                         x = Atp1x, dims = c(7, 7))

    Aepsx <- numeric(1)
    Aepsx[1] = exp(pf[3] * log(v[5]))
    Aeps <- sparseMatrix(i = c(1),
                         j = c(1),
                         x = Aepsx, dims = c(7, 1))

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
