# Generated on 2023-11-10 14:31:13 by gEcon ver. 1.2.1 (2023-01-18)
# http://gecon.r-forge.r-project.org/

# Model name: RW_commit

# info
info__ <- c("RW_commit", "C:/Users/fm007/Documents/GitHub/OptimalRSW/RW_commit.gcn", "2023-11-10 14:31:13", "false")

# index sets
index_sets__ <- list()

# variables
variables__ <- c("etapi",
                 "lambda__RW_1",
                 "pi",
                 "y",
                 "U")

variables_tex__ <- c("{e\\!t\\!a\\!p\\!i}",
                     "\\lambda^{\\mathrm{RW}^{\\mathrm{1}}}",
                     "\\pi",
                     "y",
                     "U")

# shocks
shocks__ <- c("epsilon_pi")

shocks_tex__ <- c("\\epsilon^{\\pi}")

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
                           0,
                           0,
                           1,
                           6)

# equations
equations__ <- c("-etapi[] + exp(epsilon_pi[] + phi * log(etapi[-1])) = 0",
                 "beta * kappa * E[][lambda__RW_1[1]] - kappa * theta^-1 * y[] = 0",
                 "pistar - pi[] + beta * lambda__RW_1[] - beta * E[][lambda__RW_1[1]] = 0",
                 "-pi[-1] + log(etapi[-1]) + beta * pi[] + kappa * y[-1] = 0",
                 "U[] + 0.5 * (-pistar + pi[])^2 - beta * E[][U[1]] + 0.5 * kappa * theta^-1 * y[]^2 = 0")

# calibrating equations
calibr_equations__ <- character(0)

# variables / equations map
vareqmap__ <- sparseMatrix(i = c(1, 2, 2, 3, 3, 4, 4, 4, 5, 5,
                                 5),
                           j = c(1, 2, 4, 2, 3, 1, 3, 4, 3, 4,
                                 5),
                           x = c(3, 4, 2, 6, 2, 1, 3, 1, 2, 2,
                                 6),
                           dims = c(5, 5))

# variables / calibrating equations map
varcalibreqmap__ <- sparseMatrix(i = NULL, j = NULL, dims = c(0, 5))

# calibrated parameters / equations map
calibrpareqmap__ <- sparseMatrix(i = NULL, j = NULL, dims = c(5, 0))

# calibrated parameters / calibrating equations map
calibrparcalibreqmap__ <- sparseMatrix(i = NULL, j = NULL, dims = c(0, 0))

# free parameters / equations map
freepareqmap__ <- sparseMatrix(i = c(1, 2, 2, 2, 3, 3, 4, 4, 5, 5,
                                     5, 5),
                               j = c(3, 1, 2, 6, 1, 4, 1, 2, 1, 2,
                                     4, 6),
                               x = rep(1, 12), dims = c(5, 6))

# free parameters / calibrating equations map
freeparcalibreqmap__ <- sparseMatrix(i = NULL, j = NULL, dims = c(0, 6))

# shocks / equations map
shockeqmap__ <- sparseMatrix(i = c(1),
                             j = c(1),
                             x = rep(1, 1), dims = c(5, 1))

# steady state equations
ss_eq__ <- function(v, pc, pf)
{
    r <- numeric(5)
    r[1] = -v[1] + exp(pf[3] * log(v[1]))
    r[2] = pf[1] * pf[2] * v[2] - pf[2] * pf[6]^-1 * v[4]
    r[3] = pf[4] - v[3]
    r[4] = -v[3] + log(v[1]) + pf[1] * v[3] + pf[2] * v[4]
    r[5] = v[5] + 0.5 * (-pf[4] + v[3])^2 - pf[1] * v[5] + 0.5 * pf[2] * pf[6]^-1 * v[4]^2

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
    jac <- numeric(10)
    jac[1] = -1 + pf[3] * v[1]^-1 * exp(pf[3] * log(v[1]))
    jac[2] = pf[1] * pf[2]
    jac[3] = -pf[2] * pf[6]^-1
    jac[4] = -1
    jac[5] = v[1]^-1
    jac[6] = -1 + pf[1]
    jac[7] = pf[2]
    jac[8] = -pf[4] + v[3]
    jac[9] = pf[2] * pf[6]^-1 * v[4]
    jac[10] = 1 - pf[1]
    jacob <- sparseMatrix(i = c(1, 2, 2, 3, 4, 4, 4, 5, 5, 5),
                          j = c(1, 2, 4, 3, 1, 3, 4, 3, 4, 5),
                          x = jac, dims = c(5, 5))

    return(jacob)
}

# 1st order perturbation
pert1__ <- function(v, pc, pf)
{
    Atm1x <- numeric(4)
    Atm1x[1] = pf[3] * v[1]^-1 * exp(pf[3] * log(v[1]))
    Atm1x[2] = v[1]^-1
    Atm1x[3] = -1
    Atm1x[4] = pf[2]
    Atm1 <- sparseMatrix(i = c(1, 4, 4, 4),
                         j = c(1, 1, 3, 4),
                         x = Atm1x, dims = c(5, 5))

    Atx <- numeric(8)
    Atx[1] = -1
    Atx[2] = -pf[2] * pf[6]^-1
    Atx[3] = pf[1]
    Atx[4] = -1
    Atx[5] = pf[1]
    Atx[6] = -pf[4] + v[3]
    Atx[7] = pf[2] * pf[6]^-1 * v[4]
    Atx[8] = 1
    At <- sparseMatrix(i = c(1, 2, 3, 3, 4, 5, 5, 5),
                       j = c(1, 4, 2, 3, 3, 3, 4, 5),
                         x = Atx, dims = c(5, 5))

    Atp1x <- numeric(3)
    Atp1x[1] = pf[1] * pf[2]
    Atp1x[2] = -pf[1]
    Atp1x[3] = -pf[1]
    Atp1 <- sparseMatrix(i = c(2, 3, 5),
                         j = c(2, 2, 5),
                         x = Atp1x, dims = c(5, 5))

    Aepsx <- numeric(1)
    Aepsx[1] = exp(pf[3] * log(v[1]))
    Aeps <- sparseMatrix(i = c(1),
                         j = c(1),
                         x = Aepsx, dims = c(5, 1))

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
