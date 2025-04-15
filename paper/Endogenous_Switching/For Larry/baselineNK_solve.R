# ############################################################################
# (c) 3Eq NK model a la Rotember Woodford model with MS Regimes              #
#                                                                            #
# Authors:  Matyas Farkas, ECB and Kai Christoffel, ECB                      #
# ############################################################################
setwd("C:/Users/fm007/Documents/GitHub/OptimalRSW/paper/Endogenous_Switching/For Larry/")
# load gEcon package
library(gEcon)

# make and load the model
nk <- make_model("baseline_NK.gcn")
# find and print steady-state values
nk <- steady_state(nk)
get_residuals(nk)
list_eq(nk, eq_idx = c(15))
nk <-steady_state(model = nk, calibration = FALSE)
get_ss_values(nk, to_tex = TRUE)

# find and print perturbation solution
nk <- solve_pert(model = nk,loglin = TRUE)
check_bk(nk)
var_info(nk,all =T)

get_pert_solution(nk, to_tex = TRUE)

# set and print the shock distribution parameters
nk <- set_shock_cov_mat(nk, cov_matrix = diag(1, 5, 5),
                         shock_order = c("epsilon_Z","eta_pi","eta_p","eta_R","eta_G"))
shock_info(nk, all = TRUE)

# compute and print correlations
nk <- compute_model_stats(nk, ref_var = "pi", n_leadlags = 5)
get_model_stats(model = nk, 
                basic_stats = TRUE, 
                corr = TRUE, 
                autocorr = TRUE, 
                var_dec = FALSE,
                to_tex = TRUE)

# compute and print the IRFs
nk_irf <- compute_irf(nk, 
                       variables = c("pi", "C"))
plot_simulation(nk_irf, to_eps = TRUE)

# print summary of the model results
summary(nk)
