# ############################################################################
# (c) 3Eq NK model a la Rotember Woodford model with MS Regimes              #
#                                                                            #
# Authors:  Matyas Farkas, ECB and Kai Christoffel, ECB                      #
# ############################################################################
setwd("C:/Users/fm007/Documents/GitHub/OptimalRSW/paper/Endogenous_Switching/For Larry/")
# load gEcon package
library(gEcon)

# make and load the model
nkrs <- make_model("NK_RS.gcn")
# find and print steady-state values
nkrs <- steady_state(nkrs)
get_residuals(nkrs)
get_ss_values(nkrs, to_tex = TRUE)

# find and print perturbation solution
nkrs <- solve_pert(model = nkrs,loglin = TRUE)
check_bk(nkrs)
var_info(nkrs,all =T)

get_pert_solution(nkrs, to_tex = TRUE)

# set and print the shock distribution parameters
nkrs <- set_shock_cov_mat(nkrs, cov_matrix = diag(1, 5, 5),
                         shock_order = c("epsilon_Z","eta_pi","eta_p","eta_R","eta_G"))
shock_info(nkrs, all = TRUE)

# compute and print correlations
nkrs <- compute_model_stats(nkrs, ref_var = "pi", n_leadlags = 5)
get_model_stats(model = nkrs, 
                basic_stats = TRUE, 
                corr = TRUE, 
                autocorr = TRUE, 
                var_dec = FALSE,
                to_tex = TRUE)

# compute and print the IRFs
nkrs_irf <- compute_irf(nkrs, 
                       variables = c("pi","R", "C","pH"))
plot_simulation(nkrs_irf, to_eps = TRUE)

# print summary of the model results
summary(nkrs)
