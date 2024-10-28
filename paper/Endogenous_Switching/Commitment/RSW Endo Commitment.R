# ############################################################################
# (c) 3Eq NK model a la Rotember Woodford model with MS Regimes              #
#                                                                            #
# Authors:  Matyas Farkas, ECB and Kai Christoffel, ECB                      #
# ############################################################################
setwd("C:/Users/fm007/Documents/GitHub/OptimalRSW/paper/Endogenous_Switching/Commitment/")
# load gEcon package
library(gEcon)

# make and load the model
rsw <- make_model("RSW_RP_ONEOBJ_ENDO_COMMITMENT.gcn")
# find and print steady-state values
rsw <- steady_state(rsw)
get_ss_values(rsw, to_tex = TRUE)

# find and print perturbation solution
rsw <- solve_pert(model = rsw,loglin = TRUE)
check_bk(rsw)
var_info(rsw,all =T)

get_pert_solution(rsw, to_tex = TRUE)

# set and print the shock distribution parameters
rsw <- set_shock_cov_mat(rsw, cov_matrix = matrix(c(0.01), 1, 1),
                         shock_order = c("epsilon_pi"))
shock_info(rsw, all = TRUE)

# compute and print correlations
rsw <- compute_model_stats(rsw, ref_var = "piH", n_leadlags = 5)
get_model_stats(model = rsw, 
                basic_stats = TRUE, 
                corr = TRUE, 
                autocorr = TRUE, 
                var_dec = FALSE,
                to_tex = TRUE)

# compute and print the IRFs
rsw_irf <- compute_irf(rsw, 
                       variables = c("piH", "piL", "yH", "yL","pH","pL"))
plot_simulation(rsw_irf, to_eps = TRUE)

# print summary of the model results
summary(rsw)
