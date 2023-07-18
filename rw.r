# ############################################################################
# (c) 3Eq NK model a la Rotember Woodford model with MS Regimes              #
#                                                                            #
# Authors:  Matyas Farkas, ECB and Kai Christoffel, ECB                      #
# ############################################################################
setwd("C:/Users/fm007/Documents/GitHub/OptimalRSW")
# load gEcon package
library(gEcon)

# make and load the model
rw <- make_model("RW.gcn")
# find and print steady-state values
rw <- steady_state(rw)
get_ss_values(rw, to_tex = TRUE)

# find and print perturbation solution
rw <- solve_pert(model = rw,loglin = TRUE)
check_bk(rw)
var_info(rw,all =T)

get_pert_solution(rw, to_tex = TRUE)

# set and print the shock distribution parameters
rw <- set_shock_cov_mat(rw, cov_matrix = matrix(c(0.01), 1, 1),
                         shock_order = c("epsilon_Z"))
shock_info(rw, all = TRUE)

# compute and print correlations
rw <- compute_model_stats(rw, ref_var = "i", n_leadlags = 5)
get_model_stats(model = rw, 
                basic_stats = TRUE, 
                corr = TRUE, 
                autocorr = TRUE, 
                var_dec = FALSE,
                to_tex = TRUE)

# compute and print the IRFs
rw_irf <- compute_irf(rw, 
                       variables = c("i", "pi", "y"))
plot_simulation(rw_irf, to_eps = TRUE)

# print summary of the model results
summary(rw)
