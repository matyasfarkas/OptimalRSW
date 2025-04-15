# ############################################################################
# (c) 3Eq NK model                                                           #
#                                                                            #
# Authors:  Matyas Farkas, IMF April 2025                                    #
# ############################################################################
setwd("C:/Users/fm007/Documents/GitHub/OptimalRSW/paper/Endogenous_Switching/For Larry/")
# load gEcon package
library(gEcon)

# make and load the model
nk <- make_model("baseline_NK.gcn")
# find and print steady-state values
steady_state_values <- c(
  epsilon_G = 1,
  g_1 = 7.3514,
  g_2 = 4.9009,
  lambda = 1.5467,
  mc = 0.6667,
  pi_obj = 1,
  nu_p = 1,
  pi = 1,
  pi_star = 1,
  q = 1.5467,
  r = 0.0351,
  B = 0,
  C = 0.3255,
  Div = 0.1601,
  G = 0.0865,
  I = 0.0684,
  K_s = 2.7374,
  L_s = 0.2279,
  Q = 1,
  R = 1.0101,
  T = 0.0865,
  U = -167.8256,
  W = 0.9837,
  Y = 0.4804,
  Y_j = 0.4804,
  Y_s = 0.4804,
  Z = 1
)
nk <- initval_var(nk, steady_state_values)

nk <- steady_state(nk)
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
nk_irf <- compute_irf(nk,variables = c("pi","C","R"))
plot_simulation(nk_irf, to_eps = TRUE)

# print summary of the model results
summary(nk)
