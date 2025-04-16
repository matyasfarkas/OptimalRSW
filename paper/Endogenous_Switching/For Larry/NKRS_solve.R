# ############################################################################
# (c) 3Eq NK model with endogenous probability                               #
#                                                                            #
# Authors:  Matyas Farkas, IMF 2025 April                                    #
# ############################################################################
setwd("C:/Users/fm007/Documents/GitHub/OptimalRSW/paper/Endogenous_Switching/For Larry/")
# load gEcon package
library(gEcon)

# make and load the model
nkrs <- make_model("NK_RS.gcn")
# find and print steady-state values

steady_state_values <- c(
  epsilon_G = 1,
  g_1 = 7.3514,
  g_2 = 4.9009,
  inflation_gap = 1,
  lambda = 1.5467,
  mc = 0.6667,
  nu_p = 1,
  perceived_pi_obj = 1,
  pi = 1,
  pi_star = 1,
  pi_obj = 1,
  pH = 0.95,
  pL = 0.05,
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
nkrs <- initval_var(nkrs, steady_state_values)

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
                var_dec = TRUE,
                to_tex = TRUE)

# compute and print the IRFs
nkrs_irf <- compute_irf(nkrs, 
                       variables = c("pi","pH","R","percieved_pi_obj","inflation_gap"))
plot_simulation(nkrs_irf, to_eps = TRUE)

# print summary of the model results
summary(nkrs)

