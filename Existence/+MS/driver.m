%
% Status : main Dynare file
%
% Warning : this file is generated automatically by Dynare
%           from model file (.mod)

if isoctave || matlab_ver_less_than('8.6')
    clear all
else
    clearvars -global
    clear_persistent_variables(fileparts(which('dynare')), false)
end
tic0 = tic;
% Define global variables.
global M_ options_ oo_ estim_params_ bayestopt_ dataset_ dataset_info estimation_info ys0_ ex0_
options_ = [];
M_.fname = 'MS';
M_.dynare_version = '5.1';
oo_.dynare_version = '5.1';
options_.dynare_version = '5.1';
%
% Some global variables initialization
%
global_initialization;
M_.exo_names = cell(2,1);
M_.exo_names_tex = cell(2,1);
M_.exo_names_long = cell(2,1);
M_.exo_names(1) = {'rnsh'};
M_.exo_names_tex(1) = {'rnsh'};
M_.exo_names_long(1) = {'rnsh'};
M_.exo_names(2) = {'mkpsh'};
M_.exo_names_tex(2) = {'mkpsh'};
M_.exo_names_long(2) = {'mkpsh'};
M_.endo_names = cell(7,1);
M_.endo_names_tex = cell(7,1);
M_.endo_names_long = cell(7,1);
M_.endo_names(1) = {'yH'};
M_.endo_names_tex(1) = {'yH'};
M_.endo_names_long(1) = {'yH'};
M_.endo_names(2) = {'yL'};
M_.endo_names_tex(2) = {'yL'};
M_.endo_names_long(2) = {'yL'};
M_.endo_names(3) = {'piH'};
M_.endo_names_tex(3) = {'piH'};
M_.endo_names_long(3) = {'piH'};
M_.endo_names(4) = {'piL'};
M_.endo_names_tex(4) = {'piL'};
M_.endo_names_long(4) = {'piL'};
M_.endo_names(5) = {'iH'};
M_.endo_names_tex(5) = {'iH'};
M_.endo_names_long(5) = {'iH'};
M_.endo_names(6) = {'iL'};
M_.endo_names_tex(6) = {'iL'};
M_.endo_names_long(6) = {'iL'};
M_.endo_names(7) = {'AUX_EXO_LAG_6_0'};
M_.endo_names_tex(7) = {'AUX\_EXO\_LAG\_6\_0'};
M_.endo_names_long(7) = {'AUX_EXO_LAG_6_0'};
M_.endo_partitions = struct();
M_.param_names = cell(9,1);
M_.param_names_tex = cell(9,1);
M_.param_names_long = cell(9,1);
M_.param_names(1) = {'kappa'};
M_.param_names_tex(1) = {'kappa'};
M_.param_names_long(1) = {'kappa'};
M_.param_names(2) = {'beta'};
M_.param_names_tex(2) = {'beta'};
M_.param_names_long(2) = {'beta'};
M_.param_names(3) = {'theta'};
M_.param_names_tex(3) = {'theta'};
M_.param_names_long(3) = {'theta'};
M_.param_names(4) = {'sigma'};
M_.param_names_tex(4) = {'sigma'};
M_.param_names_long(4) = {'sigma'};
M_.param_names(5) = {'pH'};
M_.param_names_tex(5) = {'pH'};
M_.param_names_long(5) = {'pH'};
M_.param_names(6) = {'pL'};
M_.param_names_tex(6) = {'pL'};
M_.param_names_long(6) = {'pL'};
M_.param_names(7) = {'pitH'};
M_.param_names_tex(7) = {'pitH'};
M_.param_names_long(7) = {'pitH'};
M_.param_names(8) = {'pitL'};
M_.param_names_tex(8) = {'pitL'};
M_.param_names_long(8) = {'pitL'};
M_.param_names(9) = {'pitCB'};
M_.param_names_tex(9) = {'pitCB'};
M_.param_names_long(9) = {'pitCB'};
M_.param_partitions = struct();
M_.exo_det_nbr = 0;
M_.exo_nbr = 2;
M_.endo_nbr = 7;
M_.param_nbr = 9;
M_.orig_endo_nbr = 6;
M_.aux_vars(1).endo_index = 7;
M_.aux_vars(1).type = 3;
M_.aux_vars(1).orig_index = 1;
M_.aux_vars(1).orig_lead_lag = 0;
M_.aux_vars(1).orig_expr = 'rnsh';
options_.varobs = cell(6, 1);
options_.varobs(1)  = {'piH'};
options_.varobs(2)  = {'piL'};
options_.varobs(3)  = {'yH'};
options_.varobs(4)  = {'yL'};
options_.varobs(5)  = {'iH'};
options_.varobs(6)  = {'iL'};
options_.varobs_id = [ 3 4 1 2 5 6  ];
M_ = setup_solvers(M_);
M_.Sigma_e = zeros(2, 2);
M_.Correlation_matrix = eye(2, 2);
M_.H = 0;
M_.Correlation_matrix_ME = 1;
M_.sigma_e_is_diagonal = true;
M_.det_shocks = [];
M_.surprise_shocks = [];
M_.heteroskedastic_shocks.Qvalue_orig = [];
M_.heteroskedastic_shocks.Qscale_orig = [];
options_.linear = false;
options_.block = false;
options_.bytecode = false;
options_.use_dll = false;
M_.nonzero_hessian_eqs = [];
M_.hessian_eq_zero = isempty(M_.nonzero_hessian_eqs);
M_.orig_eq_nbr = 6;
M_.eq_nbr = 7;
M_.ramsey_eq_nbr = 0;
M_.set_auxiliary_variables = exist(['./+' M_.fname '/set_auxiliary_variables.m'], 'file') == 2;
M_.epilogue_names = {};
M_.epilogue_var_list_ = {};
M_.orig_maximum_endo_lag = 0;
M_.orig_maximum_endo_lead = 0;
M_.orig_maximum_exo_lag = 1;
M_.orig_maximum_exo_lead = 0;
M_.orig_maximum_exo_det_lag = 0;
M_.orig_maximum_exo_det_lead = 0;
M_.orig_maximum_lag = 1;
M_.orig_maximum_lead = 0;
M_.orig_maximum_lag_with_diffs_expanded = 1;
M_.lead_lag_incidence = [
 0 2;
 0 3;
 0 4;
 0 5;
 0 6;
 0 7;
 1 8;]';
M_.nstatic = 6;
M_.nfwrd   = 0;
M_.npred   = 1;
M_.nboth   = 0;
M_.nsfwrd   = 0;
M_.nspred   = 1;
M_.ndynamic   = 1;
M_.dynamic_tmp_nbr = [0; 0; 0; 0; ];
M_.model_local_variables_dynamic_tt_idxs = {
};
M_.equations_tags = {
  1 , 'name' , 'High Regime - IS' ;
  2 , 'name' , 'High Regime - PC' ;
  3 , 'name' , 'High Regime - Optimal MP' ;
  4 , 'name' , 'Low Regime - IS' ;
  5 , 'name' , 'Low Regime - PC' ;
  6 , 'name' , 'Low Regime - Optimal MP' ;
};
M_.mapping.yH.eqidx = [1 2 3 4 ];
M_.mapping.yL.eqidx = [1 4 5 6 ];
M_.mapping.piH.eqidx = [1 2 3 4 5 ];
M_.mapping.piL.eqidx = [1 2 4 5 6 ];
M_.mapping.iH.eqidx = [1 ];
M_.mapping.iL.eqidx = [4 ];
M_.mapping.rnsh.eqidx = [1 4 ];
M_.mapping.mkpsh.eqidx = [2 5 ];
M_.static_and_dynamic_models_differ = false;
M_.has_external_function = false;
M_.state_var = [7 ];
M_.exo_names_orig_ord = [1:2];
M_.maximum_lag = 1;
M_.maximum_lead = 0;
M_.maximum_endo_lag = 1;
M_.maximum_endo_lead = 0;
oo_.steady_state = zeros(7, 1);
M_.maximum_exo_lag = 0;
M_.maximum_exo_lead = 0;
oo_.exo_steady_state = zeros(2, 1);
M_.params = NaN(9, 1);
M_.endo_trends = struct('deflator', cell(7, 1), 'log_deflator', cell(7, 1), 'growth_factor', cell(7, 1), 'log_growth_factor', cell(7, 1));
M_.NNZDerivatives = [28; 0; -1; ];
M_.static_tmp_nbr = [0; 0; 0; 0; ];
M_.model_local_variables_static_tt_idxs = {
};
M_.params(1) = 0.2465;
kappa = M_.params(1);
M_.params(2) = 0.99;
beta = M_.params(2);
M_.params(3) = 6;
theta = M_.params(3);
M_.params(4) = 1;
sigma = M_.params(4);
M_.params(5) = 0.99;
pH = M_.params(5);
M_.params(6) = 0.99;
pL = M_.params(6);
M_.params(8) = (-2);
pitL = M_.params(8);
M_.params(7) = 2;
pitH = M_.params(7);
M_.params(9) = 0;
pitCB = M_.params(9);
%
% SHOCKS instructions
%
M_.exo_det_length = 0;
M_.Sigma_e(1, 1) = 1;
M_.Sigma_e(2, 2) = 1;
steady;
resid;                  
oo_.dr.eigval = check(M_,options_,oo_);
model_diagnostics(M_,options_,oo_);
options_.drop = 50;
options_.irf = 0;
options_.order = 1;
options_.periods = 150;
var_list_ = {};
[info, oo_, options_, M_] = stoch_simul(M_, options_, oo_, var_list_);
options_.rplottype=2;
options_.TeX = 1;
var_list_ = {'piH';'piL';'yH';'yL';'iH';'iL'};
rplot(var_list_);
estim_params_.var_exo = zeros(0, 10);
estim_params_.var_endo = zeros(0, 10);
estim_params_.corrx = zeros(0, 11);
estim_params_.corrn = zeros(0, 11);
estim_params_.param_vals = zeros(0, 10);
estim_params_.param_vals = [estim_params_.param_vals; 1, NaN, (-Inf), Inf, 5, 0.25, 0.1443375672974065, NaN, NaN, NaN ];
estim_params_.param_vals = [estim_params_.param_vals; 4, NaN, (-Inf), Inf, 5, 5.5, 2.598076211353316, NaN, NaN, NaN ];
options_.nograph=0; 
options_gsa = struct();
options_gsa.Nsam = 5000;
options_gsa.prior_range = 0;
options_gsa.stab = 1;
dynare_sensitivity(options_gsa);


oo_.time = toc(tic0);
disp(['Total computing time : ' dynsec2hms(oo_.time) ]);
if ~exist([M_.dname filesep 'Output'],'dir')
    mkdir(M_.dname,'Output');
end
save([M_.dname filesep 'Output' filesep 'MS_results.mat'], 'oo_', 'M_', 'options_');
if exist('estim_params_', 'var') == 1
  save([M_.dname filesep 'Output' filesep 'MS_results.mat'], 'estim_params_', '-append');
end
if exist('bayestopt_', 'var') == 1
  save([M_.dname filesep 'Output' filesep 'MS_results.mat'], 'bayestopt_', '-append');
end
if exist('dataset_', 'var') == 1
  save([M_.dname filesep 'Output' filesep 'MS_results.mat'], 'dataset_', '-append');
end
if exist('estimation_info', 'var') == 1
  save([M_.dname filesep 'Output' filesep 'MS_results.mat'], 'estimation_info', '-append');
end
if exist('dataset_info', 'var') == 1
  save([M_.dname filesep 'Output' filesep 'MS_results.mat'], 'dataset_info', '-append');
end
if exist('oo_recursive_', 'var') == 1
  save([M_.dname filesep 'Output' filesep 'MS_results.mat'], 'oo_recursive_', '-append');
end
if ~isempty(lastwarn)
  disp('Note: warning(s) encountered in MATLAB/Octave code')
end
