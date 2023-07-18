% =========================================================================
% Declare endogenous variables
% =========================================================================

var 
yH 
yL 
piH 
piL 
iH 
iL 
;

varexo rnsh mkpsh;

% =========================================================================
% Declare parameters
% =========================================================================

parameters
kappa 
beta 
theta
sigma 
pH
pL
pitH
pitL
pitCB;

% =========================================================================
% Calibrate parameter values
% =========================================================================

kappa =  0.2465; %(1-alpha)(1-alpha*beta)/(alpha*(1+eta*theta)) * ((sigma^-1) + eta), where alpha =0.33, beta =0.99, sigma =2, theta=6
beta = 0.99;
theta = 6;
sigma = 1;
pH = 0.99;
pL = 0.99;
pitL = -2;
pitH = 2;
pitCB = 0;
% =========================================================================
% Model equations
% =========================================================================
model;
#rn = 1/beta -1;

[name='High Regime - IS']
yH = pH * yH + (1-pH)* yL  + sigma * ( pH*(piH) + (1-pH)*(piL) - iH + rn + rnsh + 0.9*rnsh(-1));

[name='High Regime - PC']
(piH) = kappa * yH + beta * ( pH*(piH) + (1-pH)*(piL) ) + mkpsh;

[name='High Regime - Optimal MP']
0 = kappa * (piH-pitCB) + kappa/theta * yH;

[name='Low Regime - IS']
yL = pL * yL + (1-pL)* yH  + sigma * ( pL*(piL) + (1-pL)*(piH) - iL + rn + rnsh + 0.9*rnsh(-1)) ;

[name='Low Regime - PC']
(piL) = kappa * yL + beta * (  pL*(piL) + (1-pL)*(piH) ) + mkpsh;

[name='Low Regime - Optimal MP']
0 = kappa * (piL-pitCB) + kappa/theta * yL;

end;

% =========================================================================
% Steady state Model
% =========================================================================
steady_state_model;
piL = pitCB;
piH = pitCB;
yH = 0;
yL = 0;
iH = 1/beta -1;
iL = 1/beta -1;
end;

shocks;
var rnsh = 1;
var mkpsh = 1;
end;


steady;                 % compute steady state given the starting values
resid;                  % check the residuals of model equations evaluated at steady state
check;                  % check Blanchard-Kahn conditions
model_diagnostics;      % check obvious model errors

stoch_simul(order=1,periods=150,drop=50,irf=0);

options_.rplottype=2;
options_.TeX = 1;
rplot piH piL yH yL iH iL;



%varobs piH piL yH yL iH iL;
%identification;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%% The Determinacy Region 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

                
// A similar figure could have been obtained using Dynare's dynare_sensitivity command. The difference
// is that the parameter would have been randomly sampled from the prior instead of uniformly on the grid
// Moreover, many draws (nsam) are required to get a mapping figure of the type obtained above.
//         
%specify parameters for which to map sensitivity
estimated_params;
kappa,uniform_pdf,(0+0.5)/2,sqrt(12)^(-1)*(0.5-0);
sigma,uniform_pdf,((1)+10)/2,sqrt(12)^(-1)*(10-(1));      
end;
 
varobs piH piL yH yL iH iL;
options_.nograph=0; 
dynare_sensitivity(prior_range=0,stab=1,nsam=5000);