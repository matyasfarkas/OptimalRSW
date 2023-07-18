function [rp, gp, rpp, gpp, hp, g3p] = dynamic_params_derivs(y, x, params, steady_state, it_, ss_param_deriv, ss_param_2nd_deriv)
%
% Compute the derivatives of the dynamic model with respect to the parameters
% Inputs :
%   y         [#dynamic variables by 1] double    vector of endogenous variables in the order stored
%                                                 in M_.lead_lag_incidence; see the Manual
%   x         [nperiods by M_.exo_nbr] double     matrix of exogenous variables (in declaration order)
%                                                 for all simulation periods
%   params    [M_.param_nbr by 1] double          vector of parameter values in declaration order
%   steady_state  [M_.endo_nbr by 1] double       vector of steady state values
%   it_       scalar double                       time period for exogenous variables for which to evaluate the model
%   ss_param_deriv     [M_.eq_nbr by #params]     Jacobian matrix of the steady states values with respect to the parameters
%   ss_param_2nd_deriv [M_.eq_nbr by #params by #params] Hessian matrix of the steady states values with respect to the parameters
%
% Outputs:
%   rp        [M_.eq_nbr by #params] double    Jacobian matrix of dynamic model equations with respect to parameters 
%                                              Dynare may prepend or append auxiliary equations, see M_.aux_vars
%   gp        [M_.endo_nbr by #dynamic variables by #params] double    Derivative of the Jacobian matrix of the dynamic model equations with respect to the parameters
%                                                           rows: equations in order of declaration
%                                                           columns: variables in order stored in M_.lead_lag_incidence
%   rpp       [#second_order_residual_terms by 4] double   Hessian matrix of second derivatives of residuals with respect to parameters;
%                                                              rows: respective derivative term
%                                                              1st column: equation number of the term appearing
%                                                              2nd column: number of the first parameter in derivative
%                                                              3rd column: number of the second parameter in derivative
%                                                              4th column: value of the Hessian term
%   gpp      [#second_order_Jacobian_terms by 5] double   Hessian matrix of second derivatives of the Jacobian with respect to the parameters;
%                                                              rows: respective derivative term
%                                                              1st column: equation number of the term appearing
%                                                              2nd column: column number of variable in Jacobian of the dynamic model
%                                                              3rd column: number of the first parameter in derivative
%                                                              4th column: number of the second parameter in derivative
%                                                              5th column: value of the Hessian term
%   hp      [#first_order_Hessian_terms by 5] double   Jacobian matrix of derivatives of the dynamic Hessian with respect to the parameters;
%                                                              rows: respective derivative term
%                                                              1st column: equation number of the term appearing
%                                                              2nd column: column number of first variable in Hessian of the dynamic model
%                                                              3rd column: column number of second variable in Hessian of the dynamic model
%                                                              4th column: number of the parameter in derivative
%                                                              5th column: value of the Hessian term
%   g3p      [#first_order_g3_terms by 6] double   Jacobian matrix of derivatives of g3 (dynamic 3rd derivs) with respect to the parameters;
%                                                              rows: respective derivative term
%                                                              1st column: equation number of the term appearing
%                                                              2nd column: column number of first variable in g3 of the dynamic model
%                                                              3rd column: column number of second variable in g3 of the dynamic model
%                                                              4th column: column number of third variable in g3 of the dynamic model
%                                                              5th column: number of the parameter in derivative
%                                                              6th column: value of the Hessian term
%
%
% Warning : this file is generated automatically by Dynare
%           from model file (.mod)

T = NaN(1,1);
T(1) = (-((-params(1))*(params(3)+params(3))))/(params(3)*params(3)*params(3)*params(3));
rp = zeros(7, 9);
rp(1, 2) = (-(params(4)*(-1)/(params(2)*params(2))));
rp(1, 4) = (-(x(it_, 1)+1/params(2)-1+params(5)*y(4)+(1-params(5))*y(5)-y(6)+0.9*y(1)));
rp(1, 5) = (-(y(2)-y(3)+params(4)*(y(4)-y(5))));
rp(2, 1) = (-y(2));
rp(2, 2) = (-(params(5)*y(4)+(1-params(5))*y(5)));
rp(2, 5) = (-(params(2)*(y(4)-y(5))));
rp(3, 1) = (-(y(4)-params(9)+y(2)*1/params(3)));
rp(3, 3) = (-(y(2)*(-params(1))/(params(3)*params(3))));
rp(3, 9) = params(1);
rp(4, 2) = (-(params(4)*(-1)/(params(2)*params(2))));
rp(4, 4) = (-(x(it_, 1)+1/params(2)-1+y(5)*params(6)+y(4)*(1-params(6))-y(7)+0.9*y(1)));
rp(4, 6) = (-(y(3)-y(2)+params(4)*(y(5)-y(4))));
rp(5, 1) = (-y(3));
rp(5, 2) = (-(y(5)*params(6)+y(4)*(1-params(6))));
rp(5, 6) = (-(params(2)*(y(5)-y(4))));
rp(6, 1) = (-(y(5)-params(9)+y(3)*1/params(3)));
rp(6, 3) = (-(y(3)*(-params(1))/(params(3)*params(3))));
rp(6, 9) = params(1);
gp = zeros(7, 10, 9);
gp(1, 2, 5) = (-1);
gp(1, 3, 5) = 1;
gp(1, 4, 4) = (-params(5));
gp(1, 4, 5) = (-params(4));
gp(1, 5, 4) = (-(1-params(5)));
gp(1, 5, 5) = params(4);
gp(1, 6, 4) = 1;
gp(1, 9, 4) = (-1);
gp(1, 1, 4) = (-0.9);
gp(2, 2, 1) = (-1);
gp(2, 4, 2) = (-params(5));
gp(2, 4, 5) = (-params(2));
gp(2, 5, 2) = (-(1-params(5)));
gp(2, 5, 5) = params(2);
gp(3, 2, 1) = (-(1/params(3)));
gp(3, 2, 3) = (-((-params(1))/(params(3)*params(3))));
gp(3, 4, 1) = (-1);
gp(4, 2, 6) = 1;
gp(4, 3, 6) = (-1);
gp(4, 4, 4) = (-(1-params(6)));
gp(4, 4, 6) = params(4);
gp(4, 5, 4) = (-params(6));
gp(4, 5, 6) = (-params(4));
gp(4, 7, 4) = 1;
gp(4, 9, 4) = (-1);
gp(4, 1, 4) = (-0.9);
gp(5, 3, 1) = (-1);
gp(5, 4, 2) = (-(1-params(6)));
gp(5, 4, 6) = params(2);
gp(5, 5, 2) = (-params(6));
gp(5, 5, 6) = (-params(2));
gp(6, 3, 1) = (-(1/params(3)));
gp(6, 3, 3) = (-((-params(1))/(params(3)*params(3))));
gp(6, 5, 1) = (-1);
if nargout >= 3
rpp = zeros(14,4);
rpp(1,1)=1;
rpp(1,2)=2;
rpp(1,3)=2;
rpp(1,4)=(-(params(4)*(params(2)+params(2))/(params(2)*params(2)*params(2)*params(2))));
rpp(2,1)=1;
rpp(2,2)=2;
rpp(2,3)=4;
rpp(2,4)=(-((-1)/(params(2)*params(2))));
rpp(3,1)=1;
rpp(3,2)=4;
rpp(3,3)=2;
rpp(3,4)=rpp(2,4);
rpp(4,1)=1;
rpp(4,2)=4;
rpp(4,3)=5;
rpp(4,4)=(-(y(4)-y(5)));
rpp(5,1)=1;
rpp(5,2)=5;
rpp(5,3)=4;
rpp(5,4)=rpp(4,4);
rpp(6,1)=2;
rpp(6,2)=2;
rpp(6,3)=5;
rpp(6,4)=(-(y(4)-y(5)));
rpp(7,1)=2;
rpp(7,2)=5;
rpp(7,3)=2;
rpp(7,4)=rpp(6,4);
rpp(8,1)=3;
rpp(8,2)=1;
rpp(8,3)=3;
rpp(8,4)=(-(y(2)*(-1)/(params(3)*params(3))));
rpp(9,1)=3;
rpp(9,2)=3;
rpp(9,3)=1;
rpp(9,4)=rpp(8,4);
rpp(10,1)=3;
rpp(10,2)=1;
rpp(10,3)=9;
rpp(10,4)=1;
rpp(11,1)=3;
rpp(11,2)=9;
rpp(11,3)=1;
rpp(11,4)=rpp(10,4);
rpp(12,1)=3;
rpp(12,2)=3;
rpp(12,3)=3;
rpp(12,4)=(-(y(2)*T(1)));
rpp(13,1)=4;
rpp(13,2)=2;
rpp(13,3)=2;
rpp(13,4)=(-(params(4)*(params(2)+params(2))/(params(2)*params(2)*params(2)*params(2))));
rpp(14,1)=4;
rpp(14,2)=2;
rpp(14,3)=4;
rpp(14,4)=(-((-1)/(params(2)*params(2))));
rpp(15,1)=4;
rpp(15,2)=4;
rpp(15,3)=2;
rpp(15,4)=rpp(14,4);
rpp(16,1)=4;
rpp(16,2)=4;
rpp(16,3)=6;
rpp(16,4)=(-(y(5)-y(4)));
rpp(17,1)=4;
rpp(17,2)=6;
rpp(17,3)=4;
rpp(17,4)=rpp(16,4);
rpp(18,1)=5;
rpp(18,2)=2;
rpp(18,3)=6;
rpp(18,4)=(-(y(5)-y(4)));
rpp(19,1)=5;
rpp(19,2)=6;
rpp(19,3)=2;
rpp(19,4)=rpp(18,4);
rpp(20,1)=6;
rpp(20,2)=1;
rpp(20,3)=3;
rpp(20,4)=(-(y(3)*(-1)/(params(3)*params(3))));
rpp(21,1)=6;
rpp(21,2)=3;
rpp(21,3)=1;
rpp(21,4)=rpp(20,4);
rpp(22,1)=6;
rpp(22,2)=1;
rpp(22,3)=9;
rpp(22,4)=1;
rpp(23,1)=6;
rpp(23,2)=9;
rpp(23,3)=1;
rpp(23,4)=rpp(22,4);
rpp(24,1)=6;
rpp(24,2)=3;
rpp(24,3)=3;
rpp(24,4)=(-(y(3)*T(1)));
gpp = zeros(12,5);
gpp(1,1)=1;
gpp(1,2)=4;
gpp(1,3)=4;
gpp(1,4)=5;
gpp(1,5)=(-1);
gpp(2,1)=1;
gpp(2,2)=4;
gpp(2,3)=5;
gpp(2,4)=4;
gpp(2,5)=gpp(1,5);
gpp(3,1)=1;
gpp(3,2)=5;
gpp(3,3)=4;
gpp(3,4)=5;
gpp(3,5)=1;
gpp(4,1)=1;
gpp(4,2)=5;
gpp(4,3)=5;
gpp(4,4)=4;
gpp(4,5)=gpp(3,5);
gpp(5,1)=2;
gpp(5,2)=4;
gpp(5,3)=2;
gpp(5,4)=5;
gpp(5,5)=(-1);
gpp(6,1)=2;
gpp(6,2)=4;
gpp(6,3)=5;
gpp(6,4)=2;
gpp(6,5)=gpp(5,5);
gpp(7,1)=2;
gpp(7,2)=5;
gpp(7,3)=2;
gpp(7,4)=5;
gpp(7,5)=1;
gpp(8,1)=2;
gpp(8,2)=5;
gpp(8,3)=5;
gpp(8,4)=2;
gpp(8,5)=gpp(7,5);
gpp(9,1)=3;
gpp(9,2)=2;
gpp(9,3)=1;
gpp(9,4)=3;
gpp(9,5)=(-((-1)/(params(3)*params(3))));
gpp(10,1)=3;
gpp(10,2)=2;
gpp(10,3)=3;
gpp(10,4)=1;
gpp(10,5)=gpp(9,5);
gpp(11,1)=3;
gpp(11,2)=2;
gpp(11,3)=3;
gpp(11,4)=3;
gpp(11,5)=(-T(1));
gpp(12,1)=4;
gpp(12,2)=4;
gpp(12,3)=4;
gpp(12,4)=6;
gpp(12,5)=1;
gpp(13,1)=4;
gpp(13,2)=4;
gpp(13,3)=6;
gpp(13,4)=4;
gpp(13,5)=gpp(12,5);
gpp(14,1)=4;
gpp(14,2)=5;
gpp(14,3)=4;
gpp(14,4)=6;
gpp(14,5)=(-1);
gpp(15,1)=4;
gpp(15,2)=5;
gpp(15,3)=6;
gpp(15,4)=4;
gpp(15,5)=gpp(14,5);
gpp(16,1)=5;
gpp(16,2)=4;
gpp(16,3)=2;
gpp(16,4)=6;
gpp(16,5)=1;
gpp(17,1)=5;
gpp(17,2)=4;
gpp(17,3)=6;
gpp(17,4)=2;
gpp(17,5)=gpp(16,5);
gpp(18,1)=5;
gpp(18,2)=5;
gpp(18,3)=2;
gpp(18,4)=6;
gpp(18,5)=(-1);
gpp(19,1)=5;
gpp(19,2)=5;
gpp(19,3)=6;
gpp(19,4)=2;
gpp(19,5)=gpp(18,5);
gpp(20,1)=6;
gpp(20,2)=3;
gpp(20,3)=1;
gpp(20,4)=3;
gpp(20,5)=(-((-1)/(params(3)*params(3))));
gpp(21,1)=6;
gpp(21,2)=3;
gpp(21,3)=3;
gpp(21,4)=1;
gpp(21,5)=gpp(20,5);
gpp(22,1)=6;
gpp(22,2)=3;
gpp(22,3)=3;
gpp(22,4)=3;
gpp(22,5)=(-T(1));
end
if nargout >= 5
hp = zeros(0,5);
end
if nargout >= 6
g3p = zeros(0,6);
end
end
