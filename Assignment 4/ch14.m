%CH14     Program for Chapter 14
%
% Computes implied volatility for a European call

%%%%%%%%%%%%%%%%%%%%%%% parameters %%%%%%%%%%%%%%%%%%%%%%
r = 0.03; S = 2; E = 2; T = 3; tau = T; sigma_true = 0.3;
[C_true, Cdelta, P, Pdelta] = ch08(S,E,r,sigma_true,tau);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%starting value
sigmahat = sqrt(2*abs( (log(S/E) + r*T)/T ) );

%%%%%%%%%%%%%% Newton's Method %%%%%%%%%%%
tol = 1e-8;
sigma = sigmahat;
sigmadiff = 1;
k = 1;
kmax = 100;
while (sigmadiff >= tol & k < kmax)
    [C, Cdelta, Cvega, P, Pdelta, Pvega] = ch10(S,E,r,sigma,tau);
    increment = (C-C_true)/Cvega;
    sigma = sigma - increment;
    k = k+1;
    sigmadiff = abs(increment);
end
sigma
