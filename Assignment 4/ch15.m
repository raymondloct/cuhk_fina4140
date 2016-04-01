%CH15      Program for Chapter 15
%
%   Monte Carlo for a European put

randn('state',100)

%%%%%%%%%%%%%%%%% Problem and method parameters %%%%%%%%%%%%%%%
S = 4; E = 5; sigma = 0.3; r = 0.04; T = 1; 
Dt = 1e-3; N = T/Dt; M = 1e4;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

V = zeros(M,1);
for i = 1:M
    Sfinal = S*exp((r-0.5*sigma^2)*T+sigma*sqrt(T)*randn);
    V(i) = exp(-r*T)*max(E-Sfinal,0);
end
aM = mean(V); bM = std(V);
conf = [aM - 1.96*bM/sqrt(M), aM + 1.96*bM/sqrt(M)]
