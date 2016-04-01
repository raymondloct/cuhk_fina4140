function [C, Cdelta, Cvega, P, Pdelta, Pvega] = ch10(S,E,r,sigma,tau)
% Program for Chapter 10
% This is a MATLAB function
%
% Input arguments: S = asset price at time t
%                  E = Exercise price
%                  r = interest rate
%                  sigma = volatility
%                  tau = time to expiry (T-t) 
%
% Output arguments: C = call value, Cdelta = delta value of call 
%                   Cvega = vega value of call
%                   P = Put value, Pdelta = delta value of put
%                   Pvega = vega value of put
%
%   function [C, Cdelta, Cvega, P, Pdelta, Pvega] = ch10(S,E,r,sigma,tau)

if tau > 0
   d1 = (log(S/E) + (r + 0.5*sigma^2)*(tau))/(sigma*sqrt(tau));
   d2 = d1 - sigma*sqrt(tau);
   N1 = 0.5*(1+erf(d1/sqrt(2)));
   N2 = 0.5*(1+erf(d2/sqrt(2)));
   C = S*N1-E*exp(-r*(tau))*N2;
   Cdelta = N1;
   Cvega = S*sqrt(tau)*exp(-0.5*d1^2)/sqrt(2*pi);
   P = C + E*exp(-r*tau) - S;
   Pdelta = Cdelta - 1;
   Pvega = Cvega;
else 
   C = max(S-E,0);
   Cdelta =  0.5*(sign(S-E) + 1);
   Cvega = 0;
   P = max(E-S,0);
   Pdelta = Cdelta - 1;
   Pvega = 0;
end
