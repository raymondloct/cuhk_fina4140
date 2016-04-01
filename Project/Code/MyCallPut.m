function [V0, V]=MyCallPut(S0,r,q,T,E,flag)

% Price of European option with underlying satisfying CEV model
% using Crank-Nicolson Method
% dS/S = (mu-q) * dt + alpha * S^(1-beta) * dW
% 
% Pre:
% S0: initial underlying price
% r: constant interest rate
% q: continuous dividend yield
% T: expiry of option
% K: strike price
% flag: 1 for call, other values for put
% 
% Post:
% V0: initial option value
% V: (optional) option values for different underlying prices at different times, as an (Nx+1)*(Nt+1) matrix
% 
% Internal setting:
% alpha: 20
% beta: 2
% Nx=1600
% Nt=1600

alpha=20;
beta=2;
Nx=1600;
Nt=1600;

if nargout>1
    [V0,V]=MyCallPutExp(S0,r,q,T,E,flag,alpha,beta,Nx,Nt);
else
    V0=MyCallPutExp(S0,r,q,T,E,flag,alpha,beta,Nx,Nt);
end