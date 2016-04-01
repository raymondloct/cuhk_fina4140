function [CallValSim,CallValAnti] = part4cMonte(N)
% problem and method parameters
L = 2;
S0 = ones(1,L).*22151.06; % From part 1.1
w = [0.5, 0.5];
K = sum(S0.*w);
sigma = [0.3, 0.3];
rho = [1 0.5; 0.5 1]; % Correlation between asset i and j
T = 1;
r = 0.05;
if ~exist('N','var')
    N = 100000;
end

covar = zeros(L,L);
% Calculate covariance matrix
for i = 1:L
  for j = 1:L
    covar(i,j) = sigma(i)*sigma(j)*rho(i,j);
  end
end

% Cholesky decomposition of covariance matrix
covarL = zeros(L,L);
for k = 1:L
  temp = 0;
  for s = 1:k-1
    temp = temp + covarL(k,s)^2;
  end
  covarL(k,k) = sqrt(covar(k,k) - temp);
  for i = k+1:L
    temp = 0;
    for s = 1:k-1
      temp = temp + covarL(i,s) * covarL(k,s);
    end
    covarL(i,k) = (covar(i,k) - temp) / covarL(k,k);
  end
end

% Monte Carlo simulation
Vsim = 0;
Vanti = 0;
temp1 = (r-0.5*sigma(1)^2)*T;
temp2 = (r-0.5*sigma(2)^2)*T;
for run = 1:N
    X = (covarL*randn(L,1))';
    S1sim = S0(1)*exp(temp1 + sqrt(T)*X(1));
    S1anti = S0(1)*exp(temp1 - sqrt(T)*X(1));
    S2sim = S0(2)*exp(temp2 + sqrt(T)*X(2));
    S2anti = S0(2)*exp(temp2 - sqrt(T)*X(2));
    Vsim = Vsim + max(S1sim*w(1)+S2sim*w(2)-K,0);
    Vanti = Vanti + max(S1anti*w(1)+S2anti*w(2)-K,0);
end
CallValSim = (Vsim/N)*exp(-r*T);
CallValAnti = (Vsim + Vanti)/(2*N)*exp(-r*T);

% Results
disp('Number of simulations for Monte Carlo simulation is');
disp(N);
disp('Estimate of basket call value by Monte Carlo simulation without antithetic variate is');
disp(CallValSim);
disp('Estimate of basket call value by Monte Carlo simulation with antithetic variate is');
disp(CallValAnti);