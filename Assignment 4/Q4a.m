clear;
randn('state',100);

% Input parameters
S = input('Input initial stock price: ');
r = input('Input interest rate: ');
sigma = input('Input volatility: ');
K = input('Input strike price: ');
T = input('Input time to maturity: ');
N = input('Input number of simulations of the payoff: ');

Z = randn(N,1);
Svals = S .* exp((r - 0.5 .* sigma.^2) .* T + sigma .* sqrt(T) .* Z);
Santi = S .* exp((r - 0.5 .* sigma.^2) .* T - sigma .* sqrt(T) .* Z);
Pvals = exp(-r*T).*max(Svals-K,0);
Panti = exp(-r*T).*max(Santi-K,0);

% Output
Pvals = (Pvals + Panti)./2;
Pmean = mean(Pvals)
width = 1.96*std(Pvals)/sqrt(N);
conf = [Pmean - width, Pmean + width]