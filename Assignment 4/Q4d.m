clear;
randn('state',100);

% Input parameters
S = 100; r = 0.03; sigma = 0.25;
N = 2.^[5:1:16];
[K,T] = meshgrid([80,100,120],[0.5,1]);
K = K(:);
T = T(:);

% Calculate analytical value
d1 = (log(S./K) + (r + 0.5*sigma^2).*T)./(sigma.*sqrt(T));
d2 = d1 - sigma.*sqrt(T);
N1 = 0.5.*(1+erf(d1./sqrt(2)));
N2 = 0.5.*(1+erf(d2./sqrt(2)));
Canalytic = S.*N1-K.*exp(-r.*T).*N2;
Panalytic = Canalytic + K.*exp(-r.*T) - S;
STanalytic = Canalytic + Panalytic;

% Monte Carlo method using antithetic variates
STmean = zeros(length(N),length(T));
for i = 1:length(N)
    for j = 1:length(T)
        Z = randn(N(i),1);
        Svals = S .* exp((r - 0.5 .* sigma.^2) .* T(j) + sigma .* sqrt(T(j)) .* Z);
        Santi = S .* exp((r - 0.5 .* sigma.^2) .* T(j) - sigma .* sqrt(T(j)) .* Z);
        STvals = exp(-r.*T(j)).*abs(Svals-K(j));
        STanti = exp(-r.*T(j)).*abs(Svals-K(j));
        STvals = (STvals + STanti)./2;
        STmean(i,j) = mean(STvals);
    end
end

% Output
disp('Analytical value for Straddle:')
disp(STanalytic')
disp('Computed price for Straddle using 2^16 steps:')
disp(STmean(length(N),:))
STerror = STmean;
for i = 1:length(N)
    for j = 1:length(T)
        STerror(i,j) = STerror(i,j) - STanalytic(j);
    end
end
STerror = abs(STerror);
for i = 1:length(T)
    subplot(2,3,i)
    semilogx(N,STerror(:,i))
    xlabel('Number of MC steps')
    ylabel('Absolute error')
    xlim([2^5,2^16])
end