clear;
randn('state',100);

% Input parameters
S = 100; r = 0.03; sigma = 0.25;
[K,T,power] = meshgrid([80,100,120],[0.5,1],2:1:5);
K = K(:);
T = T(:);
power = power(:);
N = 2.^[5:1:16];

% Calculate analytical value
APCanalytic = zeros(1,length(T));
for j = 1:length(T)
    Z = randn(100000,1);
    Svals = S .* exp((r - 0.5 .* sigma.^2) .* T(j) + sigma .* sqrt(T(j)) .* Z);
    Santi = S .* exp((r - 0.5 .* sigma.^2) .* T(j) - sigma .* sqrt(T(j)) .* Z);
    APCvals = exp(-r.*T(j)).*max(Svals.^power(j)-K(j).^power(j),0);
    APCanti = exp(-r.*T(j)).*max(Santi.^power(j)-K(j).^power(j),0);
    APCvals = (APCvals + APCanti)./2;
    APCanalytic(j) = mean(APCvals);
end

% Monte Carlo method using antithetic variates
APCmean = zeros(length(N),length(T));
for i = 1:length(N)
    for j = 1:length(T)
        Z = randn(N(i),1);
        Svals = S .* exp((r - 0.5 .* sigma.^2) .* T(j) + sigma .* sqrt(T(j)) .* Z);
        Santi = S .* exp((r - 0.5 .* sigma.^2) .* T(j) - sigma .* sqrt(T(j)) .* Z);
        APCvals = exp(-r.*T(j)).*max(Svals.^power(j)-K(j).^power(j),0);
        APCanti = exp(-r.*T(j)).*max(Santi.^power(j)-K(j).^power(j),0);
        APCvals = (APCvals + APCanti)./2;
        APCmean(i,j) = mean(APCvals);
    end
end

% Output
disp('Analytical value for Asymmetric Power Call (using 100000 steps):')
disp(APCanalytic)
disp('Computed price for Asymmetric Power Call using 2^16 steps:')
disp(APCmean(length(N),:))
APCerror = APCmean;
for i = 1:length(N)
    for j = 1:length(T)
        APCerror(i,j) = APCerror(i,j) - APCanalytic(j);
    end
end
APCerror = abs(APCerror);
for i = 1:length(T)
    subplot(4,6,i)
    semilogx(N,APCerror(:,i))
    xlabel('Number of MC steps')
    ylabel('Absolute error')
    xlim([2^5,2^16])
end