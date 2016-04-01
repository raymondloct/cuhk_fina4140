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

% Monte Carlo method using antithetic variates
Pmean = zeros(length(N),length(T));
Cmean = zeros(length(N),length(T));
for i = 1:length(N)
    for j = 1:length(T)
        Z = randn(N(i),1);
        Svals = S .* exp((r - 0.5 .* sigma.^2) .* T(j) + sigma .* sqrt(T(j)) .* Z);
        Santi = S .* exp((r - 0.5 .* sigma.^2) .* T(j) - sigma .* sqrt(T(j)) .* Z);
        Pvals = exp(-r.*T(j)).*max(K(j)-Svals,0);
        Panti = exp(-r.*T(j)).*max(K(j)-Santi,0);
        Cvals = exp(-r.*T(j)).*max(Svals-K(j),0);
        Canti = exp(-r.*T(j)).*max(Santi-K(j),0);
        Pvals = (Pvals + Panti)./2;
        Pmean(i,j) = mean(Pvals);
        Cvals = (Cvals + Canti)./2;
        Cmean(i,j) = mean(Cvals);
    end
end

% Output
disp('Analytical value for European call:')
disp(Canalytic')
disp('Computed price for European call using 2^16 steps:')
disp(Cmean(length(N),:))
disp('Analytical value for European put:')
disp(Panalytic')
disp('Computed price for European put using 2^16 steps:')
disp(Pmean(length(N),:))
Cerror = Cmean;
Perror = Pmean;
for i = 1:length(N)
    for j = 1:length(T)
        Cerror(i,j) = Cerror(i,j) - Canalytic(j);
        Perror(i,j) = Perror(i,j) - Panalytic(j);
    end
end
Cerror = abs(Cerror);
Perror = abs(Perror);
for i = 1:length(T)
    subplot(2,3,i)
    semilogx(N,Cerror(:,i),N,Perror(:,i))
    legend('European Call','European Put')
    xlabel('Number of MC steps')
    ylabel('Absolute error')
    xlim([2^5,2^16])
end