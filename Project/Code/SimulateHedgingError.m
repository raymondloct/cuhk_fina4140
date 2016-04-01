clear;

% problem and method parameters
mu = 0.15; q = 0; r = 0.04; S0 = 100;
alpha = 20; beta = 2;
T = 1; M = 10000; dt = T/M;
trials = 10000;
Nx=1600; Nt=Nx;

[V0, D, x]=MyCallPutDelta(S0,r,q,T,S0,-1,alpha,beta,Nx,Nt);

% simulate stock prices
S=zeros(M+1,trials);
S(1,:)=S0*ones(1, trials);
delta = zeros(M+1,trials);
sigma = min(100, alpha*S0^(1-beta));
delta(1,:) = getDelta(S0,T,T,x,D)*ones(1,trials); % Initial value of delta
HedgeError = zeros(1,trials);

for run = 1:trials
    for i = 1:M
        tau = (M-i)*dt;
        sigma = min(100, alpha*S(i,run)^(1-beta));
        dS = S(i,run) * ((mu-q)*dt + sigma * sqrt(dt) * randn);
        S(i+1,run) = max(S(i,run) + dS,0);
        delta(i+1,run)=getDelta(S(i+1,run),tau,T,x,D);
        HedgeError(run) = HedgeError(run) - delta(i,run)*(exp(-r*i*dt)*S(i+1,run)-exp(-r*(i-1)*dt)*S(i,run));
    end
    HedgeError(run) = HedgeError(run) + exp(-r*T)*max(0,S0-S(M+1,run)) - V0;    
end

% output
meanE=mean(HedgeError);
stdE=std(HedgeError);
upE=quantile(HedgeError,0.95);
downE=quantile(HedgeError,0.05);

disp('Mean of hedging error is')
disp(meanE)
disp('Standard deviation of hedging error is')
disp(stdE)
disp('95% quantile of hedging error is')
disp(upE)
disp('95% shortfall of hedging error is')
disp(downE)

hist(HedgeError,20)
