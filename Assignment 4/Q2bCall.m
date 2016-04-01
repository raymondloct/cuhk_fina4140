randn('state',100);

%%%%%%%%%%%%%%%%% Problem and method parameters %%%%%%%%%%%%%%%
S = 4; E = 5; sigma = 0.3; r = 0.04; T = 1; 
Dt = 1e-3; N = T/Dt;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
AnalyticVal = ch08(S,E,r,sigma,T);

start = 5;
finish = 17;
M = (2.^[start:finish])';
aM = zeros(length(M),1);
bM = zeros(length(M),1);
for i = 1:length(M)
    V = zeros(M(i),1);
    for j = 1:M(i)
        Sfinal = S*exp((r-0.5*sigma^2)*T+sigma*sqrt(T)*randn);
        V(j) = exp(-r*T)*max(Sfinal-E,0);
    end
    aM(i) = mean(V);
    bM(i) = 1.96*std(V)/sqrt(M(i));
end
fig = figure;
errorbar(M,aM,bM,'x');
ax = get(fig,'CurrentAxes');
set(ax,'XScale','log','YScale','log');
xlabel('Num samples')
ylabel('Option value approximation')
axis([10^1,10^6,10^-2,10^0])
hold on
plot([10^1 10^6],[AnalyticVal AnalyticVal],':') % Add the dashed line for analytical value
hold off
