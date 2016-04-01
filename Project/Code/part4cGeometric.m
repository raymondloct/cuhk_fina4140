S0 = 22151.06; % From part 1
L = 2;
w1 = 0.5;
w2 = 0.5;
sigma1 = 0.3;
sigma2 = 0.3;
rho = 0.5;
T = 1;
r = 0.05;

var = L.^-2 * (sigma1^2+sigma2^2+2*sigma1*sigma2*rho);
G = S0;
K = S0;
d1 = (log(G/K)+(r+var-0.5/L*(sigma1^2+sigma2^2))*T)/sqrt(var*T);
d2 = d1 - sqrt(var*T);
Expec = G*exp((r-0.5/L*(sigma1^2+sigma2^2))*T+0.5*var*T);
V = exp(-r*T)*Expec*normcdf(d1)-K*exp(-r*T)*normcdf(d2)