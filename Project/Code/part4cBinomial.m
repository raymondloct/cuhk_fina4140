function Vresult = part4cBinomial(M)
% problem and method parameters
S0 = 22151.06; % From part 1.1
w = [0.5,0.5];
K = S0;
sigma = [0.3 0.3];
rho = 0.5;
T = 1;
r = 0.05;
if ~exist('M','var')
    M = 100;
end
dt = T/M;

% Binomial model
V = @(S1,S2) max(0,S1*w(1)+S2*w(2)-K);
Vbin = zeros(M+1,M+1);
u1 = exp(sigma(1)*sqrt(dt));
d1 = 1/u1;
u2 = exp(sigma(2)*sqrt(dt));
d2 = 1/u2;
p1 = 0.25*(1+rho+((r-0.5*sigma(1)^2)/sigma(1)+(r-0.5*sigma(2)^2)/sigma(2))*sqrt(dt));
p2 = 0.25*(1-rho+((r-0.5*sigma(1)^2)/sigma(1)-(r-0.5*sigma(2)^2)/sigma(2))*sqrt(dt));
p3 = 0.25*(1-rho+(-(r-0.5*sigma(1)^2)/sigma(1)+(r-0.5*sigma(2)^2)/sigma(2))*sqrt(dt));
p4 = 0.25*(1+rho-((r-0.5*sigma(1)^2)/sigma(1)+(r-0.5*sigma(2)^2)/sigma(2))*sqrt(dt));
for i = 1:M+1
    for j = 1:M+1
        Vbin(i,j) = V(S0*u1^(i-1)*d1^(M-i+1),S0*u2^(j-1)*d2^(M-j+1));
    end
end

for time = M:-1:1
    for i = 1:time
        for j = 1:time
            Vbin(i,j) = exp(-r*dt)*(p1*Vbin(i+1,j+1)+p2*Vbin(i+1,j)+p3*Vbin(i,j+1)+p4*Vbin(i,j));
        end
    end
end

Vresult = Vbin(1,1);
disp('Number of lattice steps is')
disp(M)
disp('Basket call value is')
disp(Vresult)