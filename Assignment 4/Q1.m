% Input parameters
otype = input('Input option type (Enter c for call, p for put): ','s');
r = input('Input interest rate: ');
S = input('Input initial stock price: ');
T = input('Input time to maturity: ');
E = input('Input strike price: ');
otrue = input('Input observed option price: ');
N = input('Input number of steps to be carried out: ');

% Initialization
sigmahat = sqrt(2*abs((log(S/E) + r*T)/T));
tol = 1e-8;
sigma = sigmahat;
k = 0;

% Newton's method
if strncmp(otype, 'c', 1) == 1 || strncmp(otype, 'C', 1) == 1
    while k < N
        [C, Cdelta, Cvega, P, Pdelta, Pvega] = ch10(S,E,r,sigma,T);
        increment = (C-otrue)/Cvega;
        sigma = sigma - increment;
        k = k+1;
    end
elseif strncmp(otype, 'p', 1) == 1 || strncmp(otype, 'P', 1) == 1
    while k < N
        [C, Cdelta, Cvega, P, Pdelta, Pvega] = ch10(S,E,r,sigma,T);
        increment = (P-otrue)/Pvega;
        sigma = sigma - increment;
        k = k+1;
    end
end
sigma