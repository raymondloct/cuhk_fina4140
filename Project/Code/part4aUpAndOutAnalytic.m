function V = part4aUpAndOutAnalytic(S)
% problem and method parameters
if ~exist('S','var')
    S = 100;
end
Su = 120; K = 90; r = 0.05; sigma = 0.2; T = 1;

% Analytic solution
if S >= Su
    V = 0
else
    d1 = (log(S./K) + (r + 0.5*sigma^2)*T)/(sigma*sqrt(T));
    d2 = (log(S./K) + (r - 0.5*sigma^2)*T)/(sigma*sqrt(T));
    d3 = (log(S./Su) + (r + 0.5*sigma^2)*T)/(sigma*sqrt(T));
    d4 = (log(S./Su) + (r - 0.5*sigma^2)*T)/(sigma*sqrt(T));
    d5 = (log(S./Su) - (r - 0.5*sigma^2)*T)/(sigma*sqrt(T));
    d6 = (log(S./Su) - (r + 0.5*sigma^2)*T)/(sigma*sqrt(T));
    d7 = (log(S.*K./Su.^2) - (r - 0.5*sigma^2)*T)/(sigma*sqrt(T));
    d8 = (log(S.*K./Su.^2) - (r + 0.5*sigma^2)*T)/(sigma*sqrt(T));
    N1 = normcdf(d1);
    N2 = normcdf(d2);
    N3 = normcdf(d3);
    N4 = normcdf(d4);
    N5 = normcdf(d5);
    N6 = normcdf(d6);
    N7 = normcdf(d7);
    N8 = normcdf(d8);
    V = S.*(N1-N3-(Su./S).^(1+2.*r./sigma.^2).*(N6-N8))-K.*exp(-r.*T).*(N2-N4-(Su./S).^(-1+2.*r/sigma.^2).*(N5-N7))
end
