A = xlsread('q7data.xls');
M = A(1:51,1:1);
error = A(1:51,2:2);
plot(M,error,'x')
xlabel('Mesh size')
ylabel('Error')
title('Binomial Model Estimation Error')