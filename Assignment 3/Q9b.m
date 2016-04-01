clf

dsig = 0.25;
dx = 0.5;
[X,LAMBDA] = meshgrid(0:dx:20,1:dsig:5);
Y = LAMBDA.*exp(-LAMBDA.*X);
waterfall(X,LAMBDA,Y)
xlabel('x')
ylabel('\lambda')
zlabel('f(x)')
title('exp(\lambda) density for various \lambda')