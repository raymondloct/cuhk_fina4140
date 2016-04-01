clf;

%%%%%%% Problem and method parameters %%%%%%%
E = 4; sigma = 0.3; r = 0.03; T = 1; bound = 2;
L = 10; Nx = 50; Nt = 50; k = T/Nt; h = (L-bound)/Nx;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

T1 = diag(ones(Nx-2,1),1) - diag(ones(Nx-2,1),-1);
T2 = -2*eye(Nx-1,Nx-1) + diag(ones(Nx-2,1),1) + diag(ones(Nx-2,1),-1);
mvec = bound / h + [1:Nx-1];
D1 = diag(mvec);
D2 = diag(mvec.^2);
F = (1-r*k)*eye(Nx-1,Nx-1) + 0.5*k*sigma^2*D2*T2 + 0.5*k*r*D1*T1;
B = (1+r*k)*eye(Nx-1,Nx-1) - 0.5*k*sigma^2*D2*T2 - 0.5*k*r*D1*T1;
A1 = 0.5*(eye(Nx-1,Nx-1) + F); % 0.5*(I+F)
A2 = 0.5*(eye(Nx-1,Nx-1) + B); % 0.5*(I+B)

U = zeros(Nx-1,Nt+1);
U(:,1) = max([bound+h:h:L-h]'-E,0); % At maturity

for i = 1:Nt
    tau = (i-1)*k;
    c1 = ch08(L,E,r,sigma,tau);
    c2 = ch08(L,E,r,sigma,tau+k);
    p = zeros(Nx-1,1);
    p(Nx-1) = (k/2*sigma^2*(bound+(Nx-1)*h)^2/h^2-r*k*(bound+(Nx-1)*h)/(2*h))*c1;
    q = zeros(Nx-1,1);
    q(Nx-1) = (k/2*sigma^2*(bound+(Nx-1)*h)^2/h^2-r*k*(bound+(Nx-1)*h)/(2*h))*c2;
    rhs = A1*U(:,i) + 0.5*(p+q);
    X = A2\rhs;
    U(:,i+1) = X;
end

bca = zeros(1,Nt+1);
bcb = (L-E)*ones(1,Nt+1);
U = [bca;U;bcb];
mesh([0:k:T],[bound:h:L],U)
xlabel('T-t'), ylabel('S'), zlabel('Down-and-out Call Value')
