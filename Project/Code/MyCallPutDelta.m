function [V0, D, x]=MyCallPutDelta(S0,r,q,T,E,flag,alpha,beta,Nx,Nt)

% Delta values of European option with underlying satisfying CEV
% using Crank-Nicolson Method
% dS/S = (mu-q) dt + alpha * S^(1-beta) dW
% 
% Pre:
% S0: initial underlying price
% r: constant interest rate
% q: continuous dividend yield
% T: expiry of option
% K: strike price
% flag: 1 for call, other values for put
% alpha: CEV model parameter
% beta: CEV model parameter, beta=1 -> BSM model
% Nx: number of spatial grids
% Nt: number of time grids
% 
% Post:
% V0: initial option value
% D: option deltas for different underlying prices at different times,
%    (Nx-1)*(Nt+1) matrix
% x: spatial grid points

k=T/Nt;
% Need to generate a non-uniform grid x=(x0, x1, ... , xN)
% generate the grid with highest density about E
x=E*(1+nthroot([0:4/Nx:4]'-1,5).^7);
% Boundary values when t=0
if flag==1
    % For European Call
    TBound=max(x(2:Nx)-E,0);
    xBound1=zeros(Nt+1,1); % S=0: C=0
    xBound2=x(Nx)-E*exp(-r*[0:k:T]'); % S>>E: C=S-Ee^{-rt}
else
    % For European Put
    TBound=max(E-x(2:Nx),0);
    xBound1=E*exp(-r*[0:k:T]'); % S=0; P=Ee^{-rt}
    xBound2=zeros(Nt+1,1); % S>>E: P~0
end
% construct matrix equations
hi=x(2:Nx+1)-x(1:Nx); % width of grids of x, h(i)=x(i+1)-x(i)
ai=2./(hi(1:Nx-1).*(hi(1:Nx-1)+hi(2:Nx)));
bi=2./(hi(2:Nx).*(hi(1:Nx-1)+hi(2:Nx)));
% F1,F2,F3 represents the three diagonals of I+F
F1=0.5*k*(min(100,alpha*(x(2:Nx).^(1-beta))).^2).*(x(2:Nx).^2).*ai-k*(r-q)*x(2:Nx)./(hi(1:Nx-1)+hi(2:Nx));
F2=(2-r*k)*ones(Nx-1,1)-0.5*k*(min(100,alpha*(x(2:Nx).^(1-beta))).^2).*(x(2:Nx).^2).*(ai+bi);
F3=0.5*k*(min(100,alpha*(x(2:Nx).^(1-beta))).^2).*(x(2:Nx).^2).*bi+k*(r-q)*x(2:Nx)./(hi(1:Nx-1)+hi(2:Nx));
% B1,B2,B3 represents the three diagonals of I+B
B1=-F1;
B2=4*ones(Nx-1,1)-F2;
B3=-F3;
% The equation become (I+B)*V(i+1) = (I+F)*V(i) + ri
% LU factorization of (I+B)
[L,U]=TriDiLU(B2,B1,B3);
% i.e. L*U*V(i+1) = (I+F)*V(i) + ri
vi=TBound; % values of the options
D=zeros(Nx-1,Nt+1);
D(:,1)=([vi(2:Nx-1);xBound2(1)]-[xBound1(1);vi(1:Nx-2)])./(hi(1:Nx-1)+hi(2:Nx));
for i=1:Nt
    ri=zeros(Nx-1,1);
    ri(1)=F1(1)*(xBound1(i)+xBound1(i+1));
    ri(Nx-1)=F3(Nx-1)*(xBound2(i)+xBound2(i+1));
    % First solve for y in L*y = (I+F)*V(i) + ri
    y=LBidiSol(L,F2.*vi+[0;F1(2:Nx-1).*vi(1:Nx-2)]+[F3(1:Nx-2).*vi(2:Nx-1);0]+ri);
    % Next solve for vi in U*vi=y
    vi=UBidiSol(U,B3,y);
    % calculate delta
    D(:,i+1)=([vi(2:Nx-1);xBound2(i+1)]-[xBound1(i+1);vi(1:Nx-2)])./(hi(1:Nx-1)+hi(2:Nx));
end

% Calculate value of option at S=S0 with interpolation
if S0<0 || S0>x(Nx)
    disp('Please enter an initial stock price in the range [0,5.66E]');
    V0=-1;
else
    vf=[xBound1(Nt+1);vi;xBound2(Nt+1)];
    i=1;
    while S0>x(i)
        i=i+1;
    end
    if S0==x(i)
        V0=vf(i);
    else
        V0=((x(i)-S0)*vf(i-1)+(S0-x(i-1))*vf(i))/(x(i)-x(i-1));
    end
end