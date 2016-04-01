ImV=zeros(4,5);
for i=1:4
    for j=1:5
        E=70+10*j;
        T=0.25*i;
        p=MyCallPutExp(100,0.04,0.02,T,E,-1,20,2,1600,1600);
        ImV(i,j)=blsimpv(100,E,0.04,T,p,10,0.02,1e-6,false);
    end
end