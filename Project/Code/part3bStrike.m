c=zeros(11,1);
p=zeros(11,1);
for i=0:10
    c(i+1)=MyCallPutExp(100,0.04,0.02,1,90+2*i,1,20,2,1600,1600);
    p(i+1)=MyCallPutExp(100,0.04,0.02,1,90+2*i,-1,20,2,1600,1600);
end