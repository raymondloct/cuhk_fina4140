M = dlmread('Q13.txt');
n = M(1:5,1:1);
box = M(1:5,2:2);
polar = M(1:5,3:3);
plot(n,box,'x-',n,polar,'x-')
xlabel('n')
ylabel('Time required (s)')
legend('Box-Muller method','Marsaglia polar method')