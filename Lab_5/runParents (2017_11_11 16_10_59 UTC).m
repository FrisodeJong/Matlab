x0=3;
tol=10^-12;
maxIt=25;

f = @(x) log(x*exp(x));
df = @(x) (1/x)+1;
a0 = -1/df(x0);
phi = @(x) staticParent(f, a0, x, 0, 1);
phiN = @(x) newtonParent(f, df, x, 0, 1);

for k=0:3
[root, flag, convHist] = aitkenExtrap(phi, x0, tol, maxIt,k);
[rootN, flagN, convHistN] = aitkenExtrap(phiN, x0, tol, maxIt,k);
for i=1:length(convHist)-2
    ConvOrder(i,k+1)=log(convHist(i+2)/convHist(i+1))/log(convHist(i+1)/convHist(i));
end
for i=1:length(convHistN)-2
    ConvOrderN(i,k+1)=log(convHistN(i+2)/convHistN(i+1))/log(convHistN(i+1)/convHistN(i));
end
figure
semilogy(convHist) 
hold on
semilogy(convHistN)
end


