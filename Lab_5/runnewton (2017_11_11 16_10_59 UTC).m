tol=10^-10;
maxIt=25;
n=3;
x0=eye(n);
B=[2 1 0; 0 2 1; 0 0 2];
dfFunc = @(x) matJac(x, 'exact', [], B);
vecF = @(x) matFunc(x, B);
[root, flag, convHist] = newton(vecF, dfFunc, x0, tol,...
    maxIt)

dfFunc = @(x) matJac(x, 'approx', 10^-9, B);
[rootapp9, flag, convHistapp9] = newton(vecF, dfFunc, x0, tol,...
    maxIt)

dfFunc = @(x) matJac(x, 'approx', 10^-1, B);
[rootapp1, flag, convHistapp1] = newton(vecF, dfFunc, x0, tol,...
    maxIt)
figure(1)
semilogy(convHist)
hold on
semilogy((convHistapp9),'o')
semilogy((convHistapp1),'.-.')