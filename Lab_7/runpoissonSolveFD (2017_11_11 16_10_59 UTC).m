f=@(x) pi^2*sin(pi*x);
g=@(x) 0;
solExact=@(x) sin(pi*x);
solExactArray=0;
tol=10^-10;
for i = 2:13
    N=2^i;
    Narray(i)=N;
    h=1/N;
    [sol, xVals]=poissonSolveFD(f, g, N, 'gs', tol, h);
    harray(i)=h;
    for k=1:N+1
        solExactArray(k)=solExact(xVals(k));
    end
    error(i)=norm(max(sol-solExactArray));
end
loglog(harray,error,harray,harray.^2)
legend('maxerror','h against h^2')