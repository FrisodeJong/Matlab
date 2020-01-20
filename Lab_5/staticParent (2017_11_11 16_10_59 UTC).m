% INPUT
% f         function of rootfinding problem
% a0        static parameter: xnew = x + a0*f(x)
% x0        initial guess
% tol       desired tolerance
% maxIt     maximum number of iterations
% OUTPUT
% root      root of f
% flag      if 0: attained desired tolerance
%           if 1: reached maxIt nr of iterations
% convHist  convergence history
function [root, flag, convHist] = staticParent(f, a0, x0, tol,...
    maxIt)

convHist=zeros(maxIt,1);
x=x0;
flag=1;
for k=1:maxIt
    xnew=x+a0*f(x);
    convHist(k)=norm(xnew-x);
    if norm(x-xnew)< tol
        flag=0;
        break
    end
    x=xnew;
end
root=x;